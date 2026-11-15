"""ChEMBL database access.

Two things live here, both shared by the inspection and extraction scripts:

1. A read-only connection. The dump is opened immutable so no code path can
   modify it -- the downloaded file with its recorded checksum is the data of
   record.
2. A temporary target -> family table. Family assignment is a prefix match on
   the flattened class path, and joining that match into every profile query
   turns a five-table join into a scan of 20M activities per family. Resolving
   it once into an indexed temp table takes seconds and makes the rest fast.

How the family mapping works, and why it is a prefix match on `short_name`
paths, is documented in config/families.yaml.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

# Activity filters shared by every query. Kept as SQL fragments rather than
# rebuilt per call site so that the inspection counts and the extraction counts
# cannot drift apart.
ACTIVITY_FILTER_SQL = """
    act.standard_type IN ('Ki', 'Kd', 'IC50')
    AND act.standard_units = 'nM'
    AND act.standard_relation = '='
    AND act.pchembl_value IS NOT NULL
    AND act.data_validity_comment IS NULL
"""

ASSAY_FILTER_SQL = """
    a.confidence_score >= :min_confidence
"""

TARGET_FILTER_SQL = """
    td.target_type = :target_type
    AND td.organism = :organism
"""


def connect(db_path: Path | str) -> sqlite3.Connection:
    """Open the dump read-only and immutable."""
    db_path = Path(db_path)
    if not db_path.exists():
        raise SystemExit(
            f"{db_path} not found.\n"
            "Run: python scripts/00_download_chembl.py --version 37"
        )
    conn = sqlite3.connect(f"file:{db_path}?immutable=1", uri=True)
    # A generous cache: the profile queries scan large tables and this is a
    # read-only workload on a 30 GB file.
    conn.execute("PRAGMA cache_size = -512000")  # 512 MB
    conn.execute("PRAGMA temp_store = FILE")
    return conn


def build_family_target_table(
    conn: sqlite3.Connection, config, table: str = "fam_target"
) -> dict[str, int]:
    """Resolve target -> family into an indexed temp table.

    Returns the number of targets per family, plus an "ambiguous" entry for
    targets that matched more than one family and were therefore dropped.

    Ambiguity is resolved by dropping, never by picking a winner: a target
    counted in two families would leave its records in training when one of
    those families is held out, which would silently invalidate R4.
    """
    filters = {
        "target_type": config["filtering"]["target_type"],
        "organism": config["filtering"]["organism"],
    }

    conn.execute(f"DROP TABLE IF EXISTS {table}")
    conn.execute(
        f"CREATE TEMP TABLE {table} (tid INTEGER PRIMARY KEY, family TEXT NOT NULL)"
    )

    # One pass per family into a staging table, then resolve ambiguity in SQL.
    conn.execute("DROP TABLE IF EXISTS fam_target_raw")
    conn.execute(
        "CREATE TEMP TABLE fam_target_raw (tid INTEGER NOT NULL, family TEXT NOT NULL)"
    )

    for family in config.active_families:
        prefixes = config.family_def(family)["class_path_prefixes"]
        for prefix in prefixes:
            conn.execute(
                f"""
                INSERT INTO fam_target_raw (tid, family)
                SELECT DISTINCT td.tid, :family
                FROM protein_classification pc
                JOIN component_class   cc ON cc.protein_class_id = pc.protein_class_id
                JOIN target_components tc ON tc.component_id     = cc.component_id
                JOIN target_dictionary td ON td.tid              = tc.tid
                WHERE LOWER(pc.protein_class_desc) LIKE :prefix || '%'
                  AND {TARGET_FILTER_SQL}
                """,
                {"family": family, "prefix": prefix.lower(), **filters},
            )

    # A target may appear several times for one family (several matching
    # classes); that is not ambiguity. Ambiguity is >1 DISTINCT family.
    conn.execute(
        f"""
        INSERT INTO {table} (tid, family)
        SELECT tid, MIN(family)
        FROM fam_target_raw
        GROUP BY tid
        HAVING COUNT(DISTINCT family) = 1
        """
    )
    n_ambiguous = conn.execute(
        """
        SELECT COUNT(*) FROM (
            SELECT tid FROM fam_target_raw
            GROUP BY tid HAVING COUNT(DISTINCT family) > 1
        )
        """
    ).fetchone()[0]
    conn.execute("DROP TABLE fam_target_raw")

    counts = {
        family: n
        for family, n in conn.execute(
            f"SELECT family, COUNT(*) FROM {table} GROUP BY family"
        )
    }
    counts["ambiguous_dropped"] = n_ambiguous
    return counts


def family_activity_counts(
    conn: sqlite3.Connection, config, table: str = "fam_target"
) -> list[tuple[str, int, int, int]]:
    """(family, activities, distinct compounds, distinct targets) after filtering."""
    params = {
        "target_type": config["filtering"]["target_type"],
        "organism": config["filtering"]["organism"],
        "min_confidence": config["filtering"]["min_assay_confidence_score"],
    }
    return list(
        conn.execute(
            f"""
            SELECT ft.family,
                   COUNT(*)                        AS n_activities,
                   COUNT(DISTINCT act.molregno)    AS n_compounds,
                   COUNT(DISTINCT td.tid)          AS n_targets
            FROM {table} ft
            JOIN target_dictionary td ON td.tid = ft.tid
            JOIN assays     a   ON a.tid        = td.tid
            JOIN activities act ON act.assay_id = a.assay_id
            WHERE {TARGET_FILTER_SQL}
              AND {ASSAY_FILTER_SQL}
              AND {ACTIVITY_FILTER_SQL}
            GROUP BY ft.family
            ORDER BY n_activities DESC
            """,
            params,
        )
    )


def family_binding_domain_coverage(
    conn: sqlite3.Connection, config, table: str = "fam_target"
) -> list[tuple[str, int, int]]:
    """(family, filtered activities, how many have a predicted binding domain).

    This is the measurement the stage-0 checkpoint turns on: the binding domain
    is the primary protein representation, and if coverage is low in the very
    families being compared, the representation shrinks the dataset unevenly --
    which is the confounder it was adopted to avoid.
    """
    params = {
        "target_type": config["filtering"]["target_type"],
        "organism": config["filtering"]["organism"],
        "min_confidence": config["filtering"]["min_assay_confidence_score"],
    }
    return list(
        conn.execute(
            f"""
            SELECT ft.family,
                   COUNT(*) AS n_activities,
                   SUM(CASE WHEN pbd.activity_id IS NOT NULL THEN 1 ELSE 0 END)
                       AS n_mapped
            FROM {table} ft
            JOIN target_dictionary td ON td.tid = ft.tid
            JOIN assays     a   ON a.tid        = td.tid
            JOIN activities act ON act.assay_id = a.assay_id
            LEFT JOIN (
                SELECT DISTINCT activity_id FROM predicted_binding_domains
            ) pbd ON pbd.activity_id = act.activity_id
            WHERE {TARGET_FILTER_SQL}
              AND {ASSAY_FILTER_SQL}
              AND {ACTIVITY_FILTER_SQL}
            GROUP BY ft.family
            ORDER BY n_activities DESC
            """,
            params,
        )
    )


def family_sequence_lengths(
    conn: sqlite3.Connection, table: str = "fam_target"
) -> list[tuple[str, int, int, int, int]]:
    """(family, n targets with a sequence, min, median, max full-sequence length).

    Reported per family because the length distribution is exactly what makes
    truncation a confounder: families differ in it, so a single cutoff hits them
    at different rates.
    """
    rows = []
    for (family,) in conn.execute(f"SELECT DISTINCT family FROM {table}"):
        lengths = [
            r[0]
            for r in conn.execute(
                f"""
                SELECT LENGTH(cs.sequence)
                FROM {table} ft
                JOIN target_components tc ON tc.tid = ft.tid
                JOIN component_sequences cs ON cs.component_id = tc.component_id
                WHERE ft.family = ?
                  AND cs.sequence IS NOT NULL
                """,
                (family,),
            )
        ]
        if not lengths:
            rows.append((family, 0, 0, 0, 0))
            continue
        lengths.sort()
        rows.append(
            (
                family,
                len(lengths),
                lengths[0],
                lengths[len(lengths) // 2],
                lengths[-1],
            )
        )
    return rows
