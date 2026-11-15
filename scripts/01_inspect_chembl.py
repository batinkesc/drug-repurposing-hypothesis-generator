"""Inspect the downloaded ChEMBL database before writing any extraction query.

Purpose: verify the schema assumptions this project rests on, and print the
real protein-classification inventory so that `config/families.yaml` can be
mapped onto actual ChEMBL class names rather than guessed ones.

Three questions this answers:

1. Do the tables and columns the pipeline needs actually exist, with the names
   we assume? (ChEMBL schema changes between releases.)
2. What are the real protein class hierarchy strings? The design maps five
   families onto ChEMBL classes, but ChEMBL stores the hierarchy as a
   `parent_id` chain plus a flattened `protein_class_desc`, not as separate
   l1/l2 columns -- so the mapping has to be written against what is there.
3. How much of the data does `predicted_binding_domains` actually cover, and
   can a domain be cut out of a sequence? The binding-domain representation is
   the primary protein representation, so if `component_domains` carries no
   usable start/end positions the whole representation plan needs rethinking.

Writes a report to reports/tables/01_chembl_inventory.md and prints a summary.
Read-only: it opens the database in immutable mode.

Usage::

    python scripts/01_inspect_chembl.py
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dti.chembl import (  # noqa: E402
    build_family_target_table,
    connect,
    family_activity_counts,
    family_binding_domain_coverage,
    family_sequence_lengths,
)
from dti.config import REPO_ROOT, load_config  # noqa: E402

# Tables the pipeline depends on, with the reason each one is needed.
REQUIRED_TABLES = {
    "target_dictionary": "target records, target_type, organism",
    "target_components": "target -> protein component link",
    "component_sequences": "UniProt accession and amino-acid sequence",
    "component_class": "component -> protein class link",
    "protein_classification": "the family hierarchy used to build splits",
    "activities": "the measurements, including pchembl_value",
    "assays": "assay type and confidence_score",
    "molecule_dictionary": "compound records, max_phase for the retrieval test",
    "compound_structures": "canonical SMILES",
    "predicted_binding_domains": "PPDMs activity -> binding site mapping",
    "binding_sites": "binding site records",
    "site_components": "binding site -> component / domain link",
    "component_domains": "domain start/end positions within a sequence",
    "domains": "Pfam domain names",
}


def table_exists(conn: sqlite3.Connection, table: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
    ).fetchone()
    return row is not None


def columns(conn: sqlite3.Connection, table: str) -> list[str]:
    return [r[1] for r in conn.execute(f"PRAGMA table_info({table})")]


def count(conn: sqlite3.Connection, table: str) -> int:
    return conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]


def scalar(conn: sqlite3.Connection, sql: str, params: tuple = ()) -> Any:
    row = conn.execute(sql, params).fetchone()
    return row[0] if row else None


def section(lines: list[str], title_en: str, title_tr: str) -> None:
    lines.append("")
    lines.append(f"## {title_en} / {title_tr}")
    lines.append("")


def check_schema(conn: sqlite3.Connection, lines: list[str]) -> list[str]:
    """Report which required tables are present. Returns the missing ones."""
    section(lines, "Schema check", "Şema denetimi")
    lines.append("| Table | Rows | Needed for |")
    lines.append("| --- | --- | --- |")

    missing = []
    for table, reason in REQUIRED_TABLES.items():
        if table_exists(conn, table):
            lines.append(f"| `{table}` | {count(conn, table):,} | {reason} |")
        else:
            missing.append(table)
            lines.append(f"| `{table}` | **MISSING** | {reason} |")

    print(f"  schema: {len(REQUIRED_TABLES) - len(missing)}/{len(REQUIRED_TABLES)} "
          f"required tables present")
    if missing:
        print(f"  MISSING: {', '.join(missing)}")
    return missing


def report_class_hierarchy(conn: sqlite3.Connection, lines: list[str]) -> None:
    """Print the real protein class hierarchy, level by level.

    This is the output that `config/families.yaml` has to be written against.
    """
    section(lines, "Protein class hierarchy", "Protein sınıf hiyerarşisi")

    cols = columns(conn, "protein_classification")
    lines.append(f"`protein_classification` columns: {', '.join(cols)}")
    lines.append("")

    levels = [
        r[0]
        for r in conn.execute(
            "SELECT DISTINCT class_level FROM protein_classification "
            "ORDER BY class_level"
        )
    ]
    lines.append(f"Class levels present / Mevcut sınıf seviyeleri: {levels}")

    # Level 1 and 2 with the number of human single-protein targets under each.
    # Target counts matter more than class counts: a class with three targets
    # cannot carry a family.
    for level in (1, 2):
        lines.append("")
        lines.append(f"### class_level = {level}")
        lines.append("")
        lines.append("| pref_name | Classes | Human single-protein targets |")
        lines.append("| --- | --- | --- |")
        rows = conn.execute(
            """
            SELECT pc.pref_name,
                   COUNT(DISTINCT pc.protein_class_id) AS n_classes,
                   COUNT(DISTINCT td.tid)              AS n_targets
            FROM protein_classification pc
            LEFT JOIN component_class    cc ON cc.protein_class_id = pc.protein_class_id
            LEFT JOIN target_components  tc ON tc.component_id     = cc.component_id
            LEFT JOIN target_dictionary  td ON td.tid              = tc.tid
                 AND td.target_type = 'SINGLE PROTEIN'
                 AND td.organism    = 'Homo sapiens'
            WHERE pc.class_level = ?
            GROUP BY pc.pref_name
            ORDER BY n_targets DESC
            """,
            (level,),
        ).fetchall()
        for pref_name, n_classes, n_targets in rows:
            lines.append(f"| {pref_name} | {n_classes} | {n_targets:,} |")

    # The flattened hierarchy string, which is often the easiest thing to match
    # a family definition against.
    lines.append("")
    lines.append("### Sample `protein_class_desc` values")
    lines.append("")
    lines.append("```")
    for (desc,) in conn.execute(
        "SELECT DISTINCT protein_class_desc FROM protein_classification "
        "WHERE class_level = 2 ORDER BY protein_class_desc LIMIT 40"
    ):
        lines.append(desc)
    lines.append("```")
    print(f"  class hierarchy: levels {levels} written to the report")


def report_family_profile(
    conn: sqlite3.Connection, lines: list[str], config
) -> None:
    """Per-family targets, activities and compounds after the stage-2 filters.

    The family mapping is a prefix match on the flattened class path; see
    config/families.yaml for why, and `dti.chembl` for the implementation. A
    family whose count comes back zero has a wrong prefix, not an empty ChEMBL.
    """
    section(lines, "Family profile", "Aile profili")

    target_counts = build_family_target_table(conn, config)
    n_ambiguous = target_counts.pop("ambiguous_dropped", 0)

    lines.append(
        "Filters applied: single human protein targets, Ki/Kd/IC50 in nM with "
        "an exact relation, assay confidence >= "
        f"{config['filtering']['min_assay_confidence_score']}, "
        "no data-validity comment, pChEMBL present."
    )
    lines.append("")
    lines.append("| Family | Targets | Activities | Distinct compounds |")
    lines.append("| --- | --- | --- | --- |")

    rows = family_activity_counts(conn, config)
    total_activities = 0
    for family, n_activities, n_compounds, n_targets in rows:
        label = config.family_label(family, "en")
        total_activities += n_activities
        lines.append(
            f"| {label} | {n_targets:,} | {n_activities:,} | {n_compounds:,} |"
        )
        print(
            f"  {family:<18} targets={n_targets:>5,}  "
            f"activities={n_activities:>9,}  compounds={n_compounds:>8,}"
        )

    lines.append(f"| **Total** | {sum(target_counts.values()):,} "
                 f"| {total_activities:,} | |")
    lines.append("")
    lines.append(
        f"- Targets matching more than one family (dropped): **{n_ambiguous}**"
    )
    lines.append(
        "  Zero is the precondition for R4 meaning anything: a target counted "
        "in two families would leave its records in training when one of those "
        "families is held out."
    )
    print(f"  ambiguous targets dropped: {n_ambiguous}")


def report_family_lengths(
    conn: sqlite3.Connection, lines: list[str], config
) -> None:
    """Full-sequence length distribution per family.

    Without this table the cross-family comparison cannot be interpreted: the
    families differ in sequence length, so a single ESM-2 cutoff truncates them
    at different rates, and a performance gap could come from truncation rather
    than from transfer.
    """
    section(lines, "Sequence lengths per family",
            "Aile başına dizi uzunlukları")
    limit = config["protein"]["esm_max_residues"]
    lines.append(f"ESM-2 limit: {limit} residues.")
    lines.append("")
    lines.append("| Family | Sequences | Min | Median | Max | Over the limit |")
    lines.append("| --- | --- | --- | --- | --- | --- |")

    for family, n, low, median, high in family_sequence_lengths(conn):
        n_over = scalar(
            conn,
            """
            SELECT COUNT(*)
            FROM fam_target ft
            JOIN target_components tc ON tc.tid = ft.tid
            JOIN component_sequences cs ON cs.component_id = tc.component_id
            WHERE ft.family = ? AND LENGTH(cs.sequence) > ?
            """,
            (family, limit),
        )
        label = config.family_label(family, "en")
        pct = 100 * n_over / n if n else 0.0
        lines.append(
            f"| {label} | {n:,} | {low:,} | {median:,} | {high:,} "
            f"| {n_over:,} ({pct:.0f}%) |"
        )
        print(
            f"  {family:<18} median={median:>5,} aa  max={high:>6,} aa  "
            f"over limit={n_over:>4,} ({pct:.0f}%)"
        )


def report_family_domain_coverage(
    conn: sqlite3.Connection, lines: list[str], config
) -> None:
    """Binding-domain mapping coverage per family -- the stage-0 checkpoint.

    The binding domain is the primary protein representation precisely to avoid
    uneven truncation across families. If its coverage is itself uneven across
    families, it trades one confounder for another, and that has to be visible
    before embeddings are extracted rather than after.
    """
    section(lines, "Binding-domain coverage per family",
            "Aile başına bağlanma alanı kapsamı")
    lines.append("| Family | Activities | With a mapped domain | Coverage |")
    lines.append("| --- | --- | --- | --- |")

    coverage = {}
    for family, n_activities, n_mapped in family_binding_domain_coverage(
        conn, config
    ):
        label = config.family_label(family, "en")
        pct = 100 * n_mapped / n_activities if n_activities else 0.0
        coverage[family] = pct
        lines.append(
            f"| {label} | {n_activities:,} | {n_mapped:,} | **{pct:.1f}%** |"
        )
        print(f"  {family:<18} domain coverage={pct:>5.1f}%  "
              f"({n_mapped:,} / {n_activities:,})")

    if coverage:
        spread = max(coverage.values()) - min(coverage.values())
        lines.append("")
        lines.append(f"- Spread between families: **{spread:.1f} percentage points**")
        lines.append(
            "  A large spread means the binding-domain representation keeps a "
            "different share of each family's data, which is a confounder of "
            "the same kind truncation would have been. Decide the primary "
            "representation on this number."
        )


def report_binding_domains(conn: sqlite3.Connection, lines: list[str]) -> None:
    """Check whether the binding-domain representation is actually buildable.

    Needs three things to hold: the PPDMs mapping covers a usable share of
    activities, the mapping reaches a protein component, and domain start/end
    positions exist so a sub-sequence can be cut.
    """
    section(lines, "Binding-domain coverage", "Bağlanma alanı kapsamı")

    if not table_exists(conn, "predicted_binding_domains"):
        lines.append("`predicted_binding_domains` is absent in this release.")
        lines.append("The binding-domain representation is not available; the")
        lines.append("full-sequence representation becomes primary.")
        print("  binding domains: TABLE ABSENT -- representation plan needs revising")
        return

    pbd_cols = columns(conn, "predicted_binding_domains")
    cd_cols = columns(conn, "component_domains")
    lines.append(f"`predicted_binding_domains` columns: {', '.join(pbd_cols)}")
    lines.append(f"`component_domains` columns: {', '.join(cd_cols)}")
    lines.append("")

    n_mapped = count(conn, "predicted_binding_domains")
    n_distinct_act = scalar(
        conn, "SELECT COUNT(DISTINCT activity_id) FROM predicted_binding_domains"
    )
    n_act_total = scalar(
        conn,
        "SELECT COUNT(*) FROM activities WHERE pchembl_value IS NOT NULL "
        "AND standard_type IN ('Ki','Kd','IC50') AND standard_relation = '='",
    )

    lines.append(f"- Mapping rows: {n_mapped:,}")
    lines.append(f"- Distinct activities mapped: {n_distinct_act:,}")
    lines.append(f"- Candidate activities (Ki/Kd/IC50, exact, pChEMBL): {n_act_total:,}")
    if n_act_total:
        lines.append(f"- **Coverage: {100 * n_distinct_act / n_act_total:.1f}%**")

    # Confidence distribution, if the release records one. PPDMs reports ~90%
    # correct domain assignment in validation; the per-row confidence is how we
    # would filter.
    if "confidence" in pbd_cols:
        lines.append("")
        lines.append("| confidence | rows |")
        lines.append("| --- | --- |")
        for conf, n in conn.execute(
            "SELECT confidence, COUNT(*) FROM predicted_binding_domains "
            "GROUP BY confidence ORDER BY confidence"
        ):
            lines.append(f"| {conf} | {n:,} |")

    # Can we actually cut a domain out of a sequence?
    has_positions = {"start_position", "end_position"} <= set(cd_cols)
    lines.append("")
    if has_positions:
        n_with_pos = scalar(
            conn,
            "SELECT COUNT(*) FROM component_domains "
            "WHERE start_position IS NOT NULL AND end_position IS NOT NULL",
        )
        median_len = scalar(
            conn,
            """
            SELECT (end_position - start_position + 1) AS len
            FROM component_domains
            WHERE start_position IS NOT NULL AND end_position IS NOT NULL
            ORDER BY len
            LIMIT 1 OFFSET (
                SELECT COUNT(*) / 2 FROM component_domains
                WHERE start_position IS NOT NULL AND end_position IS NOT NULL
            )
            """,
        )
        over_limit = scalar(
            conn,
            "SELECT COUNT(*) FROM component_domains "
            "WHERE end_position - start_position + 1 > 1022",
        )
        lines.append(f"- Domains with usable positions: {n_with_pos:,}")
        lines.append(f"- Median domain length: {median_len} residues")
        lines.append(f"- Domains above the ESM-2 1022 limit: {over_limit:,}")
        print(f"  binding domains: coverage "
              f"{100 * n_distinct_act / n_act_total:.1f}%, "
              f"median length {median_len} aa, {over_limit:,} above the ESM-2 limit")
    else:
        lines.append("- **`component_domains` carries no start/end positions.**")
        lines.append("  A domain sub-sequence cannot be cut from the full sequence;")
        lines.append("  the binding-domain representation would have to come from")
        lines.append("  another source (e.g. InterPro) or be dropped.")
        print("  binding domains: NO POSITION COLUMNS -- cannot cut sub-sequences")


def report_activity_profile(conn: sqlite3.Connection, lines: list[str]) -> None:
    """Assay-type and pChEMBL distribution over the candidate records."""
    section(lines, "Activity profile", "Aktivite profili")

    lines.append("| standard_type | Records (exact, nM, with pChEMBL) |")
    lines.append("| --- | --- |")
    for stype, n in conn.execute(
        """
        SELECT standard_type, COUNT(*)
        FROM activities
        WHERE standard_type IN ('Ki','Kd','IC50')
          AND standard_units = 'nM'
          AND standard_relation = '='
          AND pchembl_value IS NOT NULL
        GROUP BY standard_type ORDER BY COUNT(*) DESC
        """
    ):
        lines.append(f"| {stype} | {n:,} |")

    # Does ChEMBL's own pchembl_value agree with our conversion? If it does we
    # use theirs (config: normalisation.use_chembl_pchembl) instead of silently
    # recomputing.
    n_disagree = scalar(
        conn,
        """
        SELECT COUNT(*) FROM activities
        WHERE standard_type IN ('Ki','Kd','IC50')
          AND standard_units = 'nM'
          AND standard_relation = '='
          AND pchembl_value IS NOT NULL
          AND standard_value > 0
          AND ABS(pchembl_value - (9.0 - LOG10(standard_value))) > 0.01
        """,
    )
    lines.append("")
    lines.append(f"- Records where `pchembl_value` disagrees with "
                 f"9 - log10(nM) by more than 0.01: **{n_disagree:,}**")
    lines.append("  (A nonzero count is expected for records whose pChEMBL was")
    lines.append("  derived from a different standard_value; they are dropped or")
    lines.append("  recomputed in stage 3.)")
    print(f"  pchembl_value disagreements: {n_disagree:,}")


def main() -> int:
    config = load_config()
    db_path = config.path("chembl_sqlite")

    if not db_path.exists():
        raise SystemExit(
            f"{db_path} not found.\n"
            "Run: python scripts/00_download_chembl.py --version 37"
        )

    print(f"Inspecting {db_path.relative_to(REPO_ROOT)} "
          f"({db_path.stat().st_size / 1e9:.1f} GB)")

    # Read-only, immutable: this script must never be able to modify the dump.
    conn = connect(db_path)

    lines: list[str] = [
        "# ChEMBL inventory / ChEMBL envanteri",
        "",
        "Generated by `scripts/01_inspect_chembl.py`. Read-only diagnostic:",
        "verifies the schema assumptions of the pipeline and reports the real",
        "protein-classification inventory that `config/families.yaml` must be",
        "written against.",
        "",
        "Bu bir teşhis çıktısıdır, tez tablosu değil. Veri hattının şema",
        "varsayımlarını doğrular ve `config/families.yaml` eşlemesinin",
        "yazılacağı gerçek sınıf envanterini raporlar.",
        "",
        f"- Database: `{db_path.relative_to(REPO_ROOT)}`",
        f"- Size: {db_path.stat().st_size / 1e9:.1f} GB",
    ]

    try:
        missing = check_schema(conn, lines)
        report_class_hierarchy(conn, lines)
        report_family_profile(conn, lines, config)
        report_family_lengths(conn, lines, config)
        report_family_domain_coverage(conn, lines, config)
        report_binding_domains(conn, lines)
        report_activity_profile(conn, lines)
    finally:
        conn.close()

    out = config.path("tables") / "01_chembl_inventory.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nReport written to {out.relative_to(REPO_ROOT)}")

    if missing:
        print("\nSchema assumptions do not hold; fix the queries before extracting.")
        return 1
    print("Next: map config/families.yaml onto the class names in the report.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
