"""Extract the raw dataset from ChEMBL into parquet files.

Stage 1 of the pipeline: pull targets, activities, compounds and binding-domain
boundaries, and record how many records each filter removes. No
standardisation happens here -- that is stage 2 (scripts/03_clean.py), so that
the expensive molecule work runs on an already-filtered set and each stage can
be re-run on its own.

Outputs under data/interim/:
    activities.parquet   one row per measurement that passed the filters
    targets.parquet      one row per target: accession, sequence
    compounds.parquet    one row per compound: raw SMILES, clinical phase
    domains.parquet      activity -> binding domain, with boundaries
    ../../reports/tables/03_drop_table.md   what each filter removed

The drop table is a thesis table, not a log: the transparency of the pipeline
is half the credibility of the results.

Usage::

    python scripts/02_extract.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dti.chembl import build_family_target_table, connect  # noqa: E402
from dti.config import REPO_ROOT, load_config  # noqa: E402

# The filters, in the order the design document applies them, each with the
# SQL predicate it adds. Counting them cumulatively is what makes the drop
# table honest: each row shows what remains after that filter, so the reader
# can see which step costs what.
FILTER_STAGES = [
    (
        "family_matched",
        "Aile sınıfı eşleşen tüm biyoaktivite kayıtları",
        "All bioactivity records on a family-matched target",
        "1 = 1",
    ),
    (
        "assay_confidence",
        "Deney güven skoru >= eşik (tek proteine atanmış deney)",
        "Assay confidence score >= threshold (single-protein assignment)",
        "a.confidence_score >= :min_confidence",
    ),
    (
        "assay_type",
        "Ölçüm türü Ki / Kd / IC50",
        "Measurement type is Ki / Kd / IC50",
        "act.standard_type IN ('Ki', 'Kd', 'IC50')",
    ),
    (
        "units_nm",
        "Standart birim nM",
        "Standard unit is nM",
        "act.standard_units = 'nM'",
    ),
    (
        "exact_relation",
        "İlişki işareti kesin (=), sansürlü kayıtlar dışarıda",
        "Exact relation (=); censored records excluded",
        "act.standard_relation = '='",
    ),
    (
        "pchembl_present",
        "pChEMBL değeri mevcut",
        "A pChEMBL value is present",
        "act.pchembl_value IS NOT NULL",
    ),
    (
        "validity_ok",
        "Veri doğrulama bayrağı yok",
        "No data-validity comment",
        "act.data_validity_comment IS NULL",
    ),
    (
        "pchembl_range",
        "pChEMBL makul aralıkta",
        "pChEMBL within a plausible range",
        "act.pchembl_value BETWEEN :pchembl_min AND :pchembl_max",
    ),
    (
        "has_smiles",
        "Bileşiğin SMILES kaydı var",
        "The compound has a SMILES record",
        "cst.canonical_smiles IS NOT NULL",
    ),
    (
        "has_sequence",
        "Hedefin dizi kaydı var",
        "The target has a sequence record",
        "EXISTS (SELECT 1 FROM target_components tc2 "
        "JOIN component_sequences cs2 ON cs2.component_id = tc2.component_id "
        "WHERE tc2.tid = td.tid AND cs2.sequence IS NOT NULL)",
    ),
]

def base_from(extra_joins: str = "") -> str:
    """The FROM/WHERE shared by every query here, plus any extra joins.

    One definition so that the drop-table counts and the extracted rows cannot
    drift apart: if the join path differed between them, the table would
    describe a dataset that was never written.
    """
    return f"""
    FROM fam_target ft
    JOIN target_dictionary td  ON td.tid        = ft.tid
    JOIN assays            a   ON a.tid         = td.tid
    JOIN activities        act ON act.assay_id  = a.assay_id
    LEFT JOIN compound_structures cst ON cst.molregno = act.molregno
    {extra_joins}
    WHERE td.target_type = :target_type
      AND td.organism    = :organism
"""


def all_filters() -> str:
    return " AND ".join(f"({p})" for _, _, _, p in FILTER_STAGES)


def query_params(config) -> dict:
    return {
        "target_type": config["filtering"]["target_type"],
        "organism": config["filtering"]["organism"],
        "min_confidence": config["filtering"]["min_assay_confidence_score"],
        "pchembl_min": config["normalisation"]["pchembl_range"][0],
        "pchembl_max": config["normalisation"]["pchembl_range"][1],
    }


def cumulative_drop_table(conn, config) -> pd.DataFrame:
    """Count what remains after each filter, applied cumulatively."""
    params = query_params(config)
    predicates: list[str] = []
    rows = []
    previous = None

    for key, label_tr, label_en, predicate in FILTER_STAGES:
        predicates.append(predicate)
        where = " AND ".join(f"({p})" for p in predicates)
        n = conn.execute(
            f"SELECT COUNT(*) {base_from()} AND {where}", params
        ).fetchone()[0]
        removed = None if previous is None else previous - n
        rows.append(
            {
                "stage": key,
                "label_tr": label_tr,
                "label_en": label_en,
                "remaining": n,
                "removed": removed,
            }
        )
        print(
            f"  {key:<18} remaining={n:>9,}"
            + (f"  removed={removed:>9,}" if removed is not None else "")
        )
        previous = n

    return pd.DataFrame(rows)


def extract_activities(conn, config) -> pd.DataFrame:
    """All measurements passing every filter, one row each."""
    return pd.read_sql_query(
        f"""
        SELECT act.activity_id,
               act.assay_id,
               act.molregno,
               td.tid,
               td.chembl_id        AS target_chembl_id,
               ft.family,
               act.standard_type,
               act.standard_value,
               act.pchembl_value,
               a.confidence_score,
               a.description       AS assay_description
        {base_from()} AND {all_filters()}
        """,
        conn,
        params=query_params(config),
    )


def extract_targets(conn, config) -> pd.DataFrame:
    """Target sequences. The longest component sequence per target is the one
    that would be embedded, so it is the one recorded."""
    params = query_params(config)
    return pd.read_sql_query(
        """
        SELECT ft.tid,
               ft.family,
               td.chembl_id  AS target_chembl_id,
               td.pref_name,
               cs.accession,
               cs.component_id,
               cs.sequence,
               LENGTH(cs.sequence) AS sequence_length
        FROM fam_target ft
        JOIN target_dictionary   td ON td.tid          = ft.tid
        JOIN target_components   tc ON tc.tid          = ft.tid
        JOIN component_sequences cs ON cs.component_id = tc.component_id
        WHERE cs.sequence IS NOT NULL
          AND td.target_type = :target_type
          AND td.organism    = :organism
        """,
        conn,
        params=params,
    )


def extract_compounds(conn, config) -> pd.DataFrame:
    """Compounds appearing in the filtered activities, with clinical phase.

    `max_phase` is what the retrieval protocol (D6) selects approved drugs on,
    so it is carried from the start rather than joined back later.
    """
    joins = "JOIN molecule_dictionary md ON md.molregno = act.molregno"
    return pd.read_sql_query(
        f"""
        SELECT DISTINCT act.molregno,
               md.chembl_id AS compound_chembl_id,
               md.max_phase,
               cst.canonical_smiles AS smiles_raw
        {base_from(joins)} AND {all_filters()}
        """,
        conn,
        params=query_params(config),
    )


def extract_domains(conn, config) -> pd.DataFrame:
    """Activity -> binding domain mapping with sequence boundaries.

    One activity can map to more than one domain (a binding site may sit at the
    interface of two). Those are not resolved here: the ambiguity is carried
    forward as a count per activity so the cleaning stage can drop them under
    `multi_domain_policy` and report how many were lost.
    """
    joins = """
    JOIN predicted_binding_domains pbd ON pbd.activity_id = act.activity_id
    JOIN site_components   sc ON sc.site_id      = pbd.site_id
    JOIN component_domains cd ON cd.domain_id    = sc.domain_id
                             AND cd.component_id = sc.component_id
    JOIN domains           d  ON d.domain_id     = cd.domain_id
    """
    frame = pd.read_sql_query(
        f"""
        SELECT act.activity_id,
               ft.family,
               sc.component_id,
               cd.domain_id,
               d.domain_name,
               d.source_domain_id,
               cd.start_position,
               cd.end_position,
               pbd.confidence AS domain_confidence
        {base_from(joins)} AND {all_filters()}
          AND cd.start_position IS NOT NULL
          AND cd.end_position   IS NOT NULL
        """,
        conn,
        params=query_params(config),
    )
    frame = frame[
        frame["start_position"].notna() & frame["end_position"].notna()
    ].copy()
    frame["domain_length"] = (
        frame["end_position"] - frame["start_position"] + 1
    ).astype(int)
    # How many distinct domains this activity maps to: 1 is unambiguous.
    frame["n_domains_for_activity"] = frame.groupby("activity_id")[
        "domain_id"
    ].transform("nunique")
    return frame


def write_drop_table(frame: pd.DataFrame, config) -> Path:
    lines = [
        "# Drop table / Düşürme tablosu",
        "",
        "Generated by `scripts/02_extract.py`. Her süzme adımının kaç kayıt",
        "düşürdüğü kümülatif olarak sayılır — veri hattının şeffaflığı,",
        "sonuçların güvenilirliğinin yarısı.",
        "",
        "Counts are cumulative: each row is what survives that filter and",
        "every filter above it.",
        "",
        "| # | Adım / Stage | Kalan / Remaining | Düşen / Removed |",
        "| --- | --- | --- | --- |",
    ]
    for i, row in enumerate(frame.itertuples(), start=1):
        removed = "—" if pd.isna(row.removed) else f"{int(row.removed):,}"
        lines.append(
            f"| {i} | {row.label_tr}<br>*{row.label_en}* "
            f"| {row.remaining:,} | {removed} |"
        )

    first, last = frame["remaining"].iloc[0], frame["remaining"].iloc[-1]
    kept = 100 * last / first if first else 0.0
    lines += [
        "",
        f"- Başlangıç / Start: **{first:,}**",
        f"- Son / Final: **{last:,}** (%{kept:.1f})",
    ]
    out = config.path("tables") / "03_drop_table.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> int:
    config = load_config()
    conn = connect(config.path("chembl_sqlite"))
    interim = config.path("interim")
    interim.mkdir(parents=True, exist_ok=True)

    print("Building family target table")
    counts = build_family_target_table(conn, config)
    print(f"  {counts}")

    print("\nCounting filter stages")
    drops = cumulative_drop_table(conn, config)
    drop_path = write_drop_table(drops, config)
    drops.to_parquet(interim / "drop_table.parquet", index=False)

    print("\nExtracting")
    for name, extractor in [
        ("activities", extract_activities),
        ("targets", extract_targets),
        ("compounds", extract_compounds),
        ("domains", extract_domains),
    ]:
        frame = extractor(conn, config)
        path = interim / f"{name}.parquet"
        frame.to_parquet(path, index=False)
        print(f"  {name:<12} {len(frame):>9,} rows  -> {path.name}")

    conn.close()
    print(f"\nDrop table: {drop_path.relative_to(REPO_ROOT)}")
    print("Next: python scripts/03_clean.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
