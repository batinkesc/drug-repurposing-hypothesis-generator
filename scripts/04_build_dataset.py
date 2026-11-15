"""Assemble the modelling dataset.

Stage 3: join the extracted tables, aggregate replicate measurements, attach a
protein representation, and write the dataset the splits and models will use.

Three decisions are implemented here, each of which changes the data:

1. **Compounds are identified by canonical SMILES, not by molregno.** Two
   ChEMBL compound records that standardise to the same molecule are one
   molecule. Keying on molregno would leave them as two, and the R4b split --
   which asks whether a test ligand is chemically new -- could then put the
   same molecule on both sides under two identifiers.

2. **Replicate measurements of one target-compound pair collapse to their
   median**, and the spread is recorded. A pair whose replicates disagree by
   more than the configured threshold is flagged, not dropped: an ablation
   drops them, so the flag has to survive into the dataset.

3. **The protein representation is the binding domain**, with no fallback to a
   truncated full sequence -- see config/config.yaml for the measurements
   behind that. An activity whose binding site maps ambiguously to the
   sequence is dropped under `multi_domain_policy` -- both when it maps to
   several domain families and when it maps to several copies of one family,
   which is the case ChEMBL cannot resolve at all.

Input:  data/interim/{activities,targets,domains,compounds_standardised}.parquet
Output: data/processed/dataset.parquet
        reports/tables/04_dataset_profile.md
        reports/tables/05_drop_table_stage2.md

Usage::

    python scripts/04_build_dataset.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dti.config import REPO_ROOT, load_config  # noqa: E402


class DropLog:
    """Accumulates how many rows each step removed, for the thesis table."""

    def __init__(self) -> None:
        self.rows: list[dict] = []

    def record(self, label_tr: str, label_en: str, before: int, after: int) -> None:
        self.rows.append(
            {
                "label_tr": label_tr,
                "label_en": label_en,
                "remaining": after,
                "removed": before - after,
            }
        )
        print(f"  {label_en:<46} {after:>9,}  (-{before - after:,})")

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


def attach_domains(
    activities: pd.DataFrame,
    domains: pd.DataFrame,
    targets: pd.DataFrame,
    config,
    log: DropLog,
) -> pd.DataFrame:
    """Attach one binding-domain sub-sequence per activity."""
    policy = config["protein"]["multi_domain_policy"]

    if policy != "drop_ambiguous":
        raise SystemExit(f"unsupported multi_domain_policy: {policy}")

    # Ambiguity has two forms, and counting only the first one is a mistake
    # that is easy to make: `n_domains_for_activity` counts DISTINCT domain
    # ids, so it sees a binding site that maps to two different Pfam families,
    # but not a site that maps to several COPIES of one family within the same
    # protein. The second form is the common one here -- a voltage-gated
    # channel carries four Ion_trans repeats, a JAK carries two kinase
    # domains -- and ChEMBL cannot resolve it: `component_domains` holds one
    # row per copy with its own boundaries, while `site_components` references
    # the domain id alone, with no position. The source genuinely does not say
    # which copy the site sits in.
    #
    # MEASURED (ChEMBL 37): copies affect 20.9% of ion-channel activities,
    # 7.0% of kinase, 6.3% of protease, 0.5% of GPCR and 0.0% of nuclear
    # receptor -- a 20.9-point spread. Picking one copy arbitrarily would feed
    # the model a quarter of an ion channel's pore while giving nuclear
    # receptors their whole domain: systematic, family-dependent information
    # loss, which is exactly the confounder the binding-domain representation
    # was chosen to avoid (config/config.yaml, `representation`). Dropping
    # these activities costs sample size instead, and sample loss costs power
    # rather than creating bias.
    #
    # The biologically motivated alternative -- concatenating the copies, so a
    # channel is represented by all four repeats -- is ablation 6, not the
    # primary path: the concatenation runs to 1,050-1,161 aa for 34% of the
    # affected activities, back over the ESM-2 limit, so it reintroduces
    # truncation for ion channels specifically.
    copies = domains.groupby("activity_id")[
        ["domain_id", "start_position", "end_position"]
    ].nunique()
    ambiguous_ids = set(copies.index[copies.max(axis=1) > 1])
    n_distinct_families = domains.loc[
        domains["n_domains_for_activity"] > 1, "activity_id"
    ].nunique()
    n_repeated_copies = len(ambiguous_ids) - n_distinct_families
    print(
        f"  ambiguous mappings dropped: {len(ambiguous_ids):,} activities "
        f"({n_distinct_families:,} map to several domain families, "
        f"{n_repeated_copies:,} to several copies of one)"
    )

    # One row per activity now that ambiguity is gone: what remains are
    # duplicate rows describing the SAME copy via different components, so any
    # of them gives the same sub-sequence. Sorted so re-runs agree.
    unambiguous = (
        domains[~domains["activity_id"].isin(ambiguous_ids)]
        .sort_values(["activity_id", "component_id"])
        .drop_duplicates("activity_id", keep="first")
    )

    merged = activities.merge(
        unambiguous[
            [
                "activity_id",
                "component_id",
                "domain_id",
                "domain_name",
                "start_position",
                "end_position",
                "domain_length",
            ]
        ],
        on="activity_id",
        how="inner",
    )
    log.record(
        "Bağlanma alanı eşlemesi olan kayıtlar",
        "Records with an unambiguous binding-domain mapping",
        len(activities),
        len(merged),
    )

    # Cut the domain out of its component's sequence. ChEMBL positions are
    # 1-based and inclusive.
    sequences = targets[["component_id", "sequence", "accession"]].drop_duplicates(
        "component_id"
    )
    merged = merged.merge(sequences, on="component_id", how="left")

    starts = merged["start_position"].astype("Int64")
    ends = merged["end_position"].astype("Int64")
    merged["domain_sequence"] = [
        seq[int(s) - 1 : int(e)] if isinstance(seq, str) and pd.notna(s) else None
        for seq, s, e in zip(merged["sequence"], starts, ends, strict=True)
    ]

    before = len(merged)
    min_residues = config["protein"]["min_residues"]
    max_residues = config["protein"]["esm_max_residues"]
    lengths = merged["domain_sequence"].str.len()
    merged = merged[
        merged["domain_sequence"].notna()
        & (lengths >= min_residues)
        & (lengths <= max_residues)
    ]
    log.record(
        f"Alan dizisi çıkarılabilir ve {min_residues}-{max_residues} aa arasında",
        f"Domain sub-sequence extractable and {min_residues}-{max_residues} aa",
        before,
        len(merged),
    )
    return merged


def aggregate_pairs(frame: pd.DataFrame, config, log: DropLog) -> pd.DataFrame:
    """Collapse replicates of one target-compound pair to their median.

    The pair key uses canonical SMILES, so two ChEMBL records of the same
    molecule are one compound here.
    """
    threshold = config["normalisation"]["inconsistency_log_units"]

    grouped = frame.groupby(["tid", "canonical_smiles"], sort=False)
    aggregated = grouped.agg(
        family=("family", "first"),
        target_chembl_id=("target_chembl_id", "first"),
        accession=("accession", "first"),
        domain_sequence=("domain_sequence", "first"),
        domain_name=("domain_name", "first"),
        domain_length=("domain_length", "first"),
        inchikey=("inchikey", "first"),
        compound_chembl_id=("compound_chembl_id", "first"),
        max_phase=("max_phase", "max"),
        pchembl=("pchembl_value", "median"),
        pchembl_min=("pchembl_value", "min"),
        pchembl_max=("pchembl_value", "max"),
        n_measurements=("pchembl_value", "size"),
        assay_types=("standard_type", lambda s: "|".join(sorted(set(s)))),
    ).reset_index()

    aggregated["pchembl_spread"] = (
        aggregated["pchembl_max"] - aggregated["pchembl_min"]
    )
    # Flagged, not dropped: ablation 4 drops these, so the flag must be in the
    # dataset rather than applied here.
    aggregated["inconsistent"] = aggregated["pchembl_spread"] > threshold
    aggregated["mixed_assay_types"] = aggregated["assay_types"].str.contains("|",
                                                                             regex=False)

    log.record(
        "Hedef–bileşik çifti bazında medyan toplama",
        "Collapsed to one row per target-compound pair (median)",
        len(frame),
        len(aggregated),
    )
    n_flagged = int(aggregated["inconsistent"].sum())
    print(
        f"  pairs flagged as inconsistent (>{threshold} log units, kept): "
        f"{n_flagged:,}"
    )
    return aggregated


def filter_sparse_targets(frame: pd.DataFrame, config, log: DropLog) -> pd.DataFrame:
    """Drop targets with too few measurements to support a per-target metric."""
    minimum = config["filtering"]["min_activities_per_target"]
    counts = frame["tid"].value_counts()
    keep = counts[counts >= minimum].index
    before = len(frame)
    filtered = frame[frame["tid"].isin(keep)]
    log.record(
        f"Hedef başına en az {minimum} ölçüm",
        f"At least {minimum} measurements per target",
        before,
        len(filtered),
    )
    return filtered


def write_profile(frame: pd.DataFrame, config) -> Path:
    """The stage-1 profile table: what the dataset actually contains."""
    lines = [
        "# Dataset profile / Veri seti profili",
        "",
        "Generated by `scripts/04_build_dataset.py`.",
        "",
        "Bir satır = bir hedef–bileşik çifti (tekrarlı ölçümler medyanla "
        "birleştirilmiş). Bileşik kimliği kanonik SMILES'tır, ChEMBL "
        "kayıt numarası değil.",
        "",
        "One row = one target-compound pair, replicates collapsed to their "
        "median. Compound identity is the canonical SMILES, not the ChEMBL "
        "record id.",
        "",
        "## Per family / Aile başına",
        "",
        "| Family | Pairs | Targets | Compounds | pChEMBL median | "
        "Pairs/target (median) |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for family, group in frame.groupby("family"):
        per_target = group["tid"].value_counts()
        lines.append(
            f"| {config.family_label(family, 'en')} | {len(group):,} "
            f"| {group['tid'].nunique():,} "
            f"| {group['canonical_smiles'].nunique():,} "
            f"| {group['pchembl'].median():.2f} "
            f"| {per_target.median():.0f} |"
        )
    lines.append(
        f"| **Total** | {len(frame):,} | {frame['tid'].nunique():,} "
        f"| {frame['canonical_smiles'].nunique():,} "
        f"| {frame['pchembl'].median():.2f} | |"
    )

    lines += [
        "",
        "## pChEMBL distribution / pChEMBL dağılımı",
        "",
        "Aileler arası ortalama farkı, R4'te üç RMSE varyantının (ham, "
        "ortalama-kaydırmalı, lineer kalibrasyonlu) neden gerekli olduğunun "
        "doğrudan gerekçesi: dışarıda bırakılan aileye geçişte model "
        "sistematik bir kayma yapar.",
        "",
        "The spread of family means is the direct justification for the three "
        "RMSE variants in R4: transferring to a held-out family produces a "
        "systematic offset.",
        "",
        "| Family | Mean | Std | Min | p25 | p50 | p75 | Max |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for family, group in frame.groupby("family"):
        values = group["pchembl"]
        lines.append(
            f"| {config.family_label(family, 'en')} | {values.mean():.2f} "
            f"| {values.std():.2f} | {values.min():.2f} "
            f"| {values.quantile(0.25):.2f} | {values.median():.2f} "
            f"| {values.quantile(0.75):.2f} | {values.max():.2f} |"
        )
    means = frame.groupby("family")["pchembl"].mean()
    lines.append("")
    lines.append(
        f"- Aileler arası ortalama farkı / spread of family means: "
        f"**{means.max() - means.min():.2f} log units**"
    )

    lines += [
        "",
        "## Measurement quality / Ölçüm kalitesi",
        "",
        "| Metric | Value |",
        "| --- | --- |",
        f"| Pairs with several measurements | "
        f"{int((frame['n_measurements'] > 1).sum()):,} |",
        f"| Pairs flagged inconsistent (ablation 4 drops these) | "
        f"{int(frame['inconsistent'].sum()):,} |",
        f"| Pairs mixing assay types (ablation 3) | "
        f"{int(frame['mixed_assay_types'].sum()):,} |",
        f"| Domain length, median | "
        f"{frame['domain_length'].median():.0f} aa |",
    ]

    out = config.path("tables") / "04_dataset_profile.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_drop_table(log: DropLog, config) -> Path:
    frame = log.to_frame()
    lines = [
        "# Drop table, stage 2 / Düşürme tablosu, 2. aşama",
        "",
        "Generated by `scripts/04_build_dataset.py`. Çekme sonrası adımların",
        "düşürdüğü kayıtlar. 1. aşama tablosu: `03_drop_table.md`.",
        "",
        "Rows removed after extraction. Stage 1 is in `03_drop_table.md`.",
        "",
        "| # | Adım / Stage | Kalan / Remaining | Düşen / Removed |",
        "| --- | --- | --- | --- |",
    ]
    for i, row in enumerate(frame.itertuples(index=False), start=1):
        lines.append(
            f"| {i} | {row.label_tr}<br>*{row.label_en}* "
            f"| {row.remaining:,} | {row.removed:,} |"
        )
    out = config.path("tables") / "05_drop_table_stage2.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> int:
    config = load_config()
    interim = config.path("interim")
    processed = config.path("processed")
    processed.mkdir(parents=True, exist_ok=True)

    required = {
        "activities": interim / "activities.parquet",
        "targets": interim / "targets.parquet",
        "domains": interim / "domains.parquet",
        "compounds": interim / "compounds_standardised.parquet",
    }
    for name, path in required.items():
        if not path.exists():
            raise SystemExit(
                f"{path.relative_to(REPO_ROOT)} missing ({name}). "
                "Run scripts/02_extract.py and 03_standardise_molecules.py first."
            )

    print("Loading")
    activities = pd.read_parquet(required["activities"])
    targets = pd.read_parquet(required["targets"])
    domains = pd.read_parquet(required["domains"])
    compounds = pd.read_parquet(required["compounds"])
    print(f"  activities {len(activities):>9,}  targets {len(targets):>6,}  "
          f"domains {len(domains):>8,}  compounds {len(compounds):>8,}")

    log = DropLog()
    print("\nBuilding")

    # Standardised compounds only.
    kept_compounds = compounds[compounds["canonical_smiles"].notna()]
    before = len(activities)
    frame = activities.merge(
        kept_compounds[
            [
                "molregno",
                "canonical_smiles",
                "inchikey",
                "compound_chembl_id",
                "max_phase",
                "n_heavy_atoms",
            ]
        ],
        on="molregno",
        how="inner",
    )
    log.record(
        "Standardizasyonu geçen bileşiklere ait kayıtlar",
        "Records whose compound passed standardisation",
        before,
        len(frame),
    )

    frame = attach_domains(frame, domains, targets, config, log)
    pairs = aggregate_pairs(frame, config, log)
    pairs = filter_sparse_targets(pairs, config, log)

    # Deterministic row order so the file is byte-stable across re-runs.
    pairs = pairs.sort_values(["family", "tid", "canonical_smiles"]).reset_index(
        drop=True
    )

    out_path = processed / "dataset.parquet"
    pairs.to_parquet(out_path, index=False)

    profile_path = write_profile(pairs, config)
    drop_path = write_drop_table(log, config)

    print(f"\nDataset: {out_path.relative_to(REPO_ROOT)}  ({len(pairs):,} pairs)")
    print(f"Profile: {profile_path.relative_to(REPO_ROOT)}")
    print(f"Drops:   {drop_path.relative_to(REPO_ROOT)}")

    means = pairs.groupby("family")["pchembl"].mean()
    print("\n  pChEMBL mean per family:")
    for family, mean in means.sort_values(ascending=False).items():
        n = int((pairs["family"] == family).sum())
        print(f"    {family:<18} {mean:.2f}   ({n:,} pairs)")
    print(f"  spread across families: {means.max() - means.min():.2f} log units")
    return 0


if __name__ == "__main__":
    sys.exit(main())
