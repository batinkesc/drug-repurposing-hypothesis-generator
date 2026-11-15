"""Measure the two candidate protein representations, side by side.

The stage-0 checkpoint: the binding-domain representation was adopted because
truncating full sequences at the ESM-2 limit hits the families at different
rates, which is a confounder for the study's main question. But the domain
mapping covers only part of the data, and if that coverage is itself uneven
across families it trades one confounder for another of the same kind.

So the decision cannot be made on coverage alone. This script measures both
paths in the same units -- what share of each family's activities survives --
and reports the spread across families for each, because the spread is the
confounder, not the average.

Two measurements the first inventory did not make:

1. Truncation weighted by activities, not targets. "21% of kinase targets are
   over the limit" says little; those targets may hold 2% or 60% of the
   family's measurements.
2. The length distribution of the domains that are actually mapped, rather
   than of every Pfam domain in the database.

Usage::

    python scripts/01b_representation_decision.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dti.chembl import (  # noqa: E402
    ACTIVITY_FILTER_SQL,
    ASSAY_FILTER_SQL,
    TARGET_FILTER_SQL,
    build_family_target_table,
    connect,
)
from dti.config import REPO_ROOT, load_config  # noqa: E402


def filter_params(config) -> dict:
    return {
        "target_type": config["filtering"]["target_type"],
        "organism": config["filtering"]["organism"],
        "min_confidence": config["filtering"]["min_assay_confidence_score"],
    }


def truncation_by_activity(conn, config) -> dict[str, tuple[int, int]]:
    """Per family: (activities, activities on targets over the ESM-2 limit).

    A target's length is taken as its longest component sequence: that is what
    would be fed to ESM-2 and therefore what would be truncated.
    """
    params = filter_params(config) | {"limit": config["protein"]["esm_max_residues"]}
    rows = conn.execute(
        f"""
        WITH target_length AS (
            SELECT ft.tid, ft.family, MAX(LENGTH(cs.sequence)) AS len
            FROM fam_target ft
            JOIN target_components   tc ON tc.tid          = ft.tid
            JOIN component_sequences cs ON cs.component_id = tc.component_id
            WHERE cs.sequence IS NOT NULL
            GROUP BY ft.tid, ft.family
        )
        SELECT tl.family,
               COUNT(*) AS n_activities,
               SUM(CASE WHEN tl.len > :limit THEN 1 ELSE 0 END) AS n_truncated
        FROM target_length tl
        JOIN target_dictionary td ON td.tid = tl.tid
        JOIN assays     a   ON a.tid        = td.tid
        JOIN activities act ON act.assay_id = a.assay_id
        WHERE {TARGET_FILTER_SQL}
          AND {ASSAY_FILTER_SQL}
          AND {ACTIVITY_FILTER_SQL}
        GROUP BY tl.family
        """,
        params,
    ).fetchall()
    return {family: (n_act, n_trunc) for family, n_act, n_trunc in rows}


def mapped_domain_lengths(conn, config) -> dict[str, list[int]]:
    """Per family: the lengths of the binding domains actually mapped.

    Joins the PPDMs activity -> site mapping through to the domain boundaries,
    which is what makes the representation buildable at all: without
    start/end positions there is no sub-sequence to embed.
    """
    params = filter_params(config)
    rows = conn.execute(
        f"""
        SELECT ft.family,
               cd.end_position - cd.start_position + 1 AS domain_length
        FROM fam_target ft
        JOIN target_dictionary td ON td.tid = ft.tid
        JOIN assays     a   ON a.tid        = td.tid
        JOIN activities act ON act.assay_id = a.assay_id
        JOIN predicted_binding_domains pbd ON pbd.activity_id = act.activity_id
        JOIN site_components   sc ON sc.site_id      = pbd.site_id
        JOIN component_domains cd ON cd.domain_id    = sc.domain_id
                                 AND cd.component_id = sc.component_id
        WHERE {TARGET_FILTER_SQL}
          AND {ASSAY_FILTER_SQL}
          AND {ACTIVITY_FILTER_SQL}
          AND cd.start_position IS NOT NULL
          AND cd.end_position   IS NOT NULL
        """,
        params,
    ).fetchall()

    lengths: dict[str, list[int]] = {}
    for family, length in rows:
        lengths.setdefault(family, []).append(length)
    return lengths


def percentile(values: list[int], q: float) -> int:
    if not values:
        return 0
    ordered = sorted(values)
    index = min(int(q * len(ordered)), len(ordered) - 1)
    return ordered[index]


def main() -> int:
    config = load_config()
    conn = connect(config.path("chembl_sqlite"))
    limit = config["protein"]["esm_max_residues"]

    print("Building family target table")
    build_family_target_table(conn, config)

    print("Measuring truncation weighted by activity")
    truncation = truncation_by_activity(conn, config)

    print("Measuring mapped binding-domain lengths")
    domains = mapped_domain_lengths(conn, config)

    lines = [
        "# Protein representation decision / Protein temsili kararı",
        "",
        "Generated by `scripts/01b_representation_decision.py`.",
        "",
        "İki aday temsil aynı birimde karşılaştırılıyor: her ailenin "
        "aktivitelerinin ne kadarı korunuyor. Karar ortalamaya değil "
        "**aileler arası yayılıma** bakılarak verilir, çünkü karıştırıcı olan "
        "yayılımın kendisi.",
        "",
        "The two candidate representations in the same units: what share of "
        "each family's activities survives. The decision rests on the spread "
        "across families, not the average, because the spread is the "
        "confounder.",
        "",
        f"ESM-2 limit: {limit} residues.",
        "",
        "## Path A — full sequence, truncated / Tam dizi, kırpılmış",
        "",
        "| Family | Activities | On truncated targets | Share truncated |",
        "| --- | --- | --- | --- |",
    ]

    trunc_share = {}
    for family in sorted(truncation, key=lambda f: -truncation[f][0]):
        n_act, n_trunc = truncation[family]
        share = 100 * n_trunc / n_act if n_act else 0.0
        trunc_share[family] = share
        lines.append(
            f"| {config.family_label(family, 'en')} | {n_act:,} "
            f"| {n_trunc:,} | **{share:.1f}%** |"
        )
        print(f"  {family:<18} truncated activities={share:>5.1f}%")

    lines.append("")
    lines.append(
        "No activity is lost on this path -- a truncated sequence is still "
        "embedded, just incompletely. The share above is the fraction of "
        "measurements whose protein is only partly seen."
    )
    lines.append("")
    lines.append(
        "Bu yolda hiçbir aktivite kaybolmaz; kırpılmış dizi yine gömülür, "
        "yalnız eksik görülür. Yukarıdaki oran, proteini kısmen görülen "
        "ölçümlerin payıdır."
    )

    lines += [
        "",
        "## Path B — binding domain / Bağlanma alanı",
        "",
        "| Family | Mapped activities | Median | p95 | Max | Over the limit |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for family in sorted(domains, key=lambda f: -len(domains[f])):
        values = domains[family]
        n_over = sum(1 for v in values if v > limit)
        lines.append(
            f"| {config.family_label(family, 'en')} | {len(values):,} "
            f"| {percentile(values, 0.5):,} | {percentile(values, 0.95):,} "
            f"| {max(values):,} | {n_over:,} |"
        )
        print(
            f"  {family:<18} mapped={len(values):>8,}  "
            f"median={percentile(values, 0.5):>5,} aa  "
            f"max={max(values):>5,} aa  over limit={n_over:,}"
        )

    out = config.path("tables") / "02_representation_decision.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    conn.close()

    spread = max(trunc_share.values()) - min(trunc_share.values())
    print(f"\n  truncation spread across families: {spread:.1f} points")
    print(f"Report written to {out.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
