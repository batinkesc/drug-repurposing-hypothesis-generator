"""Standardise every compound, in parallel.

Stage 2 of the pipeline, kept as its own script because it is the expensive
step: 497,867 compounds at ~37.5 ms each is about 5.2 hours on one core. The
work is per-molecule independent, so it is split across processes; the result
is identical either way, which is what lets us parallelise a step that is part
of the reproducibility contract.

The tautomer canonicalisation inside it cannot be skipped even though it is
what costs the time: it changes the canonical SMILES of 19% of molecules, and
those 19% would otherwise survive as two records of one compound -- able to
land on both sides of the R4b split and silently leak the "chemically new
ligand" definition.

Input:  data/interim/compounds.parquet
Output: data/interim/compounds_standardised.parquet

Usage::

    python scripts/03_standardise_molecules.py
    python scripts/03_standardise_molecules.py --force     # redo existing work
    python scripts/03_standardise_molecules.py --limit 5000  # quick trial run
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dti.config import REPO_ROOT, load_config  # noqa: E402
from dti.molecules import Standardiser  # noqa: E402

# One Standardiser per worker process, built once in the initialiser: the RDKit
# standardisation objects are expensive to construct and hold internal caches,
# so building one per molecule would dominate the runtime.
_standardiser: Standardiser | None = None


def _init_worker(settings: dict) -> None:
    global _standardiser
    _standardiser = Standardiser(**settings)


def _standardise_chunk(smiles_batch: list[str]) -> list[tuple]:
    """Standardise one batch. Returns a row tuple per input SMILES."""
    assert _standardiser is not None, "worker initialiser did not run"
    rows = []
    for smiles in smiles_batch:
        result = _standardiser.standardise(smiles)
        rows.append(
            (
                smiles,
                result.canonical_smiles,
                result.inchikey,
                result.n_heavy_atoms,
                result.reason,
            )
        )
    return rows


def chunked(items: list, size: int):
    for start in range(0, len(items), size):
        yield items[start : start + size]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="redo existing output")
    parser.add_argument("--limit", type=int, help="only the first N compounds")
    parser.add_argument("--workers", type=int, help="override worker count")
    parser.add_argument("--chunk-size", type=int, default=500)
    args = parser.parse_args()

    config = load_config()
    interim = config.path("interim")
    in_path = interim / "compounds.parquet"
    out_path = interim / "compounds_standardised.parquet"

    if not in_path.exists():
        raise SystemExit(f"{in_path} missing. Run scripts/02_extract.py first.")
    if out_path.exists() and not args.force and not args.limit:
        print(f"{out_path.relative_to(REPO_ROOT)} already exists; nothing to do.")
        print("Pass --force to redo it.")
        return 0

    compounds = pd.read_parquet(in_path)
    if args.limit:
        compounds = compounds.head(args.limit)

    settings = config["molecule"]
    # Leave one core free so the machine stays usable during the ~20 minutes.
    workers = (
        args.workers
        or settings.get("n_workers")
        or max(1, (os.cpu_count() or 2) - 1)
    )
    # Deduplicate the input SMILES: the same raw string appears under several
    # molregnos, and standardising it twice is wasted work that cannot change
    # the answer.
    unique_smiles = compounds["smiles_raw"].dropna().unique().tolist()

    print(f"{len(compounds):,} compounds, {len(unique_smiles):,} distinct raw SMILES")
    print(f"Standardising over {workers} workers "
          f"(tautomer canonicalisation: {settings['canonical_tautomer']})")

    worker_settings = {
        "strip_salts": settings["strip_salts"],
        "neutralise_charge": settings["neutralise_charge"],
        "canonical_tautomer": settings["canonical_tautomer"],
        "min_heavy_atoms": settings["min_heavy_atoms"],
        "max_heavy_atoms": settings["max_heavy_atoms"],
    }

    batches = list(chunked(unique_smiles, args.chunk_size))
    rows: list[tuple] = []
    started = time.perf_counter()

    with ProcessPoolExecutor(
        max_workers=workers, initializer=_init_worker, initargs=(worker_settings,)
    ) as pool:
        for i, batch_rows in enumerate(
            pool.map(_standardise_chunk, batches), start=1
        ):
            rows.extend(batch_rows)
            if i % 20 == 0 or i == len(batches):
                done = len(rows)
                elapsed = time.perf_counter() - started
                rate = done / elapsed if elapsed else 0
                remaining = (len(unique_smiles) - done) / rate if rate else 0
                print(
                    f"\r  {done:,} / {len(unique_smiles):,} "
                    f"({100 * done / len(unique_smiles):.1f}%)  "
                    f"{rate:.0f} mol/s  ~{remaining / 60:.1f} min left",
                    end="",
                    flush=True,
                )
    print()

    standardised = pd.DataFrame(
        rows,
        columns=[
            "smiles_raw",
            "canonical_smiles",
            "inchikey",
            "n_heavy_atoms",
            "drop_reason",
        ],
    )

    # Join back onto the full compound table: one row per molregno, carrying
    # the standardisation outcome.
    merged = compounds.merge(standardised, on="smiles_raw", how="left")
    merged.to_parquet(out_path, index=False)

    elapsed = time.perf_counter() - started
    kept = merged["canonical_smiles"].notna().sum()
    n_distinct = merged["canonical_smiles"].nunique()

    print(f"\nDone in {elapsed / 60:.1f} min -> {out_path.relative_to(REPO_ROOT)}")
    print(f"  kept:                {kept:,} / {len(merged):,}")
    print(f"  distinct molecules:  {n_distinct:,}")
    # Collapsing is the point, not a side effect: salt, charge and tautomer
    # variants of one compound become one molecule, which is what keeps the
    # R4b split honest.
    print(f"  collapsed by standardisation: {kept - n_distinct:,} records")
    print("\n  drop reasons:")
    for reason, n in merged["drop_reason"].value_counts().items():
        print(f"    {reason:<28} {n:>8,}")
    print("\nNext: python scripts/04_build_dataset.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
