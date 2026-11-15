"""Extract frozen ChemBERTa-2 and ESM-2 embeddings for the dataset.

Stage 6. Runs each encoder once over the distinct inputs and writes the result
to disk, so that training a head -- 13 regimes x 5 seeds -- never touches a
language model again. That is also what lets the product run on CPU: the 565
target embeddings are a 1.4 MB file.

Keyed on the dataset's own identities, not on row order: molecules by canonical
SMILES, targets by `tid`. Each array ships with an index parquet giving the key
for every row, and a JSON manifest recording the model id, the pinned revision,
the pooling and the input counts. Anything that reads an embedding must go
through the index -- an array alone is a matrix of numbers with no way to tell
whether it still matches the dataset.

Outputs under data/embeddings/:
    molecule_{mean,cls}.npy        (n_molecules, 384)   float16
    molecule_index.parquet         canonical_smiles per row
    protein_{mean,cls}.npy         (n_targets, 1280)    float16
    protein_index.parquet          tid, target_chembl_id, accession per row
    manifest.json                  what produced all of the above

Usage::

    python scripts/06_extract_embeddings.py
    python scripts/06_extract_embeddings.py --force
    python scripts/06_extract_embeddings.py --limit 2000   # trial run
    python scripts/06_extract_embeddings.py --cpu
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dti.config import REPO_ROOT, load_config  # noqa: E402
from dti.embeddings import POOLINGS, embed_texts, load_encoder  # noqa: E402


def write_side(
    out_dir: Path,
    name: str,
    arrays: dict[str, np.ndarray],
    index: pd.DataFrame,
) -> None:
    for pooling in POOLINGS:
        np.save(out_dir / f"{name}_{pooling}.npy", arrays[pooling])
    index.to_parquet(out_dir / f"{name}_index.parquet", index=False)


def verify_against_dataset(out_dir: Path, data: pd.DataFrame) -> dict:
    """Re-read what was written and check it still describes this dataset.

    Written as a check rather than an assumption because the failure it catches
    is silent: an embedding file left over from an earlier dataset has the right
    dtype and a plausible shape, and every downstream number would be computed
    from the wrong rows.
    """
    molecules = pd.read_parquet(out_dir / "molecule_index.parquet")
    proteins = pd.read_parquet(out_dir / "protein_index.parquet")
    checks = {}

    for pooling in POOLINGS:
        array = np.load(out_dir / f"molecule_{pooling}.npy")
        checks[f"molecule_{pooling}_shape"] = list(array.shape)
        if len(array) != len(molecules):
            raise AssertionError(
                f"molecule_{pooling}.npy has {len(array)} rows but the index "
                f"has {len(molecules)}"
            )
        if not np.isfinite(array.astype(np.float32)).all():
            raise AssertionError(f"molecule_{pooling}.npy contains non-finite values")
        array = np.load(out_dir / f"protein_{pooling}.npy")
        checks[f"protein_{pooling}_shape"] = list(array.shape)
        if len(array) != len(proteins):
            raise AssertionError(
                f"protein_{pooling}.npy has {len(array)} rows but the index "
                f"has {len(proteins)}"
            )
        if not np.isfinite(array.astype(np.float32)).all():
            raise AssertionError(f"protein_{pooling}.npy contains non-finite values")

    missing_molecules = set(data["canonical_smiles"]) - set(molecules["canonical_smiles"])
    missing_targets = set(data["tid"]) - set(proteins["tid"])
    if missing_molecules or missing_targets:
        raise AssertionError(
            f"{len(missing_molecules)} molecules and {len(missing_targets)} "
            f"targets in the dataset have no embedding"
        )
    checks["covers_dataset"] = True
    checks["n_molecules"] = len(molecules)
    checks["n_targets"] = len(proteins)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="redo existing output")
    parser.add_argument("--limit", type=int, help="only the first N molecules")
    parser.add_argument("--cpu", action="store_true", help="disable the GPU")
    parser.add_argument("--chem-batch", type=int, help="override molecule batch size")
    parser.add_argument("--prot-batch", type=int, help="override protein batch size")
    args = parser.parse_args()

    config = load_config()
    out_dir = config.path("embeddings")
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "manifest.json"

    if manifest_path.exists() and not args.force and not args.limit:
        print(f"{manifest_path.relative_to(REPO_ROOT)} already exists; nothing to do.")
        print("Pass --force to redo it.")
        return 0

    dataset_path = config.path("processed") / "dataset.parquet"
    if not dataset_path.exists():
        raise SystemExit(f"{dataset_path} missing. Run scripts/04_build_dataset.py.")

    settings = config["models"]["m1_frozen_lm"]
    device = torch.device(
        "cuda" if (torch.cuda.is_available() and not args.cpu) else "cpu"
    )

    print("Loading dataset")
    data = pd.read_parquet(dataset_path)
    molecules = sorted(data["canonical_smiles"].unique())
    if args.limit:
        molecules = molecules[: args.limit]
    targets = (
        data.drop_duplicates("tid")
        .sort_values("tid")[
            ["tid", "target_chembl_id", "accession", "family", "domain_sequence"]
        ]
        .reset_index(drop=True)
    )
    print(
        f"  {len(data):,} pairs -> {len(molecules):,} distinct molecules, "
        f"{len(targets):,} targets  device={device}"
    )
    if device.type == "cpu":
        print("  WARNING: running on CPU; ESM-2 650M will be slow")

    reports = {}

    # --- Protein tower -----------------------------------------------------
    # First because it is the one that can run out of VRAM: 33 layers at 1280
    # hidden over up to 586 residues. If it fits, the molecule tower will.
    print(f"\nProtein tower: {settings['prot_model']}")
    protein_encoder = load_encoder(
        settings["prot_model"],
        settings["prot_revision"],
        # +2 for <cls> and <eos>; the residue budget itself is the config's.
        max_length=config["protein"]["esm_max_residues"] + 2,
        dtype=settings["dtype"],
        device=device,
    )
    print(f"  revision {protein_encoder.revision[:12]}  hidden {protein_encoder.hidden_size}")
    protein_arrays, reports["protein"] = embed_texts(
        targets["domain_sequence"].tolist(),
        protein_encoder,
        batch_size=args.prot_batch or settings["batch_size_embedding"],
        device=device,
        progress_every=10,
        label="proteins ",
    )
    write_side(
        out_dir,
        "protein",
        protein_arrays,
        targets[["tid", "target_chembl_id", "accession", "family"]],
    )
    del protein_encoder, protein_arrays
    if device.type == "cuda":
        torch.cuda.empty_cache()

    # --- Molecule tower ----------------------------------------------------
    print(f"\nMolecule tower: {settings['chem_model']}")
    chem_encoder = load_encoder(
        settings["chem_model"],
        settings["chem_revision"],
        max_length=512,
        dtype=settings["dtype"],
        device=device,
    )
    print(f"  revision {chem_encoder.revision[:12]}  hidden {chem_encoder.hidden_size}")
    molecule_arrays, reports["molecule"] = embed_texts(
        molecules,
        chem_encoder,
        batch_size=args.chem_batch or settings["batch_size_embedding_chem"],
        device=device,
        progress_every=100,
        label="molecules ",
    )
    write_side(
        out_dir,
        "molecule",
        molecule_arrays,
        pd.DataFrame({"canonical_smiles": molecules}),
    )
    del chem_encoder, molecule_arrays
    if device.type == "cuda":
        torch.cuda.empty_cache()

    # --- Manifest and verification ----------------------------------------
    print("\nVerifying what was written")
    checks = verify_against_dataset(out_dir, data if not args.limit else data.head(0))
    manifest = {
        "created": date.today().isoformat(),
        "dataset_rows": len(data),
        "pooling_written": list(POOLINGS),
        "pooling_used": settings["pooling"],
        "dtype": "float16",
        "device": device.type,
        "limited_run": bool(args.limit),
        "towers": reports,
        "checks": checks,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    for tower, report in reports.items():
        print(
            f"  {tower:<9} {report['n_inputs']:>8,} inputs  "
            f"{report['seconds'] / 60:>5.1f} min  "
            f"max {report['max_tokens_seen']:>4} tokens  "
            f"truncated {report['n_truncated']}"
        )
    truncated = sum(r["n_truncated"] for r in reports.values())
    if truncated:
        print(
            f"  WARNING: {truncated} inputs hit the token limit and were "
            f"truncated; see the manifest"
        )
    print(f"\nEmbeddings: {out_dir.relative_to(REPO_ROOT)}")
    print(f"Manifest:   {manifest_path.relative_to(REPO_ROOT)}")
    print("\nNext: python scripts/07_baselines.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
