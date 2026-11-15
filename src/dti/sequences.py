"""Pairwise sequence identity between targets, computed exactly.

The design document specified MMseqs2 for the R3 split. MMseqs2 is the right
tool for millions of sequences; here there are **565**, and that changes the
answer. 565 targets is 159,330 pairs, which exact Smith-Waterman alignment
settles in about six minutes on one core and under a minute across the
machine. At that size MMseqs2's value -- a k-mer prefilter that avoids most
alignments -- buys nothing, while its cost is real: the prefilter is a
heuristic with a sensitivity setting, and pairs just above the 0.40 identity
cutoff are exactly the ones a prefilter can miss. A missed pair in R3 is a
test target left similar to a training target, i.e. the leakage the split
exists to prevent.

So this computes every pair exactly, which also removes a Windows-only binary
install from the pipeline. The cost of being exact is 159,330 alignments; the
benefit is that the R3 guarantee holds by construction rather than up to a
heuristic's recall.

Identity is computed on the **domain** sequence, not the full protein, because
the domain is what the model is given (see `protein.representation` in the
config). Two proteins can differ across their full length while sharing a
near-identical catalytic domain; clustering on the full sequence would call
them distinct targets and then hand the model the same input on both sides of
the split.

Three numbers are returned per pair, because the obvious one is a trap:

* `identity`  = identical aligned columns / alignment length (gaps included)
* `coverage`  = aligned residue pairs / length of the SHORTER sequence
* `identity_over_shorter` = identical residues / length of the SHORTER sequence

`identity` on its own is close to useless here and the first run of the audit
proved it: local alignment finds the best-scoring segment, so two unrelated
proteins sharing a ten-residue motif come back at 70-100% "identity" over a
window covering 3% of the sequence. Every regime's audit therefore reported a
maximum identity of 1.000, including R3, where no target is on both sides at
all. The artefact is in the metric, not in the split.

`identity_over_shorter` has no such loophole: it asks what fraction of the
smaller protein is identical to the larger one, so a short motif can never
exceed the length ratio. It is the number the leakage audit reports.

For clustering, `identity` and `coverage` are used as a pair and both must
clear their threshold -- equivalent in spirit, and kept because it is what the
R3 groups were built with.
"""

from __future__ import annotations

import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

# Set once per worker process: the aligner object is expensive to build and the
# sequence list would otherwise be pickled with every chunk.
_aligner = None
_sequences: list[str] = []


def _init_worker(sequences: list[str]) -> None:
    global _aligner, _sequences
    from Bio import Align

    # blastp defaults: BLOSUM62, gap open -12, extend -1. Named rather than
    # hand-set so the scoring is a documented standard, not an arbitrary choice.
    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "local"
    _aligner = aligner
    _sequences = sequences


def _align_pairs(pairs: np.ndarray) -> list[tuple[int, int, float, float, float]]:
    assert _aligner is not None, "worker initialiser did not run"
    out = []
    for i, j in pairs:
        a, b = _sequences[i], _sequences[j]
        alignment = _aligner.align(a, b)[0]
        counts = alignment.counts()
        aligned_columns = counts.identities + counts.mismatches + counts.gaps
        aligned_pairs = counts.identities + counts.mismatches
        shorter = min(len(a), len(b))
        identity = counts.identities / aligned_columns if aligned_columns else 0.0
        coverage = aligned_pairs / shorter
        identity_over_shorter = counts.identities / shorter
        out.append((int(i), int(j), identity, coverage, identity_over_shorter))
    return out


def identity_matrix(
    sequences: list[str],
    *,
    workers: int | None = None,
    chunk_size: int = 400,
    progress: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Symmetric (n, n) identity, coverage and identity-over-shorter matrices.

    Diagonals are 1.0: a sequence is wholly identical to itself.
    """
    n = len(sequences)
    identity = np.eye(n, dtype=np.float32)
    coverage = np.eye(n, dtype=np.float32)
    over_shorter = np.eye(n, dtype=np.float32)
    if n < 2:
        return identity, coverage, over_shorter

    pairs = np.array(
        [(i, j) for i in range(n) for j in range(i + 1, n)], dtype=np.int64
    )
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    chunks = [pairs[k : k + chunk_size] for k in range(0, len(pairs), chunk_size)]

    started = time.perf_counter()
    done = 0
    with ProcessPoolExecutor(
        max_workers=workers, initializer=_init_worker, initargs=(sequences,)
    ) as pool:
        for results in pool.map(_align_pairs, chunks):
            for i, j, ident, cov, over in results:
                identity[i, j] = identity[j, i] = ident
                coverage[i, j] = coverage[j, i] = cov
                over_shorter[i, j] = over_shorter[j, i] = over
            done += len(results)
            elapsed = time.perf_counter() - started
            if progress and elapsed > 0:
                rate = done / elapsed
                print(
                    f"\r  aligning {done:,} / {len(pairs):,} pairs  "
                    f"{rate:.0f}/s  ~{(len(pairs) - done) / rate / 60:.1f} min left",
                    end="",
                    flush=True,
                )
    if progress:
        print()
    return identity, coverage, over_shorter


def load_or_compute(
    cache: Path,
    target_ids: list[str],
    sequences: list[str],
    *,
    workers: int | None = None,
    force: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cached identity/coverage matrices.

    The cache stores the target ids it was computed for and is rejected if they
    differ, so a dataset rebuild cannot silently reuse a matrix whose rows no
    longer line up with the targets. That failure would be invisible: the
    matrix would still be the right shape if the target count happened to match.
    """
    ids = np.asarray(target_ids, dtype=object)
    if cache.exists() and not force:
        stored = np.load(cache, allow_pickle=True)
        matches = (
            len(stored["target_ids"]) == len(ids)
            and (stored["target_ids"] == ids).all()
            # A cache written before `identity_over_shorter` existed is not
            # usable; the audit depends on it.
            and "identity_over_shorter" in stored
        )
        if matches:
            return (
                stored["identity"],
                stored["coverage"],
                stored["identity_over_shorter"],
            )
        print(f"  cache {cache.name} is stale or incomplete; recomputing")

    identity, coverage, over_shorter = identity_matrix(sequences, workers=workers)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache,
        identity=identity,
        coverage=coverage,
        identity_over_shorter=over_shorter,
        target_ids=ids,
    )
    return identity, coverage, over_shorter
