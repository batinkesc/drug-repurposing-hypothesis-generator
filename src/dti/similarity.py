"""Large-scale Tanimoto similarity, tiled onto the GPU.

The split code has to answer two questions over 190,339 compounds:

* how many neighbours above a cutoff does each compound have (Butina needs
  this to order candidates), and
* what is the highest similarity between each test compound and the whole
  training set (the R4b filter and the leakage audit need this).

Both are all-pairs questions: 190,339^2 / 2 is 1.8e10 pairs. A dense matrix is
out of the question (3.6e10 float32 = 145 GB), and a Python loop over RDKit's
bulk similarity would take days. Expressed as a matrix product, though, the
intersection of two bit sets is just a dot product, so the whole pass is one
big fp16 matmul -- 7.4e13 multiply-accumulates, which the RTX 4050 does in a
couple of minutes. `molecules.tanimoto_matrix` stays as the straightforward
dense reference implementation; the tests check this module against it.

Two details make the fp16 arithmetic safe:

* ECFP4 fingerprints have 30-80 bits set out of 2048, and fp16 represents
  integers below 2048 exactly, so the intersection and popcount values are
  exact -- only the final division is approximate.
* Threshold tests avoid the division entirely. `sim >= c` is rearranged to
  `inter * (1 + c) >= c * (pq + pr)`, which compares exact integers. A
  similarity sitting exactly on the cutoff therefore resolves the same way
  here as it would in exact arithmetic, which matters because the cutoff is an
  ablation parameter and a borderline pair must not flip between runs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

# Query rows per tile, and reference columns per chunk inside a tile. The
# product bounds the largest temporary: 256 x 32,768 float32 is 34 MB, which
# leaves the 6 GB card room for the reference matrix itself (190,339 x 2048
# fp16 is 780 MB).
QUERY_TILE = 256
REF_CHUNK = 32_768


def pick_device(prefer_gpu: bool = True) -> torch.device:
    if prefer_gpu and torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


@dataclass
class FingerprintSet:
    """A fingerprint matrix resident on a device, with its popcounts.

    Holding this explicitly matters: the reference set is uploaded once and
    reused across every query tile, every regime and every seed. Re-uploading
    780 MB per call would dominate the runtime.
    """

    bits: torch.Tensor      # (n, n_bits) float16
    popcount: torch.Tensor  # (n,) float32
    device: torch.device

    @property
    def n(self) -> int:
        return int(self.bits.shape[0])


def to_device(fp: np.ndarray, device: torch.device) -> FingerprintSet:
    """Upload a dense uint8 bit matrix as fp16, with popcounts in fp32."""
    if fp.ndim != 2:
        raise ValueError(f"expected a 2-D bit matrix, got shape {fp.shape}")
    bits = torch.from_numpy(np.ascontiguousarray(fp)).to(
        device=device, dtype=torch.float16
    )
    popcount = bits.sum(dim=1, dtype=torch.float32)
    return FingerprintSet(bits=bits, popcount=popcount, device=device)


def _as_set(fp: np.ndarray | FingerprintSet, device: torch.device) -> FingerprintSet:
    if isinstance(fp, FingerprintSet):
        return fp
    return to_device(fp, device)


def _similarity_block(
    query: FingerprintSet,
    ref: FingerprintSet,
    q0: int,
    q1: int,
    r0: int,
    r1: int,
) -> torch.Tensor:
    """Exact-as-fp32 Tanimoto for one (query tile x reference chunk) block."""
    q_bits = query.bits[q0:q1]
    r_bits = ref.bits[r0:r1]
    inter = (q_bits @ r_bits.T).to(torch.float32)
    union = (
        query.popcount[q0:q1, None] + ref.popcount[None, r0:r1] - inter
    )
    # An all-zero fingerprint has no bits in common with anything; defining its
    # similarity as 0 rather than nan keeps the reductions below meaningful.
    return torch.where(union > 0, inter / union, torch.zeros_like(inter))


def _above_cutoff_block(
    query: FingerprintSet,
    ref: FingerprintSet,
    q0: int,
    q1: int,
    r0: int,
    r1: int,
    cutoff: float,
) -> torch.Tensor:
    """Boolean `similarity >= cutoff`, computed without dividing.

    See the module docstring: the rearrangement keeps every quantity an exact
    small integer, so a pair sitting on the cutoff does not depend on rounding.
    """
    q_bits = query.bits[q0:q1]
    r_bits = ref.bits[r0:r1]
    inter = (q_bits @ r_bits.T).to(torch.float32)
    total = query.popcount[q0:q1, None] + ref.popcount[None, r0:r1]
    return (inter * (1.0 + cutoff) >= cutoff * total) & (inter > 0)


def max_similarity(
    query_fp: np.ndarray | FingerprintSet,
    ref_fp: np.ndarray | FingerprintSet,
    *,
    device: torch.device | None = None,
    same_set: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Highest similarity from each query fingerprint to the reference set.

    Returns `(max_similarity, argmax_index)`, the index being a position in the
    reference set -- the audit reports which training molecule a test molecule
    is closest to, not only how close.

    `same_set=True` excludes the diagonal, for the case where query and
    reference are the same set and a molecule's similarity to itself (1.0)
    would otherwise win.
    """
    device = device or pick_device()
    query = _as_set(query_fp, device)
    ref = _as_set(ref_fp, device)

    best = np.zeros(query.n, dtype=np.float32)
    best_idx = np.full(query.n, -1, dtype=np.int64)

    for q0 in range(0, query.n, QUERY_TILE):
        q1 = min(q0 + QUERY_TILE, query.n)
        tile_best = torch.full((q1 - q0,), -1.0, device=device)
        tile_idx = torch.full((q1 - q0,), -1, dtype=torch.int64, device=device)
        for r0 in range(0, ref.n, REF_CHUNK):
            r1 = min(r0 + REF_CHUNK, ref.n)
            sim = _similarity_block(query, ref, q0, q1, r0, r1)
            if same_set:
                rows = torch.arange(q0, q1, device=device)[:, None]
                cols = torch.arange(r0, r1, device=device)[None, :]
                sim = torch.where(rows == cols, torch.full_like(sim, -1.0), sim)
            chunk_best, chunk_arg = sim.max(dim=1)
            improved = chunk_best > tile_best
            tile_best = torch.where(improved, chunk_best, tile_best)
            tile_idx = torch.where(improved, chunk_arg + r0, tile_idx)
        best[q0:q1] = tile_best.clamp(min=0.0).cpu().numpy()
        best_idx[q0:q1] = tile_idx.cpu().numpy()

    return best, best_idx


def neighbour_counts(
    fp: np.ndarray | FingerprintSet,
    cutoff: float,
    *,
    device: torch.device | None = None,
) -> np.ndarray:
    """How many other fingerprints each one is within `cutoff` of.

    Excludes self. This is the ordering Butina sphere exclusion needs: the
    densest compound becomes a cluster centre first.
    """
    device = device or pick_device()
    fps = _as_set(fp, device)
    counts = np.zeros(fps.n, dtype=np.int64)

    for q0 in range(0, fps.n, QUERY_TILE):
        q1 = min(q0 + QUERY_TILE, fps.n)
        tile = torch.zeros(q1 - q0, dtype=torch.int64, device=device)
        for r0 in range(0, fps.n, REF_CHUNK):
            r1 = min(r0 + REF_CHUNK, fps.n)
            mask = _above_cutoff_block(fps, fps, q0, q1, r0, r1, cutoff)
            rows = torch.arange(q0, q1, device=device)[:, None]
            cols = torch.arange(r0, r1, device=device)[None, :]
            mask = mask & (rows != cols)
            tile += mask.sum(dim=1)
        counts[q0:q1] = tile.cpu().numpy()

    return counts


def neighbours_above(
    query_fp: np.ndarray | FingerprintSet,
    ref_fp: np.ndarray | FingerprintSet,
    cutoff: float,
    *,
    device: torch.device | None = None,
) -> list[np.ndarray]:
    """Indices in the reference set within `cutoff` of each query fingerprint."""
    device = device or pick_device()
    query = _as_set(query_fp, device)
    ref = _as_set(ref_fp, device)

    out: list[np.ndarray] = []
    for q0 in range(0, query.n, QUERY_TILE):
        q1 = min(q0 + QUERY_TILE, query.n)
        tile_rows: list[np.ndarray] = []
        tile_cols: list[np.ndarray] = []
        for r0 in range(0, ref.n, REF_CHUNK):
            r1 = min(r0 + REF_CHUNK, ref.n)
            mask = _above_cutoff_block(query, ref, q0, q1, r0, r1, cutoff)
            hits = mask.nonzero()
            if hits.numel():
                hits_np = hits.cpu().numpy()
                tile_rows.append(hits_np[:, 0])
                tile_cols.append(hits_np[:, 1] + r0)

        if not tile_rows:
            out.extend(np.empty(0, np.int64) for _ in range(q1 - q0))
            continue

        # Grouped with one sort rather than a Python loop over hits: at a 0.4
        # cutoff a dense compound has thousands of neighbours, and iterating
        # them in Python costs more than the GPU pass that found them.
        rows_all = np.concatenate(tile_rows)
        cols_all = np.concatenate(tile_cols)
        order = np.lexsort((cols_all, rows_all))
        rows_all, cols_all = rows_all[order], cols_all[order]
        boundaries = np.searchsorted(rows_all, np.arange(q1 - q0 + 1))
        out.extend(
            cols_all[boundaries[i] : boundaries[i + 1]].astype(np.int64)
            for i in range(q1 - q0)
        )
    return out
