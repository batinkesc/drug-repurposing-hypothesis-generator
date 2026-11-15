"""Pooling correctness.

These embeddings are written to disk once and then read by every model, every
regime and every seed. A pooling bug would not crash anything -- it would just
make every number slightly wrong in a way correlated with sequence length, and
nothing downstream could detect it. So the invariants are pinned here rather
than trusted.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from dti.embeddings import (
    POOLINGS,
    _length_sorted_batches,
    load_encoder,
    pool_hidden_states,
)


def test_mean_ignores_padding():
    """The headline invariant: padding must not reach the average.

    Two identical content rows, different amounts of padding. If the mask were
    ignored, the padded row's mean would be pulled towards the pad vector.
    """
    hidden = torch.tensor(
        [
            [[1.0, 1.0], [3.0, 3.0], [99.0, 99.0]],   # third position is padding
            [[1.0, 1.0], [3.0, 3.0], [2.0, 2.0]],     # no padding
        ]
    )
    attention = torch.tensor([[1, 1, 0], [1, 1, 1]])
    special = torch.zeros_like(attention)
    pooled = pool_hidden_states(hidden, attention, special)
    assert pooled["mean"][0].tolist() == [2.0, 2.0]
    assert pooled["mean"][1].tolist() == [2.0, 2.0]


def test_mean_excludes_special_tokens():
    """A `<cls>`/`<eos>` wrapper must not enter a residue average.

    Including them would dilute short domains more than long ones -- a
    length-correlated bias, which is the worst kind here because domain length
    differs systematically between protein families.
    """
    hidden = torch.tensor(
        [[[100.0], [2.0], [4.0], [100.0]]]   # special, content, content, special
    )
    attention = torch.tensor([[1, 1, 1, 1]])
    special = torch.tensor([[1, 0, 0, 1]])
    pooled = pool_hidden_states(hidden, attention, special)
    assert pooled["mean"][0].item() == pytest.approx(3.0)


def test_cls_is_the_first_token_untouched():
    hidden = torch.tensor([[[7.0, 8.0], [1.0, 1.0], [1.0, 1.0]]])
    attention = torch.tensor([[1, 1, 1]])
    special = torch.tensor([[1, 0, 0]])
    pooled = pool_hidden_states(hidden, attention, special)
    assert pooled["cls"][0].tolist() == [7.0, 8.0]


def test_mean_is_invariant_to_batch_padding():
    """Same molecule, two different batch companions, same embedding.

    This is the property that makes length-sorted batching safe: if it failed,
    an embedding would depend on which batch it happened to land in, and the
    whole dataset would be irreproducible under a different batch size.
    """
    content = [[1.0], [2.0], [3.0]]
    alone = torch.tensor([content])
    padded = torch.tensor([content + [[0.0], [0.0]]])
    mask_alone = torch.tensor([[1, 1, 1]])
    mask_padded = torch.tensor([[1, 1, 1, 0, 0]])
    special_alone = torch.zeros_like(mask_alone)
    special_padded = torch.zeros_like(mask_padded)

    a = pool_hidden_states(alone, mask_alone, special_alone)["mean"]
    b = pool_hidden_states(padded, mask_padded, special_padded)["mean"]
    assert a.item() == pytest.approx(b.item())


def test_all_special_input_is_rejected_not_divided_by_zero():
    hidden = torch.zeros((1, 2, 4))
    attention = torch.tensor([[1, 1]])
    special = torch.tensor([[1, 1]])
    with pytest.raises(ValueError, match="no content"):
        pool_hidden_states(hidden, attention, special)


def test_pooling_is_done_in_float32():
    """float16 summation over a few hundred residues loses real precision, and
    this value is what gets written to disk and reused everywhere."""
    hidden = torch.full((1, 512, 8), 0.01, dtype=torch.float16)
    attention = torch.ones((1, 512), dtype=torch.long)
    special = torch.zeros((1, 512), dtype=torch.long)
    pooled = pool_hidden_states(hidden, attention, special)
    assert pooled["mean"].dtype is torch.float32
    assert pooled["mean"][0, 0].item() == pytest.approx(0.01, abs=1e-5)


def test_length_sorted_batches_cover_everything_once():
    texts = ["a" * n for n in (5, 1, 9, 3, 7, 2)]
    batches = _length_sorted_batches(texts, 2)
    flat = [i for batch in batches for i in batch]
    assert sorted(flat) == list(range(len(texts)))


def test_length_sorted_batches_are_deterministic():
    """Ties broken by index, so the floating-point summation order is fixed."""
    texts = ["ab", "cd", "ef", "gh"]   # all the same length
    assert _length_sorted_batches(texts, 2) == _length_sorted_batches(texts, 2)
    assert _length_sorted_batches(texts, 2) == [[0, 1], [2, 3]]


def test_load_encoder_refuses_an_unpinned_revision():
    """K37 as an assertion: no code path may resolve to the branch head."""
    with pytest.raises(ValueError, match="pinned revision is required"):
        load_encoder("DeepChem/ChemBERTa-77M-MTR", "", max_length=512)


def test_poolings_are_the_two_we_write():
    assert POOLINGS == ("mean", "cls")


@pytest.mark.parametrize("name", POOLINGS)
def test_pooling_outputs_have_the_hidden_dimension(name):
    hidden = torch.randn(3, 6, 11)
    attention = torch.ones((3, 6), dtype=torch.long)
    special = torch.zeros((3, 6), dtype=torch.long)
    special[:, 0] = 1
    pooled = pool_hidden_states(hidden, attention, special)
    assert pooled[name].shape == (3, 11)
    assert np.isfinite(pooled[name].numpy()).all()
