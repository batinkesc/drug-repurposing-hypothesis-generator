"""Frozen language-model embeddings for the two towers.

This module exists because pooling is where embedding pipelines go wrong
quietly. Three specific traps, each of which produces plausible-looking
numbers:

1. **Unmasked mean pooling.** Averaging over padding positions makes a
   molecule's embedding depend on the longest molecule in its batch. Nothing
   crashes, nothing looks odd, and short inputs are corrupted worst -- so the
   damage is correlated with a chemical property. Every mean here is masked.

2. **Special tokens inside the mean.** ESM-2 wraps a sequence in `<cls>` and
   `<eos>`; ChemBERTa in `<s>` and `</s>`. Including them in a residue average
   dilutes short domains more than long ones, which is again a
   length-correlated bias. `mean` pooling averages *content* tokens only,
   identified through the tokenizer's own `special_tokens_mask` rather than by
   hard-coded positions, so it stays correct across both tokenizers.

3. **An unpinned revision.** `from_pretrained` without `revision=` follows the
   branch head. The config pins the commit hashes (K37) and `load_encoder`
   *requires* the revision as a positional argument -- there is no default, so
   a caller cannot forget it.

4. **`pooler_output`.** `AutoModel` attaches a pooler head to both checkpoints,
   and for both it is RANDOMLY INITIALISED -- transformers says so on load:
   "Some weights of EsmModel were not initialized ... ['esm.pooler.dense.weight']".
   Reading `pooler_output` would therefore return a random projection of the
   sequence, which looks like a perfectly ordinary embedding: finite, the right
   shape, different for different inputs. Everything here reads
   `last_hidden_state` and pools it explicitly; `pooler_output` is never used.

Pooling runs in float32 even though the model runs in float16: summing a few
hundred float16 values accumulates visible error, and this is the one number
that gets written to disk and reused by every model downstream.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

POOLINGS = ("mean", "cls")


@dataclass(frozen=True)
class EncoderBundle:
    """A loaded encoder with the identity of exactly what was loaded."""

    model: object
    tokenizer: object
    model_id: str
    revision: str
    hidden_size: int
    max_length: int

    def provenance(self) -> dict:
        return {
            "model_id": self.model_id,
            "revision": self.revision,
            "hidden_size": self.hidden_size,
            "max_length": self.max_length,
        }


def load_encoder(
    model_id: str,
    revision: str,
    *,
    max_length: int,
    dtype: str = "float16",
    device: torch.device | None = None,
) -> EncoderBundle:
    """Load a frozen encoder at an exact revision.

    `revision` is required, not optional: the whole point of pinning it in the
    config is that no code path can quietly resolve to the branch head.
    """
    from transformers import AutoModel, AutoTokenizer

    if not revision:
        raise ValueError(
            f"{model_id}: a pinned revision is required (see K37 and "
            f"config/config.yaml models.m1_frozen_lm.*_revision)"
        )
    device = device or (
        torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    )
    torch_dtype = getattr(torch, dtype)
    # float16 is a GPU optimisation; on CPU it is slower than float32 and some
    # kernels are unimplemented, so fall back rather than fail obscurely.
    if device.type == "cpu" and torch_dtype is torch.float16:
        torch_dtype = torch.float32

    tokenizer = AutoTokenizer.from_pretrained(model_id, revision=revision)
    model = AutoModel.from_pretrained(
        model_id, revision=revision, torch_dtype=torch_dtype
    )
    model.eval().to(device)
    for parameter in model.parameters():
        parameter.requires_grad_(False)

    return EncoderBundle(
        model=model,
        tokenizer=tokenizer,
        model_id=model_id,
        revision=revision,
        hidden_size=int(model.config.hidden_size),
        max_length=max_length,
    )


def pool_hidden_states(
    hidden: torch.Tensor,
    attention_mask: torch.Tensor,
    special_tokens_mask: torch.Tensor,
) -> dict[str, torch.Tensor]:
    """Mean over content tokens, plus the first token, both in float32.

    `mean` excludes padding AND special tokens; `cls` is position 0, which is
    `<cls>` for ESM-2 and `<s>` for ChemBERTa -- the dedicated sequence-level
    token, kept separately rather than mixed into the average.
    """
    hidden = hidden.to(torch.float32)
    content = (attention_mask.bool() & ~special_tokens_mask.bool()).to(torch.float32)
    counts = content.sum(dim=1, keepdim=True)
    if (counts == 0).any():
        raise ValueError(
            "an input tokenised to special tokens only; it has no content to "
            "average and would produce a division by zero"
        )
    mean = (hidden * content.unsqueeze(-1)).sum(dim=1) / counts
    return {"mean": mean, "cls": hidden[:, 0, :]}


def _length_sorted_batches(texts: list[str], batch_size: int) -> list[list[int]]:
    """Group indices into batches of similar length.

    Padding is wasted computation, and with masked pooling it changes nothing
    about the result, so inputs are grouped by length. The order is fully
    determined by (length, index), so the batching -- and therefore the
    floating-point summation order -- is reproducible across runs.
    """
    order = sorted(range(len(texts)), key=lambda i: (len(texts[i]), i))
    return [order[i : i + batch_size] for i in range(0, len(order), batch_size)]


def embed_texts(
    texts: list[str],
    bundle: EncoderBundle,
    *,
    batch_size: int,
    device: torch.device | None = None,
    progress_every: int = 50,
    label: str = "",
) -> tuple[dict[str, np.ndarray], dict]:
    """Embed every text, returning one array per pooling plus a report.

    Rows of the returned arrays line up with `texts`, not with the internal
    length-sorted order.
    """
    import time

    device = device or next(bundle.model.parameters()).device  # type: ignore[attr-defined]
    outputs = {
        name: np.zeros((len(texts), bundle.hidden_size), dtype=np.float16)
        for name in POOLINGS
    }
    batches = _length_sorted_batches(texts, batch_size)
    truncated = 0
    max_tokens = 0
    started = time.perf_counter()
    done = 0

    for step, indices in enumerate(batches, start=1):
        batch = [texts[i] for i in indices]
        encoded = bundle.tokenizer(  # type: ignore[operator]
            batch,
            padding=True,
            truncation=True,
            max_length=bundle.max_length,
            return_special_tokens_mask=True,
            return_tensors="pt",
        )
        special = encoded.pop("special_tokens_mask")
        lengths = encoded["attention_mask"].sum(dim=1)
        max_tokens = max(max_tokens, int(lengths.max()))
        # A sequence filling max_length exactly is the signature of truncation;
        # counted and reported rather than passing silently.
        truncated += int((lengths >= bundle.max_length).sum())

        encoded = {k: v.to(device) for k, v in encoded.items()}
        with torch.inference_mode():
            hidden = bundle.model(**encoded).last_hidden_state  # type: ignore[operator]
        pooled = pool_hidden_states(
            hidden, encoded["attention_mask"], special.to(device)
        )
        for name, values in pooled.items():
            outputs[name][indices] = values.cpu().numpy().astype(np.float16)

        done += len(indices)
        if step % progress_every == 0 or step == len(batches):
            elapsed = time.perf_counter() - started
            rate = done / elapsed if elapsed else 0.0
            remaining = (len(texts) - done) / rate if rate else 0.0
            print(
                f"\r  {label}{done:,} / {len(texts):,} "
                f"({100 * done / len(texts):.1f}%)  {rate:.0f}/s  "
                f"~{remaining / 60:.1f} min left",
                end="",
                flush=True,
            )
    print()

    report = {
        "n_inputs": len(texts),
        "n_batches": len(batches),
        "batch_size": batch_size,
        "max_tokens_seen": max_tokens,
        "n_truncated": truncated,
        "seconds": round(time.perf_counter() - started, 1),
        **bundle.provenance(),
    }
    return outputs, report
