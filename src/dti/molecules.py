"""Molecule standardisation and fingerprints.

The order of operations is a reproducibility contract, not an implementation
detail: strip salts, neutralise charge, pick a canonical tautomer, and only then
write canonical SMILES. Deduplication and the R4b chemical split both key on
that final string, so changing the order changes the splits and therefore the
results. The RDKit version is pinned for the same reason.

Everything configurable lives in `config/config.yaml` under `molecule`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

# RDKit logs a warning per unparseable record. At ChEMBL scale that is tens of
# thousands of lines that say nothing the drop counts do not already say.
RDLogger.DisableLog("rdApp.*")

# Tautomer canonicalisation is the expensive step and its cost grows with the
# number of enumerable tautomers. The default cap (1000) lets a handful of
# molecules dominate the runtime of the whole pipeline, so it is lowered. Hitting
# the cap is recorded rather than silently accepted: Canonicalize still returns
# a deterministic result for a given RDKit version and cap, which is what
# reproducibility requires.
MAX_TAUTOMERS = 250


@dataclass(frozen=True)
class StandardisationResult:
    """One molecule's outcome, including why it was dropped if it was."""

    canonical_smiles: str | None
    inchikey: str | None
    n_heavy_atoms: int | None
    reason: str | None = None  # None means kept

    @property
    def ok(self) -> bool:
        return self.canonical_smiles is not None


class Standardiser:
    """Applies the fixed standardisation sequence to SMILES strings.

    One instance is reused across a run: the RDKit standardisation objects are
    expensive to construct and hold internal caches.
    """

    def __init__(
        self,
        *,
        strip_salts: bool = True,
        neutralise_charge: bool = True,
        canonical_tautomer: bool = True,
        min_heavy_atoms: int = 6,
        max_heavy_atoms: int = 100,
        max_tautomers: int = MAX_TAUTOMERS,
    ) -> None:
        self.strip_salts = strip_salts
        self.neutralise_charge = neutralise_charge
        self.canonical_tautomer = canonical_tautomer
        self.min_heavy_atoms = min_heavy_atoms
        self.max_heavy_atoms = max_heavy_atoms

        self._fragment_chooser = rdMolStandardize.LargestFragmentChooser()
        self._uncharger = rdMolStandardize.Uncharger()

        params = rdMolStandardize.CleanupParameters()
        params.maxTautomers = max_tautomers
        self._tautomer_enumerator = rdMolStandardize.TautomerEnumerator(params)

    @classmethod
    def from_config(cls, config) -> Standardiser:
        settings = config["molecule"]
        return cls(
            strip_salts=settings["strip_salts"],
            neutralise_charge=settings["neutralise_charge"],
            canonical_tautomer=settings["canonical_tautomer"],
            min_heavy_atoms=settings["min_heavy_atoms"],
            max_heavy_atoms=settings["max_heavy_atoms"],
        )

    def standardise(self, smiles: str | None) -> StandardisationResult:
        """Run the sequence on one SMILES string.

        Returns a result carrying the canonical SMILES, or a reason the record
        was dropped. Reasons are aggregated into the drop table, so they are
        coarse on purpose -- one reason per stage, not per RDKit exception.
        """
        if not smiles or not isinstance(smiles, str):
            return StandardisationResult(None, None, None, "empty_smiles")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return StandardisationResult(None, None, None, "unparseable")

        try:
            # Cleanup does sanitisation, metal disconnection, normalisation and
            # reionisation in one pass. It must run before fragment choice:
            # otherwise a metal still bonded to the ligand looks like one
            # fragment and the "largest fragment" is the wrong molecule.
            mol = rdMolStandardize.Cleanup(mol)

            if self.strip_salts:
                mol = self._fragment_chooser.choose(mol)
            if self.neutralise_charge:
                mol = self._uncharger.uncharge(mol)
            if self.canonical_tautomer:
                mol = self._tautomer_enumerator.Canonicalize(mol)
        except Exception:
            # RDKit raises a variety of C++-backed exceptions here. Which one
            # it was does not change what we do: the record is dropped and
            # counted.
            return StandardisationResult(None, None, None, "standardisation_failed")

        if mol is None or mol.GetNumAtoms() == 0:
            return StandardisationResult(None, None, None, "empty_after_cleanup")

        n_heavy = mol.GetNumHeavyAtoms()
        if n_heavy < self.min_heavy_atoms:
            return StandardisationResult(None, None, n_heavy, "too_small")
        if n_heavy > self.max_heavy_atoms:
            # Peptide-like outliers: ChemBERTa-2 was not pretrained on them, so
            # their embeddings are not comparable to the rest.
            return StandardisationResult(None, None, n_heavy, "too_large")

        try:
            canonical = Chem.MolToSmiles(mol)
            inchikey = Chem.MolToInchiKey(mol) or None
        except Exception:
            return StandardisationResult(None, None, n_heavy, "serialisation_failed")

        if not canonical:
            return StandardisationResult(None, None, n_heavy, "empty_canonical_smiles")

        return StandardisationResult(canonical, inchikey, n_heavy, None)


def make_fingerprint_generator(radius: int = 2, n_bits: int = 2048):
    """ECFP generator used for R2/R4b clustering and as M0's molecule features.

    ECFP4 means radius 2. Uses the current generator API rather than the
    deprecated `GetMorganFingerprintAsBitVect`, which emits warnings and is
    scheduled for removal.
    """
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=n_bits
    )


def fingerprint(mol_or_smiles, generator) -> object | None:
    """Return the ECFP bit vector for a molecule, or None if it cannot be read."""
    mol = (
        Chem.MolFromSmiles(mol_or_smiles)
        if isinstance(mol_or_smiles, str)
        else mol_or_smiles
    )
    if mol is None:
        return None
    return generator.GetFingerprint(mol)


def fingerprint_matrix(smiles_list: list[str], generator) -> np.ndarray:
    """Dense uint8 matrix of fingerprints, one row per input SMILES.

    Dense on purpose: 200k x 2048 bits is ~400 MB as uint8, which fits in RAM
    and makes the pairwise Tanimoto work in the split code straightforward.
    Rows for unreadable SMILES are all zero; callers standardise first, so this
    should not happen in the pipeline.
    """
    rows = np.zeros((len(smiles_list), generator.GetOptions().fpSize), dtype=np.uint8)
    for i, smiles in enumerate(smiles_list):
        fp = fingerprint(smiles, generator)
        if fp is not None:
            rows[i] = np.frombuffer(
                bytes(fp.ToBitString(), "ascii"), dtype=np.uint8
            ) - ord("0")
    return rows


def tanimoto_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pairwise Tanimoto similarity between two dense fingerprint matrices.

    Used by the R4b post-split filter, which has to know the maximum similarity
    between every test molecule and the whole training set -- Butina clustering
    alone does not guarantee that bound.
    """
    a = a.astype(np.float32)
    b = b.astype(np.float32)
    intersection = a @ b.T
    popcount_a = a.sum(axis=1)[:, None]
    popcount_b = b.sum(axis=1)[None, :]
    union = popcount_a + popcount_b - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        similarity = np.where(union > 0, intersection / union, 0.0)
    return similarity
