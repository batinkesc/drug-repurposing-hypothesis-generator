# Pinned versions / Sabitlenen sürümler

Bu dosya tezde birebir alıntılanır. / Quoted verbatim in the thesis.

Neden: RDKit'in kanonik atom sıralaması ve kanonik SMILES üretimi sürümler
arasında değişebiliyor; tautomer numaralandırma kodu da değişti. Bu projede
kanonik SMILES yalnızca bir gösterim değil, tekilleştirmenin ve R4b
kümelemesinin temeli — sürüm değişirse bölmeler değişir ve sonuçlar
tekrarlanamaz hale gelir. Aynı gerekçe veri kaynağı sürümleri için de geçerli:
ChEMBL sürümleri arasında aktivite kayıtları ve `predicted_binding_domains`
kapsamı farklılaşıyor.

Why: RDKit's canonical atom ranking and canonical SMILES output can differ
between releases, and the tautomer enumeration code has changed. Here canonical
SMILES is not a display format but the basis of deduplication and of the R4b
chemical split, so a version change changes the splits. The same argument
applies to data-source releases.

## Software

Authoritative pins live in [`requirements.txt`](../requirements.txt). Summary:

| Component | Version | Note |
| --- | --- | --- |
| Python | 3.12.10 | 3.14 lacks reliable RDKit / PyTorch wheels; 3.11 not installed |
| RDKit | 2024.03.5 | canonical SMILES + `MolStandardize` tautomer canonicalisation |
| PyTorch | 2.4.1 (cu124) | single RTX 4050 Laptop, 6 GB VRAM |
| transformers | 4.44.2 | |
| MMseqs2 | 15-6f452 | R3 clustering; installed separately, not a pip package |
| LightGBM | 4.5.0 | M0 regressor |

## Model checkpoints

Revision hashes are pinned in `config/config.yaml` under
`models.m1_frozen_lm.*_revision`. They are `null` until the first embedding run;
**fill them in before extracting embeddings**, otherwise a silently updated Hub
checkpoint would change every downstream number.

| Model | Hub id | Revision |
| --- | --- | --- |
| ChemBERTa-2 | `DeepChem/ChemBERTa-77M-MTR` | _to be pinned_ |
| ESM-2 650M | `facebook/esm2_t33_650M_UR50D` | _to be pinned_ |
| ESM-2 35M (ablation) | `facebook/esm2_t12_35M_UR50D` | _to be pinned_ |

## Data sources

Filled in automatically by `scripts/00_download_chembl.py`. ChEMBL 37 is the
most recent release as of 2026-10-26.

## ChEMBL

- Release / Sürüm: **ChEMBL 37**
- Archive / Arşiv: `chembl_37_sqlite.tar.gz`
- SHA-256: `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`
- Downloaded / İndirilme tarihi: 2026-10-27
- Source / Kaynak: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/
