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
| biopython | 1.84 | exact Smith-Waterman for the R3 target split |
| ~~MMseqs2~~ | — | **artık kullanılmıyor / no longer used** (K22): 565 hedef için tam hizalama hem hızlı hem sezgisel içermiyor |
| LightGBM | 4.5.0 | M0 regressor |

## Model checkpoints

Pinned **2026-11-15** in `config/config.yaml` under
`models.m1_frozen_lm.*_revision`. Every loader must pass the revision
explicitly: `from_pretrained` without one follows the branch head, so a Hub
update would change every embedding and every number downstream while the code
and config still looked correct.

| Model | Hub id | Revision | Hub son değişiklik | Lisans |
| --- | --- | --- | --- | --- |
| ChemBERTa-2 | `DeepChem/ChemBERTa-77M-MTR` | `66b895cab8adebea0cb59a8effa66b2020f204ca` | 2022-01-20 | **beyan edilmemiş** |
| ESM-2 650M | `facebook/esm2_t33_650M_UR50D` | `08e4846e537177426273712802403f7ba8261b6c` | 2023-03-21 | MIT |
| ESM-2 35M (ablasyon) | `facebook/esm2_t12_35M_UR50D` | `6fbf070e65b0b7291e7bbcd451118c216cff79d8` | 2023-03-21 | MIT |

**Lisans uyarısı / Licence warning.** `DeepChem/ChemBERTa-77M-MTR` deposunda
beyan edilmiş bir lisans yok (Hub'da `license:` etiketi bulunmuyor). ESM-2'nin
ikisi MIT. K38 kontrol listesi kontrol noktası lisanslarının model kartına
yazılmasını ve **kullanım öncesi doğrulanmasını** şart koşuyor; ürün bu modelle
çıkarım yapıp türetilmiş veri dağıttığı için bu açık bir sorudur ve model
kartına "lisans beyan edilmemiş" olarak yazılmalı, MIT varsayılmamalı.

## Data sources

Filled in automatically by `scripts/00_download_chembl.py`. ChEMBL 37 is the
most recent release as of 2026-10-26.

## ChEMBL

- Release / Sürüm: **ChEMBL 37**
- Archive / Arşiv: `chembl_37_sqlite.tar.gz`
- SHA-256: `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`
- Downloaded / İndirilme tarihi: 2026-10-27
- Source / Kaynak: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/
