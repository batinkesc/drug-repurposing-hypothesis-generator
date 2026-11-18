# Üçüncü taraf lisanslar ve atıf / Third-party licences and attribution

[`LICENSE`](LICENSE) yalnız bu depodaki **kodu** kapsar (MIT). Kodun indirdiği
veriyi ve o veriden türetilenleri kapsamaz. Bu dosya onları listeler.

`LICENSE` (MIT) covers the **code** in this repository only. It does not cover
the data the code downloads, nor anything derived from it. This file lists
those terms.

> Bu dosya `LICENSE`'tan ayrı tutuluyor çünkü GitHub'ın lisans algılayıcısı
> `LICENSE` dosyasında MIT metninden başka bir şey bulduğunda lisansı
> tanımıyor ve depoyu `NOASSERTION` olarak işaretliyor.

---

## Veri / Data

Bu depo **veri içermez.** Veri `scripts/00_download_chembl.py` ile kaynağından
indirilir; sabitlenen sürüm ve SHA-256 değerleri
[`docs/versions.md`](docs/versions.md) içinde.

| Kaynak | Sürüm | Lisans |
| --- | --- | --- |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) (EMBL-EBI) | 37 | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) |

**Önemli:** ChEMBL'in lisansı **atıf** gerektirir ve **share-alike** maddesi
taşır. ChEMBL'den türetilmiş veri dağıtırsanız — ölçüm tablosu, komşuluk
indeksi, gömme matrisi gibi — onu alanlara da aynı yeniden dağıtma hakkını
vermeniz gerekir. Bu koşul MIT'ten devralınmaz; MIT yalnız kodu kapsar.

Atıf metni, türetilmiş veri dağıtılan her yerde:

> Veri kaynağı: ChEMBL 37, EMBL-EBI, CC BY-SA 3.0.
> https://www.ebi.ac.uk/chembl/

---

## Model kontrol noktaları / Model checkpoints

Commit hash'leriyle sabitlenmiştir; değerler ve gerekçe
[`docs/versions.md`](docs/versions.md) içinde.

| Model | Hub id | Lisans |
| --- | --- | --- |
| ESM-2 650M | `facebook/esm2_t33_650M_UR50D` | MIT |
| ESM-2 35M | `facebook/esm2_t12_35M_UR50D` | MIT |
| ChemBERTa-2 | `DeepChem/ChemBERTa-77M-MTR` | **beyan edilmemiş / none declared** |

**ChemBERTa-2 açık sorusu:** Hugging Face deposunda beyan edilmiş bir lisans
yok (`license:` etiketi bulunmuyor). Bu, "izin verildiği" anlamına gelmez —
beyan yokluğu, hakların belirsiz olduğu anlamına gelir. MIT varsayılmamalı.
Model kartında ve tezde "beyan edilmemiş" olarak yazılır, ve dağıtım
planlanıyorsa kaynağa sorulması gerekir.

---

## Yazılım bağımlılıkları / Software dependencies

Sürümler [`requirements.txt`](requirements.txt) içinde sabitlenmiştir. Hepsi
kendi lisanslarıyla gelir; başlıcaları: RDKit (BSD-3-Clause), PyTorch
(BSD-3-Clause), transformers (Apache-2.0), biopython (Biopython License),
NumPy / pandas / SciPy / scikit-learn (BSD-3-Clause), LightGBM (MIT).
