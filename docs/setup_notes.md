# Setup notes / Kurulum notları

## 1. Python environment

```bash
py -3.12 -m venv .venv
.venv/Scripts/activate          # Windows; on Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

The system default Python (3.14.6) is **not** usable: RDKit and PyTorch do not
ship reliable wheels for it. Python 3.12.10 is installed and every pinned
package publishes a cp312 wheel, so pip alone is enough — conda is not needed.

Sistemin varsayılan Python'u (3.14.6) kullanılamaz; RDKit ve PyTorch için
güvenilir tekerlek yok. Kurulu olan 3.12.10 kullanılır ve sabitlenen her paketin
cp312 tekerleği var, yani pip yeterli — conda'ya gerek yok.

PyTorch indirmesi tek başına ~2,5 GB; ChEMBL indirmesiyle aynı anda yapmayın,
ikisi de yavaşlar.

## 2. MMseqs2 — artık gerekmiyor / no longer required

**Kurulması gerekmiyor.** R3 bölmesi 565 hedef alanı arasındaki ikili dizi
kimliğini istiyor; bu 159.330 çift demek ve tam Smith-Waterman hizalamasıyla
makinenin tamamında bir dakikanın altında hesaplanıyor (`biopython`,
`src/dti/sequences.py`). MMseqs2'nin k-mer ön süzgeci bu ölçekte hiçbir şey
kazandırmıyor, karşılığında bir sezgisel getiriyor: ön süzgecin duyarlılığı
tam da %40 kimlik eşiğinin yakınında önem taşır ve kaçırılan bir çift, R3'te
eğitim hedefine benzer kalmış bir test hedefi demektir. Gerekçe: K22.

**No installation needed.** R3 needs pairwise identity between 565 target
domains -- 159,330 pairs, which exact alignment settles in under a minute
across the machine. MMseqs2's prefilter buys nothing at this size and its
recall matters exactly at the cutoff that defines the split.

`config/config.yaml` içindeki `tools.mmseqs_bin` anahtarı, ileride bir ablasyon
MMseqs2 ile karşılaştırma yapmak isterse diye bırakıldı; boru hattı onu
kullanmıyor.

## 3. GPU budget

Single NVIDIA RTX 4050 Laptop, **6 GB VRAM**. What this allows:

| Task | Feasible | Note |
| --- | --- | --- |
| ESM-2 650M frozen inference, fp16 | yes | batch 4 at ~400 residues; batch 1 at the 1022 limit |
| ESM-2 35M frozen inference | yes | ablation 5 |
| ESM-2 3B | **no** | out of scope on this hardware |
| ESM-2 fine-tuning | **no** | hence `freeze_protein: true` in M2 |
| ChemBERTa-2 partial fine-tune (last 2 layers) | yes | M2 |
| Shallow fusion head training | yes | runs on cached embeddings, CPU-feasible |

This is why the architecture is built on frozen embeddings written to disk: not
only a design preference but a hardware constraint. Worth stating as such in the
thesis.

Mimarinin donuk gömme üzerine kurulması bu donanım kısıtının da sonucu; tezde
böyle gerekçelendirilmeli.

## 4. Disk budget

| Item | Size |
| --- | --- |
| ChEMBL 37 SQLite archive | 5.8 GB (deleted after extraction by default) |
| ChEMBL 37 SQLite database | ~25 GB |
| PyTorch + CUDA wheels in `.venv` | ~3 GB |
| ESM-2 650M embeddings (1280-dim, ~5k targets) | < 100 MB |
| ChemBERTa-2 embeddings (384-dim, ~200k compounds) | ~300 MB |

Embeddings are small because they are extracted once per unique target/compound,
not per interaction record.

## 5. ChEMBL download

```bash
python scripts/00_download_chembl.py --list
python scripts/00_download_chembl.py --version 37
```

**Why the download is parallel.** EBI throttles each connection to roughly
1.5 MB/s regardless of local bandwidth, so a single stream makes the 5.8 GB
archive a ~1 hour affair. The server does honour byte ranges
(`Accept-Ranges: bytes`), so the script splits the remaining bytes across 8
connections by default. Use `--connections 1` to force a single stream.

An interrupted download resumes: only whole, size-verified ranges are appended
to the `.part` file, so re-running the script continues from where it stopped
rather than starting over.

**İndirme neden paralel.** EBI her bağlantıyı ~1,5 MB/s ile sınırlıyor, yerel
bant genişliği ne olursa olsun. Sunucu bayt aralıklarını desteklediği için
script kalan baytları varsayılan olarak 8 bağlantıya bölüyor. Kesilen indirme
kaldığı yerden devam eder.

The resolved release number and SHA-256 are appended to `docs/versions.md`
automatically.
