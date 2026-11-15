# Aileler Arası Genellenebilir İlaç–Hedef Etkileşimi Modeli

## Cross-Family Generalisable Drug–Target Interaction Model

**Araştırma sorusu:** Önceden eğitilmiş dil modelleriyle kurulan bir ilaç–hedef
etkileşimi modeli, eğitimde hiç görmediği bir protein ailesine aktarılabilir mi
— ve aktarılamıyorsa bu nerede, ne kadar kırılıyor?

**Research question:** Can a drug–target interaction model built on pretrained
language-model representations transfer to a protein family it never saw during
supervised training — and if not, where and how badly does it break?

Tam tasarım: [`docs/00_project_design_tr.md`](docs/00_project_design_tr.md)

---

## Durum / Status

| Aşama / Stage | Durum |
| --- | --- |
| 0. Proje iskeleti, ortam, yapılandırma | ✅ tamam |
| 1. Veri hattı (çekme → süzme → temizlik → profil) | ⏳ sırada |
| 2. Bölmeler ve sızıntı denetimi (R1–R4b) | ⬜ |
| 3. Gömme çıkarma (ChemBERTa-2, ESM-2) | ⬜ |
| 4. Alt sınırlar (D0) ve M0 | ⬜ |
| 5. M1, M2 — rejim taraması (D1) | ⬜ |
| 6. Transfer matrisi (D2a, D2b) | ⬜ |
| 7. Teşhis ve yanlılık (D3–D5) | ⬜ |
| 8. Geri alma protokolü (D6) | ⬜ |
| 9. Yazım | ⬜ |

---

## Kurulum / Setup

```bash
py -3.12 -m venv .venv
.venv/Scripts/activate
pip install -r requirements.txt
```

MMseqs2 bir pip paketi değil, ayrı kurulur — yalnız R3 bölmesi için gerekli, o
yüzden 2. aşamaya kadar bekleyebilir. GPU ve disk bütçesiyle birlikte
[`docs/setup_notes.md`](docs/setup_notes.md) içinde.

Verinin indirilmesi (ChEMBL 37, ~6 GB arşiv, ~25 GB açılmış):

```bash
python scripts/00_download_chembl.py --list
python scripts/00_download_chembl.py --version 37
```

---

## Dizin yapısı / Layout

```
config/          config.yaml (tüm eşikler), families.yaml (aile tanımları)
src/dti/         kütüphane kodu
scripts/         numaralı hat aşamaları, sırayla çalıştırılır
data/raw/        indirilen kaynaklar — depoya girmez
data/interim/    aşama çıktıları — depoya girmez
data/processed/  temizlenmiş veri seti
data/splits/     R1–R4b bölmeleri
reports/tables/  tezdeki her tablo, CSV + Markdown (TR/EN)
reports/figures/ figürler, vektörel
docs/            tasarım, karar defteri, sürüm sabitleri, kurulum notları
```

### Tasarım ilkeleri / Design invariants

Bu kurallar kodda da yazılı, ama burada bir kez toplanıyor:

1. **Aile etiketi asla model girdisi değildir.** Yalnız değerlendirme bölmeleri
   kurmak için kullanılır. Girdi olsaydı model etkileşimi değil aile
   sınıflandırmasını öğrenirdi.
   *The family label is never a model input.*

2. **Her eşik `config/config.yaml`'de, gerekçesiyle birlikte durur.** `src/`
   içinde elle yazılmış hiçbir kesme değeri yok — tez bu sayıları alıntılamak
   zorunda ve bir kısmı ablasyon parametresi.
   *Every cutoff lives in config, with its reason next to it.*

3. **Test kümesi bir kez açılır**, tüm modeller ve rejimler hazır olduğunda.
   Baştan yazılı olması pazarlık konusu değil.
   *The test set is opened exactly once.*

4. **Alt sınır modelleri her rejimde yeniden kurulur**, yalnız o rejimin izin
   verdiği eğitim verisiyle. B1'i tüm veriden hesaplayıp aile-dışı rejimde
   kullanmak sızıntıdır.
   *Baselines are refitted per regime.*

5. **Sızıntı denetimi birincil sonuçtur, destekleyici istatistik değil.**
   Bölmeler arası özdeş molekül sayısı, en yüksek Tanimoto benzerliği, en
   yüksek dizi kimliği — bu rakamlar olmadan hiçbir skor yorumlanamaz.
   *Leakage audit numbers are a primary result.*

6. **Her aşama kendi çıktısını diske yazar** ve tekrar çalıştırıldığında aynı
   sonucu verir. Düşürülen kayıt sayısı her aşamada tabloya yazılır.
   *Every stage is idempotent and reports how many records it dropped.*

---

## Donanım / Hardware

Tek NVIDIA RTX 4050 Laptop, 6 GB VRAM. ESM-2 650M yalnızca donuk çıkarım için
kullanılabilir; ESM-2 3B ve ESM-2 ince ayarı bu donanımda kapsam dışı. Mimarinin
donuk gömme + sığ başlık olması bu kısıtın da sonucu.

Single RTX 4050 Laptop, 6 GB VRAM. ESM-2 650M is usable for frozen inference
only; ESM-2 3B and ESM-2 fine-tuning are out of scope on this hardware, which
is part of why the architecture is built on cached frozen embeddings.
