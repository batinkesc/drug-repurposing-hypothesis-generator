# İlaç Yeniden Konumlandırma Hipotezi Üreteci

## Drug Repurposing Hypothesis Generator

ChEMBL biyoaktivite verisi üzerine kurulu bir ilaç–hedef aktivite tahmin modeli
ve onun üzerine, tahmin edilen hedefleri hastalık ilişkileriyle birleştirip
**sıralanmış ve belirsizliği işaretli bir hipotez listesi** üreten bir araç.
Lisans bitirme projesi.

A drug–target activity model built on ChEMBL bioactivity data, and on top of it
a tool that combines predicted targets with disease associations to produce a
**ranked, uncertainty-flagged list of hypotheses**. Undergraduate final-year
project.

---

## ⚠️ Ne iddia edilmiyor / What this does NOT claim

Bu bölüm projenin en önemli kısmı ve kodda da bir karar kaydına bağlı (K27).

> **Araç hiçbir çıktıda terapötik fayda iddia etmez. Çıktı birimi, yönü
> bilinmeyen bir bağlanma hipotezidir.**

- **Terapötik fayda değil.** Bir hedefe bağlanmanın tahmin edilmesi, bir ilacın
  bir hastalıkta işe yaradığı anlamına gelmez.
- **Etki yönü değil.** Afinite "ne kadar sıkı bağlanıyor"u söyler; aktive mi
  inhibe mi ettiğini söylemez. Yön yalnız *bilinen* çiftler için, ChEMBL
  mekanizma tablosundan gösterilir.
- **Nedensellik değil.** Hedef–hastalık ilişkisi "bu hedefi vurmak hastalığı
  iyileştirir" demek değildir.
- **"Bağlanmaz" yargısı değil.** Veride gerçek inaktif kayıt yok; ölçülmemiş
  bir çift negatif değil, bilinmiyordur.
- **Yeni kullanım keşfi değil.** Bilimsel iddia **geri bulma**dır: eğitimden
  çıkarılmış, *bilinen* ilaç–hedef çiftlerinin ne kadarı aday listesinin üstüne
  taşınıyor ve bu oran kimyasal yenilik arttıkça nereye düşüyor.

Araç bir **araştırma ve eğitim amaçlı triyaj aracıdır**, klinik karar desteği
değildir.

**In English:** the tool never claims therapeutic benefit. Its output unit is a
binding hypothesis of unknown direction. The scientific claim is *retrieval* of
known drug–target pairs held out from training, not discovery of new uses. It is
a research and teaching triage tool, not clinical decision support.

---

## Durum / Status

| Aşama / Stage | Durum |
| --- | --- |
| 0 · Ortam ve veri kaynağı (ChEMBL 37, sürüm + SHA-256) | ✅ |
| 1 · ChEMBL envanteri ve aile eşlemesi | ✅ |
| 1b · Protein temsili kararı (bağlanma alanı) | ✅ |
| 2 · Veri çekme ve süzme (739.699 aktivite) | ✅ |
| 3 · Molekül standardizasyonu (368.307 benzersiz molekül) | ✅ |
| 4 · Veri seti kurulumu (**283.353 hedef–bileşik çifti**) | ✅ |
| 5 · Bölmeler ve kaçak denetimi (13 rejim × 5 tohum) | ✅ |
| 6 · Gömme çıkarımı (ChemBERTa-2 + ESM-2 650M) | ✅ |
| 7 · Taban çizgileri (B1–B3), M0, M1 | ⏳ sırada |
| 8 · Deneyler, teşhisler, geri alma protokolü | ⬜ |
| 9 · Uygulama ve tez | ⬜ |

**Veri seti:** 283.353 hedef–bileşik çifti, 565 hedef, 190.339 molekül, 5
protein ailesi, pChEMBL ortanca 6,92.
**Test:** 98 test geçiyor.

### Şimdiye kadarki bulgular / Findings so far

Model eğitilmeden, yalnız veri hattından çıkan sonuçlar:

- **"Görülmemiş bileşik" bölmesi iddia ettiği kadar zor değil.** Kümeleme
  tabanlı bölmede (R2) test bileşiklerinin **%81,7'si** eğitim kümesinde ≥0,40
  Tanimoto'da bir komşuya sahip, ortanca benzerlik 0,548 (beş tohumda
  0,805–0,828). Butina kümelemesi kümeler *arası* benzerliğe sınır koymaz, ve
  kaçak denetimi olmadan bu görünmez olurdu.
- **Aile dışı test verisinin yaklaşık yarısı kimyasal olarak tanıdık.**
  Kimyasal yenilik süzgeci aile başına kayıtların %40–56'sını atıyor.
- **Protein gömmeleri aile yapısını, aile hiç girdi olmadan kurtarıyor.**
  Aile içi kosinüs ortanca 0,871, aileler arası 0,782.
- **Güven skoru eşiği, bağlanma alanı temsili yüzünden zaten fiilen 9'du:**
  ChEMBL `predicted_binding_domains`'i yalnız doğrudan atanmış deneyler için
  hesaplıyor, dolayısıyla homolog (skor 8) kayıtların hiçbiri alan eşlemesi
  almıyor.

---

## Kurulum / Setup

```bash
py -3.12 -m venv .venv
.venv/Scripts/activate
pip install -r requirements.txt
```

MMseqs2 **gerekmiyor**: dizi kimliği 565 hedef için tam Smith-Waterman
hizalamasıyla hesaplanıyor (`biopython`). Donanım ve disk bütçesi
[`docs/setup_notes.md`](docs/setup_notes.md) içinde.

Verinin indirilmesi (ChEMBL 37, ~6 GB arşiv, ~30 GB açılmış):

```bash
python scripts/00_download_chembl.py --list
python scripts/00_download_chembl.py --version 37
```

Hattın sırayla koşturulması:

```bash
python scripts/01_inspect_chembl.py
python scripts/02_extract.py
python scripts/03_standardise_molecules.py
python scripts/04_build_dataset.py
python scripts/05_build_splits.py
python scripts/06_extract_embeddings.py
```

---

## Dizin yapısı / Layout

```
config/          config.yaml (tüm eşikler, gerekçeleriyle), families.yaml
src/dti/         kütüphane kodu
scripts/         numaralı hat aşamaları, sırayla çalıştırılır
tests/           98 test — havuzlama, benzerlik, bölme, metrik değişmezleri
data/            indirilen ve türetilen veri — depoya girmez (.gitignore)
reports/tables/  tezdeki tablolar (TR/EN)
docs/            tasarım, karar defteri (K1–K57), sürüm sabitleri
```

---

## Tasarım ilkeleri / Design invariants

Kodda da yazılı, burada bir kez toplanıyor:

1. **Her eşik `config/config.yaml`'de, gerekçesiyle durur.** `src/` içinde elle
   yazılmış kesme değeri yok — tez bu sayıları alıntılamak zorunda.
2. **Test kümesi bir kez açılır**, tüm modeller ve rejimler hazır olduğunda.
3. **Taban çizgileri her rejimde yeniden fit edilir**, yalnız o rejimin izin
   verdiği eğitim verisiyle. B1'i tüm veriden hesaplayıp başka rejimde
   kullanmak kaçaktır.
4. **Kaçak denetimi birincil sonuçtur**, destekleyici istatistik değil.
   Bölmeler arası özdeş molekül sayısı, en yüksek Tanimoto ve en yüksek dizi
   kimliği olmadan hiçbir skor yorumlanamaz.
5. **Eşikler veriye bakılmadan sabitlenir.** Sonradan seçilen eşik, katılık
   iddiasını yok eder.
6. **Her aşama kendi çıktısını diske yazar**, tekrarlanabilir, ve kaç kayıt
   düşürdüğünü tabloya yazar.
7. **Aile etiketi model girdisi değildir** — yalnız gruplandırma ve
   değerlendirme için kullanılır.

Her kararın gerekçesi [`docs/decisions.md`](docs/decisions.md) içinde numaralı
kayıt olarak duruyor (K1–K57), ölçümleriyle birlikte — ikisi kendi hatamı
düzeltme kaydı.

---

## Donanım / Hardware

Tek NVIDIA RTX 4050 Laptop, **6 GB VRAM**. ESM-2 650M yalnız donuk çıkarım için
kullanılabilir; ESM-2 3B ve ESM-2 ince ayarı kapsam dışı. Mimarinin donuk gömme
+ sığ başlık olması bu kısıtın sonucu — ve aynı zamanda ürün tarafının CPU'da
çalışabilmesinin sebebi, çünkü 565 hedef gömmesi 1,4 MB'lık bir dosya.

---

## Veri kaynağı, lisans ve atıf / Data, licence and attribution

Bu depo **veri içermez**; veri `scripts/00_download_chembl.py` ile kaynağından
indirilir. Sabitlenen sürümler ve SHA-256 değerleri
[`docs/versions.md`](docs/versions.md) içinde.

- **ChEMBL 37** — EMBL-EBI, [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/).
  ChEMBL'den türetilmiş veri dağıtılırsa aynı lisans koşulları geçerlidir
  (share-alike).
- **ESM-2** (`facebook/esm2_t33_650M_UR50D`) — MIT.
- **ChemBERTa-2** (`DeepChem/ChemBERTa-77M-MTR`) — **depoda beyan edilmiş
  lisans yok**; MIT varsayılmamalı, açık soru olarak izleniyor.

Model kontrol noktaları commit hash'leriyle sabitlenmiştir; gerekçesi ve
değerleri `docs/versions.md` içinde.

---

## Dondurulan çalışma / Frozen earlier work

Proje önce **aileler arası genellenebilirlik** sorusu olarak tasarlandı: bir
DTI modeli eğitimde hiç görmediği bir protein ailesine aktarılabilir mi? O soru
iptal edilmedi, **donduruldu** — altyapısı (aile dışı R4a/R4b bölmeleri, 65
bölme dosyası, transfer deney tasarımı) olduğu gibi duruyor.

Durduğu nokta `v0-transfer-frozen` etiketiyle işaretli. Tasarımı
[`docs/00_project_design_tr_rev2.md`](docs/00_project_design_tr_rev2.md),
ürüne geçiş gerekçesi
[`docs/00_proje_tasarimi_v2_urun.md`](docs/00_proje_tasarimi_v2_urun.md)
içinde.
