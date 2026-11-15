# Kod Tarafına Devir Notu

**Tarih:** 15 Kasım 2026
**Konum:** `docs/02_kod_devir_notu.md`
**Okuyucu:** bu depoda çalışmaya başlayan oturum (insan ya da ajan)
**Önce okunacaklar:** `docs/00_project_design_tr_rev2.md` (veri hattı),
`docs/00_proje_tasarimi_v2_urun.md` (ürün tasarımı), `docs/decisions.md`
(K1–K55), `docs/01_durum_ve_kapsam.md` (sayılar)

> **Düzeltme notu (15 Kasım 2026).** Bu belge depoya alınırken §4'teki iki
> madde depodaki gerçekle çelişiyordu ve düzeltildi: ortam Python 3.11 değil
> **3.12.10**, ve MMseqs2 artık **gerekli değil** (K22 ile tam Smith-Waterman
> hizalamasıyla değiştirildi). Orijinal ifadeler aşağıda `~~üstü çizili~~`
> bırakıldı. Ayrıca §6'daki tablonun son satırı K12'ye atıf yapıyordu; K12 bu
> depoda "sıfır aile çakışması" kaydıdır, test kümesi kuralı ise `training.
> test_set_locked` yapılandırmasına bağlıdır — satır ona göre düzeltildi.

---

## 0. Proje tek paragrafta

Bu depo artık "aileler arası genellenebilir ilaç–hedef etkileşimi" çalışması
değil. Aynı veri seti ve aynı hattın üzerine, **ilaç yeniden konumlandırma
hipotezi üreten bir araç** kuruluyor: bir molekül gir, tahmin edilen protein
hedeflerini güven işaretleriyle ve o hedeflerin hastalık ilişkileriyle birlikte
al. Aileler arası transfer çalışması **donduruldu** — dosyaları duruyor,
dokunulmuyor, sonraki bir çalışmanın giriş noktası.

Bilimsel iddia **geri bulma**dır, keşif değil: eğitimden çıkarılmış bilinen
ilaç–hedef çiftlerinin ne kadarı aday listesinin üstüne taşınıyor, ve bu oran
yenilik arttıkça nereye düşüyor.

---

## 1. Şu anki durum

| Aşama | Durum |
| --- | --- |
| 0 · Ortam ve veri kaynağı (ChEMBL 37, sürüm + SHA-256) | Bitti |
| 1 · ChEMBL envanteri ve aile eşlemesi (5 aile, 1.078 hedef) | Bitti |
| 1b · Protein temsili kararı (bağlanma alanı) | Bitti |
| 2 · Veri çekme ve süzme (739.699 aktivite, eşik 9 — K56) | Bitti |
| 3 · Molekül standardizasyonu (368.307 benzersiz molekül) | Bitti |
| 4 · Veri seti kurulumu (**283.353 hedef–bileşik çifti**) | Bitti |
| 5 · Bölmeler ve kaçak denetimi (13 rejim × 5 tohum = 65 dosya) | Bitti |
| 6 · Gömme çıkarımı (565×1280 protein, 190.339×384 molekül — K57) | Bitti |
| 7 · Taban çizgileri ve modeller | **Sıradaki** |
| 8 · Deneyler ve teşhisler | Başlamadı |
| 9 · Uygulama ve tez | Başlamadı |

**Veri seti:** 565 hedef, 190.339 molekül, 5 aile, pChEMBL ortanca 6,92.
**Kod:** ~2.050 satır kütüphane (`src/dti/`), ~2.730 satır script, ~1.000
satır test — **98 test geçiyor**.
**Disk:** 29 GB ham ChEMBL, 390 MB ara dosya, 12 MB veri seti, 87 MB bölmeler,
285 MB gömme.
**Git:** `master`, son commit `1f67354`, `v0-transfer-frozen` etiketi
`4434ea5`'te. Çalışma ağacı temiz.

---

## 2. Nerede kaldık (16 Kasım 2026)

Devir notunun ilk üç işi **kapandı**:

1. ✅ İlk commit `4434ea5` + `v0-transfer-frozen` etiketi.
2. ✅ HF kontrol noktası hash'leri sabitlendi (K37). Lisans uyarısı:
   `ChemBERTa-77M-MTR` deposunda **beyan edilmiş lisans yok**; model kartına
   "MIT" değil "beyan edilmemiş" yazılmalı ve K38 kontrol listesinde açık
   soru olarak duruyor.
3. ✅ Güven skoru doğrulandı ve **9'a çekildi** (K56). ChEMBL'in kendi
   `confidence_score_lookup` tablosundan: 8 homolog, 9 doğrudan atama.
   Beklenmedik bulgu — eşik zaten fiilen 9'du, çünkü ChEMBL
   `predicted_binding_domains`'i yalnız doğrudan atanmış deneyler için
   hesaplıyor. Veri hattı yeniden koşturuldu, veri seti **birebir aynı**
   kaldı, bölmeler geçerli.

Ardından 6. aşama da bitti (K57): gömmeler diskte, 1,1 dakikada çıktı, hiçbir
girdi kırpılmadı.

### Sıradaki iş: 7. aşama, taban çizgileri

`scripts/07_baselines.py` yok, yazılacak. Sıra:

```
B1 hedef ortalaması  →  B2 yalnız ligand  →  B3 yalnız protein  →  M0  →  M1
                                                        → 12 ARALIK KAPISI
```

**Başlamadan önce okunacak iki not:**

- **Taban çizgileri her rejimde yeniden fit edilir.** B1'i tüm veri üzerinde
  hesaplayıp rejimler arasında yeniden kullanmak kaçaktır
  (`config/config.yaml`, `models.baselines` yorumu).
- **Girdi ölçeklerini standartlaştır.** Norm ortancaları kollar arasında
  2,67–10,02 arasında değişiyor (K57). 384 boyutlu molekül vektörüyle 1280
  boyutlu protein vektörü yan yana eklenecekse, standartlaştırma olmazsa
  büyük normlu kol baskın gelir — bu bir mimari tercih değil, kaza olur.
- **E6 havuzlama ablasyonu geri geldi** (K57): `mean` ve `cls` dizileri de
  diskte. CLS, ECFP4 Tanimoto ile mean'in iki katı korele (0,437 vs 0,217),
  ama bu "daha iyi" demek değil — seçim M1'in doğrulama kümesinde yapılacak.

Gömmeler `.gitignore`'da ve 1,1 dakikada yeniden üretiliyor:
`python scripts/06_extract_embeddings.py --force`.

## 3. Değişmez kurallar

Bunlar üslup tercihi değil; ihlal edilirse üretilen sayı geçersiz olur. Her
biri `decisions.md`'de bir K kaydı.

1. **K27 — dil kuralı.** Ürün hiçbir çıktıda terapötik fayda iddia etmez. Çıktı
   birimi, yönü bilinmeyen bir bağlanma hipotezidir. "Tedavi eder",
   "etkilidir", "önerilir" kelimeleri arayüzde, API alan adlarında, model
   kartında, tezde ve sunumda geçmez.
2. **K40 — çıktı adı.** Model "tahmin edilen aktivite (pChEMBL)" üretir.
   "Bağlanma afinitesi" yalnız Ki/Kd alt kümesinden söz edilirken kullanılır
   (veride IC50 var).
3. **K42 — RT bölmesi tekilleştirme öncesi kayıtlardan kurulur.** Normalizasyon
   adımı aynı çiftin ölçümlerinin medyanını alıyor; tekilleştirilmiş veri
   setinde çiftin tarihi yok. Zamansal bölme o setin üzerine kurulamaz.
4. **K53 — bölme ölçütü "ilk kamuya açık raporlama tarihi"** (ChEMBL doküman
   tarihi), "ilk ölçüm tarihi" değil. Tarihi olmayan kayıtlar eğitim tarafına
   yazılır, test evrenine hiç alınmaz.
5. **K43/K44 — iki belirsizlik kaynağı birleştirilmez.** Tohum yayılımı ve
   split-conformal aralık ayrı gösterilir. Üretilen şey bir *tahmin
   aralığı*dır (prediction interval), güven aralığı değil. Gerçekleşen kapsama
   oranı her yenilik kovası için ayrı ölçülür; kalibre edilemeyen kovada aralık
   gösterilmez.
6. **K45/K54 — hedef güvenilirlik işareti kat dışı hesaplanır** ve minimum
   örnek sayısı şartına bağlıdır. Altında kalan hedefte gri işaret değil, ayrı
   bir "değerlendirme yetersiz" işareti çıkar.
7. **K46 — birleşik skor yok.** Hastalık ilişki skoru ile tahmin edilen
   aktivite çarpılmaz ya da toplanmaz; iki ayrı sıralama gösterilir.
8. **K50 — D6' ile RT'nin sırası sabit:** önce RT uygunluğu, sonra pozitif
   seçimi. Bir çift tek kez sayılır. Sonuçlar ayrı raporlanır.
9. **K52 — birincil metrik `Recall@%1`.** "Zenginleşme" dili kullanılmaz;
   veride gerçek inaktif yok.
10. **Test kümesi bir kez açılır**, tüm modeller ve rejimler hazır olduğunda.
11. **R4a / R4b kullanılmaz.** Bölme dosyaları duruyor ama bu çalışmada
    okunmuyor.
12. **Eşikler veriye bakılmadan sabitlenir** (E3'ün mantığı): RT kesim yılı,
    minimum pozitif sayısı, minimum kat dışı örnek sayısı, yenilik kovası
    sınırları. Sonradan seçilen eşik, katılık iddiasını yok eder.

---

## 4. Bu depoya özgü tuzaklar

- **Medyan alma tarihi yok ediyor.** (3) numaralı kuralın sebebi. RT için
  `activities` seviyesine geri dönmek gerekiyor.
- **pChEMBL yeniden hesaplanmaz, doğrulanır.** ChEMBL alanı hazır veriyor; 0,01
  toleransla karşılaştırılıyor ve ChEMBL'in değeri kullanılıyor (E2). Uyuşmayan
  kayıt sayısı düşürme tablosunda ayrı satır. *(Ölçüldü: 0 uyuşmazlık — K13.)*
- **RDKit işlem sırası tekrarlanabilirlik sözleşmesinin parçası:** tuz ayıklama
  → yük nötrleştirme → kanonik tautomer → kanonik SMILES yazımı. Tekilleştirme
  son dizgiye göre yapıldığı için sıra değişirse bölmeler de değişir.
- **Ortam:** **Python 3.12.10** + `venv`/`pip` (conda kurulu değil).
  ~~Python 3.11 (3.14 için RDKit/PyTorch tekerleği yok).~~
- **MMseqs2 gerekmiyor.** R3'ün dizi kimliği, 565 hedef için tam
  Smith-Waterman hizalamasıyla hesaplanıyor (`biopython==1.84`,
  `src/dti/sequences.py`) — K22. ~~MMseqs2'nin conda paketi Windows'ta yok —
  yalnız R3 için gerekli, statik ikili ya da WSL2.~~
- **Dizi kimliği tanımı tek:** `identity_over_shorter` (aynı kalıntı / kısa
  alanın uzunluğu). Yerel hizalama kimliği kullanılmaz — kullanıldığında
  denetim her rejim için 1,000 veriyordu (K25).
- **Donanım:** tek RTX 4050, 6 GB VRAM. ESM-2 650M donuk çıkarım evet, ince
  ayarı hayır. ESM-2 3B kapsam dışı. Ürün tarafı CPU'da çalışacak şekilde
  kurgulandı (hedef gömmeleri önbellekte).
- **Aile ataması tek kaynaktan** (ChEMBL protein sınıflandırma ağacı),
  `protein_class_desc` ön eki üzerinden; ayırıcı **iki boşluk** ve yol
  `short_name` alanından kuruluyor (K11). Yapı sürümler arasında
  değişebiliyor, envanter adımı atlanmaz.
- **Büyük Tanimoto geçişleri GPU'da** (`src/dti/similarity.py`): 1,8e10 çift
  bir fp16 matris çarpımı olarak, 0,5 dakika. `molecules.tanimoto_matrix`
  yoğun referans uygulaması olarak kalıyor ve testler GPU yolunu ona karşı
  doğruluyor (K24).

---

## 5. Teknik iş sırası

```
gömme çıkarımı              ✅ BİTTİ (tahmin 2–3 saat, gerçek 1,1 dakika)
  → B1/B2/B3 taban çizgileri (1–2 saat)        ← SIRADAKİ
  → M0: ECFP4 + LightGBM (3–5 saat)
  → M1: donuk gömmeler + sığ ağ (4–8 saat)
  → 12 ARALIK KAPISI
```

Gömme çıkarımının tahminden iki mertebe hızlı çıkması, kalan tahminleri de
yukarı yönlü şüpheli kılıyor: ChemBERTa-2 beklediğimden çok küçük (3 katman,
384 boyut; "77M" eğitim kümesinin boyutu, modelin değil) ve M1'in başlığı
önceden hesaplanmış gömmeler üzerinde çalıştığı için dil modeline hiç
dokunmuyor. M1'in 65 koşusu saatler değil dakikalar sürebilir.

**Kapı sorusu:** M1, B1–B3'ü ve M0'ı geçiyor mu?

- Geçiyorsa plan olduğu gibi devam eder.
- Geçmiyorsa ürünün vaadi "model tahmini" yerine "kanıt paneli + erişim aracı"
  olarak yeniden yazılır. Mimari ve kod aynı kalır; sıralamanın kaynağı ve
  arayüz dili değişir. Bu tezi zayıflatmaz, ölçülmüş bir bulgu olarak yazılır.

Kararı o tarihte ver; uygulama geliştirmesi bu karardan sonra başlıyor.

---

## 6. Her karar bir test

Yeni karar eklemeye gerek yok; eksik olan kararların kodda uygulanması.
Önerilen testler, karar numarası test adında geçecek şekilde:

| Karar | Test | Ne doğruluyor |
| --- | --- | --- |
| K42 | `test_k42_rt_no_pre_cutoff_pairs` | RT test kümesinde kesim öncesi raporlanmış çift yok |
| K42 | `test_k42_rt_built_from_raw_activities` | RT bölmesi tekilleştirilmiş veri setinden üretilmemiş |
| K53 | `test_k53_undated_records_excluded_from_test` | Tarihsiz kayıtlar test evrenine girmemiş |
| K45 | `test_k45_reliability_is_out_of_fold` | Güvenilirlik tablosu kat dışı tahminlerden üretilmiş |
| K54 | `test_k54_reliability_min_samples` | Minimum altı hedeflerde "değerlendirme yetersiz" işareti var |
| K44 | `test_k44_coverage_reported_per_bucket` | Her kova için gerçekleşen kapsama hesaplanmış |
| K44 | `test_k44_no_interval_when_uncalibrated` | Kalibre edilemeyen kovada aralık üretilmiyor |
| K50 | `test_k50_d6_positives_within_rt_universe` | D6' pozitifleri RT-uygun evrenden seçilmiş |
| K50 | `test_k50_pair_counted_once` | Aynı çift birden çok kayıt üzerinden tekrar sayılmıyor |
| K46 | `test_k46_no_composite_score` | Hastalık skoru ile aktivite tahmini tek sayıda birleşmiyor |
| §3·10 (`training.test_set_locked`) | `test_test_set_opened_once` | Test kümesi erişim sayacı 1 |

Mevcut 86 testin yanına eklenirler.

---

## 7. Dokunulmayacaklar

- `R4a` / `R4b` bölme dosyaları — dondurulmuş çalışmanın girdisi.
- `docs/00_project_design_tr_rev2.md`'nin veri hattı bölümleri — geçerli,
  yeniden türetilmiyor.
- Yedi hazır tablo (`reports/tables/01_chembl_inventory.md` …
  `07_leakage_audit.md`) — script'lerden yeniden üretilebilir, elle
  düzenlenmez.
- `docs/versions.md` — indirme scripti yazıyor, elle değiştirilmez.

---

## 8. Karara bağlanmadan ilgili kod yazılmayacak açık sorular

- [x] ~~Güven skoru eşiği 8 mi 9 mu~~ — **kapandı (K56): 9'a çekildi.** Veri
      hattı yeniden koşturuldu; veri seti birebir aynı kaldı, 65 bölme dosyası
      ve kaçak denetimi geçerli. `03_drop_table.md` ve
      `05_drop_table_stage2.md` yeniden üretildi
- [ ] Ölçüm envanteri: Ki/Kd ile IC50 aynı model hedefinde mi kalıyor (envanter
      deneysel heterojenliği de ölçmeli — K55)
- [ ] RT kesim yılı, ve tekilleştirme öncesi kayıtlardan kurulan testte kaç
      çift kalıyor
- [ ] Geri alma havuzunun ve gerçek pozitiflerin tam tanımı; D6'nın iki kipi
      için minimum pozitif sayısı
- [ ] Konformal aralık hangi örnekler üzerinde kalibre edilecek, hangi kovalarda
      geçerli sayılacak
- [ ] Bire-çok ve çok bileşenli hedef eşlemesinde seçim kuralı
- [ ] SMILES tabanlı sistemde biyolojik ilaçlar, tuzlar, karışımlar ve ana
      molekül eşlemesi
- [ ] Aşama 2 kapsamına (tüm insan tek-protein hedefleri) geçilecek mi — karar
      1. aşama sonuçlarından sonra

---

## 9. Hatırlatma

Tasarım tarafı kapandı ve açılış işleri de bitti: commit, hash'ler, güven
skoru kararı, gömme çıkarımı. Kalan tek şey 12 Aralık kapısına kadar olan
modelleme zinciri — **taban çizgileri, M0, M1** — ve orada yeni karar
eklemek değil, mevcut kararları gerçekten uygulamak gerekiyor.

Kapıya kadarki en kritik tek madde: **taban çizgilerini atlamamak.** B1–B3
en ucuz ve en belirleyici parça; onlar olmadan M1'in skoru yorumlanamaz ve
12 Aralık kapısının sorusu ("M1, B1–B3'ü ve M0'ı geçiyor mu?") sorulamaz.
