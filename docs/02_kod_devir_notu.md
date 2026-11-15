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
| 2 · Veri çekme ve süzme (1.080.040 aktivite) | Bitti |
| 3 · Molekül standardizasyonu (487.220 benzersiz molekül) | Bitti |
| 4 · Veri seti kurulumu (**283.353 hedef–bileşik çifti**) | Bitti |
| 5 · Bölmeler ve kaçak denetimi (13 rejim × 5 tohum = 65 dosya) | Bitti |
| 6 · Gömme çıkarımı | **Başlamadı** |
| 7 · Taban çizgileri ve modeller | Başlamadı |
| 8 · Deneyler ve teşhisler | Başlamadı |
| 9 · Uygulama ve tez | Başlamadı |

**Veri seti:** 565 hedef, 190.339 molekül, 5 aile, pChEMBL ortanca 6,92.
**Kod:** ~1.806 satır kütüphane (`src/dti/`), ~2.497 satır script, ~881 satır
test — **86 test geçiyor**.
**Disk:** 29 GB ham ChEMBL, 478 MB ara dosya, 12 MB veri seti, 87 MB bölmeler.

---

## 2. İlk üç iş — başka hiçbir şeyden önce

**1. Commit at.** Depoda hiç commit yok. Kod, scriptler, testler, `docs/`
içeri; veri `.gitignore`'a. Transfer çalışmasının durduğu noktayı
`v0-transfer-frozen` etiketiyle işaretle. Git geçmişi aynı zamanda tezin "ne
zaman ne yaptım" kaydı ve sonradan üretilemez.

**2. Hugging Face kontrol noktası commit hash'lerini sabitle.**
`config/config.yaml` içinde şu an `null`. Gömme çıkarımı bundan önce
koşturulmaz — sabitlenmeden çıkarılan gömme bütün aşağı akışı tekrarlanamaz
kılar (E1/K37).

**3. Güven skoru 8 mi 9 mu, doğrula.** Yerel ChEMBL 37'de `confidence_score`
dağılımını çıkar ve 8'in homolog tek-protein atamalarını içerip içermediğini
kontrol et. Mevcut süzgeç ≥8 kullanıyor ve gerekçesi "deneyin tek bir proteine
atandığı" diye yazılmış; ikincil kaynaklara göre 9 doğrudan, 8 homolog atama.
Sonucu ve kararı `decisions.md`'ye yaz. **Veri hattı yeniden koşturulmadan önce
karara bağlanmalı.**

---

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
gömme çıkarımı (2–3 saat)
  → B1/B2/B3 taban çizgileri (1–2 saat)
  → M0: ECFP4 + LightGBM (3–5 saat)
  → M1: donuk gömmeler + sığ ağ (4–8 saat)
  → 12 ARALIK KAPISI
```

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

Tasarım tarafı kapandı. Buradan sonrası kağıt değil kod: commit, hash'ler,
gömme çıkarımı, 12 Aralık kapısı. Yeni karar eklemek değil, mevcut kararları
gerçekten uygulamak.
