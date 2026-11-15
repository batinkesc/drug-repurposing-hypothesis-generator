# Karar defteri / Decision log

Her karar: ne, neden, nerede uygulandı. Tasarım dokümanından (`00_project_design_tr.md`)
sapmalar ve ona yapılan eklemeler burada izlenir.

Each entry: what, why, where it is implemented. Deviations from and additions
to the design document are tracked here.

---

## 2026-10-26 — Kurulum kararları / Setup decisions

### K1. ChEMBL tam SQLite dağıtımı, web API değil
**TR:** Veri ChEMBL'in numaralı sürüm SQLite dağıtımından çekilir.
**EN:** Data comes from a numbered ChEMBL SQLite release, not the web API.

**Neden / Why:** Sürüm sabitlemesi tekrarlanabilirlik iddiasının parçası; web
API sorguları sürüm damgası taşımaz ve `predicted_binding_domains` tablosuna
erişim kısıtlı. Ayrıca offline çalışma on binlerce sorguyu saatler yerine
dakikalara indirir.

**Nerede / Where:** `scripts/00_download_chembl.py`, `docs/versions.md`

### K2. Veri kaynağı sürümleri de sabitlenir
**TR:** ChEMBL sürüm numarası ve arşiv SHA-256'sı, UniProt release tarihi
`docs/versions.md`'ye yazılır.
**EN:** The ChEMBL release number with the archive SHA-256, and the UniProt
release date, are recorded in `docs/versions.md`.

**Neden / Why:** Tasarım dokümanı RDKit/transformers/PyTorch sürümlerini
sabitliyor ama veri kaynağının sürümünü değil. ChEMBL sürümleri arasında
aktivite kayıtları ve PPDMs alan eşleme kapsamı değişiyor — bu boşluk
tekrarlanabilirlik iddiasını eksik bırakıyordu. *(Tasarım dokümanına ek.)*

**Nerede / Where:** `scripts/00_download_chembl.py:record_version`

### K3. Python 3.12.10 + venv/pip, conda değil
**TR:** Ortam `py -3.12 -m venv` ile kurulur, paketler `requirements.txt`'ten
pip ile gelir. conda kullanılmaz.
**EN:** A `venv` on Python 3.12.10 with pip-installed pinned requirements;
conda is not used.

**Neden / Why:** Sistemin varsayılan Python'u 3.14.6 ve RDKit/PyTorch için
güvenilir tekerlek yok. İlk tasarımda conda öngörülmüştü ama makinede conda
kurulu değil, kurulu olan 3.12.10 ise sabitlenen her paketin cp312 tekerleğini
destekliyor — RDKit, PyTorch cu124, LightGBM dahil. conda'nın tek gerçek
avantajı MMseqs2'yi paket olarak vermesiydi; o da Windows'ta zaten ayrı binary
olarak kurulacağı için avantaj kalmadı. İki ortam dosyası tutmak yerine tek
yol: pip.

**Nerede / Where:** `requirements.txt`, `docs/setup_notes.md`
*(Bu karar `environment.yml`'in yerini aldı; o dosya kaldırıldı.)*

### K3b. İndirme paralel ve kaldığı yerden devam eder
**TR:** ChEMBL arşivi varsayılan olarak 8 bağlantıyla, bayt aralıkları
bölüştürülerek indirilir; kesilen indirme `.part` dosyasından devam eder.
**EN:** The ChEMBL archive is fetched over 8 byte-range connections by default
and an interrupted download resumes from the `.part` file.

**Neden / Why:** Ölçüldü: EBI her bağlantıyı ~1,5 MB/s ile sınırlıyor, yerel
bant genişliğinden (1 Gbit) bağımsız olarak. Tek akışla 5,8 GB arşiv ~1 saat
sürüyordu ve indirme kesilirse baştan başlamak gerekiyordu. Sunucu
`Accept-Ranges: bytes` döndürdüğü için kalan baytlar bölüştürülebiliyor.
Resume kuralı: `.part` dosyasına yalnız tam ve boyutu doğrulanmış aralıklar
eklenir, yani yarım bir parça asla diske yazılmış sayılmaz.

**Nerede / Where:** `scripts/00_download_chembl.py:download`

### K4. Kod İngilizce, raporlar iki dilli
**TR:** Değişken/fonksiyon isimleri ve kod yorumları İngilizce; tablolar,
ölçüm çıktıları, loglar ve dokümanlar hem Türkçe hem İngilizce.
**EN:** Code and code comments in English; tables, measurement outputs, logs
and docs in both Turkish and English.

**Neden / Why:** Tez iki dilde yazılacak. Çıktıları baştan iki dilli üretmek,
sonradan tüm tabloları çevirmekten ucuz.

**Nerede / Where:** `config.yaml: project.report_languages`,
`config.py:family_label`, `families.yaml: name_tr / name_en`

### K5. Rapor çıktıları tez şablonundan bağımsız
**TR:** Tüm tablolar `reports/tables/` içinde ayrı CSV + Markdown; figürler
`reports/figures/` içinde vektörel. Tez metnine gömülü üretilmez.
**EN:** Tables as standalone CSV + Markdown, figures as vector files; never
generated inline into the thesis text.

**Neden / Why:** Bölümün zorunlu tez formatı henüz verilmedi, sonradan
gelebilir. Şablon geldiğinde metin taşınır, tablo ve figürler yeniden
üretilmez.

**Nerede / Where:** `config.yaml: paths.tables / paths.figures`

---

## Tasarım dokümanına eklenen ölçüm kararları / Additions to the design

Bunlar tasarım dokümanında yoktu; gerekçeleri ilgili yapılandırma satırının
yanında da duruyor.

### K6. R4'te RMSE üç biçimde raporlanır
**TR:** Ham, ortalama-kaydırmalı, ve lineer kalibrasyonlu.
**EN:** Raw, mean-shifted, and linearly calibrated.

**Neden / Why:** Aileler pChEMBL dağılımının ortalaması bakımından farklı.
Dışarıda bırakılan aileye geçişte model sistematik bir kayma yapar: RMSE
patlar ama Spearman korunabilir. Tek bir ham RMSE "transfer çöktü" der, oysa
sıralama bilgisi aktarılmış, yalnız ölçek kaymış olabilir. Üç varyant arasındaki
fark, transferin ne kadarının sıralama ne kadarının ölçek olduğunu doğrudan
verir — ve yeniden konumlandırma iddiası (D6) için asıl önemli olan sıralama.

**Nerede / Where:** `config.yaml: evaluation.r4_rmse_variants`

### K7. Karıştırma kontrolü (D3) iki varyantlı
**TR:** Global karıştırma + aile içi karıştırma.
**EN:** Global shuffling plus within-family shuffling.

**Neden / Why:** Protein gömmesini rastgele *herhangi* bir proteinle
değiştirmek, model yalnızca aileyi ayırt ediyorsa bile skoru düşürür — yani
testi geçmek kolay. Aile içi karıştırmada performans düşmüyorsa model proteini
değil aileyi okuyordur, ki bu doğrudan çalışmanın ana sorusu. Ek maliyeti yok.

**Nerede / Where:** `config.yaml: evaluation.diagnostics.shuffle_control.variants`

### K8. R4b için minimum örneklem eşiği baştan sabit
**TR:** 0,4 Tanimoto süzmesinden sonra 200 kayıttan az kalırsa o aile için R4b
skorlanmaz, yalnız kalan kayıt sayısı bildirilir.
**EN:** If fewer than 200 records survive the similarity filter, R4b is not
scored for that family; only the surviving count is reported.

**Neden / Why:** Süzme bazı ailelerde test kümesini birkaç düzineye indirebilir;
o skorun güven aralığı anlamsız olur. Eşik veriye bakılmadan sabitlenmezse
sonradan seçilmiş eşik olur.

**Nerede / Where:** `config.yaml: splits.r4_family.r4b_min_test_records`

### K9. Havuzlama yöntemi ablasyon listesine eklendi
**TR:** Ortalama havuzlama vs CLS token (hem ESM-2 hem ChemBERTa-2).
**EN:** Mean pooling vs CLS token, both towers.

**Neden / Why:** Donuk gömmelerle bedava bir ablasyon ve etkisi bazen ESM-2
650M vs 35M farkından büyük.

**Nerede / Where:** `config.yaml: models.m1_frozen_lm.pooling`

### K10. pChEMBL yeniden hesaplanmaz, doğrulanır
**TR:** ChEMBL'in `pchembl_value` alanı kullanılır, kendi dönüşümümüzle 0,01
toleransla karşılaştırılır.
**EN:** ChEMBL's own `pchembl_value` is used and cross-checked against our
conversion within 0.01.

**Neden / Why:** Sessizce yeniden hesaplamak yerine doğrulamak hem hızlı hem
denetlenebilir; uyuşmayan kayıtlar düşürme tablosunda görünür.

**Nerede / Where:** `config.yaml: normalisation.use_chembl_pchembl`

---

---

## 2026-10-27 — ChEMBL 37 envanteri / ChEMBL 37 inventory

### K11. Aile eşlemesi düzleştirilmiş sınıf yolunun ön ek eşleşmesiyle yapılır
**TR:** Aile ataması, `protein_class_desc` dizgisinin küçük harfe çevrilmiş
halinde ön ek eşleşmesiyle yapılır — l1/l2 kolonlarıyla değil.
**EN:** Family assignment is a prefix match on the lowercased
`protein_class_desc` path, not on assumed l1/l2 columns.

**Neden / Why:** İlk eşleme varsayıma dayanıyordu ve yanlıştı. ChEMBL 37'de
ölçülen iki özellik eşlemeyi belirliyor:

1. Seviyeler **iki boşlukla** ayrılıyor: `enzyme  kinase  protein kinase`.
2. Yol her sınıfın `short_name`'inden kuruluyor, `pref_name`'inden değil —
   Family A GPCR'ler yolda `membrane receptor  7tm1` olarak geçiyor.

Yalnız l1 adını eşleştirmek yanlış: `enzyme%` bütün enzimleri topluyor, bu
yüzden kinaz ve proteaz ilk denemede birebir aynı sayıları (2.425 hedef,
1.083.520 aktivite) döndürdü. l2 `pref_name`'i eşleştirmek de yanlış: hedefler
ait oldukları en spesifik sınıfa bağlanıyor, bu yüzden l2 satırları neredeyse
boş görünüyor — ChEMBL 37'de "Kinase"a doğrudan 2 insan hedefi, "Protease"a 0
bağlıyken, kinaz ve proteaz **alt ağaçları** 446 ve 216 hedef taşıyor.

Doğrulama: nükleer reseptör sayısı 47 çıktı, literatürdeki ~48 ile uyumlu. Bu
uyum, ön eklere güvenmenin gerekçesi.

**Nerede / Where:** `config/families.yaml: class_path_prefixes`,
`src/dti/chembl.py:build_family_target_table`

### K12. Aile çakışması yok — ölçüldü, varsayılmadı
**TR:** On aile çifti üzerinde ölçüldü: hiçbir hedef birden fazla aileye
girmiyor.
**EN:** Measured across all ten family pairs: zero targets match more than one
family.

**Neden / Why:** Tasarım dokümanı bunu bir risk olarak listelemişti ("bir hedef
birden fazla aileye giriyor"). Risk gerçekleşmedi. Bu, R4'ün anlamlı olmasının
önkoşulu: hem kinaz hem proteaz sayılan bir hedef olsaydı, proteaz ailesi
dışarıda bırakıldığında o hedefin kayıtları eğitimde kalırdı. Çakışma
politikası (`ambiguous_policy: drop`) yine de kodda duruyor — ileride aile
listesi değişirse sessizce bozulmaması için.

**Nerede / Where:** `src/dti/chembl.py:build_family_target_table`

### K13. pChEMBL doğrulaması geçti: 0 uyuşmazlık
**TR:** Süzme ölçütlerine uyan kayıtların tamamında ChEMBL'in `pchembl_value`
alanı, bizim 9 − log10(nM) dönüşümümüzle 0,01 toleransı içinde uyuşuyor.
**EN:** Zero records disagree between ChEMBL's `pchembl_value` and our own
conversion within 0.01.

**Neden / Why:** E2'nin doğrulama adımının sonucu. Uyuşmazlık sıfır olduğu için
ChEMBL'in alanı olduğu gibi kullanılır ve yeniden hesaplama yapılmaz. Bu, bir
varsayımın ölçülerek kapatılmasıdır: sonuç farklı çıksaydı normalizasyon adımı
değişecekti.

**Nerede / Where:** `scripts/01_inspect_chembl.py:report_activity_profile`

### K14. Veri ölçeği beklenenden büyük — beş aile de kapsamda kalıyor
**TR:** Süzme sonrası beş aile toplam 1.078 hedef ve 1.081.125 aktivite
taşıyor. Dokümandaki "on binlerce ölçüm" beklentisi düşük bir tahmindi.
**EN:** After filtering the five families hold 1,078 targets and 1,081,125
activities; the design document's "tens of thousands" was an underestimate.

| Aile | Hedef | Aktivite | Benzersiz molekül |
| --- | --- | --- | --- |
| Kinazlar | 437 | 499.401 | 202.248 |
| GPCR'ler | 238 | 298.071 | 144.910 |
| Proteazlar | 205 | 146.441 | 77.032 |
| İyon kanalları | 156 | 86.774 | 54.028 |
| Nükleer reseptörler | 42 | 50.438 | 27.780 |

Bu tablo, aynı ölçümün ilk turundaki sayıların (1.259 hedef, 1.136.700
aktivite) yerini alıyor. İlk sorgu aktiviteleri şişiriyordu: hedefi
`target_components` üzerinden sınıfa bağlıyordu ve birden fazla bileşeni olan
bir hedefin aktiviteleri her bileşen için bir kez sayılıyordu. Geçici
`fam_target` tablosu hedef başına tek satır tuttuğu için bu çift sayım ortadan
kalktı. Buradaki sayılar yalnız aktivitesi olan hedefleri içerir; sınıf
eşleşmesi olup süzmeden geçen hiç ölçümü olmayan hedefler (örneğin 124 GPCR)
burada görünmez.

**Neden / Why:** Dokümanın 1. aşama kontrol noktası "veri profili beklenenden
küçük çıkarsa kapsamı üç aileye indir" diyordu. Tersi oldu: en küçük aile olan
nükleer reseptörlerde bile 50 binden fazla ölçüm var, yani beş ailenin hepsi
makul bir eşiği geçiyor. Kapsam daraltma kararı gerekmedi.

Not: iyon kanalları hedef sayısı bakımından nükleer reseptörlerden büyük ama
dokümanda "en zor transfer testi" olarak işaretlenmişti. Zorluk veri
miktarından değil, bu ailenin diğerlerinden yapısal farklılığından beklenmişti;
bu beklenti deneyle sınanacak.

**Nerede / Where:** `config/families.yaml` (başlıktaki yorum bloğu),
`reports/tables/01_chembl_inventory.md`

### K15. Birincil protein temsili bağlanma alanı olarak kalıyor — ölçüldü
**TR:** Tasarımın kararı doğrulandı. Ek olarak tam diziye geri düşme
(`fallback_to_full_sequence`) birincil temsilde **kapatıldı**.
**EN:** The design's choice is confirmed by measurement. In addition, the
full-sequence fallback is turned **off** for the primary representation.

**Ölçüm / Measurement** (`scripts/01b_representation_decision.py`, ChEMBL 37).
İki yol aynı birimde karşılaştırıldı — her ailenin aktivitelerinin etkilenen
payı:

| | Kinaz | GPCR | Proteaz | İyon kanalı | Nükleer reseptör | Yayılım |
| --- | --- | --- | --- | --- | --- | --- |
| A: Tam dizi — kırpılan ölçüm payı | %39,0 | %1,3 | %3,4 | %43,9 | %0,0 | **43,9 puan** |
| B: Bağlanma alanı — korunan ölçüm payı | %23,6 | %54,4 | %35,2 | %35,0 | %45,0 | **30,8 puan** |

Eşlenmiş alan uzunlukları: medyan 178–265 aa, en uzun 586 aa, **sınır üstü
alan sayısı sıfır**.

**Neden / Why:** Yol B iki ölçütte de önde — yayılım daha küçük (30,8 < 43,9)
ve kırpma tamamen ortadan kalkıyor. Ama asıl gerekçe sayıların büyüklüğü değil,
**zararların türünün farklı olması**:

- Kırpma **sistematik bilgi kaybıdır**: protein kısmen görülür ve bu, aileye
  bağlı oranlarda olur. Kinaz (birincil eğitim ailesi) ölçümlerinin %39'u,
  iyon kanalı (en zor transfer hedefi) ölçümlerinin %43,9'u etkilenirken GPCR
  %1,3 ve nükleer reseptör %0. Yani aileler arası performans farkı transferden
  değil kırpmadan gelebilir — dokümanın "ana soruyu tehdit eden karıştırıcı"
  dediği şey tam olarak budur ve ölçüm bu korkunun haklı olduğunu gösteriyor.
- Kapsam sınırı **örneklem kaybıdır**: daha az ölçüm kalır, kalanlar eksiksiz.

Örneklem kaybı istatistiksel gücü azaltır; sistematik kayıp yanlılık üretir.
Karşılaştırmanın kendisi aileler arası olduğu için yanlılık öldürücüdür, güç
kaybı değil. Üstelik Yol B bol veri bırakıyor: beş ailede toplam ~415 bin
aktivite (kinaz 125.601, GPCR 162.861, proteaz 54.778, iyon kanalı 49.152,
nükleer reseptör 22.689).

**Geri düşme neden kapatıldı:** eşlemesi olmayan hedefler için kırpılmış tam
diziye düşmek, iki temsili aileye bağlı oranlarda karıştırır — yani iki
karıştırıcıdan birini seçmek yerine ikisini birleştirir. Birincil temsil saf
alan; geri düşme ayrı bir ablasyon olarak ölçülür.

**Nerede / Where:** `config/config.yaml: protein`,
`reports/tables/02_representation_decision.md`

### K16. Tautomer adımı korunuyor, paralelleştirilerek
**TR:** Tautomer kanonikleştirme 37,5 ms/molekül (kapalıyken 5,0 ms) ve
moleküllerin %19'unun kanonik SMILES'ini değiştiriyor. Adım korunur,
çok çekirdeğe dağıtılır.
**EN:** Tautomer canonicalisation costs 37.5 ms/molecule against 5.0 without
it, and changes the canonical SMILES of 19% of molecules. Kept, parallelised.

**Neden / Why:** 200 bin molekül için ~2 saat, tautomer olmadan ~17 dakika.
Ama %19 atlanamaz: o moleküller aynı bileşiğin iki kaydı olarak kalır ve R4b
bölmesinin iki tarafına düşebilir — "kimyasal olarak yeni ligand" tanımı
sızdırır. Standardizasyon molekül başına bağımsız olduğu için paralelleştirme
sonucu değiştirmez; maliyet bir kez ödenir ve çıktı diske yazılır.

**Nerede / Where:** `config/config.yaml: molecule.n_workers`

### K17. Deney güven skoru: 8 eşiği etkisiz, 9 bir ablasyon olmalı
**TR:** `confidence_score >= 8` filtresi hiçbir kayıt düşürmüyor. Eşiğin 9'a
çıkarılması yeni bir ablasyon olarak eklendi.
**EN:** The `confidence_score >= 8` filter removes nothing; raising it to 9 is
added as an ablation.

**Ölçüm / Measurement:** Aile sınıfı eşleşen insan tek-protein hedeflerine ait
assay'lerin tamamı 8 veya 9 skoruna sahip (ChEMBL 37: 1.323.461 adet 8,
2.654.404 adet 9). Veritabanının genelinde 0–7 skorlu assay'ler var ama bunlar
tek proteinli insan hedeflerine bağlı değil. Yani `target_type = SINGLE
PROTEIN` kısıtı güven skoru eşiğini zaten içeriyor.

**Neden önemli / Why it matters:** Filtre etkisiz olduğu için kaldırılabilirdi,
ama altındaki ayrım ölçüm kalitesiyle doğrudan ilgili. ChEMBL tanımları:

- **9** — doğrudan tek protein ataması
- **8** — **homolog** tek protein ataması, yani ölçüm hedefin başka bir türdeki
  homoloğu üzerinde yapılmış olabilir

Hedefin organizması insan olarak süzülüyor, ama deneyin homolog üzerinde
yapılmış olması bu süzmeyle engellenmiyor. Bu durumda model, gömdüğü insan
dizisiyle ölçümün yapıldığı proteinin aynı olmadığı kayıtlarla eğitilir —
aileler arası transfer sorusunda gürültü kaynağı. Eşiği 9'a çıkarmak kayıtların
yaklaşık üçte birine mal oluyor; bu bedelin değip değmediği ablasyonla ölçülür.

**Nerede / Where:** `config/config.yaml: filtering.min_assay_confidence_score`

### K18. Çekme tamamlandı — düşürme tablosu
**TR/EN:** `scripts/02_extract.py`, ChEMBL 37.

| # | Adım | Kalan | Düşen |
| --- | --- | --- | --- |
| 1 | Aile eşleşen tüm kayıtlar | 3.977.865 | — |
| 2 | Güven skoru ≥ 8 | 3.977.865 | 0 |
| 3 | Ölçüm türü Ki/Kd/IC50 | 1.643.368 | 2.334.497 |
| 4 | Birim nM | 1.502.340 | 141.028 |
| 5 | Kesin ilişki (=) | 1.095.200 | 407.140 |
| 6 | pChEMBL mevcut | 1.081.431 | 13.769 |
| 7 | Doğrulama bayrağı yok | 1.081.125 | 306 |
| 8 | pChEMBL makul aralıkta | 1.081.069 | 56 |
| 9 | SMILES var | 1.080.040 | 1.029 |
| 10 | Dizi var | 1.080.040 | 0 |

Başlangıcın **%27,2'si** kalıyor. İki adım kaybın neredeyse tamamını
açıklıyor: ölçüm türü süzmesi (2,33 milyon — EC50, Potency, Inhibition gibi
türler) ve sansürlü kayıtların atılması (407 bin). İkisi de tasarımda
gerekçelendirilmiş kararlar; sansürlü bir değer regresyon etiketi değildir.

Çıktılar: 1.080.040 aktivite, 1.259 hedef, **497.867 benzersiz bileşik**,
414.858 alan eşleme satırı.

**Pratik sonuç:** Bileşik sayısı 200 bin tahminimin iki buçuk katı. Tautomer
kanonikleştirmesi tek çekirdekte 497.867 × 37,5 ms ≈ **5,2 saat** demek, yani
K16'daki paralelleştirme bir iyileştirme değil zorunluluk.

**Nerede / Where:** `reports/tables/03_drop_table.md`,
`data/interim/*.parquet`

---

---

## K19 — Molekül standardizasyonu tamamlandı / Standardisation complete

**Karar / Decision:** Kanonik SMILES, tekilleştirmenin anahtarı olarak
doğrulandı.

**Ölçüm / Measurement:** 497.867 bileşik, 15 işlemcide **67,4 dakika**
(ortalama 123 molekül/s). 493.803 kayıt korundu; 4.000 kayıt `too_large`
(>100 ağır atom), 57 kayıt standardizasyon hatası, 7 kayıt `too_small`.

Korunan 493.803 kayıt **487.220 benzersiz moleküle** çöktü: 6.583 kayıt,
tuz/yük/tautomer varyantı olduğu için başka bir kaydın aynısı çıktı. Bu sayı
K16'nın gerekçesinin doğrulanmasıdır — bu 6.583 kayıt `molregno` üzerinden
tekilleştirilseydi R4b bölmesinde aynı molekül iki kimlikle iki yakaya
düşebilirdi.

**Nerede / Where:** `data/interim/compounds_standardised.parquet`

---

## K20 — Tekrarlı alan kopyaları: belirsizliğin kaçırdığım ikinci türü

**Sorun / Problem:** `02_extract.py` içindeki `n_domains_for_activity`
sütununu belirsizlik ölçütü olarak kullanmıştım. Bu sütun **farklı** Pfam
alanlarını sayıyor (`nunique(domain_id)`), dolayısıyla bir bağlanma bölgesinin
iki ayrı alan ailesine eşlenmesini yakalıyor. Ama aynı proteinde **tek alanın
birden fazla kopyasına** eşlenmesini yakalamıyor: gerilim kapılı bir kanalda
dört `Ion_trans` tekrarı, bir JAK'ta iki kinaz alanı vardır ve hepsinin
`domain_id`'si aynıdır.

İlk kurulumda bu kayıtlar "belirsiz değil" sayılıp dört kopyadan biri
`drop_duplicates` ile keyfî biçimde seçilmişti — yani model bir iyon kanalının
gözeneğinin dörtte birini görüyordu.

**Ölçüm / Measurement:** ChEMBL 37'de bu iki türün dağılımı:

| Tür / Form | Aktivite |
| --- | --- |
| Birden fazla alan **ailesine** eşlenen | **0** |
| Tek alanın birden fazla **kopyasına** eşlenen | **18.553** |

Yani fiilen var olan tek belirsizlik türü, kontrolümün kaçırdığı türdü.
Etkisi aile bağımlı:

| Aile | Etkilenen aktivite | Oran |
| --- | --- | --- |
| İyon kanalı | 6.329 / 30.329 | **%20,9** |
| Kinaz | 8.184 / 117.349 | %7,0 |
| Proteaz | 3.239 / 51.512 | %6,3 |
| GPCR | 801 / 161.922 | %0,5 |
| Nükleer reseptör | 0 / 22.677 | %0,0 |

Yayılım **20,9 puan**.

**ChEMBL bunu çözemiyor.** `component_domains` her kopya için kendi
sınırlarıyla bir satır tutuyor; `site_components` ise yalnızca `domain_id`'ye
atıfta bulunuyor, konum bilgisi yok. Kaynak, bölgenin hangi kopyada olduğunu
söylemiyor — belirsizlik indirgenemez.

**Karar / Decision:** Bu aktiviteler düşürülüyor (`drop_ambiguous`). Gerekçe
K15'in argümanının aynısı: keyfî kopya seçimi *sistematik, aile bağımlı bilgi
kaybı* (sapma) üretir; düşürmek *örneklem kaybı* (güç) demektir ve
karşılaştırmanın kendisi aileler arası olduğunda sapma öldürücüdür.

Biyolojik olarak doğru alternatif — kopyaları dizi sırasında birleştirmek, ki
bir kanalı dört tekrarının tamamı temsil eder — **ablasyon 6** olarak
kaydedildi, birincil yol değil: birleşim, etkilenen aktivitelerin %34'ünde
1.050–1.161 aa'ya çıkıyor, yani `esm_max_residues` sınırını yine aşıyor ve
kesmeyi tam da iyon kanalları için geri getiriyor.

**Bedeli dürüstçe / The cost, stated plainly:** Düzeltme, zaten en ince olan
aileyi daraltıyor. İyon kanalı **50 → 35 hedef**, 24.580 → 19.545 çift.
35 hedefle R3 (görülmemiş hedef) iyon kanalları içinde ince kalacak; bu, bölme
aşamasında aile başına rapor edilmeli.

**Nerede / Where:** `scripts/04_build_dataset.py` (`attach_domains`),
`config/config.yaml` (`multi_domain_policy`)

---

## K21 — Veri seti kuruldu / Dataset assembled

**Sonuç / Result:** **283.353 hedef–bileşik çifti**, 565 hedef, 190.339
benzersiz molekül.

| Aile | Çift | Hedef | Bileşik | pChEMBL ortanca | Çift/hedef (ortanca) |
| --- | --- | --- | --- | --- | --- |
| GPCR | 122.142 | 175 | 79.880 | 7,04 | 318 |
| Kinaz | 85.163 | 235 | 55.287 | 6,99 | 130 |
| Proteaz | 39.639 | 89 | 28.304 | 6,80 | 182 |
| İyon kanalı | 19.545 | 35 | 18.756 | 6,20 | 69 |
| Nükleer reseptör | 16.864 | 31 | 12.391 | 6,75 | 238 |

**E4'ün gerekçesi, ölçülmüş hâli:** Aile ortalamaları 6,36 (iyon kanalı) ile
7,07 (GPCR) arasında, **yayılım 0,72 log birim**. Üç RMSE varyantı bu yüzden
gerekli: dışarıda bırakılan aileye geçişte sistematik kayma bu büyüklükte
olacak.

**Ama bu sayı tek başına yanıltıcı ve tezde öyle sunulmamalı.** Aile *içi*
standart sapma 1,21–1,41 log birim, yani aileler arası kaymanın yaklaşık iki
katı. Pratik sonucu şu: ortalama-kaydırma düzeltmesinin kurtarabileceği hata
en fazla 0,72 log birimlik bir bileşendir. R4'te ham RMSE 2,0 çıkarsa, kaymayı
düzeltmek bunu 1,9'a indirir — yani "transfer başarısız" sonucu kaymayla
açıklanamaz. Üç varyantın işi transferi *kurtarmak* değil, başarısızlığın
ölçek kaynaklı mı sıralama kaynaklı mı olduğunu **ayırt etmek**.

**Ölçüm kalitesi / Measurement quality:** 49.424 çiftin birden fazla ölçümü
var; 13.935'i tutarsız olarak işaretlendi (>1,0 log birim yayılım — ablasyon 4
bunları düşürür), 11.691'i farklı ölçüm türlerini karıştırıyor (ablasyon 3).
Alan uzunluğu ortancası 260 aa.

**Nerede / Where:** `data/processed/dataset.parquet`,
`reports/tables/04_dataset_profile.md`,
`reports/tables/05_drop_table_stage2.md`


---

## K22 — MMseqs2 yerine tam hizalama / Exact alignment instead of MMseqs2

**Karar / Decision:** R3 bölmesindeki dizi kimliği, MMseqs2 ile değil **tam
Smith-Waterman** hizalamasıyla hesaplanıyor.

**Gerekçe / Reasoning:** Tasarım belgesi MMseqs2'yi şart koşmuştu. MMseqs2
milyonlarca dizi için doğru araç; burada **565** hedef var ve bu cevabı
değiştiriyor. 565 hedef 159.330 çift demek; ölçtüm, çift başına 2,41 ms, yani
tek çekirdekte 385 saniye, makinenin tamamında bir dakikanın altı.

Bu ölçekte MMseqs2'nin sunduğu şey — hizalamaların çoğunu atlayan k-mer ön
süzgeci — hiçbir şey kazandırmıyor. Karşılığında getirdiği şey ise gerçek: ön
süzgeç, duyarlılık ayarı olan bir sezgisel ve kaçırabileceği çiftler tam da
%40 kimlik eşiğinin hemen üstündeki çiftler. R3'te kaçırılan bir çift, eğitim
hedefine benzer kalmış bir test hedefi demek — yani bölmenin var olma sebebi
olan kaçağın kendisi.

Yan kazanç: boru hattından Windows'a özgü bir ikili kurulumu kalkıyor.

**İki tanım açıkça yazıldı** (evrensel olmadıkları için):
- `identity` = aynı hizalanmış kolon / hizalama uzunluğu (boşluklar dahil)
- `coverage` = hizalanmış kalıntı çifti / **kısa** dizinin uzunluğu

Bir çiftin bağlı sayılması için **ikisinin de** eşiği geçmesi gerekiyor;
yoksa proteinin onda birini kaplayan %95 kimlikli kısa bir motif, iki alakasız
aileyi birbirine zincirler. Bu testle birlikte kayıtlı.

**Kimlik, tam dizi üzerinde değil ALAN üzerinde hesaplanıyor** — çünkü modele
verilen şey alan. İki protein tüm uzunlukları boyunca farklı olup neredeyse
aynı katalitik alanı paylaşabilir; tam dizi üzerinden kümelemek bunları ayrı
hedef sayar ve sonra modele bölmenin iki yakasında aynı girdiyi verir.

**Nerede / Where:** `src/dti/sequences.py`, `config/config.yaml`
(`splits.r3_target`), `requirements.txt` (biopython==1.84)

---

## K23 — Bağlantı türü asimetrisi: hedeflerde tek bağlantı, bileşiklerde Butina

**Karar / Decision:** Hedefler **tek bağlantı** (bağlı bileşenler), bileşikler
**Butina küre dışlama** ile gruplanıyor. Bu asimetri bilinçli ve tezde
belirtilmesi gereken bir sonuç, bir uygulama ayrıntısı değil.

**Neden hedeflerde tek bağlantı:** R3, hiçbir test hedefinin **hiçbir** eğitim
hedefine eşik üstü benzememesini istiyor. Bunu *garanti eden* tek gruplama,
"eşik üstü" grafiğinin bağlı bileşenleridir: iki hedef bağlıysa aynı gruptadır,
dolayısıyla aynı yakaya düşer. Zincirleme (a–b–c'de a ile c benzemez) tek
bağlantıya yapılan olağan itirazdır, ama 565 hedefle gruplar incelenebilir
kalıyor ve burada garanti, derli toplu kümelerden daha değerli.

**Neden bileşiklerde Butina:** Aynı garanti 190.339 ilaç benzeri bileşik için
mevcut değil. 0,4 Tanimoto eşiğinde komşuluk grafiği, bağlı bileşenlerin
kimyasal uzayın büyük kısmını kapsayan tek bir dev bileşene çökeceği kadar
yoğun — o durumda hiçbir bölme mevcut olmaz. Küre dışlama kümeleri derli toplu
tutuyor ama **farklı** kümelerdeki bileşikler arası benzerliğe sınır koymuyor.

**Sonuç, dürüstçe:** R2 ve R4a yalnızca "test bileşiği, hiçbir eğitim
bileşiğinin gelmediği bir kümeden geliyor" garantisini veriyor — R3'ün
garantisinden zayıf. Bunun yeterince sıkı olup olmadığı ampirik bir soru; bu
yüzden kaçak denetimi her rejim için eğitim–test azami Tanimoto'sunu ölçüyor ve
R4b açık bir bölme sonrası süzgeç olarak var. Zayıf garantiyi dürüst kılan şey,
denetim sayısını rapor etmek.

**Nerede / Where:** `src/dti/clustering.py` (modül başlığı bu argümanı
taşıyor), `tests/test_clustering.py`,
`config/config.yaml` (`splits.r2_compound.clustering`)

---

## K24 — 1,8e10 çiftlik Tanimoto geçişi: matris çarpımı olarak GPU'da

**Sorun / Problem:** Bölme kodu 190.339 bileşik üzerinde iki tüm-çiftler
sorusunu cevaplamak zorunda: her bileşiğin eşik üstü kaç komşusu var (Butina'nın
sıralaması için) ve her test bileşiğinin tüm eğitim kümesine azami benzerliği ne
(R4b süzgeci ve kaçak denetimi için). Bu 190.339² / 2 = **1,8e10 çift**. Yoğun
matris imkânsız (3,6e10 float32 = 145 GB) ve RDKit'in toplu benzerliği üzerinden
Python döngüsü günler sürer.

**Karar / Decision:** İki bit kümesinin kesişimi bir iç çarpımdır, dolayısıyla
tüm geçiş tek bir büyük fp16 matris çarpımı — 7,4e13 çarp-topla, RTX 4050'nin
dakikalar içinde yaptığı bir iş. `molecules.tanimoto_matrix` yoğun ve apaçık
doğru referans uygulaması olarak kalıyor; `src/dti/similarity.py` hızlı olan.

**fp16 aritmetiğini güvenli kılan iki ayrıntı:**

1. ECFP4 parmak izlerinde 2048 bitten 30–80'i kurulu ve fp16, 2048'in altındaki
   tam sayıları **tam** temsil ediyor. Yani kesişim ve popcount değerleri tam;
   yalnızca son bölme yaklaşık.
2. Eşik testleri bölmeyi tamamen atlıyor. `sim >= c` şuna dönüştürüldü:
   `inter * (1 + c) >= c * (pq + pr)` — yalnızca tam sayıları karşılaştırıyor.
   Tam eşiğin üstünde oturan bir benzerlik bu yüzden tam aritmetikteki gibi
   çözülüyor. Bu önemli: eşik bir ablasyon parametresi ve sınırdaki bir çift
   koşular arasında yakasını değiştirmemeli.

**Doğrulama / Verification:** `tests/test_similarity.py` tiled GPU yolunu yoğun
referansa karşı sınıyor — yedi farklı eşikte komşu sayıları birebir eşit, azami
benzerlikler 1e-3 toleransında, karo sınırları (QUERY_TILE=1, REF_CHUNK=2'ye
indirilerek) cevabı değiştirmiyor, tam eşik üstündeki çift doğru çözülüyor,
sıfır parmak izi 0/0 birleşimi üzerinden her şeye benzemiyor, ve CPU ile GPU
yolları aynı sonucu veriyor. 17 test.

**Nerede / Where:** `src/dti/similarity.py`, `tests/test_similarity.py`


---

## K25 — Kaçak denetimi metriği bozuktu / The audit metric was broken

**Sorun / Problem:** Kaçak denetiminin ilk koşusu **her rejim için** azami
eğitim–test dizi kimliğini 1,000 raporladı — hiçbir hedefin iki yakada
bulunmadığı R3 dahil (`shared_targets = 0` ile doğrulandı). Yani sorun
bölmede değil, metrikte.

**Sebep / Cause:** Kimliği yerel hizalama üzerinden tanımlamıştım
(`aynı kolon / hizalama uzunluğu`). Smith-Waterman en yüksek puanlı parçayı
bulur; iki alakasız protein kısa bir motif paylaşıyorsa, dizinin %2–10'unu
kaplayan bir pencerede %70–100 "kimlik" döner. Ölçtüm: farklı hedef
gruplarında olup kimliği ≥0,40 çıkan 39.755 çift vardı ve kapsamaları
0,024–0,096 aralığındaydı. Bu bir motif eşleşmesi, protein benzerliği değil.

Kümeleme bunu doğru ele alıyordu (kapsama ≥0,80 şartı vardı); hatalı olan
denetimdi, çünkü kimliği koşulsuz raporluyordu.

**Karar / Decision:** Üçüncü bir ölçü eklendi ve birincil hâline geldi:

    identity_over_shorter = aynı kalıntı sayısı / KISA alanın uzunluğu

Bu ölçünün motif açığı yapısal olarak yok: 12 kalıntılık bir motif, 250
kalıntılık bir alanın %40'ını oluşturamaz. `tests/test_sequences.py` tam bu
iddiayı sınıyor — 120'mer iki alakasız dizi, paylaştıkları 12 kalıntılık motif
sayesinde yerel hizalama kimliğinde >0,80 alırken kısa-dizi kimliğinde <0,20
alıyor.

**İkinci hata, ilkini düzeltirken ortaya çıktı.** Kümelemeyi eski tanımla
(kimlik + kapsama), denetimi yeni tanımla bıraktığımda R3 test hedeflerinin
**%10,2'si** eşik üstü göründü. Bu bir kaçak değil, iki tanımın "aynı protein"
konusunda anlaşmazlığıydı. Tek tanıma indirdim: kümeleme de denetim de
`identity_over_shorter` kullanıyor. Sonuç, R3'te azami kimlik **0,398** ve
eşik üstü hedef oranı **%0** — garanti artık inşa gereği sağlanıyor ve denetim
bunu doğruluyor.

Eski iki eşikli tanım ablasyon olarak kaldı: 206 grup yerine 183 grup veriyor.

**Nerede / Where:** `src/dti/sequences.py`, `src/dti/leakage.py`,
`src/dti/clustering.py`, `tests/test_sequences.py`,
`config/config.yaml` (`splits.r3_target.identity_definition`)

---

## K26 — Bölmeler kuruldu ve R2'nin zayıflığı ölçüldü

**Gruplama / Grouping:**

| | Grup | En büyük | Ortanca | Tekil |
| --- | --- | --- | --- | --- |
| Bileşik kümeleri (Butina, 0,4) | 13.056 | 1.378 | 3 | 4.354 (%33,3) |
| Hedef grupları (tek bağlantı, 0,4) | 183 | 48 | 1 | 98 (%53,6) |

Kümeleme korkumu boşa çıkardı: 0,4 eşiğinde kimyasal uzay tek dev bileşene
çökmedi, 13.056 küme çıktı ve en büyüğü bileşiklerin yalnızca %0,7'si.
(Koşunun başında %14'te 81 küme görünmesi yanıltıcıydı — küre dışlama en yoğun
bileşiklerden başladığı için ilk kümeler doğal olarak en büyükleri.) 1,8e10
çiftlik geçiş GPU'da **0,5 dakika** sürdü.

**Asıl bulgu — R2 iddia ettiği kadar zor değil:**

| Rejim | Eğitim–test Tanimoto ortanca | p95 | Eşik üstü test bileşiği |
| --- | --- | --- | --- |
| R1 (rastgele) | 0,886 | 1,000 | %99,4 |
| **R2 (görülmemiş bileşik)** | **0,548** | **0,823** | **%81,7** |
| R3 (görülmemiş hedef) | 0,443 | 1,000 | %64,0 |
| R4a (aile dışı) | 0,377–0,417 | 0,670–1,000 | %39,1–56,6 |
| R4b (aile dışı + kimyasal yeni) | 0,341–0,345 | 0,392–0,393 | **%0,0** |

R2'de test bileşiklerinin **%81,7'si** eğitim kümesinde ≥0,40 Tanimoto'da bir
komşuya sahip, ortanca benzerlik 0,548. Beş tohumda oran 0,805–0,828, yani
bu bir tohum tesadüfü değil. K23'te teorik olarak öngördüğüm zayıflık —
Butina'nın kümeler *arası* benzerliğe sınır koymaması — nicel hâliyle burada.

**Bunun tez için anlamı:** Literatürde "zor bölme" olarak sunulan kümeleme
tabanlı ayırma, bu ölçüyle neredeyse hiç zor değil. R2 skoru "görülmemiş
kimyada başarı" olarak okunamaz; okunabilecek şey "görülmemiş *kümede*
başarı"dır ve bu çok daha zayıf bir iddia. Denetim tablosu olmadan bu ayrım
görünmez olurdu. R4b'nin %0'ı ise bölme sonrası açık süzgecin işe yaradığını
gösteriyor — azami benzerlik 0,398–0,399, yani eşiğin hemen altında, tam
beklendiği gibi.

**R4b ne kadarını atıyor:** Aile dışı test verisinin yaklaşık yarısı kimyasal
olarak tanıdık çıkıyor.

| Aile | R4a test | R4b test | Atılan |
| --- | --- | --- | --- |
| Kinaz | 85.163 | 50.691 | %40 |
| GPCR | 122.142 | 61.746 | %49 |
| Proteaz | 39.639 | 21.912 | %45 |
| Nükleer reseptör | 16.864 | 9.685 | %43 |
| İyon kanalı | 19.545 | 8.645 | %56 |

Beş ailenin hepsi `r4b_min_test_records: 200` eşiğinin çok üstünde kalıyor
(en küçüğü 8.645), yani hiçbir aile R4b'de puanlanamaz duruma düşmedi.

**Hedef grubu başına eşik aşımı:** Gruplar bölündüğü için oranlar yaklaşık.
R3'ün test oranı tohuma göre %15,4–24,3 arasında oynuyor (en büyük hedef grubu
48 hedef taşıyor, bu da tek bir grubun oranı aşırtmasına yetiyor). R4'te test
oranı zaten ailenin boyutuna eşit: GPCR'da %43,1, nükleer reseptörde %6,0.
Bunlar seçilmiş değil, verinin kendisi.

**Aynı molekül sayısı — hangi rejimde sıfır olmalı:** R2 ve R4b'de sıfır,
ölçüldü ve sıfır. Diğerlerinde sıfır **olmaması doğru**: R1 iyimser taban
çizgisi; R3 hedefleri dışarıda bıraktığı için aynı bileşik farklı proteinlere
karşı iki yakada görünür (4.937 molekül); R4a bir aileyi dışarıda bıraktığı
için hem kinaza hem GPCR'a ölçülmüş promiskü bileşikler iki yakada yer alır
(186–3.547 molekül, aileye göre). Denetim tablosunun ilk sürümü bunu "hata"
diye etiketliyordu; düzeltildi.

**Nerede / Where:** `data/splits/` (65 dosya: 13 rejim × 5 tohum),
`reports/tables/06_splits.md`, `reports/tables/07_leakage_audit.md`,
`data/interim/leakage_audit.parquet`


---

## K56 — Güven skoru 8 vs 9: soru çözüldü, eşik fiilen zaten 9

**Soru / Question:** Devir notu §2.3 ve v2 tasarımının açık soruları, mevcut
`min_assay_confidence_score: 8` eşiğinin homolog atamaları içerip içermediğini
yerel ChEMBL 37 kopyasında doğrulamayı istiyordu. Rev2'deki gerekçe "deneyin
tek bir proteine atandığı" diye yazılmıştı; ikincil kaynaklar 9'un doğrudan,
8'in homolog atama olduğunu söylüyordu.

**ChEMBL'in kendi tanımları** (`confidence_score_lookup` tablosu, ikincil
kaynak değil):

| Skor | target_mapping | Açıklama |
| --- | --- | --- |
| 7 | Protein complex | Direct protein complex subunits assigned |
| **8** | **Homologous protein** | **Homologous single protein target assigned** |
| **9** | **Protein** | **Direct single protein target assigned** |

İkincil kaynak doğruydu ve rev2'nin gerekçesi eksikti: 8 gerçekten tek bir
proteine atama demek, ama atanan protein insan hedefinin **başka bir türdeki
homoloğu** olabilir. Yani modele gömdüğümüz insan dizisiyle ölçümün yapıldığı
protein aynı olmayabilir.

**Süzülmüş aktivite düzeyinde dağılım:**

| Skor | Aktivite | Hedef | Bileşik |
| --- | --- | --- | --- |
| 8 (homolog) | 340.341 (%31,5) | 843 | 167.720 |
| 9 (doğrudan) | 739.699 (%68,5) | 1.050 | 375.976 |

**Beklenmedik bulgu — eşik zaten 9.** Nihai veri setine etkisini ölçtüğümde
eşiği 9'a çekmenin **sıfır** çift kaybettirdiği çıktı: 283.353 çiftin %100'ü,
565 hedefin tamamı korunuyor. Sebebini aradım ve mekanizma şu:

| Skor | Aktivite | Alan eşlemesi olan | Oran |
| --- | --- | --- | --- |
| 8 | 340.341 | **0** | **%0,0** |
| 9 | 739.699 | 383.789 | %51,9 |

**ChEMBL `predicted_binding_domains`'i yalnızca doğrudan (skor 9) atanmış
deneyler için hesaplıyor.** Beş ailenin hepsinde skor-8 aktivitelerinin alan
eşlemesi oranı tam olarak %0. Dolayısıyla K15'te alınan bağlanma alanı temsili
kararı, farkında olmadan güven eşiğini de 9'a çekmiş: homolog kayıtların
tamamı güven süzgeci tarafından değil, **alan süzgeci** tarafından elenmiş.

Bu, K17'de kaydettiğim "güven süzgeci hiçbir şey düşürmüyor" bulgusunun daha
derin sebebi. K17 doğruydu ama eksikti: süzgeç `target_type = SINGLE PROTEIN`
kısıtı yüzünden değil, alan eşlemesi yüzünden etkisiz.

**KARAR (15 Kasım 2026): eşik 9'a çekildi.** Aşağıdaki öneri kabul edildi ve
`config/config.yaml` güncellendi; veri hattı yeniden koşturuldu.

### Yeniden koşturmanın sonucu

Nihai veri seti **birebir aynı** kaldı — bu sadece satır sayısı eşitliği değil,
`DataFrame.equals()` düzeyinde aynılık ve aynı satır sırası (sıra hash'i
`600a5799642cbfa0532d6f33`, öncesi ve sonrası). Bu kontrol zorunluydu: 65 bölme
dosyası veri setine **satır indeksiyle** atıf yapıyor, dolayısıyla satır sırası
değişseydi bütün bölmeler sessizce başka çiftlere işaret ediyor olacaktı.

Sonuç olarak **yeniden üretilmesi gereken hiçbir aşağı akış çıktısı yok**:
65 bölme dosyası, `06_splits.md` ve `07_leakage_audit.md` geçerli;
`04_dataset_profile.md` bit düzeyinde değişmedi.

Değişen iki tez tablosu:

| Tablo | Ne değişti |
| --- | --- |
| `03_drop_table.md` | Güven süzgeci artık **1.323.461** kayıt düşürüyor (önce 0). Son kalan 1.080.040 → **739.699**; başlangıcın %27,2'si değil **%18,6'sı** |
| `05_drop_table_stage2.md` | İlk iki satır girdi sayısına göre kaydı; "alan eşlemesi olan kayıtlar" satırı yine **360.616** |

O son sayı mekanizmanın bağımsız üçüncü doğrulaması: girdi 1.080.040'tan
739.699'a inerken alan eşlemesi olan kayıt sayısı değişmedi, çünkü düşen
340.341 kaydın hiçbirinde alan eşlemesi yoktu.

**Yan kazanç:** düşürme tablosu artık daha bilgilendirici. Önceki hâlinde güven
satırı 0 düşürüyor görünüyordu ve okuyucuya hiçbir şey söylemiyordu; şimdi
homolog kayıtların açıkça elendiği görülüyor.

**Ara dosya hizalaması — yapıldı.** Yeniden koşturmanın ardından
`compounds_standardised.parquet` eski çalıştırmadan 497.867 satır taşıyordu,
`compounds.parquet` ise 375.976'ya inmişti. Birleştirme aktiviteler üzerinden
yürüdüğü için fazlalık zararsızdı, ama dosya girdisiyle uyuşmuyordu;
`03_standardise_molecules.py --force` ile yeniden üretildi (48,7 dakika,
15 işlemci):

| | Eşik 8 (önce) | Eşik 9 (sonra) |
| --- | --- | --- |
| Girdi bileşiği | 497.867 | 375.976 |
| Korunan | 493.803 | 372.406 |
| Benzersiz molekül | 487.220 | 368.307 |
| `too_large` | 4.000 | 3.522 |

Veri seti bu üçüncü koşudan sonra da **birebir aynı** (aynı sıra hash'i,
`equals()` True) ve yeniden üretilen `04_dataset_profile.md` ile
`05_drop_table_stage2.md` **bayt düzeyinde** değişmedi — yani bu adım izlenen
hiçbir dosyayı değiştirmedi, sadece ara dosyaları girdileriyle hizaladı.
Boru hattı artık baştan sona tek bir eşikle tutarlı.

### Özgün öneri ve gerekçesi

`min_assay_confidence_score` **9'a çekilsin.**
Gerekçe, bugün hiçbir şeyi değiştirmemesi değil — asıl gerekçe ileriye dönük:

1. Mevcut davranış bir *kaza*; config 8 diyor, veri 9. Config'i 9 yapmak,
   davranışı kasıtlı hale getirir ve yanlış okunmasını engeller.
2. `fallback_to_full_sequence` ablasyonu açılırsa ya da v2'nin **Aşama 2**
   kapsamına (tüm insan tek-protein hedefleri) geçilirse, alan süzgeci
   gevşeyecek ve eşik 8 kaldığı sürece homolog kayıtlar **sessizce** veri
   setine girecek. Ürün bir hedefe bağlanma tahmini dağıttığı için bu, kalite
   beyanını fark edilmeden bozar.
3. Maliyeti ölçüldü: bugün sıfır. Veri hattı yeniden koşturulmasını
   gerektirmiyor.

Eşik 8'de bırakılırsa gereken şey: homolog/doğrudan farkının tek tablo olarak
raporlanması ve her kaydın `confidence_score` alanının veri setinde taşınması.

**Nerede / Where:** `config/config.yaml` (`filtering.min_assay_confidence_score`),
ölçüm script'leri oturum çalışma alanında; sayılar bu kayıtta.


## Açık sorular / Open questions

- **Tez formatı.** Bölümün zorunlu şablonu henüz verilmedi. Geldiğinde K5
  uyarınca yalnız metin taşınacak.
- ~~**Bağlanma alanı temsili.**~~ **Kapandı (K15).** Ölçüldü, temsil
  korunuyor, geri düşme kapatıldı.
- ~~**Aktivite başına çoklu alan eşlemesi.**~~ **Kapandı (K20).** Ölçüldü:
  birden fazla alan *ailesine* eşlenen aktivite yok (0); gerçek belirsizlik
  tek alanın tekrarlı *kopyalarında* ve 18.553 aktiviteyi etkiliyor. Sayı
  düşürme tablosuna girdi.
- ~~**Temsil birimi: hedef mi aktivite mi.**~~ **Kapandı (ölçüldü).** 804
  hedefin **803'ü** bütün aktivitelerinde tek bir alana eşleniyor; tek
  istisna, aktivitelerin %0,6'sı. Yani aktivite başına eşleme ile hedef başına
  eşleme fiilen çakışıyor: PPDMs aktivite düzeyinde ek bilgi vermiyor ve
  R3/R4'ün `tid` üzerinden gruplaması temsil birimiyle tutarlı. Bölme kodunda
  ayrı bir doğrulama gerekmiyor.
- **İyon kanallarında yardımcı altbirimler.** `ion channel` ön eki, gözenek
  oluşturmayan yardımcı altbirim ailelerini de kapsıyor. K20'den sonra aile
  35 hedefe indi; hangi hedeflerin yardımcı altbirim olduğu ve çıkarılmaları
  gerekip gerekmediği bu daralmış liste üzerinde karara bağlanacak.
- **R2'nin yeniden tanımlanması.** K26 R2'de test bileşiklerinin %81,7'sinin
  eğitimde ≥0,40 komşusu olduğunu gösterdi. Seçenekler: (a) olduğu gibi bırakıp
  denetim sayısıyla birlikte raporlamak, (b) R4b'deki gibi bölme sonrası süzgeç
  eklemek, (c) kümeleme eşiğini yükseltmek. Tez açısından en bilgilendirici
  olanı muhtemelen (a) + (b)'yi ayrı bir rejim olarak sunmak: ikisinin skor
  farkı, "görülmemiş küme" ile "görülmemiş kimya" arasındaki mesafeyi doğrudan
  ölçer. Karar kullanıcıya ait.
- **İyon kanalı ailesinin inceliği.** 35 hedef / 19.545 çift ile bu aile R3
  (görülmemiş hedef) için sınırda. R4'te dışarıda bırakılan aile olarak hâlâ
  kullanılabilir, ama R3'ün aile içi kırılımında güven aralıkları geniş
  olacak; rapor edilmesi gerekir.
