# Proje Tasarımı v2 — İlaç Yeniden Konumlandırma Hipotezi Üreteci

Nov 15, 2026 · @secret secret

## Karar özeti

Aileler arası genellenebilirlik sorusu sonraya bırakıldı; bitirme çalışması aynı veri üzerine kurulan bir **ilaç yeniden konumlandırma hipotezi üreteci** oldu. Gerekçe: aile dışı transfer deneylerinin sonucu ne çıkarsa çıksın bir ölçüm tablosudur, ortaya konulabilecek bir ürün değil — ve çalışmanın ilk niyeti "bir ilacın etiketindeki hastalık dışında başka bir hastalıkta faydası olabilir mi" sorusuydu.

Korunan varlık (hiçbiri yeniden yapılmıyor):

| Varlık | Durum |
| --- | --- |
| ChEMBL 37 yerel kopya, sürüm + SHA-256 kayıtlı | Hazır |
| 283.353 hedef–bileşik çifti, 565 hedef, 190.339 molekül, 5 aile | Hazır |
| 13 rejim × 5 tohum = 65 bölme dosyası + kaçak denetimi | Hazır |
| 1.806 satır kütüphane, 2.497 satır script, 881 satır test (86 geçiyor) | Hazır |
| Tezde yeri olan 7 tablo + 26 karar kaydı (K1–K26) | Hazır |

Revizyon 2'ye göre değişenler:

1. **Araştırma sorusu yerini ölçülebilir bir ürün iddiasına bıraktı.** "Dil modeli temsilleri görülmemiş bir aileye aktarılabilir mi" sorusu dondu; yerine "bilinen ilaç–hedef çiftlerini geri bulabiliyor muyum, ve nerede bulamıyorum" geldi.
2. **R4a/R4b (aile dışı transfer) kapsam dışı.** Bölme dosyaları duruyor, kullanılmıyor.
3. **Hastalık katmanı eklendi.** Hedef → hastalık eşlemesi, zincirin eksik son halkası.
4. **Zamansal bölme eklendi.** Ürünün kalite beyanı için en dürüst bölme.
5. **M2 (kısmi ince ayar), D2b transfer matrisi, D4 seçicilik testi ve bütün ablasyonlar kesildi.**
6. **Güven katmanı ürünün çekirdeği oldu** — süs değil, ana işlev.

Bu belge rev2'nin yerine geçer; rev2 ve `decisions.md` K1–K26 geçerliliğini korur ve veri hattı bölümleri aynen kullanılır.

## Ürün tanımı

Ürün, bir molekülün tahmin edilen protein hedeflerini o hedeflerin hastalık ilişkileriyle birleştirip **sıralanmış, belirsizliği işaretli bir hipotez listesi** üreten bir web aracıdır. İki görünüm, tek motor.

### Görünüm A — İlaçtan hastalığa

Girdi: bir SMILES ya da onaylanmış ilaç adı.

Çıktı, sırayla:

1. Hedef profili — kapsamdaki hedefler üzerinde tahmin edilen pChEMBL, aile bazında gruplandırılmış
2. Her satırda: kalibre edilmiş aralık, kimyasal yenilik işareti, hedef güvenilirlik işareti
3. Kanıt paneli — o hedef için **ölçülmüş** en yakın analoglar, gerçek değerleriyle ve ChEMBL bağlantısıyla
4. Hastalık katmanı — üst hedeflerin hastalık ilişkileri, kanıt türü ve skoruyla
5. İlacın **bilinen** hedefleri ve etki yönü ayrıca işaretli; etiket endikasyonları listelenir

İki eksenli sıralama zorunlu: ham tahmin ile hedefin eğitim dağılımına göre standartlaştırılmış skor birlikte gösterilir. Aileler arası pChEMBL ortalama farkı 0,72 log birim olduğu için ham değerler hedefler arası doğrudan karşılaştırılamaz.

### Görünüm B — Hastalıktan ilaca

Girdi: bir hastalık terimi. Çıktı: o hastalıkla ilişkilendirilmiş hedefler, ve her hedef için modelin bağlanmasını beklediği ChEMBL faz 4 bileşikleri — aynı güven işaretleriyle. Savunmada anlatılması en kolay görünüm bu.

**Sıralama kuralı:** hastalık ilişki skoru ile tahmin edilen aktivite **çarpılmaz ya da toplanmaz** — ikisi farklı türden kanıt ve birleşik bir skorun bilimsel gerekçesi yok. İki ayrı sıralama gösterilir: önce hastalıkla ilişkilendirilmiş hedefler (Open Targets skoruna göre), sonra seçilen hedef için tahmin edilen ilaçlar (aktivite tahminine göre). Böylece kullanıcı hangi kanıtın hangi sıralamayı ürettiğini görür.

**Sıra:** Görünüm A ilk sürümde çıkar, Görünüm B aynı motorun ikinci arayüzü olarak sonra eklenir. İlk çalışan sürümün yapması gereken üç şey: SMILES ya da desteklenen ilaç adını almak, hedefleri güven sinyalleriyle sıralamak, hastalık ilişkilerini ve deneysel kanıtı göstermek.

### Hedef kullanıcı ve karar

Kullanıcı: hangi deneyi önce yapacağına karar vermeye çalışan bir araştırmacı. Araç o kararı vermiyor, sıralama öneriyor ve her sıranın ne kadar güvenilir olduğunu söylüyor.

### Ürünün asla iddia etmediği şeyler

- Terapötik fayda. Çıktı birimi, yönü bilinmeyen bir **bağlanma** hipotezidir.
- Etki yönü. Afinite "ne kadar sıkı bağlanıyor"u söyler, aktive mi inhibe mi ettiğini söylemez. Yön yalnız **bilinen** çiftler için, ChEMBL mekanizma tablosundan gösterilir.
- Nedensellik. Hedef–hastalık ilişkisi "bu hedefi vurmak hastalığı iyileştirir" demek değildir.
- "Bağlanmaz" yargısı. Veride gerçek inaktif kayıt yok; araç yalnız sıralar.
- Doz, dokuda erişim, kan–beyin bariyeri, güvenlilik — hiçbiri modelde yok.
- Yeni kullanım keşfi. Bkz. bir sonraki bölüm.

## Ölçülebilir iddia

Tezin iddiası keşif değil **geri bulma**: eğitimden çıkarılmış, bilinen ilaç–hedef çiftlerinin ne kadarını aday listesinin üstüne taşıyabiliyorum, ve bu oran yenilik arttıkça nereye düşüyor?

Bu formülün seçilme sebebi: ChEMBL zaten geniş bir off-target ölçüm kütüphanesi olduğu için model büyük ölçüde *bilinen* ikincil hedefleri yeniden bulacak. "Yeni bir kullanım buldum" iddiası bu veriyle doğrulanamaz; "bilineni ne kadar iyi buluyorum" iddiası doğrulanabilir.

### İki düzey, iki ad

Projenin adı yeniden konumlandırmayı, ölçülebilir iddiası geri bulmayı anlatıyor. Bu ikisi aynı şey değil ve karıştırılmamaları için ayrı adlandırılır:

- **Bilimsel iddia:** bilinen ilaç–hedef ilişkilerini ne kadar iyi geri buluyoruz ve hipotezleri ne kadar iyi önceliklendiriyoruz.
- **Kullanım bağlamı:** bu bilgi hastalık ilişkileriyle birleştirildiğinde hangi deneylerin araştırılmaya değer olduğu.

Tezin iddia cümleleri birinci düzeyden, ürünün tanıtım cümleleri ikinci düzeyden kurulur. **Bir ilacın bilinen hedefini geri bulmak, o ilacın yeni bir hastalıkta işe yarayacağını göstermez** — bu cümle tezin girişinde ve ürünün sınırlar sayfasında birebir geçer.

Birincil metrikler, hepsi yenilik kovası başına ayrı raporlanır:

| Metrik | Ne söyler |
| --- | --- |
| Recall@%1 | Bilinen pozitiflerin ne kadarı aday havuzunun üst yüzde birine giriyor |
| MRR | Bilinen pozitifin ortalama karşılıklı sırası |
| NDCG@20 | Üst sıralamanın kalitesi |
| Aralık geçerliliği | %80 aralığın gerçekten vakaların %80'ini kapsayıp kapsamadığı |
| Aralık genişliği | Aralığın kullanışlı kalacak kadar dar olup olmadığı |

**Neden "zenginleşme" denmiyor.** Veride gerçek inaktif kayıt yok: ölçülmemiş bir ilaç–hedef çifti negatif değildir, bilinmiyordur. "Zenginleşme" dili aday havuzunun geri kalanının negatif olduğu izlenimini verir ve metriği olduğundan güçlü gösterir. Bu yüzden birincil metrik **Recall@%1** olarak adlandırılır ve üç metriğin de yalnız **anotasyonlu pozitifleri geri bulma** performansını ölçtüğü hem tezde hem arayüzde yazılır. Hesap aynı kalıyor; değişen adı ve yorumu.

Karıştırılmayacak iki şey: bilinen bir pozitifi geri bulamamak başka şey, bir çiftin gerçekten etkileşimsiz olduğunu göstermek başka şey. İkincisi bu veriyle yapılamaz ve ürün hiçbir yerde ima etmez.

Son iki metrik ürünün dürüstlük ölçüsüdür: bir güven aralığı geçerli ama aşırı genişse model "bilmiyorum" demeyi öğrenmiş olur, ki bu da bir sonuçtur. Beklenti çıpası literatürden: biyoaktivite modellemesinde %80 güven düzeyinde ortalama tahmin aralığı 1,65 pIC50 birimi civarında ölçülmüş ([Svensson ve ark.](https://discovery-pp.ucl.ac.uk/id/eprint/10060460/1/Svensson_etal_r1.pdf)). Yani ±1,6 log birimlik aralıklar normaldir ve ürünün dili buna göre kurulmalı: sıralama aracı, nokta tahmin aracı değil.

Karşılaştırma zemini: rastgele sıralama, B1 hedef ortalaması ve M0 (ECFP4 + gradyan artırma). M0'ı geçemeyen bir dil modeli hattı tezin bulgusu olur, başarısızlığı değil.

## Veri katmanları

Dört katman var; birincisi hazır, ikincisi ChEMBL'in içinde zaten duruyor, üçüncüsü yeni, dördüncüsü opsiyonel.

| # | Katman | Kaynak | Durum | Rolü |
| --- | --- | --- | --- | --- |
| 1 | Afinite çekirdeği | ChEMBL 37 | Hazır | Modelin eğitim verisi: 283.353 çift |
| 2 | İlaç ve mekanizma | ChEMBL `molecule_dictionary.max_phase`, `drug_mechanism`, `drug_indication` | Çekilecek | ChEMBL faz 4 kümesi, bilinen hedefler, **etki yönü**, etiket endikasyonları |
| 3 | Hastalık ilişkileri | Open Targets Platform | Yeni | Hedef → hastalık, kanıt türü ve skoruyla |
| 4 | Dış doğrulama | Papyrus | Opsiyonel | Bağımsız kürasyonda aynı ölçüm tekrarı |

### Katman 2 — neden kritik

Etki yönü sorununun tek dürüst çözümü bu tablo. Model yön tahmin etmiyor; ama bilinen çiftler için ChEMBL'in `action_type` alanı (inhibitor, agonist, antagonist …) var ve ürün bunu gösteriyor. Yön tahmini gelecek iş olarak [DTIAM](https://arxiv.org/pdf/2312.15252) ile adlandırılır — o çalışma etkileşim, afinite ve aktivasyon/inhibisyon mekanizmasını tek çerçevede tahmin ediyor, yani problem bilinen ve kısmen çözülmüş — ama bu projenin kapsamında değil.

**Veri evreni ve adlandırma.** `max_phase = 4` kümesi metinde ve arayüzde **ChEMBL faz 4 bileşikleri** olarak anılır, "onaylı ilaçlar" diye genellenmez: küme geri çekilmiş ilaçları da içerir ve ülke bazında güncel onay durumunu yansıtmaz. Aynı şekilde `drug_indication` kayıtları resmî ürün etiketiyle eşitlenmez; arayüzde kaynağı ve kapsamı yazılır. Geri alma protokolünün evreni bu kümedir.

### Katman 3 — eşleme ve yönü

Eşleme zinciri: ChEMBL hedef → UniProt → Ensembl gen → Open Targets. Her adımda kayıp olur; **eşleme kapsamı oranı ölçülür ve düşürme tablosuna satır olarak yazılır.**

Open Targets'ın kendi çerçevesi ürünün diline örnek alınır: ilişki skorları kanıtı sıralamak ve düzenlemek içindir, bir hedefin klinik doğrulaması değildir. Ayrıca platform sekiz kanıt kaynağı için **etki yönü değerlendirmesi** sunuyor — bir hedefin işlevindeki genetik değişimin bir özellik için risk mi koruma mı sağladığını gösteriyor. Bu, zincirin son halkasında kısmi bir yön bilgisi demek ve doğrudan arayüze girer.

### Sürüm sabitleme

K1 ve E1 kuralları aynen geçerli, iki ekle:

- **Hugging Face kontrol noktası commit hash'leri gömme çıkarımından önce sabitlenir.** `config/config.yaml` içinde şu an `null`; bu haliyle çalışılan hiçbir sonuç tekrarlanabilir değil. Bu, hattaki sonraki adımın ön koşulu.
- Open Targets sürüm etiketi ve indirme tarihi `docs/versions.md` dosyasına yazılır.

## Protein temsili ve hedef kapsamı

Beş aile kısıtı kalktı; yerine iki aşamalı bir kapsam kararı geldi. Gerekçe: beş aile, leave-one-family-out bölmelerini kurabilmek için vardı. O soru dondurulduğuna göre aile artık modelin girdisi ya da bölme ölçütü değil, yalnız bir **gruplandırma etiketi** — ve ürün için kapsam ne kadar genişse o kadar iyi.

| Aşama | Kapsam | Gerekçe |
| --- | --- | --- |
| 1 (zorunlu) | Mevcut 5 aile, 565 hedef | Veri hazır, titiz, küratlü; ürünün çekirdeği ve tezin ana ölçümleri bunun üzerinde |
| 2 (zaman kalırsa) | Tüm insan tek-protein hedefleri, aynı kalite süzgeciyle | Erişilebilir hastalık yelpazesini genişletir; seyrek hedefleri hedef güvenilirlik işareti soğurur |

Aşama 2'ye geçmeden ölçülmesi gerekenler — üçü de yeni aileler için bozulabilir:

1. **Ölçüm türü kapsamı.** Mevcut süzgeç Ki/Kd/IC50 + kesin `=` + güven skoru eşiği. Taşıyıcı ve kanal sınıflarında fonksiyonel EC50 baskın; bu süzgeç o aileleri boşaltabilir. Ölçülmeden aile eklenmez.
2. **Bağlanma alanı kapsamı (E1b tekrarı).** Çok alanlı, uzun proteinlerde PPDMs eşlemesi daha kötü. Kapsam düşerse o aile için birincil temsil tam diziye döner ve kayıt işaretlenir.
3. **Sınıf envanteri.** ChEMBL protein sınıflandırması elle bakılan özel bir sınıflandırma ve yeni sınıflar ihtiyaç oldukça ekleniyor ([ChEMBL SSS](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/target-questions)), o yüzden envanter adımı her yeni aile için yeniden koşturulur.

Aile çakışması serbest bırakılır. "Sıfır çakışma" özelliği R4'ün ön koşuluydu; üründe bir hedefin iki sınıfa girmesi sorun değil, sadece gruplamada ilk atama gösterilir.

### Açık teknik soru — güven skoru 8 mi 9 mu

Rev2'de eşik 8 ve gerekçe "deneyin tek bir proteine atandığı anlamına gelir" olarak yazılmış. ChEMBL ölçeğinde ise 8 "homolog tek protein hedefi atandı", 9 "doğrudan tek protein hedefi atandı" anlamına geliyor ([ölçek](https://targetmine.mizuguchilab.org/?p=722), ikincil kaynak). Yani ≥ 8 eşiği homolog atamaları da içeriyor olabilir.

**Yapılacak:** yerel ChEMBL 37 kopyasında doğrulanır. Doğruysa iki seçenek var — eşik 9'a çekilir, ya da 8 ile 9 arasındaki fark tek bir tablo olarak raporlanıp eşik 8'de bırakılır. Kararı veri hattı yeniden koşmadan önce ver; sonrasında maliyeti büyüyor.

## Model ve taban çizgileri

Motor rev2'nin M1'i: donuk ESM-2 ve ChemBERTa-2 gömmeleri üzerine sığ bir birleştirme başlığı, pChEMBL regresyonu. Amacı en iyi skoru almak değil, ürünün arkasında savunulabilir bir motor olmak.

| Model | Durum | Rolü |
| --- | --- | --- |
| B1 — hedef ortalaması | Zorunlu | Gerçek alt sınır. Hedef güvenilirlik işaretinin temeli |
| B2 — yalnız ligand | Zorunlu | Ligand yanlılığını ölçer |
| B3 — yalnız protein | Zorunlu | Hedefin genel ölçüm düzeyini ölçer |
| M0 — ECFP4 + gradyan artırma | Zorunlu | Gerçek rakip; dil modelini geçmesi gerçek bir olasılık |
| M1 — donuk gömmeler + sığ başlık | Zorunlu | Ürünün motoru |
| M2 — kısmi ince ayar | **Kesildi** | 15–30 saat hesap, ürüne katkısı belirsiz |

Taban çizgileri kesilemez — ve bunun literatürde açık bir gerekçesi var: DTI'yı bağımsız varlıklar gibi ele almak bilgi kaçağına yol açıyor, çünkü ilacın ya da hedefin tek başına çağrışımı bile öngörü gücü taşıyor ([Bai ve ark.](https://www.biorxiv.org/content/10.1101/2022.04.20.488898v1.full)). B2 ve B3 tam olarak bu gücü ölçer; onlar olmadan M1'in skoru yorumlanamaz.

Eğitim kurgusu:

- Kayıp: pChEMBL üzerinde ortalama kare hata
- Hiperparametreler eğitim içinden ayrılan doğrulama kümesiyle seçilir; test kümesine bakılmaz
- **3 tohum** (rev2'de 5); raporlanan değer ortalama ve yayılım
- Test kümesi bir kez açılır, tüm modeller ve rejimler hazır olduğunda
- Donanım: tek RTX 4050, 6 GB VRAM. Donuk gömme mimarisi bu donanımda zorunluluk, aynı zamanda ürünün CPU'da çalışmasının sebebi

### Ne tahmin ediliyor: ölçüm türü meselesi

Süzgeç Ki, Kd ve IC50'yi birlikte alıyor ve üçü pChEMBL ekseninde birleşiyor. Üçü aynı şeyi ölçmüyor: Kd bağlanma afinitesi, Ki inhibisyon sabiti, IC50 ise belirli bir deney koşulundaki yarı maksimum inhibisyon derişimi. Modelin çıktısını "tahmin edilen bağlanma afinitesi" diye adlandırmak, veride IC50 varken yanlış olur.

İki karar:

1. **Çıktı adı değişiyor:** her yerde "tahmin edilen aktivite (pChEMBL)". "Bağlanma afinitesi" ifadesi yalnız Ki/Kd alt kümesinden söz edilirken kullanılır.
2. **Ölçüm envanteri çıkarılır:** `standard_type` bazında aile başına dağılım, ve aynı çiftte hem Ki/Kd hem IC50 ölçülmüşse aradaki sistematik fark. Mevcut veriden bir sorgu; yeni iş değil.

Envanter IC50'nin baskın ya da sistematik olarak kaymış olduğunu gösterirse, Ki/Kd ile eğitilmiş ikinci bir model **tek ablasyon olarak** geri gelir. Göstermezse tek model kalır ve ölçüm karışımı sınırlılık olarak yazılır. Karar envanterden sonra verilir, baştan değil.

Kapsam dışı bırakılan deneyler: D2a (LOFO), D2b (transfer matrisi), D4 (seçicilik testi), ve bütün ablasyonlar. D3 karıştırma kontrolünün **aile içi varyantı kalır** — modelin proteini gerçekten okuyup okumadığını söyleyen tek test ve ucuz. D5 (hedef ön yanlılığı) kalır, bir saat.

## Bölme rejimleri

Dört rejim kullanılır; biri yeni ve ürünün kalite beyanı ona dayanır.

| Rejim | Test kümesinde ne yeni | Rolü |
| --- | --- | --- |
| R1 — Rastgele | Hiçbir şey | Yalnız üst sınır ve literatürle karşılaştırma zemini |
| R2 — Görülmemiş bileşik | Bileşik (kümelenmiş) | Yeni moleküle genelleme |
| R3 — Görülmemiş hedef | Protein (dizi kimliği kümelenmiş) | Yeni hedefe genelleme |
| **RT — Zamansal** | **İlk kamuya açık raporlama tarihi** | **Dağıtıma en yakın ölçüm; ürünün ilan edilen isabet oranı buradan çıkar** |

### RT neden birincil

Alanın bölme taksonomisinde zamansal bölme dağıtıma en yakın olanı ([derleme](https://pmc.ncbi.nlm.nih.gov/articles/PMC12818121/)). Mantık şu: gerçek kullanımda model geçmişle eğitilir, gelecekte sorgulanır. ChEMBL, biyoaktivite ölçümleri için yayın tarihi bilgisi sağlayan az sayıdaki açık veri tabanından biri, yani ek veri gerekmiyor.

Kurulum: bir kesim yılı seçilir, o tarihe kadar yayınlanmış ölçümlerle eğitilir, sonrasında yayınlanmışlarla test edilir. Kesim yılı **veriye bakılmadan**, test kümesinde yeterli kayıt bırakacak şekilde önceden sabitlenir ve gerekçesi K kaydına yazılır.

### RT'nin zorunlu sızıntı kontrolü

Zamansal bölmenin tek ciddi açığı: **aynı molekül–hedef çifti hem eğitim hem test döneminde ölçülmüş olabilir.** 2022'de ölçülen bir çift 2025'te başka bir çalışmada tekrar ölçülürse, modelin 2025 kaydını bilmesi yeni bir etkileşim bulmak değil, gördüğü bir değeri tekrarlamaktır.

Bu, veri hattının mevcut yapısıyla doğrudan çakışıyor: normalizasyon adımı aynı çiftin birden fazla ölçümünün **medyanını** alıyor, yani tekilleştirilmiş veri setinde çiftin tarihi yok. RT o veri setinin üzerine kurulamaz.

Kural:

1. RT bölmesi tekilleştirme **öncesindeki** aktivite kayıtlarından kurulur; her çift için ilk kamuya açık raporlama tarihi hesaplanır.
2. İlk kamuya açık raporlaması kesim tarihinden önce olan çiftler test kümesinden çıkarılır — kesim sonrası tekrar ölçülmüş olsalar bile.
3. Testte kalan çift sayısı ve çıkarılan çift sayısı birincil sonuç olarak raporlanır; R4b'deki süzme raporlamasının aynısı.

Bu kontrol olmadan RT skoru rastgele bölmeden **daha iyimser** bile çıkabilir — yani düzeltmeden kullanılması, bölmeyi hiç eklememekten kötü.

### Tarih alanının ne olduğu

ChEMBL'in verdiği tarih deneyin yapıldığı tarih değil, kaydın bağlı olduğu dokümanın yayın tarihidir. Bu yüzden bölme ölçütü belgede **ilk kamuya açık raporlama tarihi** olarak adlandırılır; "ilk ölçüm tarihi" denmez.

İki sonuç sınırlılıklara yazılır:

- Ölçümün yapılmasıyla kamuya açılması arasında gecikme var; kesim tarihi deneysel değil **yayınsal** bir sınırdır ve RT gerçek zaman sırasını yaklaşık olarak temsil eder.
- Kullanılabilir tarihi olmayan kayıtların ne olacağı önceden sabitlenir: eğitim tarafına yazılır, test evrenine hiç alınmaz. Sayıları düşürme tablosunda raporlanır.

### R2 hakkında bilinen şey

Kendi ölçümümüz: R2 test bileşiklerinin **%81,7'si** eğitim kümesinde ≥0,40 Tanimoto'da bir komşuya sahip, ortanca benzerlik 0,548, beş tohumda 0,805–0,828. Bu bir keşif değil, literatürdeki bulgunun kendi veri setimizde doğrulanması: 2026'daki bir çalışma iskelet ve parmak izi bölmelerinde en yakın eğitim komşusuna göre medyan 3-mer Jaccard değerini 1,00 bulurken protein ve çifte-OOD bölmelerinde 0,139 ve 0,142'ye düştüğünü gösteriyor ([Beyond Random Splits](https://arxiv.org/html/2610.03456)).

**Sunum kuralı:** bu bulgu tezde "literatürdeki bulguyu kendi veri setimde ölçtüm ve bir ürün kararına çevirdim" olarak yazılır. Kimyasal yenilik işaretinin gerekçesi budur.

### DataSAIL karşılaştırması

Bölmeler elle kuruldu ve kaçak için standart bir metrik kullanılmadı. [DataSAIL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11978981/) tam bu işin aracı: iki boyutlu etkileşim verisinde her iki boyutta benzerlik-duyarlı bölme (S2) yapıyor ve bir bölmenin kaçağını ölçmek için ölçeklenmiş L(π) metriğini tanımlıyor.

Karar: **kendi bölmeleri korunur, DataSAIL ile de üretilir, L(π) ile karşılaştırılır.** Bir hafta iş; "neden kendi bölmecini yazdın" sorusunu ölçümle kapatır ve tezde bir tablo olur.

### Dondurulan rejimler

R4a ve R4b (aile dışı transfer) kullanılmaz. 65 bölme dosyası olduğu yerde kalır; sonraki çalışmanın giriş noktası bunlar.

## Güven katmanı

Bu katman ürünün ana işlevidir, süsü değil. Alanın kendi öncelik listesinde üç baskın sorundan ikisi burada: veri kaçağı/fazlalık ve **uygulanabilirlik alanı boşlukları** ([derleme](https://pmc.ncbi.nlm.nih.gov/articles/PMC12818121/)). Üçüncüsü — pozitif/negatif dengesizliği, negatiflerin uydurma olması — pChEMBL regresyonu kurgulandığı için bu projede zaten yok; bu bir tasarım tercihi olarak tezde yazılır.

Her tahminin yanında dört şey durur:

| # | Sinyal | Nasıl hesaplanır | Eşikler |
| --- | --- | --- | --- |
| 1 | Kimyasal yenilik | Eğitim molekülleriyle maksimum ECFP4/Tanimoto | ≥0,7 yakın analog · 0,4–0,7 tanıdık kimya · <0,4 yeni kimya |
| 2 | Hedef güvenilirliği | Modelin o hedefte B1'i geçme oranı, çevrimdışı hesaplanıp tabloya yazılır | Geçmiyorsa gri işaret: "bu hedefte model ortalamadan iyi değil" |
| 3 | Tahmin aralığı | Split-conformal, yenilik kovası başına ayrı kalibre edilir; tohum yayılımı ayrı gösterilir | %80 güven düzeyi varsayılan |
| 4 | Kanıt paneli | O hedef için en yakın k ölçülmüş analog, gerçek değerleriyle | — |

### 1. sinyalin eşikleri neden bu değerler

0,40 eşiği keyfi değil: R2 bölmesinde kullanılan eşik bu ve test bileşiklerinin %81,7'sinin o eşikte bir eğitim komşusu olduğu ölçüldü. Yani kullanıcıya gösterilen "tanıdık kimya" etiketi, kendi veri setimizde ölçülmüş bir gerçeğe dayanıyor.

### 3. sinyal — öncülü ve sınırı

Yerel konformal yaklaşımın birebir öncülü var: dinamik uygulanabilirlik alanı (dAD) çalışması bileşik–hedef bağlanma afinitesi tahmininde kalibrasyon çiftlerini bileşik–hedef komşuluğundan alıyor ([dAD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10457664/)). Bu projede **basit sürüm uygulanır** (kova başına split-conformal), dAD atıf verilir ve "yerel kalibrasyon uygulanmadı" sınırlılıklara yazılır.

**İki belirsizlik kaynağı birleştirilmez.** Üç tohumun tahmin yayılımı model varyansını ölçer; split-conformal aralık tahmin hatasının dağılımını ölçer. Bunlar farklı şeyler ve tek bir "kalibre edilmiş aralık" diye sunulamazlar. Arayüzde ve raporda ayrı dururlar: tohum yayılımı küçük bir yan gösterge, tahmin aralığı birincil sayı.

Terminoloji: üretilen şey bir **tahmin aralığı** (prediction interval), güven aralığı değil. Arayüzde bu ad kullanılır.

### Aralığın geçerlilik sınırı

Konformal garanti, kalibrasyon ve test örneklerinin değiştirilebilir (exchangeable) olmasına dayanır. Eğitim verisi yalnız **ölçülmüş** aktivite kayıtlarından oluştuğu için %80 kapsama garantisi, hiç ölçülmemiş ilaç–hedef kombinasyonlarına otomatik olarak taşınmaz — ve ürünün en çok kullanılacağı yer tam olarak orası.

Dört kural:

1. **Gerçekleşen kapsama oranı** her yenilik kovası için ayrı ölçülür ve raporlanır — nominal %80 değil, ölçülen değer.
2. **Yeterli kalibrasyon örneği olmayan kovalar işaretlenir** ve o kovalarda aralık gösterilmez; yerine "kalibre edilemedi" yazar.
3. Kova sınırları kalibrasyon artıklarına bakılmadan önce sabitlenir.
4. Hedef güvenilirlik işareti (sinyal 2) **kat dışı** hesaplanır. Modelin B1'i geçip geçmediğine kendi eğitim verisiyle karar vermek başlı başına bir sızıntıdır ve işareti olduğundan iyimser gösterir.

### Hedef güvenilirlik işaretinin kendi örneklemi

Kat dışı hesaplama tek başına yetmiyor: bir hedefte yalnız sekiz kat dışı tahmin varsa B1'i geçme oranı hesaplanabilir ama güvenilirlik göstergesi olarak sunulamaz. Dört ek kural:

1. Hedef başına kat dışı değerlendirme örnek sayısı tabloda tutulur ve arayüzde görünür.
2. Minimum örnek sayısı veriye bakılmadan sabitlenir.
3. Altında kalan hedefte gri işaret değil ayrı bir **"değerlendirme yetersiz"** işareti çıkar. "Model burada iyi değil" ile "burada bilmiyoruz" aynı şey değil ve kullanıcıyı farklı yöne götürür.
4. Geçme oranının yanında bir belirsizlik aralığı verilir; binom aralığı yeterli.

İşaret **hedef düzeyinde** kalır, yenilik kovası düzeyine inmez: kovalara bölmek örneklemi parçalar ve neredeyse hiçbir hedefte yeterli sayı kalmaz.

### 4. sinyalin stratejik değeri

Kanıt paneli, modelin başarısız olması halinde **ürünü ayakta tutan şeydir**: geriye bir erişim/arama aracı kalır, kullanıcı modele güvenmek zorunda kalmadan en yakın ölçülmüş analogları görür. M1'in zayıf çıkması riski bu panelle tasarım düzeyinde kapatılmış durumda.

### Kullanıcının yaptığı şey ne değişir

Dört sinyal bir "güvenilirlik rozeti"ne indirgenmez. Üçü teknik nedenlerle farklı şeyler söyler: kimyasal yenilik molekülün, hedef güvenilirliği hedefin, aralık ikisinin birleşiminin durumunu anlatır. Tek bir skora toplamak hangisinin bozuk olduğunu saklar.

## Doğrulama protokolü

Dört test var; ilki ürünün kalite beyanı, diğerleri onu yorumlanabilir kılan şeyler.

### D6' — Geri alma (birincil)

ChEMBL faz 4 bileşiklerinin bilinen hedef ilişkilerinden bir bölüm eğitimden çıkarılır, model o çiftleri geri bulabiliyor mu diye ölçülür. Üç şey kesin olarak sabitlenir:

1. **Aday havuzunun bileşimi.** Havuz neyi içeriyorsa metrik onu ölçer; tanımı önceden yazılır ve tezde birebir geçer.
2. **Yenilik kovası başına raporlama.** Tek bir ortalama sayı bu ayrımı gizler.
3. **Özellikle etiket endikasyonundan farklı alana düşen çiftler ayrı raporlanır** — yeniden konumlandırmaya en yakın alt küme bu.

### D6' ile RT'nin ilişkisi — sıra sabittir

İki protokol üst üste biniyor; sıra belirsiz kalırsa sonuç yorumlanamaz. Sıra şudur:

1. Önce RT uygunluğu belirlenir: her çiftin ilk kamuya açık raporlama tarihi hesaplanır.
2. Raporlaması kesim tarihinden önce olan çiftler RT test evreninden çıkarılır.
3. D6' pozitifleri **bu uygun evrenin içinden** seçilir.
4. Bir çift tek kez sayılır; aynı çift farklı deney kayıtları üzerinden tekrar tekrar başarı olarak sayılmaz.

**Sonuçlar ayrı raporlanır.** D6' bağımsız bir geri bulma ölçümü, RT ise dağıtıma benzer bir rejim; tek bir başarı oranında birleştirilmezler.

**Güç sorunu, önceden kabul edilmiş.** RT kısıtı D6' pozitiflerini ciddi biçimde azaltabilir — faz 4 bileşiklerinin bilinen hedeflerinin çoğu eski kayıtlardır. Bu yüzden D6' iki kipte koşulur:

| Kip | Pozitif seçimi | İddia gücü |
| --- | --- | --- |
| D6'-RT | Yalnız RT-uygun çiftler | Güçlü, ama örneklem küçük olabilir |
| D6'-rastgele | Zamandan bağımsız tutma kümesi | Örneklem büyük, iddia zayıf |

E3 ile aynı mantıkla minimum pozitif sayısı **veriye bakılmadan** şimdi sabitlenir: altında kalan kipte skor verilmez, yalnız kalan pozitif sayısı bildirilir.

### Özgüllük kontrolü (ucuz, zorunlu)

Rastgele moleküller verildiğinde profilin düz kalması gerekir. Kalmıyorsa model ligand yanlılığı öğrenmiş demektir. Yorum sınırı: ölçümün yokluğu bağlanmamanın kanıtı değil, o yüzden bu bir **duman testi** olarak adlandırılır, özgüllük kanıtı olarak değil.

### D3 — Aile içi karıştırma

Protein gömmesi aynı aileden rastgele başka bir proteininkiyle değiştirilir. Performans düşmüyorsa model aile düzeyinde bir imza öğrenmiş, hedefe özgü bilgi değil. Maliyeti sıfır — gömmeler diskte.

### Kullanıcı testi (yeni)

Üniversitede eczacılık ya da biyoloji tarafından **bir kişiye, yirmi dakika**, hiçbir şey açıklamadan aracı kullandır. Neyi yanlış anladığı yazılır ve arayüz düzeltilir.

Çıktı tezde bir bulgu olur: "kullanıcı tahmin aralığını X olarak okudu, uyarı metnini şöyle değiştirdim." Maliyeti sıfır, "ürün" iddiasını kanıta çevirir, ve açıklanabilir/itiraz edilebilir yapay zekâ çerçevesine oturur. Tek kişi istatistik değil — "tek katılımcılı kullanılabilirlik gözlemi" olarak adlandırılır.

## Mimari ve teslim paketi

&#91;embedded content: mimari · üç katman, iki görünüm\]

Mimarinin tek kritik özelliği: 565 hedef gömmesi bir kez çıkarılıp diske yazıldığı için çıkarım anında yalnız ligand kolu çalışıyor. 6 GB VRAM kısıtı eğitim tarafında kalıyor, ürün tarafına hiç bulaşmıyor.

### Yığın

| Katman | Seçim | Gerekçe |
| --- | --- | --- |
| Servis | FastAPI | Python ekosistemi, model kodu aynı dilde |
| Çıkarım | ONNX'e çevrilmiş ChemBERTa-2, CPU | GPU gereksinimini kaldırır |
| Depolama | SQLite / parquet | Veritabanı sunucusu gerekmiyor |
| Arayüz | Angular | Mevcut yığına oturuyor |
| Paketleme | Docker Compose | Savunmada yerelde açılabilir |

### Dağıtılan varlıklar

Hedef gömme matrisi \~3 MB, başlık ağırlıkları küçük, 190k ECFP4 paketlenmiş \~50 MB, ölçüm tablosu 12 MB, hastalık eşlemesi. Toplam \~65 MB — tek imajda taşınabilir.

### Demo dayanıklılığı

Savunma günü demo çalışmıyorsa ürün yok. Üç yedek zorunlu:

- [ ] Yerelde `docker compose up` ile açılan sürüm — internet gerektirmez
- [ ] 3 dakikalık kayıtlı demo videosu
- [ ] Ekran görüntülü yedek slayt

### Teslim paketi

Çalışan uygulama, GitHub deposu, model kartı ve veri kartı, API dokümantasyonu, tez, demo videosu.

## Dil kuralları ve etik sınırlar

Bu projenin riski teknik değil **retorik**. "Bu ilaç şu hastalıkta işe yarayabilir" çıktısı, dili bir kez gevşerse tıbbi bir iddiaya dönüşür. Bu yüzden dil kuralı bir üslup tercihi değil, bir karar kaydıdır.

> **K27 — Ürün hiçbir çıktıda terapötik fayda iddia etmez; çıktı birimi, yönü bilinmeyen bir bağlanma hipotezidir.**

Sonraki her ekran, her metin ve her API yanıtı bu kayda bağlanır.

### Kullanılmayan ve kullanılan kelimeler

| Asla | Yerine |
| --- | --- |
| "tedavi eder", "etkilidir", "işe yarar" | "bağlanması bekleniyor" |
| "önerilir", "aday ilacınız" | "test edilmeye değer hipotez" |
| "yeni kullanım bulundu" | "bilinen çifti geri buldu" / "doğrulanmamış hipotez" |
| "kesin", "yüksek güvenle" | aralığın kendisi, sayıyla |

### Nerede zorunlu

Dil kuralı tek bir yerde değil, beş yerde aynı olmak zorunda: arayüz metinleri, API yanıt alanlarının adları, model kartı, tez metni, savunma sunumu. Birinde gevşerse diğerleri de geçersizleşir.

### Yön uyarısı — kalıcı ve görünür

Hastalık katmanının göründüğü her yerde, kapatılamayan bir satır: *bir hedefe bağlanmak o hedefi aktive mi inhibe mi ettiğini söylemez; bir hastalıkla ilişkili hedefe bağlanan bir ilaç o hastalığı kötüleştirebilir.*

### Kapsam beyanı

Araç "araştırma ve eğitim amaçlı triyaj aracı" olarak etiketlenir. Klinik karar desteği değildir ve bu ayrım ana ekranda yazılıdır — bir alt sayfada değil.

## Lisans ve atıf

Ürün türetilmiş veri **dağıttığı** için lisans artık bir formalite değil. ChEMBL her sürümüne CC BY-SA 3.0 uyguluyor: atıf gerektiriyor, ve share-alike maddesi veriyi başkalarına dağıttığında onlara da aynı yeniden dağıtma hakkını vermeni şart koşuyor ([ChEMBL lisans açıklaması](https://chembl.blogspot.com/2012/07/public-popular-chemistry-databases-and.html)).

Bu üründe dağıtılan türetilmiş veri: 12 MB ölçüm tablosu, komşuluk indeksi, hedef güvenilirlik tablosu. Dolayısıyla:

- [ ] Arayüzde ve depoda atıf + kullanılan ChEMBL sürümü (37) açıkça yazılı
- [ ] Dağıtılan türetilmiş veri share-alike altında; kod ayrı lisansla verilebilir, veri kısmı kapatılamaz
- [ ] Open Targets eklenince kaynak kırılımı saklanır — platform CC0 olsa da kaynak başına ayrı lisanslar işliyor
- [ ] ESM-2 ve ChemBERTa-2 kontrol noktalarının lisansları model kartına yazılır (kullanım öncesi doğrulanır)
- [ ] Papyrus eklenirse kendi lisans dosyası kontrol edilir

Atıf metni, her çıktıda aynı biçimde: veri kaynağı, sürüm etiketi, lisans bağlantısı.

Bu bir engel değil, bir satır iş — ama sonradan fark edilirse yayınlanmış bir aracı geri çekmek gerekir.

## Karar kayıtları (K27–K38)

`decisions.md` dosyasının devamı. K1–K26 geçerliliğini koruyor.

| # | Karar | Gerekçe |
| --- | --- | --- |
| K27 | Ürün hiçbir çıktıda terapötik fayda iddia etmez; çıktı birimi yönü bilinmeyen bir bağlanma hipotezidir | Retorik risk teknik riskten büyük; sonraki her ekran buna bağlanır |
| K28 | Aile dışı transfer (R4a/R4b, D2a, D2b) kapsam dışı; bölme dosyaları korunur | Sonucu teorik; lisans bitirme ölçeğini aşıyor. Dosyaların korunması sonraki çalışmanın giriş noktası |
| K29 | Aile, modelin girdisi ve bölme ölçütü olmaktan çıkıp gruplandırma etiketine indi | LOFO kalktığına göre aile kısıtının gerekçesi kalmadı; kapsam genişliği ürünün lehine |
| K30 | Zamansal bölme (RT) eklendi ve ilan edilen isabet oranının kaynağı oldu | Dağıtıma en yakın ölçüm; ChEMBL yayın tarihi zaten elde |
| K31 | Kesim yılı veriye bakılmadan sabitlenir | Sonradan seçilen eşik, katılık iddiasını zayıflatır (E3 ile aynı mantık) |
| K32 | Hastalık katmanı Open Targets'tan gelir; eşleme kapsam oranı düşürme tablosuna satır olarak yazılır | Eşleme zinciri kayıplı; ölçülmeyen kayıp raporlanamaz |
| K33 | Etki yönü tahmin edilmez; yalnız bilinen çiftler için ChEMBL mekanizma tablosundan gösterilir | Yön tahmini ayrı bir problem (DTIAM); kapsam dışı, ama boşluk açık bırakılamaz |
| K34 | Sıralama iki eksende gösterilir: ham tahmin + hedef içi standartlaştırılmış skor | Aileler arası pChEMBL ortalama farkı 0,72 log birim; ham değerler hedefler arası karşılaştırılamaz |
| K35 | Konformal aralık basit sürümde (kova başına split-conformal); dAD atıf verilir, uygulanmaz | Tam yerel kalibrasyon kapsam dışı; sınırlılık olarak yazılır |
| K36 | DataSAIL ile paralel bölme üretilir ve L(π) ile karşılaştırılır | Elle kurulmuş bölmenin kalitesi standart bir metrikle gösterilmeli |
| K37 | HF kontrol noktası commit hash'leri gömme çıkarımından önce sabitlenir | E1'in kendi kuralı; sabitlenmeden çıkarılan gömme bütün aşağı akışı tekrarlanamaz kılıyor |
| K38 | Dağıtılan türetilmiş veri CC BY-SA atıf ve share-alike koşullarıyla yayınlanır | ChEMBL lisansının gereği; ürün veri dağıttığı için artık bağlayıcı |

### Dış inceleme sonrası eklenen kararlar (K39–K49)

| # | Karar | Gerekçe |
| --- | --- | --- |
| K39 | Bilimsel iddia (geri bulma) ile kullanım bağlamı (yeniden konumlandırma) ayrı adlandırılır | İkisi aynı şey değil; karıştırılırsa iddia savunulamaz hale gelir |
| K40 | Model çıktısı "tahmin edilen aktivite (pChEMBL)"; "bağlanma afinitesi" yalnız Ki/Kd alt kümesi için | IC50 bağlanma afinitesiyle eş anlamlı değil |
| K41 | Ölçüm envanteri `standard_type` bazında çıkarılır; Ki/Kd ayrı modeli ancak envanter gerektirirse | Karar veriden sonra verilir, baştan değil |
| K42 | RT bölmesi tekilleştirme öncesi kayıtlardan kurulur; ilk ölçümü kesim öncesi olan çiftler testten çıkarılır | Medyan alma çiftin tarihini yok ediyor; kontrolsüz RT rastgele bölmeden iyimser çıkabilir |
| K43 | Tohum yayılımı ile konformal aralık ayrı gösterilir; üretilen şey tahmin aralığıdır | Farklı belirsizlik türleri; birleştirmek yanıltıcı |
| K44 | Gerçekleşen kapsama oranı kova başına raporlanır; kalibre edilemeyen kovada aralık gösterilmez | Konformal garanti ölçülmemiş kombinasyonlara taşınmaz |
| K45 | Hedef güvenilirlik işareti kat dışı hesaplanır | Kendi eğitim verisiyle B1 karşılaştırması sızıntıdır |
| K46 | Hastalık skoru ile aktivite tahmini birleşik skora çevrilmez; iki ayrı sıralama | Birleştirmenin bilimsel gerekçesi yok |
| K47 | Görünüm B ikinci sürüme bırakılır | İlk çalışan sürüm üç işi yapsın; kapsam genişlemesi burada frenlenir |
| K48 | Bire-çok hedef eşlemesinde hangi kaydın seçileceği kuralı önceden yazılır | Ölçülen kayıp yetmez; seçim kuralı da tanımlı olmalı |
| K49 | `max_phase = 4` kümesi "ChEMBL faz 4" olarak adlandırılır, "onaylı ilaçlar" diye genellenmez; `drug_indication` resmî ürün etiketiyle eşitlenmez | Faz 4 kümesi geri çekilmiş ilaçları da içerir; endikasyon kaydı etiketin kendisi değil |

### İkinci inceleme turundan gelen kararlar (K50–K55)

| # | Karar | Gerekçe |
| --- | --- | --- |
| K50 | D6' ile RT'nin sırası sabit: önce RT uygunluğu, sonra pozitif seçimi; sonuçlar ayrı raporlanır; bir çift tek kez sayılır | Sıra belirsizse sonuç yorumlanamaz; birleşik oran iki farklı iddiayı gizler |
| K51 | D6' iki kipte koşulur (RT-uygun ve zamandan bağımsız); minimum pozitif sayısı veriye bakılmadan sabitlenir | RT kısıtı pozitifleri tüketebilir; sonradan seçilen eşik iddiayı zayıflatır |
| K52 | Birincil metrik **Recall@%1** olarak adlandırılır; "zenginleşme" dili kullanılmaz ve metriklerin yalnız anotasyonlu pozitifleri ölçtüğü yazılır | Gerçek inaktif yok; zenginleşme dili havuzun geri kalanını negatif gösterir |
| K53 | Bölme ölçütü "ilk kamuya açık raporlama tarihi"; yayın gecikmesi ve tarihsiz kayıtların işlenişi sınırlılıklara yazılır | ChEMBL tarihi deney tarihi değil doküman tarihi |
| K54 | Hedef güvenilirlik işareti için minimum kat dışı örnek sayısı sabitlenir; altında "değerlendirme yetersiz" işareti ve belirsizlik aralığı | Sekiz örnekle hesaplanan oran güvenilirlik göstergesi değil |
| K55 | Ölçüm envanteri deneysel heterojenliği de ölçer (deney protokolü, güven skoru dağılımı), yalnız sayıyı değil | Aynı pChEMBL değeri farklı deney bağlamında aynı biyolojik anlamı taşımayabilir |

### Kararlardan koda

Bu noktadan sonra eklenecek yeni karar yok; eksik olan, mevcut kararların kodda gerçekten uygulanması. Her K kaydı bir teste ya da bir `assert`'e bağlanır — örneğin K42 için "RT test kümesinde kesim öncesi raporlanmış çift yok" testi, K45 için "güvenilirlik tablosu kat dışı tahminlerden üretildi" kontrolü. Mevcut 86 testin yanına eklenirler ve karar kaydı numarası test adında geçer.

### Açık sorular (karara dönüşmesi gereken)

- [ ] Güven skoru eşiği 8 mi 9 mu — yerel ChEMBL 37'de doğrulanır, **ertelenmez**
- [ ] Ölçüm envanteri: Ki/Kd ile IC50 aynı model hedefinde mi kalıyor
- [ ] RT kesim yılı, ve tekilleştirme öncesi kayıtlardan kurulan testte kaç çift kalıyor
- [ ] Geri alma havuzunun ve gerçek pozitiflerin tam tanımı
- [ ] Konformal aralık hangi örnekler üzerinde kalibre edilecek, hangi kovalarda geçerli sayılacak
- [ ] Bire-çok ve çok bileşenli hedef eşlemesinde seçim kuralı
- [ ] SMILES tabanlı sistemde biyolojik ilaçlar, tuzlar, karışımlar ve ana molekül eşlemesi nasıl ele alınacak
- [ ] Aşama 2 kapsamına geçilecek mi — karar, 1. aşama sonuçları geldikten sonra
- [ ] `decisions.md` sonundaki üç eski açık soru: R2'nin ele alınışı (çözüldü — kimyasal yenilik işareti), iyon kanalı yardımcı altbirimleri ve iyon kanalı ailesinin inceliği (Aşama 2'ye taşındı)

## Riskler

| Risk | Belirtisi | Önlem |
| --- | --- | --- |
| M1, M0'ı ve taban çizgilerini geçmiyor | Hedef ortalamasına yakın RMSE | Başarısızlık değil bulgu. Ürün kanıt paneli sayesinde erişim aracı olarak ayakta kalır; hedef güvenilirlik işareti bunu kullanıcıya söyler |
| Aralıklar kullanışsız geniş | %80 aralığı ±2 log birimi aşıyor | Beklenen: literatürde çıpa 1,65. Rapor edilir, ürünün dili "sıralama aracı"nda kalır |
| Open Targets eşlemesi düşük kapsamlı çıkıyor | Hedeflerin büyük kısmı eşlenemiyor | Kapsam oranı önce ölçülür (K32). Çok düşükse hastalık katmanı yalnız eşlenebilen hedeflerde gösterilir ve bu sayı arayüzde yazılı |
| Model yalnız bilineni buluyor | Yeni kimyada isabet oranı çöküyor | İddia zaten "geri bulma". Kova başına raporlama bunu gizlemiyor, ölçüyor |
| Dil gevşiyor | Bir ekranda "önerilir" yazıyor | K27 kontrol listesi; beş yüzey tek tek gezilir (arayüz, API, model kartı, tez, sunum) |
| Savunma günü demo çalışmıyor | Barındırma uyuyor, internet yok | Üç yedek: yerel Docker, 3 dakikalık kayıt, ekran görüntülü yedek slayt |
| Kapsam yine şişiyor | Aşama 2 + ablasyonlar + LOFO geri geliyor | Hastalık katmanı **eklenen** değil, aile dışı transferin **yerine geçen** şey. K28 bunu sabitliyor |
| Tez şablonu kapsamı daraltıyor | Sayfa/bölüm sınırı geliyor | Tablolar script'lerden üretiliyor; yalnız metin taşınır. Şablon gelene kadar metin yazımı bekletilir |
| 2209-A başvurusuyla uyumsuzluk | Başvuruda aile dışı transfer tasarımı var | **Kontrol edilmeli.** Teslim edilen işin başvurulan işle örtüşmesi beklenebilir; tek e-postayla çözülür, ama kapsamı değiştirmenin en pahalı yeri burası |

### En gerçekçi risk

Sondan ikincisi değil, sonuncusu ve dil olanı. Teknik riskler ölçümle yönetiliyor; idari uyum ve retorik tutarlılık ise yalnız dikkatle yönetiliyor.

## İş-zaman planı

&#91;embedded content: iş-zaman planı · 10 kalem, 1 kapı\]

Tarihler 17 Kasım başlangıcı varsayımıyla çizildi; bölüm teslim takvimi geldiğinde kayabilir. Süreler kaba ve bu makine (RTX 4050, 6 GB) içindir.

### Tek kapı

**12 Aralık — M1 kararı.** Taban çizgileri ve M0 koştuktan sonra sorulacak soru: M1, B1–B3'ü ve M0'ı geçiyor mu?

- **Geçiyorsa:** plan olduğu gibi devam eder.
- **Geçmiyorsa:** ürünün vaadi "model tahmini" yerine "kanıt paneli + erişim aracı" olarak yeniden yazılır; mimari ve kod aynı kalır, yalnız sıralamanın kaynağı ve arayüz dili değişir. Bu, tezi zayıflatmaz — ölçülmüş bir bulgu olarak yazılır.

Kararı bu tarihte ver, sonraya bırakma: uygulama geliştirmesi bu karardan sonra başlıyor.

### Paralel yürüyen kalemler

Güven katmanı, hastalık katmanı ve geri alma doğrulaması kısmen üst üste biniyor — üçü de bağımsız iş ve üçünün çıktısı da uygulama tarafına besleniyor. Hastalık katmanının ontoloji eşlemesi bir haftadan iki haftaya çıkabilir; plan bunu soğuruyor.

### Tez yazımı

Şablon henüz gelmedi. Metin yazımı şablon gelene kadar bekletilir; tablolar ve şekiller script'lerden üretildiği için şablondan bağımsız ve şimdiden kesinleştirilebilir.

## Sonraya bırakılan çalışma

Aileler arası genellenebilirlik sorusu iptal edilmedi, **donduruldu**. Bu bölüm onu altı ay sonra açacak kişi için yazıldı — muhtemelen sen.

### Dondurma işlemleri (bu hafta)

- [ ] Depoya ilk commit atılır. Şu an hiç commit yok: 1.806 + 2.497 + 881 satır kod ve 7 tablo tek diskte, sürümsüz
- [ ] Git geçmişi aynı zamanda tezin "ne zaman ne yaptım" kaydı; sonradan üretilemez
- [ ] Veri git'e girmez; `versions.md` + SHA-256 zaten yeniden üretilebilir kılıyor
- [ ] Transfer çalışmasının durduğu nokta bir etiketle işaretlenir (ör. `v0-transfer-frozen`)
- [ ] HF kontrol noktası hash'leri sabitlenir (K37) — hem ileriki iş hem donan iş için

### Devam notu — ne yazılmalı

Tek sayfa, bugün, akıldayken:

1. Nerede kaldı: 0–5. aşamalar bitti, 6. aşama (gömme çıkarımı) başlamadı
2. Hangi dosya hangi deneyi besliyor — 65 bölme dosyasının hangisi R4a, hangisi R4b
3. Sıradaki adım: gömme çıkarımı → D0 taban çizgileri → D1 rejim taraması
4. Açık üç soru ve o günkü düşünceleri
5. Bu ürün çalışmasından çıkan ve transfer çalışmasını etkileyen şeyler — örneğin güven skoru 8/9 kararı, DataSAIL karşılaştırması, zamansal bölme kodu

Altı ay sonra bu sayfa iki hafta kazandırır; olmazsa aynı bağlam sıfırdan kurulur.

### Transfer çalışması nerede devam eder

Yeni ürün çalışması transfer sorusunu **harcamıyor**: aile dışı deneyler el değmemiş kalıyor ve aynı altyapıyla sonradan koşulabiliyor. Üstelik ürün çalışması ona iki şey kazandırıyor: kalibre edilmiş bir belirsizlik katmanı ve DataSAIL ile doğrulanmış bölmeler. Yüksek lisans ya da makale ölçeğinde devam etmesi için uygun bir çalışma.

## Kaynaklar

**Değerlendirme ve kaçak**

- [Machine learning for drug-target interaction prediction: a comprehensive review](https://pmc.ncbi.nlm.nih.gov/articles/PMC12818121/) — Comput Struct Biotechnol J, 2026. Üç baskın sorun ve bölme taksonomisi
- [Improving the Assessment of Deep Learning Models in the Context of DTI Prediction](https://www.biorxiv.org/content/10.1101/2022.04.20.488898v1.full) — ilacın/hedefin tek başına öngörü gücü; B2/B3'ün gerekçesi
- [Beyond Random Splits](https://arxiv.org/html/2610.03456) — bölme şiddeti denetimi, en yakın eğitim komşusu ölçümleri
- [Auditing widely used biomolecular benchmarks](https://www.researchgate.net/publication/410687302_Auditing_widely_used_biomolecular_benchmarks_reveals_systematic_data_inconsistencies) — DTI veri setlerinde protein örtüşmesi
- [Data splitting to avoid information leakage with DataSAIL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11978981/) — Nat Commun, 2025 · [dokümantasyon](https://datasail.readthedocs.io/en/latest/workflow/splits.html)

**Belirsizlik ve uygulanabilirlik alanı**

- [Dynamic applicability domain (dAD)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10457664/) — Bioinformatics, 2023. Bileşik–hedef komşuluğunda yerel konformal tahmin
- [Conformal Regression for QSAR Modelling](https://discovery-pp.ucl.ac.uk/id/eprint/10060460/1/Svensson_etal_r1.pdf) — %80 güvende \~1,65 pIC50 aralık çıpası

**Yön ve hastalık katmanı**

- [DTIAM](https://arxiv.org/pdf/2312.15252) — etkileşim + afinite + aktivasyon/inhibisyon mekanizması, tek çerçeve
- [Open Targets Platform](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11701534/) — NAR, 2024. Hedef–hastalık ilişkileri ve sekiz kaynak için etki yönü değerlendirmesi · [kullanım sınırları ve lisans](https://casrai.org/guides/open-targets)

**Veri ve kapsam**

- [Papyrus](https://pmc.ncbi.nlm.nih.gov/articles/PMC9824924/) — J Cheminform, 2023. \~60M veri noktası, standardize, ChEMBL sınıf seviyeleri
- [ChEMBL hedef soruları](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/target-questions) — protein sınıflandırması ve bağlanma alanı tahmini
- [ChEMBL güven skoru ölçeği](https://targetmine.mizuguchilab.org/?p=722) — ikincil kaynak, doğrulanacak
- [ChEMBL lisans açıklaması](https://chembl.blogspot.com/2012/07/public-popular-chemistry-databases-and.html) — CC BY-SA, share-alike

**Ürün karşılaştırması**

- [SwissTargetPrediction](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4086140/) — NAR, 2014. Benzerlik tabanlı hedef tahmini
- [RepurposeDrugs](https://academic.oup.com/bib/article/25/4/bbae328/7709763) — Brief Bioinform, 2024. İlaç → endikasyon
- [OTRec](https://www.biorxiv.org/content/10.64898/2025.12.21.695803v2.full) — hastalık ↔ hedef, iki kollu öneri sistemi

**Proje içi**

- `docs/00_project_design_tr_rev2.md` — veri hattı bölümleri aynen geçerli
- `docs/decisions.md` — K1–K26
- `docs/01_durum_ve_kapsam.md` — 15 Kasım 2026 durum ölçümü
- `docs/01_chembl_inventory.md` … `07_leakage_audit.md` — tezde yeri olan yedi tablo
