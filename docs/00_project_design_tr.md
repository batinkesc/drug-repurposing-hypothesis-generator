# Proje Tasarımı: Aileler Arası Genellenebilir İlaç–Hedef Etkileşimi Modeli

Oct 16, 2026 · @secret secret

Amaç, tek bir protein ailesine bağlı kalmayan, eğitimde hiç görülmemiş bir aileye aktarılabilen bir ilaç–hedef etkileşimi modeli kurmak ve bu aktarımın nerede çöktüğünü dürüst biçimde ölçmek.

## Araştırma sorusu

**Önceden eğitilmiş dil modelleriyle kurulan bir ilaç–hedef etkileşimi modeli, eğitimde hiç görmediği bir protein ailesine aktarılabilir mi — ve aktarılamıyorsa bu nerede, ne kadar kırılıyor?**

Aileler arası aktarımın kendisi literatürde ele alınmıştır: Dalkıran ve ark. (2023) beş kaynak aileden hedef aileye transferi ölçmüş ve veri kıtlığı çekilen hedeflerde transferin avantajlı olduğunu göstermiştir. Ancak oradaki kurgu hedef ailede bir miktar veriyle ince ayar yapmayı varsayar ve klasik özniteliklere dayanır. Burada sorulan, önceden eğitilmiş dil modeli temsillerinin hedef aileden hiçbir denetimli etiket görmeden aynı sınırı aşıp aşmadığıdır — ve bunun sızıntı denetimi ile tanısal testlerle birlikte ölçülmesidir. Çalışmaların çoğu ise ya tek bir aile içinde kalıyor ya da aileleri karıştırıp rastgele bölme yapıyor; ikinci durumda test kümesindeki bir proteinin akrabası eğitimde bulunduğu için model ailenin genel imzasını ezberleyerek yüksek skor alabiliyor.

Sorunun pratik karşılığı şudur: aynı aileden gelen hastalıklar için aynı aileden bir ilacın etkileşimini çıkarabilmek, ancak model ailenin kendisini bir girdi olarak değil, protein dizisinin taşıdığı bilgiyi genel bir temsil olarak kullanıyorsa mümkün.

Beklenen üç sonuçtan hangisi çıkarsa çıksın raporlanabilir bir bulgu:

1. Transfer korunuyor — dil modeli temsilinin aile-ötesi bilgi taşıdığının kanıtı.
2. Transfer kısmen çöküyor — hangi aile çiftleri arasında aktarımın çalıştığı haritası çıkar.
3. Transfer tamamen çöküyor — alandaki "genel DTI modeli" iddialarına karşı somut bir uyarı.

## Tasarım ilkesi: neden iki kollu mimari

Hedeflenen genelleme biçimi, mimariyi doğrudan belirliyor. Hedef başına ayrı sınıflandırıcı kuran yaklaşımlar (DEEPScreen bu şekilde 704 hedef için ayrı model eğitti) yeni bir protein geldiğinde hiçbir şey söyleyemez — o protein için eğitilmiş model yoktur. Aynı şekilde yalnız SMILES okuyan bir model de hedefi hiç görmediği için aileler arası bir soru soramaz.

Bu yüzden protein **girdi olarak** verilmek zorunda. İki kollu (two-tower) kurgu tam olarak bunu yapar:

&#91;embedded content: iki kollu mimari · 2 girdi, 2 gömme, 1 çıktı\]

Protein kolu dizinin kendisini okuduğu için, eğitimde hiç görülmemiş bir hedef için de bir temsil üretilebilir. Modelin o hedefte başarılı olup olmayacağı ayrı bir soru — zaten ölçmek istediğimiz şey bu.

## Kapsam: beş protein ailesi

Seçilen beş dal, alanın klasik hedef sınıflandırmasıyla örtüşüyor — bu kategoriler DTI literatüründe uzun süredir ayrı değerlendirme kümeleri olarak kullanılıyor, yani seçim keyfi değil ve karşılaştırma zemini hazır.

| Aile | Ne yapar | İlaç örneği | Veri bolluğu | Rolü |
| --- | --- | --- | --- | --- |
| Kinazlar | Fosfat transferiyle hücre içi sinyali yönetir | imatinib | Çok yüksek | Birincil eğitim ailesi |
| GPCR'ler | Hücre dışı sinyali hücre içine taşır | propranolol | Yüksek | Transfer hedefi |
| Proteazlar | Protein keser | ritonavir | Orta | Transfer hedefi |
| Nükleer reseptörler | Gen ifadesini doğrudan düzenler | tamoksifen | Düşük–orta | Zor transfer testi |
| İyon kanalları | Zar boyunca iyon akışını denetler | amlodipin | Düşük | En zor transfer testi |

Kinazlar birincil eğitim ailesi çünkü veri en bol orada ve aile içi dizi benzerliği yüksek — bu da "görülmemiş hedef" testini aile içinde bile anlamlı kılıyor. Diğer dördü sırayla dışarıda bırakılıp transfer ölçülecek.

Önemli bir tasarım kararı: bu aileler **eğitimde bir etiket olarak kullanılmayacak.** Model aile bilgisini girdi olarak almaz; aile yalnızca değerlendirme bölmelerini kurmak için kullanılır. Aksi halde modelin öğrendiği şey etkileşim değil, aile sınıflandırması olur.

Ölçek beklentisi: bu beş aile birlikte ChEMBL'de on binlerce hedef–bileşik ölçümü anlamına gelir. Bu, tek GPU'da gömme çıkarıp sığ bir başlık eğitmek için fazlasıyla yeterli, sıfırdan model eğitmek için ise yetersiz — mimarinin donuk gömme üzerine kurulmasının pratik sebebi de bu.

## Veri hattı

**Ortam sürümleri baştan sabitlenir.** RDKit'in kendi geriye dönük uyumsuz değişiklikler belgesi, kanonik atom sıralamasının ve üretilen kanonik SMILES'ın sürümler arasında farklılaşabildiğini açıkça belirtiyor; tautomer numaralandırma kodunda da sürümler arası değişiklikler yapılmış. Kanonik SMILES bu projede yalnız bir gösterim değil, tekilleştirmenin ve R4b kümelemesinin temelidir — sürüm değişirse bölmeler de değişir ve sonuçlar tekrarlanamaz hale gelir. Bu yüzden RDKit, transformers, PyTorch sürümleri ve kullanılan model kontrol noktalarının revizyon kimlikleri bir `environment.yml` dosyasında sabitlenir ve tezde birebir yazılır.

Hat tekrarlanabilirlik ilkesiyle kurulacak: her aşama kendi çıktısını diske yazar, tekrar çalıştırıldığında aynı sonucu verir, ve her kararın gerekçesi kodda durur.

1. **Çekme.** ChEMBL'den seçilen beş ailenin hedefleri ve bunlara ait biyoaktivite kayıtları. Hedefler UniProt kimliğiyle eşlenir, protein dizileri UniProt'tan alınır.
2. **Süzme.** Yalnız tek proteinli insan hedefleri; ölçüm türü Ki, Kd veya IC50; standart birim nM; ilişki işareti kesin (`=`), yani "büyüktür/küçüktür" kayıtları dışarıda; veri doğrulama bayrağı sorunlu olanlar atılır.
3. **Normalizasyon.** Afinite pChEMBL ölçeğine çevrilir (negatif logaritma), böylece Ki/Kd/IC50 tek eksende karşılaştırılabilir. Aynı hedef–bileşik çiftinin birden fazla ölçümü varsa medyan alınır; ölçümler arası fark büyükse (örneğin 1 log birimden fazla) çift işaretlenir ve bir ablasyonda dışarıda bırakılır.
4. **Molekül temizliği.** RDKit ile SMILES kanonikleştirme, tuz/karşı iyon ayıklama, yük nötrleştirme. Ayrıştırılamayan moleküller atılır. Kanonik SMILES üzerinden mükerrer kayıtlar birleştirilir.
5. **Protein temizliği.** Çok kısa veya çok uzun diziler elenir (ESM-2 bağlam sınırı nedeniyle); aynı dizinin farklı hedef kaydı olarak geçtiği durumlar birleştirilir.
6. **Profil çıktısı.** Aile başına hedef sayısı, hedef başına bileşik sayısı, pChEMBL dağılımı, ölçüm türü dağılımı. Bölme stratejisi bu profile bakılarak kesinleştirilecek.

### Protein temsili: tam dizi değil, bağlanma alanı

ESM-2'nin sert bir sınırı var: en fazla 1.022 amino asit (iki özel tokenle birlikte 1.024 token). Bu sınır bizim için yalnızca teknik bir engel değil, **çalışmanın ana sorusunu tehdit eden bir karıştırıcı**: aileler ortalama dizi uzunluğu bakımından birbirinden çok farklı. Kinazların katalitik alanı birkaç yüz kalıntıyken, voltaj kapılı iyon kanalları iki bini aşabiliyor. Tam diziyi kırpmak, tam da karşılaştırmak istediğimiz aileleri farklı oranlarda vurur; aileler arası performans farkı transferden değil kırpmadan geliyor olabilir.

Çözüm, ChEMBL'in kendi ek açıklamasını kullanmak: **`predicted_binding_domains` tablosu**, küçük molekülün bağlandığı Pfam-A alanını tahmin ediyor ve yaklaşık 650 bin aktivite ölçümü için mevcut. Arkasındaki PPDMs çalışması, doğrulamalarda bağlanmayı doğru alana vakaların yaklaşık %90'ında eşlediğini bildiriyor.

Tam dizi yerine bağlanma alanını gömmek üç sorunu birden çözüyor:

- Uzunluklar aileler arasında benzer hale gelir, kırpma yanlılığı ortadan kalkar
- Neredeyse tüm diziler ESM-2 sınırının altına iner
- Biyolojik olarak doğru birimdir: model, ligandın hiç değmediği yüzlerce kalıntıyla uğraşmaz

**Bedeli ve karşı önlem.** PPDMs eşlemesi aktivitelerin yarısı civarını kapsıyor ve ChEMBL'deki insan hedeflerinin yaklaşık yarısı birden fazla Pfam alanına sahip; bazı bağlanmalar iki alanın arayüzünde gerçekleşiyor. Yani alan temelli temsil veriyi küçültür ve bir kısmını belirsiz bırakır. Bu yüzden:

1. **Birincil temsil:** bağlanma alanı (eşlemesi olan kayıtlar)
2. **Ablasyon:** tam dizi (1.022 sınırında kırpılmış) — iki temsilin karşılaştırılması başlı başına bir bulgu
3. Hangi yol seçilirse seçilsin, **aile başına dizi/alan uzunluğu dağılımı ve kaç hedefin kırpıldığı raporlanır**. Bu tablo olmadan aileler arası karşılaştırma yorumlanamaz.

Alan eşlemesi olmayan hedefler için tam dizi kullanılır ve bu kayıtlar işaretlenir; bir ablasyonda dışarıda bırakılarak eşleme kapsamının sonucu etkileyip etkilemediği ölçülür.

Her aşamanın kaç kayıt düşürdüğü bir tabloya yazılır. Bu tablo tezde yer alır — veri hattının şeffaflığı, sonuçların güvenilirliğinin yarısı.

## Bölme rejimleri

Projenin asıl deneysel kurgusu burada. Aynı veriden dört farklı bölme üretilir ve aynı model her birinde ayrı ayrı değerlendirilir. Aradaki fark, sonucun kendisidir.

| Rejim | Test kümesinde ne yeni | Neyi ölçer | Beklenen zorluk |
| --- | --- | --- | --- |
| R1 — Rastgele | Hiçbir şey | Literatürle karşılaştırma zemini | En kolay |
| R2 — Görülmemiş ilaç | Bileşik | Yeni moleküle genelleme | Orta |
| R3 — Görülmemiş hedef | Protein (aile içi) | Yeni hedefe genelleme | Zor |
| R4 — Aile-dışı transfer | Tüm aile | Aileler arası aktarım | En zor |

Uygulama ayrıntıları:

- **R2** bileşik kümeleme ile kurulur: ECFP parmak izi üzerinden Tanimoto benzerliğine göre kümelenip kümeler bölmelere dağıtılır. Sadece kimliğe göre ayırmak yetmez — analog moleküller iki tarafa düşerse sızıntı sürer.
- **R3** protein dizi kimliğine göre kümeleme ile kurulur (MMseqs2, örneğin %40 eşik). Eşik değeri bir ablasyon parametresi.
- **R4** en basit kuralla kurulur ve en sert testtir: bir aile tamamen dışarıda bırakılır, kalan dördüyle eğitilir. Beş aile için beş kez tekrarlanır (leave-one-family-out). Çıktı, hangi ailenin hangi ailelerden öğrenilebildiğini gösteren bir transfer matrisidir.

**R4 tek başına yetmez: aile yeniliği ile kimyasal yenilik ayrılmalı.** Bir aileyi dışarıda bırakmak, test kümesindeki moleküllerin de yeni olduğu anlamına gelmez — aynı bileşik ya da yakın analogları başka ailelerin kayıtlarında bulunabilir. Bu durumda model yeni aileye geçerken tanıdık ligandlardan yararlanıyor olabilir ve ölçtüğümüz şey saf aile transferi olmaz. Bu yüzden R4, ligand yeniliğine göre ikiye ayrılır:

| Alt rejim | Protein ailesi | Ligand | Ne ölçer |
| --- | --- | --- | --- |
| R4a | Görülmemiş | Eğitimde var ya da yakın analog | Protein tarafının tek başına aktarımı |
| R4b | Görülmemiş | Kimyasal olarak da yeni (küme dışı) | Gerçek çift yenilik |

R4a ile R4b arasındaki fark, modelin aile transferinin ne kadarının ligand tanıdıklığından geldiğini doğrudan verir. İkisinin skorları ayrı raporlanır; tek bir R4 sayısı bu ayrımı gizler.

R4b'nin kimyasal ayrım kuralı açıkça sabitlenir, aksi halde "yeni ligand" tanımı muğlak kalır: moleküller kanonik SMILES üzerinden tekilleştirilir, ECFP4 parmak izi üzerinden Butina kümelemesiyle 0,4 Tanimoto eşiğinde kümelenir, ve kümeler bölmelere bütün olarak dağıtılır — bir küme asla iki tarafa bölünmez. Tuz ve tautomer varyantları veri hattında aynı moleküle indirgendiği için farklı kayıt biçimleri aynı kümeye düşer. Eşik bir ablasyon parametresidir.

Kümeleme tek başına yetmez. Butina kümelemesi küme merkezlerini eşikle ayırır, ancak farklı kümelerdeki her molekül çiftinin eşik altında kaldığını garanti etmez — kenar moleküller komşu kümelerin üyelerine eşiğin üstünde benzer olabilir. Bu yüzden bölmeden sonra açık bir süzme adımı uygulanır: eğitim kümesindeki herhangi bir moleküle 0,4'ün üzerinde Tanimoto benzerliği gösteren test molekülleri R4b'den çıkarılır. Süzmeden sonra kalan en yüksek eğitim–test benzerliği ve çıkarılan molekül sayısı, destekleyici istatistik olarak değil **birincil sonuç** olarak raporlanır: R4b'nin ne kadar katı bir test olduğu bu iki rakamla ölçülür.

Tekilleştirme, veri hattındaki şu adımlara dayanır: tuz ve karşı iyon ayıklama, yük nötrleştirme, ve RDKit `MolStandardize` modülünün `TautomerEnumerator` bileşeniyle kanonik tautomer seçimi. Kanonik SMILES bu üç işlemden sonra üretilir, böylece aynı molekülün farklı kayıt biçimleri tek bir girdiye iner.

**Doğrulama kümesi de aynı kurala uyar.** R4'te doğrulama kümesi eğitim ailelerinden seçilir — dışarıdaki aileden bir doğrulama kümesi kullanmak, o aileyi dolaylı olarak eğitime sokar ve sonucu geçersiz kılar.

Her rejim için sızıntı denetimi raporlanır: bölmeler arası özdeş molekül sayısı, en yüksek Tanimoto benzerliği, en yüksek dizi kimliği. Bu rakamlar olmadan hiçbir skor yorumlanamaz.

## Modeller

Üç model karşılaştırılır. Amaç en iyi skoru bulmak değil, **avantajın hangi rejimde kaybolduğunu** görmek.

Bunların altında, parametresiz ya da tek kollu üç **alt sınır modeli** yer alır. Bunlar olmadan hiçbir skor yorumlanamaz:

- **B1 — Hedef ortalaması.** Her hedef için eğitim kümesindeki ortalama pChEMBL değerini tahmin eder, molekülü hiç görmez. DTI verisinde bu model şaşırtıcı derecede iyi skor alır; geçilemiyorsa model etkileşim hakkında hiçbir şey öğrenmemiştir. Aile-dışı rejimde hedefe özgü ortalama bulunmayacağı için eğitim ailelerinin genel ortalamasına düşer — bu da o rejimin gerçek alt sınırıdır.
- **B2 — Yalnız ligand.** Protein girdisi yok; yalnız molekülden tahmin. Moleküllerin genel "etkinlik eğilimi" ne kadar açıklıyor?
- **B3 — Yalnız protein.** Molekül girdisi yok. Hedefin genel ölçüm düzeyi ne kadar açıklıyor?

Her alt sınır modeli, her bölme rejiminde ayrı kurulur ve yalnız o rejimin izin verdiği eğitim verisini kullanır. B1'i tüm veriden hesaplayıp aile-dışı rejimde kullanmak sızıntıdır.

**M0 — Temel yöntem.** ECFP4 parmak izi (molekül) + amino asit bileşimi ve k-mer sayımları (protein), birleştirilip gradyan artırma ile regresyon. Literatürdeki aktivite uçurumu bulguları göz önüne alındığında bu model atılacak bir karşılaştırma değil; dil modellerini geçmesi gerçek bir olasılık.

**M1 — Dil modeli gömmeleri, donuk.** ChemBERTa-2 ve ESM-2'den gömmeler bir kez çıkarılıp diske yazılır; üzerine sığ bir birleştirme ağı eğitilir (birleştirme → iki gizli katman → tek çıktı). Birleştirme yöntemi bir ablasyon: basit yan yana ekleme, eleman bazlı çarpım, ve ikisinin birlikte kullanımı.

**M2 — Kısmi ince ayar.** M1'in üstüne, ChemBERTa-2'nin son birkaç katmanının açılması. ESM-2 donuk kalır (tek GPU'da ince ayarı pratik değil). Bu modelin M1'i geçip geçmediği, ince ayarın maliyetinin değip değmediğini söyler.

Eğitim kurgusu — her rejimde aynı, aksi halde karşılaştırma geçersiz:

- Kayıp: pChEMBL üzerinde ortalama kare hata
- Hiperparametreler: eğitim kümesinin içinden ayrılan doğrulama kümesiyle seçilir, test kümesine asla bakılmaz
- Her yapılandırma birkaç farklı rastgele tohumla tekrarlanır; raporlanan değer ortalama ve yayılımdır, tek koşu değil
- Test kümesi yalnız bir kez açılır, tüm modeller ve rejimler hazır olduğunda

Son madde pazarlık konusu değil: test kümesinin tek seferde açılacağının baştan yazılı olması, sonuçların güvenilirliğinin teminatıdır.

## Değerlendirme

**Metrikler.** Birincil: RMSE ve Pearson korelasyonu (pChEMBL üzerinde). İkincil: Spearman ve uyum indeksi (CI) — sıralama kalitesini ölçer, yeniden konumlandırma senaryosunda asıl önemli olan budur. Hedef başına ayrı ayrı da raporlanır, çünkü toplam skor birkaç veri-zengini hedefin arkasına saklanabilir.

**Yanlılık testleri.** İki tanesi literatürden alınıyor ve yeni mimariye uygulanmamış olmaları başlı başına katkı:

- *Seçicilik testi* (etiket ters çevirmenin regresyona uyarlanmış hali) — aynı molekülün iki farklı hedefte ölçülmüş, aralarında büyük aktivite farkı bulunan kayıtları ayrılır (başlangıç eşiği 2 log birim; eşik keyfi olduğundan 1,5 ve 3 log birimle de tekrarlanır ve sonucun eşiğe duyarlılığı raporlanır). Model bu farkı yakalayabiliyorsa gerçek etkileşim özelliği öğrenmiştir; iki hedef için benzer değer veriyorsa yalnız ligandın kendisine bakıyordur. Not: özgün etiket ters çevirme deneyi ikili sınıflandırma için tasarlanmış olup açıkça tanımlanmış inaktif kayıtlar gerektirir; nicel veri hattında kayıt yokluğu inaktiflik sayılamayacağı için test bu biçimde uyarlanmıştır.
- *Hedef ön yanlılığı* — eğitimde bir hedefin etiket dağılımı ne kadar dengesizse, modelin o hedefte çoğunluk etiketini tahmin etme eğilimi ölçülür. Tahminler ile hedefin eğitim ortalaması arasındaki ilişki, yanlılığın doğrudan göstergesi.
- *Karıştırma kontrolü* — protein gömmesi rastgele başka bir proteininkiyle değiştirildiğinde performans ne kadar düşüyor? Düşmüyorsa model proteini hiç kullanmıyordur. Bu, iki kollu mimarinin gerçekten iki kollu çalıştığının en basit kanıtı.

Seçicilik testiyle ilgili bir sınır: iki hedefteki ölçümler farklı deney türlerinden geliyorsa (biri Ki, diğeri IC50) aradaki fark biyolojik seçicilikten değil ölçüm yönteminden kaynaklanıyor olabilir. Bu yüzden test, yalnız aynı ölçüm türünden gelen çiftler üzerinde yürütülür. Ancak aynı ölçüm türü gerekli koşuldur, yeterli değildir: aynı türden iki kayıt bile farklı deney protokolleri ya da protein yapıları nedeniyle tam karşılaştırılabilir olmayabilir. Mümkün olduğunca deney koşulu uyumu da kontrol edilir (ChEMBL deney kimliği ve güven skoru üzerinden); sağlanamadığı durumda sonuç, biyolojik seçiciliğin kesin kanıtı olarak değil \*\*ölçüm türü eşleştirilmiş seçicilik testi\*\* olarak yorumlanır ve metinde bu adla anılır.

**Geri alma protokolü.** Yeniden konumlandırma iddiasını sınanabilir hale getiren kısım: ChEMBL'de klinik faz bilgisiyle işaretlenmiş onaylı ilaçların bilinen hedef ilişkilerinden bir bölümü eğitimden çıkarılır. Model, o ilacı o hedefin aday listesinde ne kadar yukarı taşıyor? Ölçüler: ilk %1'lik dilimdeki zenginleşme, ortalama karşılıklı sıra (MRR) ve NDCG@k. Karşılaştırma hem rastgele sıralamaya hem de M0'a karşı yapılır. Değerlendirme havuzunun hangi adaylardan oluştuğu açıkça tanımlanır; havuzun bileşimi bu metrikleri doğrudan belirler.

Bu protokol R4 ile birleştirildiğinde projenin en güçlü sonucunu verir: *dışarıda bırakılmış bir aileden, hiç görülmemiş bir hedef için, bilinen bir ilacı bulabiliyor muyuz?*

**İddianın sınırı: ön eğitim düzeyinde sıfır-atış değildir.** ESM-2 UniProt üzerinde, ChemBERTa-2 PubChem üzerinde ön eğitilmiştir. Dolayısıyla "dışarıda bırakılan" ailenin protein dizileri ve o ailenin ligandlarına benzer moleküller, ön eğitim verisinde büyük olasılıkla zaten bulunmaktadır. Bu çalışmanın sıfır-atış iddiası **denetimli öğrenme düzeyindedir**: model, hedef aileye ait hiçbir etkileşim etiketi görmemiştir. Ön eğitimin sağladığı temsil bilgisi mevcuttur ve zaten yöntemin dayanağıdır.

Bu ayrım sonucu geçersiz kılmaz ama yorumunu belirler: transfer başarılı olursa bu, "model aileyi hiç bilmiyordu" değil, "ön eğitimli temsil, etiketsiz olarak da aile-ötesi bilgi taşıyor" anlamına gelir. Sonuçlar bu çerçevede raporlanacaktır.

## Deney planı

| # | Deney | Kurgu | Cevapladığı soru |
| --- | --- | --- | --- |
| D0 | Alt sınır modelleri | B1, B2, B3 × tüm rejimler | Skorların gerçek tabanı nerede? |
| D1 | Rejim taraması | M0, M1, M2 × R1, R2, R3, R4a, R4b | Avantaj hangi rejimde kayboluyor? |
| D2a | Aile aktarımı (LOFO) | M1, leave-one-family-out | Bir aile, kalan dördünün birleşiminden öğrenilebiliyor mu? |
| D2b | Transfer matrisi (ikili) | M1, tek aileden eğit → tek aileye test | Hangi aile hangisinden öğreniliyor? |
| D3 | Karıştırma kontrolü | M1, protein gömmesi karıştırılmış | Model proteini gerçekten kullanıyor mu? |
| D4 | Seçicilik testi | M0, M1 | Gerçek etkileşim mi, ligand yanlılığı mı? |
| D5 | Hedef ön yanlılığı | M0, M1 | Tahminler hedef ortalamasına mı çakılıyor? |
| D6 | Geri alma | M1, R3 ve R4 | Bilinen ilaç–hedef çiftleri bulunabiliyor mu? |

Ablasyonlar (zaman kalırsa, öncelik sırasıyla):

1. Dizi kimliği eşiği (R3'te %30 / %40 / %50) — "görülmemiş hedef" tanımı sonucu ne kadar değiştiriyor?
2. Birleştirme yöntemi (ekleme / çarpım / ikisi birden)
3. Ölçüm türü karışımı — yalnız Ki/Kd ile eğitim, IC50 dahil edildiğinde ne değişiyor?
4. Tutarsız ölçümlerin (1 log üstü yayılım) dışarıda bırakılması
5. ESM-2 model boyutu (35M vs 650M) — küçük model yetiyor mu?

D0, D1, D2a ve D2b çekirdektir. D2a ile D2b aynı soruyu sormaz: LOFO bir ailenin öğrenilebilir olup olmadığını söyler, ikili matris ise bu bilginin hangi aileden geldiğini. Donuk gömmelerle ikili matrisin yirmi ek koşusu ucuzdur, çünkü dil modelleri yeniden çalıştırılmaz. D0 ucuzdur ve önce yapılır: alt sınırlar bilinmeden D1'in skorları yorumlanamaz. D3 de ucuz ve zorunludur, mimarinin gerçekten iki kollu çalıştığını kanıtlar. D4–D6, çalışmayı yayımlanabilir hale getiren kısımdır.

## Riskler

| Risk | Belirtisi | Önlem |
| --- | --- | --- |
| Nadir ailelerde veri yetersiz | İyon kanalları / nükleer reseptörlerde hedef başına birkaç bileşik | Eşik koy (hedef başına minimum ölçüm sayısı); yetmezse aileyi yalnız test tarafında kullan |
| R4'te performans tamamen çöküyor | Korelasyon sıfıra yakın | Bu bir başarısızlık değil, bulgudur — raporla ve nedenini analiz et |
| Gömme çıkarma yavaş | Uzun protein dizilerinde GPU belleği yetmiyor | Gömmeleri bir kez çıkarıp diske yaz; uzun dizileri parçalayıp ortalama al; gerekirse ESM-2'nin küçük sürümü |
| Aile tanımı bulanık | Bir hedef birden fazla aileye giriyor | Tek bir kaynaktan (ChEMBL hedef sınıflandırması) tek atama; belirsizler dışarıda |
| M0 her yerde kazanıyor | Dil modeli hiçbir rejimde öne geçmiyor | Beklenen sonuçlardan biri; literatürdeki aktivite uçurumu bulgularıyla tutarlı, tezin tezi bu olur |
| Kapsam şişiyor | Beş aile + altı deney + beş ablasyon yetişmiyor | D1–D3 çekirdek, gerisi opsiyonel; aile sayısı üçe indirilebilir |

En gerçekçi risk sonuncusu. Beş aile hedefi iyi bir ufuk ama ilk turda üç aileyle (kinaz, GPCR, proteaz) başlayıp genişletmek daha güvenli.

## İş-zaman planı

Hedef: dönem sonunda tez.

| Aşama | Süre | Çıktı |
| --- | --- | --- |
| 1. Veri hattı | 2–3 hafta | Temizlenmiş veri seti + profil tablosu + düşürme tablosu |
| 2. Bölmeler ve sızıntı denetimi | 1 hafta | R1, R2, R3, R4a, R4b bölmeleri + benzerlik raporu |
| 3. Gömme çıkarma | 1 hafta | ChemBERTa-2 ve ESM-2 gömmeleri diskte |
| 4. Alt sınırlar (D0) ve M0 | 1 hafta | B1–B3 ve M0'ın tüm rejimlerdeki skorları |
| 5. M1 ve M2 | 2–3 hafta | D1 tamamlanmış |
| 6. Transfer matrisi (D2) | 2 hafta | Aileler arası aktarım haritası |
| 7. Teşhis ve yanlılık (D3–D5) | 2 hafta | Teşhis sonuçları |
| 8. M2 ince ayar ve geri alma (D6) | 1–2 hafta | Yeniden konumlandırma sonuçları |
| 9. Yazım | 3 hafta | Tez |

İki kontrol noktası:

- **1. aşama sonunda** — veri profili beklenenden küçük çıkarsa kapsamı üç aileye indir. Bu kararı erken ver, sonra değil.
- **4. aşama sonunda** — M0 alt sınır modellerini geçiyor ve tüm rejimlerde çalışıyorsa hat sağlamdır. Buraya kadar olan kısım tek başına bir çalışma gövdesi oluşturur; gerisi üstüne eklenir.
