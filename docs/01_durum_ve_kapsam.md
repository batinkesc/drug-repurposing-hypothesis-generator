# Çalışmanın Durumu ve Kapsam Değerlendirmesi

**Tarih:** 15 Kasım 2026
**Amaç:** Şu ana kadar ne yapıldığını tek yerde toplamak ve kapsamın bitirme
projesi için uygun olup olmadığını değerlendirebilmen için gereken sayıları
vermek.

---

## 1. Kısa cevap: haklısın, ama şişme sandığın yerde değil

Çalışma fazla detaylı geliyorsa sebebi büyük ölçüde **henüz yapılmamış** olan
deney planı — yapılmış olan kısım değil.

Şu ana kadar bitenler bir *altyapı*: ChEMBL'den veri çekme, molekül
standardizasyonu, veri setinin kurulması, bölmeler ve kaçak denetimi. Bunların
hepsi, projenin **en küçük** sürümü için de gerekli. Yani kapsamı şimdi
daraltırsan yapılmış işin neredeyse hiçbiri çöpe gitmez.

Şişme şurada: tasarım belgesinde **7 deney (D0–D6)**, **5+ ablasyon**, **3 RMSE
varyantı**, **2 karıştırma kontrolü**, **5 aile**, **5 tohum** ve bir de
**geri alma protokolü** var. Bunların çarpımı bitirme projesi ölçeğini aşıyor.
Kesilecek yer burası ve kesmenin maliyeti bugün çok düşük.

---

## 2. Ne bitti

| Aşama | Durum | Çıktı |
| --- | --- | --- |
| 0 · Ortam ve veri kaynağı | ✅ | ChEMBL 37 (29 GB), sürüm + SHA-256 kayıtlı |
| 1 · ChEMBL envanteri ve aile eşlemesi | ✅ | 5 aile, 1.078 hedef, sıfır aile çakışması |
| 1b · Protein temsili kararı | ✅ | Bağlanma alanı, ölçümle gerekçelendirildi |
| 2 · Veri çekme ve süzme | ✅ | 1.080.040 aktivite, düşürme tablosu |
| 3 · Molekül standardizasyonu | ✅ | 487.220 benzersiz molekül (67 dk, 15 çekirdek) |
| 4 · Veri seti kurulumu | ✅ | **283.353 hedef–bileşik çifti** |
| 5 · Bölmeler ve kaçak denetimi | ✅ | 13 rejim × 5 tohum = 65 bölme dosyası |
| 6 · Gömme çıkarımı | ⬜ | — |
| 7 · Taban çizgileri ve modeller | ⬜ | — |
| 8 · Deneyler ve teşhisler | ⬜ | — |
| 9 · Tez yazımı | ⬜ | — |

**Kod:** 1.806 satır kütüphane (`src/dti/`), 2.497 satır script, 881 satır test.
**86 test geçiyor.**

**Veri:** 29 GB ham ChEMBL, 478 MB ara dosya, 12 MB veri seti, 87 MB bölmeler.

### Veri setinin özeti

| Aile | Çift | Hedef | Bileşik | pChEMBL ortanca |
| --- | --- | --- | --- | --- |
| GPCR | 122.142 | 175 | 79.880 | 7,04 |
| Kinaz | 85.163 | 235 | 55.287 | 6,99 |
| Proteaz | 39.639 | 89 | 28.304 | 6,80 |
| İyon kanalı | 19.545 | 35 | 18.756 | 6,20 |
| Nükleer reseptör | 16.864 | 31 | 12.391 | 6,75 |
| **Toplam** | **283.353** | **565** | **190.339** | 6,92 |

---

## 3. Tezde doğrudan yeri olan tablolar (7 tane, hepsi hazır)

| Dosya | Ne gösteriyor | Tezde nerede |
| --- | --- | --- |
| `01_chembl_inventory.md` | Veri kaynağının yapısı, aile eşlemesi | Yöntem |
| `02_representation_decision.md` | Tam dizi vs bağlanma alanı ölçümü | Yöntem (gerekçe) |
| `03_drop_table.md` | Her süzme adımının kaç kayıt düşürdüğü | Yöntem (şeffaflık) |
| `04_dataset_profile.md` | Aile başına dağılımlar, pChEMBL yayılımı | Veri bölümü |
| `05_drop_table_stage2.md` | İkinci aşama düşürmeleri | Veri bölümü |
| `06_splits.md` | Bölme profili, gerçekleşen oranlar | Yöntem |
| `07_leakage_audit.md` | **Kaçak denetimi** | Bulgular |

Ayrıca `docs/decisions.md` içinde 26 numaralı karar kaydı var (K1–K26); her biri
"neden böyle yaptık" sorusunun cevabı ve savunmada doğrudan kullanılabilir.

---

## 4. Şimdiye kadarki asıl bulgular

Bunlar model eğitilmeden elde edildi ve tezde **bulgu** olarak yer alabilir:

1. **R2 ("görülmemiş bileşik") iddia ettiği kadar zor değil.** Test
   bileşiklerinin **%81,7'si** eğitim kümesinde ≥0,40 Tanimoto'da bir komşuya
   sahip, ortanca benzerlik 0,548. Beş tohumda 0,805–0,828 — tesadüf değil.
   Literatürde "zor bölme" diye sunulan kümeleme tabanlı ayırmanın bu ölçüyle
   neredeyse hiç zor olmadığını söyleyebiliyoruz.
2. **Aile dışı test verisinin yarısı kimyasal olarak tanıdık.** R4b süzgeci
   aile başına kayıtların %40–56'sını atıyor. Yani "aile dışı transfer"
   ölçümleri, kimyasal yenilik şartı konmazsa büyük ölçüde bilinen ligandların
   yeniden tanınması oluyor.
3. **Aileler arası pChEMBL ortalama farkı 0,72 log birim**, ama aile *içi*
   standart sapma 1,21–1,41. Yani ortalama-kaydırma düzeltmesi R4'teki hatayı
   en fazla 0,72 birim azaltabilir — transfer başarısızlığı kaymayla
   açıklanamaz.
4. **Beş aile arasında sıfır hedef çakışması** — R4'ün anlamlı olmasının ön
   koşulu, ölçüldü.

---

## 5. Kalan iş ve kaba süre tahmini

Tahminler bu makine (RTX 4050, 6 GB) içindir ve **kabadır**.

| İş | Tahmini süre | Zorunlu mu? |
| --- | --- | --- |
| Gömme çıkarımı (ESM-2 565 alan + ChemBERTa 190k molekül) | 2–3 saat | **Evet** |
| B1/B2/B3 taban çizgileri × 13 rejim × 5 tohum | 1–2 saat | **Evet** |
| M0 (ECFP4 + LightGBM) × 13 rejim × 5 tohum | 3–5 saat | **Evet** |
| M1 (donuk gömmeler + sığ ağ) × 13 rejim × 5 tohum | 4–8 saat | **Evet** |
| M2 (kısmi ince ayar) | 15–30 saat | Hayır |
| D2b transfer matrisi (5×5 ikili) × 5 tohum | 10–20 saat | Hayır |
| D3 karıştırma kontrolü (2 varyant) | 2–3 saat | Biri evet |
| D4 seçicilik testi | 3–5 saat | Hayır |
| D5 hedef ön yanlılığı | 1 saat | Ucuz, kalsın |
| D6 geri alma protokolü | 8–15 saat | Hayır |
| Ablasyonlar (5+) | her biri 2–8 saat | Hayır |
| Tez yazımı (TR + EN) | 3–5 hafta | **Evet** |

---

## 6. Kapsam senaryoları

### A · Minimal (bitirme projesi için yeterli)

- **3 aile** (kinaz, GPCR, proteaz) — veri hazır, aile çıkarmak sadece süzme
- Rejimler: **R1, R2, R3, R4a, R4b** (hepsi hazır, ek maliyet yok)
- Modeller: **B1–B3, M0, M1** (M2 yok)
- Deneyler: **D0, D1, D2a, D3 (aile içi karıştırma), D5**
- Tohum: **3**
- RMSE varyantı: **2** (ham + ortalama-kaydırmalı)
- Ablasyon: **yok** (ya da sadece 1)

**Hesap yükü:** ~10–15 saat. **Tez:** tam bir hikâye anlatır, savunulabilir.

### B · Orta (zamanın varsa)

A + **5 aile** + D2b transfer matrisi + 2 ablasyon + 5 tohum.
**Hesap yükü:** ~35–50 saat.

### C · Tam (tasarım belgesindeki hâli)

Hepsi. **Hesap yükü:** 80–150 saat + yazım. Bu, bitirme projesi değil, bir
yüksek lisans tezi ya da makale ölçeği.

---

## 7. Neyi kesersen ne kaybedersin

### Kesilebilir — hikâye bozulmaz

| Kesilen | Kaybedilen |
| --- | --- |
| **M2 (kısmi ince ayar)** | "Donuk gömme yeterli mi?" sorusu cevapsız kalır. M1 ile M0 farkı zaten bu sorunun büyük kısmını cevaplıyor. |
| **D6 geri alma** | Pratik kullanım senaryosu gösterilemez. Güzel bir ek, zorunlu değil. |
| **D4 seçicilik testi** | Ligand yanlılığı sorusu B2 taban çizgisiyle kısmen cevaplanıyor zaten. |
| **D2b transfer matrisi** | "Hangi aile hangisinden öğreniliyor" detayı gider; D2a (LOFO) ana soruyu yine cevaplıyor. |
| **Ablasyonlar** | Eşik seçimlerinin sonuca duyarlılığı gösterilemez. Tek bir cümleyle "ölçülmedi" denip sınırlılıklara yazılabilir. |
| **İki aile (iyon kanalı + nükleer reseptör)** | En ince iki aile. Çıkarılırsa R4 beş değil üç aile üzerinden döner; aileler arası iddia zayıflar ama ayakta kalır. |
| **5 → 3 tohum** | Hata çubukları biraz genişler. Kabul edilebilir. |

### Kesilmemeli — kesersen tez savunulamaz hâle gelir

| Kesilmemeli | Neden |
| --- | --- |
| **B1–B3 taban çizgileri** | Bunlar olmadan hiçbir skor yorumlanamaz. "RMSE 1,2" tek başına bir şey söylemez; hedef ortalamasının RMSE'si 1,3 ise model neredeyse hiçbir şey öğrenmemiştir. **En ucuz, en kritik kısım.** |
| **Kaçak denetimi** | Zaten yapıldı ve projenin en güçlü kozu. R2 hakkındaki %81,7 bulgusu buradan geliyor. |
| **R4a + R4b ayrımı** | Projenin araştırma sorusu bu. R4b olmadan "aile dışı transfer" ölçümü yarısı tanıdık ligandlardan oluşur. |
| **D3 (en az aile içi varyant)** | Modelin proteini gerçekten okuyup okumadığını söyleyen tek test. Kesilirse "model sadece aileyi öğrenmiş olabilir" itirazına cevabın olmaz. |
| **Test setinin bir kez açılması** | Metodolojik dürüstlük. Bedava. |

---

## 8. Önerim

**Senaryo A ile başla, bitince karar ver.** Gerekçe:

1. Senaryo A'nın çıktısı tam bir tez. Eksik olan şeyler "sınırlılıklar"
   bölümünde bir paragrafla dürüstçe kapatılabilir.
2. A bittikten sonra B'ye geçmek için hiçbir şeyi yeniden yapmak gerekmiyor —
   aynı bölmeler, aynı gömmeler, sadece daha fazla koşu.
3. Bölüm şablonu **henüz gelmedi**. Şablon sayfa/bölüm sınırı getirirse kapsamı
   zaten daraltmak gerekebilir; o zaman A'da olmak avantaj.
4. İki aileyi (iyon kanalı, nükleer reseptör) şimdi atmamanı öneririm: veri
   hazır ve onları kullanmak ek hesap maliyeti dışında bir iş çıkarmıyor. Eğer
   süre sıkışırsa son anda süzmek bir satırlık değişiklik.

Bu bir öneri; kapsam kararı senin.

---

## 9. Bu çalışmanın "fazla detaylı" kısmı aslında neydi

Dürüst olmak gerekirse detayın bir kısmı benim eklemelerimden geldi (tasarım
belgesine E1–E8 olarak işlendi): üç RMSE varyantı, iki karıştırma kontrolü,
ek risk satırları, aşama-0 kontrol noktası. Bunların gerekçeleri sağlam ama
hepsi **bitirme projesi için zorunlu değil.**

Buna karşılık, detay gibi görünüp aslında zorunlu olan iki şey var ve ikisi de
ölçümle ortaya çıktı:

- **Kaçak denetimi.** Olmasa R2'yi "zor bölme" diye sunacaktık. %81,7 bulgusu
  bunu imkânsız kılıyor. Bu bir süs değil, sonucun kendisi.
- **Tekrarlı alan kopyaları meselesi (K20).** Kendi belirsizlik kontrolümün
  kaçırdığı bir vakaydı; düzeltilmese model bir iyon kanalının gözeneğinin
  dörtte birini görüyor olacaktı. Bedeli iyon kanalı ailesinin 50'den 35
  hedefe inmesi oldu.

Yani: planın büyüklüğü abartılı, ama yapılan işin titizliği abartılı değil —
o titizlik iki kez gerçek hata yakaladı.

---

## 10. Durum notları

- **Git:** depo başlatıldı ama **hiçbir commit yok.** Bir commit atılması
  iyi olur (geri dönülebilirlik için), ama senin onayını bekliyorum.
- **Tez şablonu:** henüz verilmedi. Geldiğinde yalnız metin taşınacak, tablolar
  script'lerden yeniden üretilebilir.
- **Model checkpoint sürümleri:** Hugging Face commit hash'leri henüz
  sabitlenmedi (`config/config.yaml` içinde `null`). Gömme çıkarımından önce
  sabitlenmeli.
- **Açık sorular:** `docs/decisions.md` sonundaki listede 3 madde var; biri
  R2'nin nasıl ele alınacağı (bkz. §4.1), biri iyon kanalı yardımcı
  altbirimleri, biri iyon kanalı ailesinin inceliği.
