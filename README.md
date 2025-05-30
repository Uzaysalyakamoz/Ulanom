# Ulanom

**Göksel Veriyi Bilgiye Dönüştüren Çok Yönlü FITS Görüntüleyici ve Spektrum Analiz Aracı**

![Ulanom Logo Tutucu](https://github.com/Uzaysalyakamoz/Ulanom/blob/main/Ulanom.png)


## Hakkında

Ulanom, astronomik FITS görüntülerini kolayca görüntülemek, kalibre etmek ve tek boyutlu (1D) spektrumları analiz etmek için tasarlanmış kapsamlı bir Python uygulamasıdır. Spektroskopik gözlemlerden anlamlı bilimsel sonuçlar çıkarmanızı sağlayacak güçlü araçlar sunar. Amacımız, göksel verileri "haberci" (Ulagut) rolünde işleyerek, evrenin sırlarını ortaya çıkaran bilgilere dönüştürmektir.

## Özellikler

* **FITS Görüntüleme:** Hızlı ve etkili FITS (Flexible Image Transport System) görüntü açma ve görüntüleme.
* **Gelişmiş Ölçekleme:** Görüntü kontrastını ve parlaklığını ayarlamak için Linear, Log, Sqrt, Squared, HistEq gibi çeşitli renk dağılımı (stretch) seçenekleri ve Min/Max, Percentile, ZScale, ZMax gibi limit algoritmaları. Özel ZScale parametre ayarları.
* **Temel Kalibrasyon:** Bias, Dark ve Flat çerçevelerini yükleyerek bilim verisini kalibre etme desteği. Poz süresine göre dark ölçeklemesi.
* **Açıklık Tanımlama ve İz Sürme:** Nesnenin izini belirlemek için interaktif bölge seçimi, piksel tabanlı iz bulma ve polinom fit desteği.
* **Spektrum Ekstraksiyonu:** Tanımlanan açıklık sınırları içinde 2D görüntüden 1D spektrum çıkarma.
* **Arkaplan Çıkarma:** Spektrum ekstraksiyonu sırasında tanımlanan bölgelerden arkaplan modellemesi ve çıkarma.
* **Dalgaboyu Kalibrasyonu:**
    * Ark lambası spektrumları yükleme ve ekstraksiyonu.
    * Referans çizgi listesi desteği.
    * Etkileşimli piksel-dalgaboyu eşleşme noktası belirleme.
    * Polinom fit ile dispersiyon çözümü oluşturma.
    * **Yeni:** Fit kovaryans matrisinden dalgaboyu hatası hesaplama ve grafikte gösterme.
    * Dalgaboyu kalibrasyon çözümünü JSON dosyası olarak kaydetme ve yükleme.
* **Spektrum Normalizasyonu:**
    * Süreklilik (continuum) fit ve çıkarma.
    * **Yeni:** Polinom veya Spline (SciPy gerektirir) fit seçenekleri.
    * Sigma clipping ile aykırı değerleri dışarıda bırakarak daha sağlam fitler.
* **Spektral Analiz Araçları:**
    * **Yeni:** Eşdeğer Genişlik (Equivalent Width - EW) hesaplama (absorpsiyon ve emisyon çizgileri için).
    * **Yeni:** Gaussian çizgi fiti ile spektral çizgilerin merkezini, genliğini, standart sapmasını ve FWHM'ini belirleme.
* **Görselleştirme:**
    * Ana FITS görüntüsü ve spektrumları (bilim/ark) için ayrı pencerede interaktif grafikler.
    * **Yeni:** Görüntü üzerinde seçilen bir X-kolonunun uzaysal profilini çizme.
    * Histogram grafiği ile piksel değer dağılımını inceleme.
* **Veri Kaydetme:**
    * Normalize edilmiş spektrumu (piksel/dalgaboyu, normal akı, hata, orijinal akı, süreklilik dahil) TXT/ASCII formatında kaydetme.
    * Ana görüntü, bilim spektrumu ve ark spektrumu grafiklerini PNG/PDF gibi popüler formatlarda kaydetme.
* **Kullanıcı Arayüzü:** Tkinter tabanlı, kullanıcı dostu ve sezgisel arayüz.

## Kurulum

Ulanom'u yerel makinenizde çalıştırmak için aşağıdaki adımları izleyin:

### Ön Gereksinimler

* Python 3.8+
* Git

### Adımlar

1.  **Depoyu Klonlayın:**
    ```bash
    git clone [https://github.com/Uzaysalyakamoz/Ulanom/](https://github.com/Uzaysalyakamoz/Ulanom/)
    cd Ulanom
    ```

2.  **Sanal Ortam Oluşturun ve Aktifleştirin:**
    Sanal ortam, projenin bağımlılıklarını sisteminizdeki diğer Python projelerinden izole eder.
    ```bash
    python -m venv venv
    ```
    * **Windows:**
        ```bash
        .\venv\Scripts\activate
        ```
    * **macOS/Linux:**
        ```bash
        source venv/bin/activate
        ```

3.  **Bağımlılıkları Yükleyin:**
    Gerekli tüm Python kütüphanelerini `requirements.txt` dosyasından yükleyin.
    ```bash
    pip install -r requirements.txt
    ```

    *Not: `scipy` kütüphanesi, süreklilik fitinde Spline seçeneği ve Gaussian çizgi fiti için gereklidir. Kurulum sırasında bir hata alırsanız veya bu özelliklere ihtiyacınız yoksa, `requirements.txt` dosyasından `scipy` satırını silebilirsiniz. Ancak tam işlevsellik için yüklenmesi tavsiye edilir.*


## 👤 Geliştirici

Emre Bilgin  
📧 emre.bilgin64@gmail.com  
🌐 [GitHub](https://github.com/Uzaysalyakamoz)

Bu yazılım, astronomi topluluğuna katkı sağlamak amacıyla açık kaynak olarak geliştirilmiştir.
## Kullanım

Sanal ortamınız aktifken, uygulamayı ana dizinden çalıştırabilirsiniz:

```bash
python ulanom.py
