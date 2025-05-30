# Ulanom

**A Versatile FITS Viewer and Spectrum Analysis Tool for Astronomical Data**

![Ulanom Logo Placeholder](https://github.com/Uzaysalyakamoz/Ulanom/blob/main/Ulanom.png)

## About

Ulanom is a comprehensive Python application designed for easily viewing, calibrating, and analyzing one-dimensional (1D) spectra from astronomical FITS images. It offers powerful tools to help you extract meaningful scientific results from spectroscopic observations. Our goal is to process celestial data, acting as a "messenger" (Ulagut), to transform it into valuable insights that unveil the universe's secrets.

## Features

* **FITS Image Viewing:** Fast and efficient opening and display of FITS (Flexible Image Transport System) images.
* **Advanced Scaling:** Various stretch options (Linear, Log, Sqrt, Squared, HistEq) and limit algorithms (Min/Max, Percentile, ZScale, ZMax) to adjust image contrast and brightness. Custom ZScale parameter settings are available.
* **Basic Calibration:** Support for loading Bias, Dark, and Flat frames to calibrate science data. Dark scaling based on exposure time.
* **Aperture Definition and Trace Finding:** Interactive region selection, pixel-based trace finding, and polynomial fitting for defining the object's trace.
* **Spectrum Extraction:** Extraction of 1D spectra from 2D images within defined aperture boundaries.
* **Background Subtraction:** Background modeling and subtraction from defined regions during spectrum extraction.
* **Wavelength Calibration:**
    * Loading and extraction of arc lamp spectra.
    * Reference line list support.
    * Interactive pixel-to-wavelength mapping point definition.
    * Generation of a dispersion solution through polynomial fitting.
    * **New:** Calculation of wavelength errors from the fit covariance matrix and display on plots.
    * Saving and loading of wavelength calibration solutions as JSON files.
* **Spectrum Normalization:**
    * Continuum fitting and subtraction.
    * **New:** Polynomial or Spline fit options (requires SciPy).
    * Robust fitting by excluding outliers using sigma clipping.
* **Spectral Analysis Tools:**
    * **New:** Equivalent Width (EW) calculation for both absorption and emission lines.
    * **New:** Gaussian line fitting to determine line center, amplitude, standard deviation, and FWHM of spectral lines.
* **Visualization:**
    * Interactive plots for the main FITS image and spectra (science/arc) in separate windows.
    * **New:** Plotting the spatial profile of a selected X-column on the image.
    * Histogram plot to examine pixel value distribution.
* **Data Saving:**
    * Saving normalized spectra (including pixel/wavelength, normalized flux, error, original flux, and continuum) in TXT/ASCII format.
    * Saving plots of the main image, science spectrum, and arc spectrum in popular formats like PNG/PDF.
* **User Interface:** Tkinter-based, user-friendly, and intuitive interface.

## Setup

To get Ulanom up and running on your local machine, follow these steps:

### Prerequisites

* Python 3.8+
* Git

### Steps

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/Uzaysalyakamoz/Ulanom](https://github.com/Uzaysalyakamoz/Ulanom)
    cd Ulanom
    ```


2.  **Create and Activate a Virtual Environment:**
    A virtual environment isolates the project's dependencies from other Python projects on your system.
    ```bash
    python -m venv venv
    ```
    * **On Windows:**
        ```bash
        .\venv\Scripts\activate
        ```
    * **On macOS/Linux:**
        ```bash
        source venv/bin/activate
        ```

3.  **Install Dependencies:**
    Install all required Python libraries from the `requirements.txt` file.
    ```bash
    pip install -r requirements.txt
    ```

    *Note: The `scipy` library is required for the Spline option in continuum fitting and for Gaussian line fitting. If you encounter an error during installation or do not need these specific features, you can remove the `scipy` line from `requirements.txt`. However, installation is recommended for full functionality.*

## 👤 Developer

Emre Bilgin  
📧 emre.bilgin64@gmail.com  
🌐 [GitHub](https://github.com/Uzaysalyakamoz)

This software is developed as an open-source contribution to the astronomy community.

## Usage

With your virtual environment activated, you can run the application from the root directory:

```bash
python ulanom.py

