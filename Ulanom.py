
# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox, scrolledtext, simpledialog, Listbox, Scrollbar
import os
import numpy as np
from functools import partial
import warnings
import json

# Astropy modülleri
from astropy.io import fits
from astropy.visualization import (
    ManualInterval, PercentileInterval, MinMaxInterval, ZScaleInterval,
    LinearStretch, LogStretch, SqrtStretch, SquaredStretch, PowerStretch, HistEqStretch,
    ImageNormalize
)
from astropy.stats import histogram as astropy_histogram, sigma_clip
from astropy.utils.exceptions import AstropyDeprecationWarning

# Matplotlib modülleri
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import RectangleSelector, SpanSelector
import matplotlib.patches as patches
import matplotlib.pyplot as plt

# SciPy importları (opsiyonel)
try:
    from scipy.interpolate import splrep, splev
    from scipy.optimize import curve_fit
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    # Program başında uyarı verilecek

# Uyarıları Yönet
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
warnings.simplefilter('ignore', category=RuntimeWarning)
warnings.simplefilter('ignore', category=np.RankWarning)

# --- Sabitler ve Varsayılanlar ---
DEFAULT_ZSCALE_CONTRAST = 0.25
DEFAULT_ZSCALE_NSAMPLES = 600
DEFAULT_WLCALIB_DEGREE = 3
DEFAULT_CONTINUUM_DEGREE = 3
PIXEL_MATCH_TOLERANCE = 0.5

# === Gaussian Fonksiyonu ===
def gaussian(x, amplitude, mean, stddev): # Çizgi merkezleme için basit Gaussian
    return amplitude * np.exp(-0.5 * ((x - mean) / stddev)**2)

def gaussian_with_offset(x, amplitude, mean, stddev, offset): # Çizgi fiti için offset'li Gaussian
    return offset + amplitude * np.exp(-0.5 * ((x - mean) / stddev)**2)

# === Dialog Sınıfları ===
class ZScaleDialog(tk.Toplevel):
    """ZScale algoritma parametrelerini ayarlamak için diyalog kutusu (Basitleştirilmiş)."""
    def __init__(self, parent, app_instance):
        super().__init__(parent); self.transient(parent); self.title("ZScale Parametreleri"); self.app = app_instance; self.geometry("320x180"); self.resizable(False, False)
        frame = ttk.Frame(self, padding="10"); frame.pack(expand=True, fill=tk.BOTH)
        params = {"contrast": ("Kontrast:", tk.DoubleVar(value=self.app.zscale_contrast)),
                  "n_samples": ("Örnek Piksel Sayısı (n_samples):", tk.IntVar(value=self.app.zscale_nsamples)),}
        self.vars = {}; self.entries = {}
        for i, (key, (label_text, var)) in enumerate(params.items()):
            self.vars[key] = var; ttk.Label(frame, text=label_text).grid(row=i, column=0, sticky=tk.W, pady=5, padx=5)
            entry = ttk.Entry(frame, textvariable=var, width=12); entry.grid(row=i, column=1, sticky=(tk.W, tk.E), pady=5, padx=5); self.entries[key] = entry
        frame.columnconfigure(1, weight=1); frame.rowconfigure(len(params), weight=1)
        button_frame = ttk.Frame(frame); button_frame.grid(row=len(params), column=0, columnspan=2, pady=15, sticky=(tk.S, tk.E))
        ttk.Button(button_frame, text="Uygula", command=self.apply).pack(side=tk.LEFT, padx=5); ttk.Button(button_frame, text="Varsayılan", command=self.set_defaults).pack(side=tk.LEFT, padx=5); ttk.Button(button_frame, text="Kapat", command=self.destroy).pack(side=tk.LEFT, padx=5)
        self.protocol("WM_DELETE_WINDOW", self.destroy); self.grab_set(); self.wait_window(self)
    def set_defaults(self): self.vars["contrast"].set(DEFAULT_ZSCALE_CONTRAST); self.vars["n_samples"].set(DEFAULT_ZSCALE_NSAMPLES)
    def apply(self):
        try:
            new_params = {}
            for key, var in self.vars.items():
                new_params[key] = var.get(); entry_widget = self.entries[key]; label = entry_widget.master.grid_slaves(row=entry_widget.grid_info()['row'], column=0)[0].cget('text').replace(":", "").strip()
                if key == "contrast" and (not isinstance(new_params[key], (float, int)) or new_params[key] <= 0): raise ValueError(f"'{label}' pozitif sayı olmalı.")
                if key == "n_samples" and (not isinstance(new_params[key], int) or new_params[key] <= 0): raise ValueError(f"'{label}' pozitif tamsayı olmalı.")
            self.app.zscale_contrast = new_params["contrast"]; self.app.zscale_nsamples = new_params["n_samples"]; print("ZScale parametreleri güncellendi.");
            if self.app.current_limit_type == 'zscale': print("ZScale aktif, görüntü güncelleniyor..."); self.app.apply_scaling()
            self.destroy()
        except ValueError as e: messagebox.showerror("Geçersiz Değer", str(e), parent=self)
        except Exception as e: messagebox.showerror("Hata", f"Parametreler uygulanırken hata: {e}", parent=self)

class ScaleDialog(tk.Toplevel):
    def __init__(self, parent, app_instance):
        super().__init__(parent); self.transient(parent); self.title("Ölçek Parametreleri"); self.app = app_instance; self.geometry("400x350"); self.resizable(False, False); frame = ttk.Frame(self, padding="10"); frame.pack(expand=True, fill=tk.BOTH); stretch_frame = ttk.LabelFrame(frame, text="Renk Dağılımı (Stretch)"); stretch_frame.pack(pady=5, fill=tk.X); self.stretch_var = tk.StringVar(value=self.app.current_stretch_type); stretches = ['linear', 'log', 'sqrt', 'squared', 'histeq'];
        self.science_file_path = None # Dosya yolunu saklamak için eklenebilir
        for s_type in stretches: ttk.Radiobutton(stretch_frame, text=s_type.capitalize(), variable=self.stretch_var, value=s_type).pack(anchor=tk.W, padx=10)
        limit_frame = ttk.LabelFrame(frame, text="Limit Algoritması"); limit_frame.pack(pady=5, fill=tk.X); initial_limit_type = self.app.current_limit_type; initial_limit_value = f"percentile_{self.app.current_percentile}" if initial_limit_type == 'percentile' else initial_limit_type; self.limit_var = tk.StringVar(value=initial_limit_value); self.percentile_var = tk.DoubleVar(value=self.app.current_percentile)
        initial_vmin = getattr(self.app, 'manual_vmin', 0.0); initial_vmax = getattr(self.app, 'manual_vmax', 1.0); self.vmin_var = tk.DoubleVar(value=initial_vmin); self.vmax_var = tk.DoubleVar(value=initial_vmax); limits_display = ['minmax', 'percentile', 'zscale', 'zmax', 'manual']; self.limit_entries = {}
        for i, l_display in enumerate(limits_display):
            rb_value = f"percentile_{self.app.current_percentile}" if l_display == 'percentile' else l_display; rb = ttk.Radiobutton(limit_frame, text=l_display.capitalize(), variable=self.limit_var, value=rb_value, command=self.toggle_limit_entries); rb.grid(row=i, column=0, sticky=tk.W, padx=10, pady=2)
            if l_display == 'percentile': self.limit_entries['percentile'] = ttk.Entry(limit_frame, textvariable=self.percentile_var, width=8, state=tk.DISABLED); self.limit_entries['percentile'].grid(row=i, column=1, padx=5)
            elif l_display == 'manual': f_manual = ttk.Frame(limit_frame); f_manual.grid(row=i, column=1, padx=5, sticky=tk.W); entry_vmin = ttk.Entry(f_manual, textvariable=self.vmin_var, width=8, state=tk.DISABLED); entry_vmax = ttk.Entry(f_manual, textvariable=self.vmax_var, width=8, state=tk.DISABLED); entry_vmin.pack(side=tk.LEFT); ttk.Label(f_manual, text="-").pack(side=tk.LEFT, padx=2); entry_vmax.pack(side=tk.LEFT); self.limit_entries['manual'] = (entry_vmin, entry_vmax)
        self.toggle_limit_entries(); button_frame = ttk.Frame(frame); button_frame.pack(pady=15, side=tk.BOTTOM, anchor=tk.E); ttk.Button(button_frame, text="Uygula", command=self.apply).pack(side=tk.LEFT, padx=5); ttk.Button(button_frame, text="Kapat", command=self.destroy).pack(side=tk.LEFT, padx=5); self.protocol("WM_DELETE_WINDOW", self.destroy); self.grab_set(); self.wait_window(self)
    def toggle_limit_entries(self):
        selected_limit_val = self.limit_var.get(); selected_limit_type = 'percentile' if selected_limit_val.startswith('percentile_') else selected_limit_val
        for key, entry_widget in self.limit_entries.items(): state = tk.NORMAL if key == selected_limit_type else tk.DISABLED
        if isinstance(entry_widget, tuple): entry_widget[0].config(state=state); entry_widget[1].config(state=state)
        else: entry_widget.config(state=state)
    def apply(self):
        try:
            new_stretch = self.stretch_var.get(); selected_limit_val = self.limit_var.get(); new_limit_type = 'percentile' if selected_limit_val.startswith('percentile_') else selected_limit_val; new_percentile = self.app.current_percentile; new_vmin = getattr(self.app, 'manual_vmin', 0.0); new_vmax = getattr(self.app, 'manual_vmax', 1.0)
            if new_limit_type == 'percentile': new_percentile = self.percentile_var.get();
            if not (0 < new_percentile <= 100): raise ValueError("Yüzdelik 0-100 arası olmalı.")
            elif new_limit_type == 'manual': new_vmin = self.vmin_var.get(); new_vmax = self.vmax_var.get();
            if new_vmin >= new_vmax: raise ValueError("Vmin < Vmax olmalı.")
            self.app.current_stretch_type = new_stretch; self.app.current_limit_type = new_limit_type; self.app.current_percentile = new_percentile; self.app.manual_vmin = new_vmin; self.app.manual_vmax = new_vmax; self.app._update_menu_selection(); print("Ölçek parametreleri diyalogdan güncellendi.")
            self.app.apply_scaling(); self.destroy()
        except ValueError as e: messagebox.showerror("Geçersiz Değer", str(e), parent=self)
        except Exception as e: messagebox.showerror("Hata", f"Parametreler uygulanırken hata: {e}", parent=self) 

class HistogramDialog(tk.Toplevel):
    """Görüntü verisinin histogramını gösteren diyalog kutusu."""
    def __init__(self, parent, app_instance, data, vmin, vmax):
        super().__init__(parent)
        self.app = app_instance
        self.transient(parent)
        self.title("Görüntü Histogramı")
        self.geometry("600x500")

        frame = ttk.Frame(self, padding="10")
        frame.pack(expand=True, fill=tk.BOTH)

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(1, 1, 1)

        valid_data = data[np.isfinite(data)].flatten()

        if valid_data.size > 0:
            # --- DEĞİŞİKLİK BURADA ---
            # Hatalı olan Astropy çağrısını kaldırıp, Matplotlib'in kendi 'freedman'
            # metodunu doğrudan kullanıyoruz. Bu daha basit ve doğru bir yöntemdir.
            self.ax.hist(valid_data, bins='freedman', color='gray', histtype='stepfilled', alpha=0.8, label='Tüm Pikseller')

            self.ax.axvline(vmin, color='cyan', linestyle='--', linewidth=2, label=f'Vmin: {vmin:.2f}')
            self.ax.axvline(vmax, color='magenta', linestyle='--', linewidth=2, label=f'Vmax: {vmax:.2f}')
            
            self.ax.set_title("Piksel Değer Dağılımı")
            self.ax.set_xlabel("Piksel Değeri (Counts/ADU)")
            self.ax.set_ylabel("Piksel Sayısı")
            self.ax.set_yscale('log') 
            self.ax.grid(True, linestyle=':', alpha=0.6)
            self.ax.legend()
        else:
            self.ax.text(0.5, 0.5, "Histogram çizmek için veri yok.", ha='center', va='center')

        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)

        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.grab_set()
        self.wait_window(self)
class SpectrumPlotWindow(tk.Toplevel):
    def __init__(self, parent, app_instance, x_axis_data, y_axis_data, title="1D Spektrum",
                 x_axis_label="Piksel/Dalga Boyu", x_error_data=None):
        super().__init__(parent)
        self.app = app_instance
        self.title(title)
        self.geometry("800x600")

        frame = ttk.Frame(self, padding="5")
        frame.pack(expand=True, fill=tk.BOTH)

        self.fig = Figure(figsize=(7, 5), dpi=100)
        self.ax = self.fig.add_subplot(1, 1, 1)

        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)

        self.line = None
        self.errorbar_container = None # Hata çubukları için container

        self.update_plot(x_axis_data, y_axis_data, x_axis_label=x_axis_label, x_error_data=x_error_data, title=title)

        self.ax.grid(True, linestyle='--', alpha=0.6)
        self.fig.tight_layout()
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def update_plot(self, x_axis_data, y_axis_data, x_axis_label="Piksel/Dalga Boyu",
                    y_axis_label="Akı / Counts", title=None, x_error_data=None):
        self.ax.cla() # Önceki tüm çizimleri temizle

        # Errorbar çizimi
        if x_error_data is not None and len(x_error_data) == len(x_axis_data) and np.any(np.isfinite(x_error_data)):
            # Errorbar çizimi için color, linewidth gibi özellikler eklendi
            self.errorbar_container = self.ax.errorbar(x_axis_data, y_axis_data, xerr=x_error_data, fmt='-',
                                                        elinewidth=0.7, capsize=2, label=y_axis_label.split('/')[0].strip(),
                                                        errorevery=max(1, len(x_axis_data) // 100), # Çok fazla nokta varsa her 100. noktaya hata çubuğu çiz
                                                        ecolor='gray', color='blue') # errorbar rengi gri, çizgi rengi mavi
            self.line = self.errorbar_container.lines[0] # errorbar'ın çizdiği çizgiyi line olarak ata
        else:
            self.line, = self.ax.plot(x_axis_data, y_axis_data, linewidth=1, color='blue', label=y_axis_label.split('/')[0].strip()) # Normal çizgi, mavi renk

        self.ax.set_xlabel(x_axis_label)
        self.ax.set_ylabel(y_axis_label)
        if title:
            self.ax.set_title(title)

        # Legend'ı güncelle
        handles, labels = self.ax.get_legend_handles_labels()
        # _nolegend_ etiketini filtrele
        filtered_handles_labels = [(h, l) for h, l in zip(handles, labels) if l != '_nolegend_']
        if filtered_handles_labels:
            filtered_handles, filtered_labels = zip(*filtered_handles_labels)
            self.ax.legend(filtered_handles, filtered_labels, fontsize='small')
        elif self.ax.get_legend() is not None:
            self.ax.get_legend().remove() # Legend yoksa kaldır

        self.ax.grid(True, linestyle='--', alpha=0.6) # Izgarayı her zaman açık tut
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw_idle()

    def on_close(self):
        if hasattr(self.app, 'science_plot_window') and self.app.science_plot_window is self: self.app.science_plot_window = None
        if hasattr(self.app, 'arc_plot_window') and self.app.arc_plot_window is self: self.app.arc_plot_window = None
        if hasattr(self.app, 'spatial_profile_window') and self.app.spatial_profile_window is self: self.app.spatial_profile_window = None
        self.destroy()

class WavelengthCalibWindow(tk.Toplevel):
    def __init__(self, parent, app_instance, pixels, flux, line_list=None):
        super().__init__(parent); self.app = app_instance; self.pixels = pixels; self.flux = flux
        self.line_list = np.array(line_list) if line_list is not None else np.array([])
        self.calib_points_internal = list(self.app.calib_points); self.fitted_coefficients = self.app.dispersion_coefficients
        self.fitted_covariance = self.app.dispersion_covariance
        self.fit_degree = tk.IntVar(value=getattr(self.app, 'wl_fit_degree', DEFAULT_WLCALIB_DEGREE)); self.applied = False; self.title("Dalga Boyu Kalibrasyon Aracı"); self.geometry("950x700")
        main_frame = ttk.Frame(self, padding=5); main_frame.pack(fill=tk.BOTH, expand=True); plot_frame = ttk.Frame(main_frame); plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)
        self.fig = Figure(figsize=(7, 6), dpi=100); self.ax_spec = self.fig.add_subplot(2, 1, 1); self.ax_resid = self.fig.add_subplot(2, 1, 2); self.fig.subplots_adjust(hspace=0.35)
        self.line_spec, = self.ax_spec.plot(self.pixels, self.flux, linewidth=1, label="Ark Spektrumu"); self.ax_spec.set_xlabel("Piksel (X)"); self.ax_spec.set_ylabel("Akı"); self.ax_spec.set_title("Ark Spektrumu - Nokta Eklemek İçin Çizgiye Tıklayın"); self.ax_spec.grid(True, alpha=0.5)

        # DEĞİŞİKLİK: İşaretli noktaların çizimi için placeholder
        self.marked_points_plot = None # Scatter plot objesi olacak
        # DEĞİŞİKLİK: Gaussian fit eğrisi için placeholder
        self.gaussian_fit_plot_line = None # Çizgi olarak eklenecek

        self.fit_curve_plot, = self.ax_spec.plot([], [], 'r--', linewidth=1, label="Fit Dalgaboyu (sağ eksen)", zorder=5)
        self.ax_spec_wl = self.ax_spec.twinx(); self.ax_spec_wl.set_ylabel("Fit Edilmiş Dalga Boyu", color='red'); self.ax_spec_wl.tick_params(axis='y', labelcolor='red'); self.ax_spec_wl.spines['right'].set_color('red'); self.ax_spec_wl.grid(False)
        self.ax_resid.set_xlabel("Piksel (X)"); self.ax_resid.set_ylabel("Artık (Girilen - Fit)"); self.ax_resid.set_title("Fit Artıkları"); self.ax_resid.grid(True, alpha=0.5); self.line_resid, = self.ax_resid.plot([], [], 'ro', markersize=4, label="Artıklar"); self.ax_resid.axhline(0, color='k', linestyle=':', linewidth=0.8)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame); self.canvas_widget = self.canvas.get_tk_widget(); self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True);
        self.cid_click = self.canvas.mpl_connect('button_press_event', self.on_plot_click)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame); self.toolbar.update(); self.toolbar.pack(side=tk.BOTTOM, fill=tk.X); control_frame = ttk.Frame(main_frame, width=300); control_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=5); control_frame.pack_propagate(False)
        self.auto_detected_peaks_pixels = [] # YENİ: Otomatik bulunan piklerin piksel konumları
        self.auto_peaks_plot = None          # YENİ: Otomatik bulunan pikleri gösteren çizim nesnesi
        points_frame = ttk.LabelFrame(control_frame, text="Tanımlanan Noktalar (P, W)"); points_frame.pack(pady=5, fill=tk.BOTH, expand=True); points_scrollbar_y = Scrollbar(points_frame, orient=tk.VERTICAL); points_scrollbar_x = Scrollbar(points_frame, orient=tk.HORIZONTAL); self.points_listbox = Listbox(points_frame, yscrollcommand=points_scrollbar_y.set, xscrollcommand=points_scrollbar_x.set, width=35, exportselection=False)
        points_scrollbar_y.config(command=self.points_listbox.yview); points_scrollbar_x.config(command=self.points_listbox.xview); points_scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y); points_scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X); self.points_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True); ttk.Button(control_frame, text="Seçili List Noktasını Sil", command=self.delete_selected_point).pack(pady=(5,10))
        auto_peak_frame = ttk.Frame(control_frame)
        auto_peak_frame.pack(pady=5, fill=tk.X)
        ttk.Button(auto_peak_frame, text="Pikleri Otomatik Bul", command=self.find_and_plot_arc_peaks).pack(side=tk.LEFT, padx=5)
        self.peak_prominence_var = tk.StringVar(value="10")
        ttk.Label(auto_peak_frame, text="Min. Prominence:").pack(side=tk.LEFT, padx=(10,2))
        ttk.Entry(auto_peak_frame, textvariable=self.peak_prominence_var, width=5).pack(side=tk.LEFT)
        fit_frame = ttk.LabelFrame(control_frame, text="Fit Ayarları"); fit_frame.pack(pady=5, fill=tk.X); ttk.Label(fit_frame, text="Polinom Derecesi:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W); ttk.Spinbox(fit_frame, from_=1, to=10, textvariable=self.fit_degree, width=5).grid(row=0, column=1, padx=5, pady=5, sticky=tk.W); ttk.Button(fit_frame, text="Fit Et ve Artıkları Göster", command=self.fit_dispersion).grid(row=1, column=0, columnspan=2, pady=10)
        # DEĞİŞİKLİK: RMS değeri için daha belirgin bir etiket
        self.rms_label_var = tk.StringVar(value="RMS: -")
        ttk.Label(fit_frame, textvariable=self.rms_label_var, font=("TkDefaultFont", 10, "bold")).grid(row=2, column=0, columnspan=2, pady=5) # Font kalınlaştırıldı

        action_frame = ttk.Frame(control_frame); action_frame.pack(side=tk.BOTTOM, pady=10, fill=tk.X); ttk.Button(action_frame, text="Uygula ve Kapat", command=self.apply_and_close).pack(side=tk.LEFT, expand=True, padx=5); ttk.Button(action_frame, text="İptal Et", command=self.on_close).pack(side=tk.LEFT, expand=True, padx=5)
        self.load_initial_points();
        if self.fitted_coefficients is not None: self.plot_fit_and_residuals()
        self.protocol("WM_DELETE_WINDOW", self.on_close); self.grab_set()

    def find_and_plot_arc_peaks(self):
        if not SCIPY_AVAILABLE:
            messagebox.showwarning("SciPy Eksik",
                                   "Otomatik pik bulma için SciPy kütüphanesi (scipy.signal.find_peaks) gereklidir.\n"
                                   "Lütfen 'pip install scipy' ile yükleyin.", parent=self)
            return

        from scipy.signal import find_peaks

        try:
            prominence_val = float(self.peak_prominence_var.get())
            if prominence_val <= 0:
                messagebox.showerror("Geçersiz Değer", "Minimum Prominence pozitif bir sayı olmalıdır.", parent=self)
                return
        except ValueError:
            messagebox.showerror("Geçersiz Değer", "Minimum Prominence sayısal bir değer olmalıdır.", parent=self)
            return

        valid_flux_indices = ~np.isnan(self.flux)
        if not np.any(valid_flux_indices):
            messagebox.showinfo("Veri Yok", "Pik aramak için geçerli akı verisi bulunamadı.", parent=self)
            return

        flux_to_search = self.flux[valid_flux_indices]
        pixels_to_search = self.pixels[valid_flux_indices]

        peak_indices_in_search, properties = find_peaks(flux_to_search,
                                                         prominence=prominence_val,
                                                         distance=5) # İki pik arası en az 5 piksel olsun (ayarlanabilir)

        self.auto_detected_peaks_pixels = pixels_to_search[peak_indices_in_search]

        print(f"Otomatik olarak {len(self.auto_detected_peaks_pixels)} pik bulundu (Prominence > {prominence_val}).")

        if self.auto_peaks_plot:
            try:
                for artist in self.auto_peaks_plot:
                    artist.remove()
            except (TypeError, ValueError):
                if hasattr(self.auto_peaks_plot, 'remove'): self.auto_peaks_plot.remove()
            self.auto_peaks_plot = None

        if len(self.auto_detected_peaks_pixels) > 0:
            self.auto_peaks_plot = []
            for peak_px in self.auto_detected_peaks_pixels:
                line = self.ax_spec.axvline(peak_px, color='orange', linestyle=':', linewidth=0.8, alpha=0.7, ymin=0.05, ymax=0.95)
                self.auto_peaks_plot.append(line)

        self.canvas.draw_idle()
        if not self.auto_detected_peaks_pixels.size:
            messagebox.showinfo("Pik Bulunamadı", "Belirtilen kriterlere uygun pik bulunamadı.", parent=self)

    def load_initial_points(self): self.calib_points_internal = list(self.app.calib_points); self.update_points_listbox(); self.update_marked_points_plot()
    def update_points_listbox(self):
        self.points_listbox.delete(0, tk.END); self.calib_points_internal.sort(key=lambda p: p[0])
        for i, (pix, wl) in enumerate(self.calib_points_internal): self.points_listbox.insert(tk.END, f"{i+1}: P={pix:.2f}, W={wl:.4f}")

    # --- update_marked_points_plot metodu (DEĞİŞİKLİK) ---
    def update_marked_points_plot(self):
        if self.marked_points_plot:
            try:
                # Eğer scatter plot ise remove metodu vardır.
                # Scatter plot collections içinde olduğu için daha karmaşık temizleme gerekebilir.
                if self.marked_points_plot in self.ax_spec.collections:
                    self.marked_points_plot.remove()
                else: # Eğer yanlışlıkla line objesi falan atanmışsa
                    if hasattr(self.marked_points_plot, 'remove'):
                        self.marked_points_plot.remove()
            except ValueError:
                pass
            self.marked_points_plot = None

        if self.calib_points_internal:
            pixels_calib = np.array([p[0] for p in self.calib_points_internal])
            # Dalga boyu değerlerini de alalım, y ekseninde göstermek için
            fluxes_at_calib_points = []

            if hasattr(self, 'pixels') and self.pixels is not None and len(self.pixels) > 1 and \
               hasattr(self, 'flux') and self.flux is not None and len(self.flux) == len(self.pixels):

                valid_flux_indices = ~np.isnan(self.flux)
                if np.any(valid_flux_indices):
                    valid_arc_pixels = self.pixels[valid_flux_indices]
                    valid_arc_flux = self.flux[valid_flux_indices]

                    if len(valid_arc_pixels) > 1:
                        # Tıklanan piksel değerlerinin, mevcut spektrumun piksel aralığı içinde olup olmadığını kontrol et
                        interp_mask = (pixels_calib >= valid_arc_pixels[0]) & (pixels_calib <= valid_arc_pixels[-1])
                        if np.any(interp_mask):
                            valid_calib_pixels_for_plot = pixels_calib[interp_mask]
                            try:
                                # Bu noktadaki akı değerlerini interpolate et
                                fluxes_at_calib = np.interp(valid_calib_pixels_for_plot, valid_arc_pixels, valid_arc_flux)
                                # DEĞİŞİKLİK: 'x' işareti yerine sadece daire (marker='o') kullanıldı
                                self.marked_points_plot = self.ax_spec.scatter(valid_calib_pixels_for_plot, fluxes_at_calib,
                                                                               c='magenta', marker='o', s=50, # Marker 'o' olarak değiştirildi
                                                                               edgecolor='black', linewidth=0.5, # Kenarlık eklendi
                                                                               label='Tanımlı Noktalar', zorder=10)
                            except Exception as interp_err:
                                print(f"Interpolation error: {interp_err}")
                                self.marked_points_plot = None
                        else:
                            self.marked_points_plot = None # Hiçbir nokta aralıkta değilse
                    else:
                        print("Warning: Not enough valid flux points for interpolation range.")
                        self.marked_points_plot = None
                else:
                    print("Uyarı: Geçerli flux verisi yok.")
                    self.marked_points_plot = None
            else:
                print("Uyarı: Geçerli piksel/flux aralığı yok.")
                self.marked_points_plot = None
        else:
            self.marked_points_plot = None # calib_points_internal boşsa

        self.update_legend()
        self.canvas.draw_idle()

    def update_legend(self):
        handles_spec, labels_spec = self.ax_spec.get_legend_handles_labels()
        handles_wl, labels_wl = self.ax_spec_wl.get_legend_handles_labels()
        all_handles = handles_spec + handles_wl
        all_labels = labels_spec + labels_wl

        unique_handles_labels = {}
        for h, l in zip(all_handles, all_labels):
            if l not in unique_handles_labels and h is not None:
                unique_handles_labels[l] = h

        # --- DEĞİŞİKLİK BAŞLANGICI: Hatalı lejant güncelleme mantığı düzeltildi ---
        # Eğer "Tanımlı Noktalar" etiketi varsa, onu sayı ile güncelleyip
        # sözlükte (dictionary) doğru anahtarla değiştiriyoruz.
        if 'Tanımlı Noktalar' in unique_handles_labels:
            handle = unique_handles_labels.pop('Tanımlı Noktalar')  # Eski kaydı al ve sil
            new_label = f"Tanımlı Noktalar ({len(self.calib_points_internal)})"
            unique_handles_labels[new_label] = handle  # Yeni etiketle tekrar ekle

        # Varsa Gauss Fit çizgisini de lejanta ekle
        if self.gaussian_fit_plot_line and self.gaussian_fit_plot_line.get_label() not in unique_handles_labels:
            label = self.gaussian_fit_plot_line.get_label()
            unique_handles_labels[label] = self.gaussian_fit_plot_line
        # --- DEĞİŞİKLİK SONU ---

        current_legend = self.ax_spec.get_legend()
        if current_legend:
            current_legend.remove()

        if unique_handles_labels:
            # Artık sözlüğün anahtarları (labels) ve değerleri (handles) doğru şekilde eşleşiyor.
            self.ax_spec.legend(handles=unique_handles_labels.values(),
                                labels=unique_handles_labels.keys(),
                                fontsize='x-small', loc='upper left')

    def delete_selected_point(self):
        selected_indices = self.points_listbox.curselection()
        if not selected_indices: messagebox.showwarning("Seçim Yok", "Silmek için listeden bir nokta seçin.", parent=self); return
        was_fitted_before = self.fitted_coefficients is not None
        del self.calib_points_internal[selected_indices[0]]; self.update_points_listbox(); self.update_marked_points_plot()
        if was_fitted_before and len(self.calib_points_internal) >= self.fit_degree.get() + 1: print("Nokta silindi, otomatik yeniden fit ediliyor..."); self.fit_dispersion()
        else: self.fitted_coefficients = None; self.fitted_covariance = None; self.fit_curve_plot.set_data([], []); self.line_resid.set_data([], []); self.rms_label_var.set("RMS: - (Fit Gerekli)"); self.ax_spec_wl.set_ylim(0,1); self.canvas.draw_idle()

    # --- on_plot_click metodu (DEĞİŞİKLİK) ---
    def on_plot_click(self, event):
        if event.inaxes not in [self.ax_spec, self.ax_spec_wl]: return
        if event.button != 1: return
        if not self.toolbar.mode == '': return # Eğer pan veya zoom modu aktifse tıklamayı engelle

        clicked_pixel_rough = event.xdata
        if clicked_pixel_rough is None: return

        # Önceki Gaussian fit çizimini temizle (varsa)
        if self.gaussian_fit_plot_line:
            try:
                self.gaussian_fit_plot_line.remove()
                self.gaussian_fit_plot_line = None
                self.canvas.draw_idle()
            except ValueError:
                pass

        try:
            # En yakın pikseli bul
            nearest_idx = np.argmin(np.abs(self.pixels - clicked_pixel_rough))
            search_radius = 5 # Yakındaki pikselleri arama penceresi
            search_start = max(0, nearest_idx - search_radius)
            search_end = min(len(self.pixels), nearest_idx + search_radius + 1)

            local_pixels = self.pixels[search_start:search_end]
            local_flux = self.flux[search_start:search_end]

            valid_indices_local = np.where(np.isfinite(local_flux))[0]

            peak_pixel_refined = None # Nihai piksel konumu

            if SCIPY_AVAILABLE and valid_indices_local.size >= 4 : # Gaussian fit için en az 4 nokta gerekir (amplitude, mean, stddev, offset)
                try:
                    x_fit = local_pixels[valid_indices_local]
                    y_fit = local_flux[valid_indices_local]

                    # Eğer veri çok düzse veya negatifse Gaussian fit zor olabilir.
                    # Bu durumda basit ağırlıklı ortalama veya peak değeri kullanılmalı.
                    if np.min(y_fit) < 0 or (np.nanmax(y_fit) - np.nanmin(y_fit) < 1e-6 and np.nanmax(y_fit) < 1e-6):
                        print("[DEBUG] Gaussian fit için data uygun değil gibi (çok düz veya negatif). Ağırlıklı ortalama deneniyor.")
                        if np.sum(y_fit) > 1e-9: # Pozitif toplam varsa ağırlıklı ortalama
                            peak_pixel_refined = np.sum(x_fit * y_fit) / np.sum(y_fit)
                        else: # Yoksa en yüksek akı noktasının pikseli
                            peak_pixel_refined = x_fit[np.argmax(y_fit)]
                        popt_gaussian = None # Gaussian fit başarısız

                    else:
                        amplitude_guess = np.max(y_fit) - np.min(y_fit)
                        mean_guess = x_fit[np.argmax(y_fit)]
                        stddev_guess = (x_fit[-1] - x_fit[0]) / 6.0 # Yaklaşık FWHM / 2.355

                        # Geniş bir aralık için p0, offsetli Gaussian
                        p0_offset_gaussian = [amplitude_guess if amplitude_guess > 0 else 1.0, mean_guess, stddev_guess if stddev_guess > 0 else 1.0, np.min(y_fit)]
                        # Bounds'u daha sıkı yapabiliriz, özellikle mean için
                        bounds_offset_gaussian = (
                            [0, x_fit[0], 0.1, -np.inf], # Amplitüd pozitif, mean aralıkta, stddev pozitif, offset herhangi bir değer
                            [np.inf, x_fit[-1], search_radius, np.inf] # Amplitüd, mean, stddev max search_radius, offset herhangi bir değer
                        )

                        popt_gaussian, pcov_gaussian = curve_fit(gaussian_with_offset, x_fit, y_fit, p0=p0_offset_gaussian, bounds=bounds_offset_gaussian, maxfev=10000)
                        fitted_mean = popt_gaussian[1] # Mean (piksel konumu)

                        # Fit edilen merkezin seçilen aralıkta olup olmadığını kontrol et
                        if np.min(local_pixels) <= fitted_mean <= np.max(local_pixels):
                            peak_pixel_refined = fitted_mean
                            print(f"[DEBUG] Peak (Gaussian Fit): {peak_pixel_refined:.4f}")

                            # Fitten sonra Gauss eğrisini çiz
                            x_gaussian_plot = np.linspace(x_fit[0], x_fit[-1], 100)
                            y_gaussian_plot = gaussian_with_offset(x_gaussian_plot, *popt_gaussian)
                            self.gaussian_fit_plot_line, = self.ax_spec.plot(x_gaussian_plot, y_gaussian_plot, color='red', linestyle='-', linewidth=1.5, label='Gauss Fit')
                            self.update_legend() # Efsaneyi güncelle
                            self.canvas.draw_idle()
                        else:
                            print("[DEBUG] Gaussian fit sonucu pencere dışında, ağırlıklı ortalama kullanılıyor.")
                            # Gauss fit başarısız olursa ağırlıklı ortalama fallback
                            flux_for_centroid_positive = y_fit - np.min(y_fit) # Pozitif akı yap
                            flux_for_centroid_positive[flux_for_centroid_positive < 0] = 0 # Negatif değerleri sıfırla
                            sum_of_weights = np.sum(flux_for_centroid_positive)
                            if sum_of_weights > 1e-9:
                                peak_pixel_refined = np.sum(x_fit * flux_for_centroid_positive) / sum_of_weights
                            else:
                                peak_pixel_refined = x_fit[np.argmax(y_fit)]
                            popt_gaussian = None # Gaussian fit başarısız

                except (RuntimeError, ValueError) as fit_err:
                    print(f"[DEBUG] Gaussian fit yakınsamadı veya hata ({fit_err}), ağırlıklı ortalama kullanılıyor.")
                    # Fallback: ağırlıklı ortalama
                    flux_for_centroid_positive = local_flux[valid_indices_local] - np.min(local_flux[valid_indices_local])
                    flux_for_centroid_positive[flux_for_centroid_positive < 0] = 0
                    sum_of_weights = np.sum(flux_for_centroid_positive)
                    if sum_of_weights > 1e-9:
                        peak_pixel_refined = np.sum(local_pixels[valid_indices_local] * flux_for_centroid_positive) / sum_of_weights
                    else:
                        peak_pixel_refined = local_pixels[np.argmax(local_flux[valid_indices_local])]
                    popt_gaussian = None # Gaussian fit başarısız
            else:
                # SciPy yoksa veya yeterli nokta yoksa, ağırlıklı ortalama veya doğrudan pikseli kullan
                print("[DEBUG] SciPy yok veya yeterli nokta yok, ağırlıklı ortalama/en yüksek piksel kullanılıyor.")
                flux_for_centroid_positive = local_flux[valid_indices_local] - np.min(local_flux[valid_indices_local])
                flux_for_centroid_positive[flux_for_centroid_positive < 0] = 0
                sum_of_weights = np.sum(flux_for_centroid_positive)
                if sum_of_weights > 1e-9:
                    peak_pixel_refined = np.sum(local_pixels[valid_indices_local] * flux_for_centroid_positive) / sum_of_weights
                else:
                    peak_pixel_refined = local_pixels[np.argmax(local_flux[valid_indices_local])]

        except Exception as e:
            print(f"[DEBUG] Error during peak finding/fitting: {e}")
            messagebox.showerror("Hata", f"Pik bulunurken/fit edilirken hata oluştu: {e}", parent=self)
            return

        if peak_pixel_refined is None: # Eğer herhangi bir nedenle piksel bulunamazsa
            messagebox.showerror("Hata", "Piksel konumu belirlenemedi.", parent=self)
            return

        prompt = f"Bulunan Piksel Merkezi: {peak_pixel_refined:.4f} için dalga boyunu girin:"
        suggested_wl = None
        if self.fitted_coefficients is not None and self.line_list.size > 0:
            try:
                predicted_wl = np.poly1d(self.fitted_coefficients)(peak_pixel_refined)
                nearest_wl_index = np.argmin(np.abs(self.line_list - predicted_wl))
                suggested_wl = self.line_list[nearest_wl_index]
                prompt += f"\n(Öneri: {suggested_wl:.4f})"
            except Exception as e:
                print(f"[DEBUG] Error suggesting wavelength: {e}")

        wavelength = simpledialog.askfloat("Dalga Boyu Girişi", prompt, parent=self, initialvalue=suggested_wl if suggested_wl is not None else None)
        if wavelength is not None:
            self.add_point(peak_pixel_refined, wavelength)

    def add_point(self, pixel, wavelength):
        was_fitted_before = self.fitted_coefficients is not None
        existing_indices = [i for i, p in enumerate(self.calib_points_internal) if abs(p[0] - pixel) < PIXEL_MATCH_TOLERANCE]
        if existing_indices:
            idx_to_replace = existing_indices[0]
            if messagebox.askyesno("Piksel Yakın", f"{pixel:.4f} pikseline çok yakın ({self.calib_points_internal[idx_to_replace][0]:.4f}) nokta zaten var.\nDalga boyu ({wavelength:.4f}) ile güncellensin mi?", parent=self): self.calib_points_internal[idx_to_replace] = (pixel, wavelength)
            else: return
        else: self.calib_points_internal.append((pixel, wavelength))
        self.update_points_listbox(); self.update_marked_points_plot()
        if was_fitted_before or len(self.calib_points_internal) >= self.fit_degree.get() + 1: print("Nokta eklendi/güncellendi, otomatik yeniden fit ediliyor..."); self.fit_dispersion()
        else: self.fitted_coefficients = None; self.fitted_covariance = None; self.fit_curve_plot.set_data([], []); self.line_resid.set_data([], []); self.rms_label_var.set("RMS: - (Fit Gerekli)"); self.ax_spec_wl.set_ylim(0,1); self.canvas.draw_idle()

    # --- plot_fit_and_residuals metodu (DEĞİŞİKLİK) ---
    def plot_fit_and_residuals(self):
        if self.fitted_coefficients is None or len(self.calib_points_internal) == 0:
            self.fit_curve_plot.set_data([], [])
            self.line_resid.set_data([], [])
            self.rms_label_var.set("RMS: -")
            self.ax_spec_wl.set_ylim(0, 1) # Sağ y eksenini sıfırla
            self.update_legend()
            self.canvas.draw_idle()
            return

        pixels_calib = np.array([p[0] for p in self.calib_points_internal])
        wavelengths_calib = np.array([p[1] for p in self.calib_points_internal])
        fit_func = np.poly1d(self.fitted_coefficients)

        # Fit eğrisini çizerken tüm spektrumun x aralığını kullan
        if self.pixels is not None and len(self.pixels) > 1:
            fit_x = self.pixels
        else:
            # Eğer spektrum verisi yoksa, kalibrasyon noktalarının min/max aralığını kullan
            fit_x = np.linspace(np.min(pixels_calib), np.max(pixels_calib), 500)

        fit_y_wl = fit_func(fit_x)
        self.fit_curve_plot.set_data(fit_x, fit_y_wl)
        self.ax_spec_wl.relim()
        self.ax_spec_wl.autoscale_view()

        # Artıkları hesapla ve çiz
        predicted_wl = fit_func(pixels_calib)
        residuals = wavelengths_calib - predicted_wl
        self.line_resid.set_data(pixels_calib, residuals)
        self.ax_resid.relim()
        self.ax_resid.autoscale_view()

        # RMS değerini hesapla ve güncelle
        rms = np.sqrt(np.mean(residuals**2))
        # DEĞİŞİKLİK: RMS'i daha iyi bir formatta göster (yüzde hata payı ile birlikte)
        mean_wl_in_calib = np.nanmean(wavelengths_calib)
        if mean_wl_in_calib != 0 and np.isfinite(mean_wl_in_calib):
            rms_percent = (rms / mean_wl_in_calib) * 100
            self.rms_label_var.set(f"RMS: {rms:.5f} ({rms_percent:.3f}%)")
        else:
            self.rms_label_var.set(f"RMS: {rms:.5f}")

        self.update_legend()
        self.canvas.draw_idle()

    def fit_dispersion(self):
        degree = self.fit_degree.get()
        if len(self.calib_points_internal) < degree + 1: messagebox.showerror("Yetersiz Nokta", f"{degree}. derece fit için en az {degree+1} nokta gerekli.", parent=self); self.fitted_coefficients = None; self.fitted_covariance = None; self.plot_fit_and_residuals(); return
        pixels = np.array([p[0] for p in self.calib_points_internal]); wavelengths = np.array([p[1] for p in self.calib_points_internal])
        try:
            coeffs, cov_matrix = np.polyfit(pixels, wavelengths, degree, cov=True)
            self.fitted_coefficients = coeffs; self.fitted_covariance = cov_matrix
            print(f"Dispersiyon Fit ({degree}. d): {self.fitted_coefficients}")
            self.plot_fit_and_residuals()
        except Exception as e: messagebox.showerror("Fit Hatası", f"Fit hatası: {e}", parent=self); self.fitted_coefficients = None; self.fitted_covariance = None; self.plot_fit_and_residuals()
    def apply_and_close(self): self.applied = True; self.on_close()
    def on_close(self):
        if hasattr(self, 'cid_click') and self.cid_click: self.canvas.mpl_disconnect(self.cid_click)
        self.destroy()

class LineFitResultsDialog(tk.Toplevel):
    def __init__(self, parent, fit_params, x_region, y_region, fitted_curve_y):
        super().__init__(parent)
        self.transient(parent)
        self.title("Gaussian Çizgi Fit Sonuçları")
        self.geometry("600x650")
        amplitude, mean, stddev, offset, fwhm = fit_params
        plot_frame = ttk.Frame(self, padding="10"); plot_frame.pack(fill=tk.BOTH, expand=True)
        fig = Figure(figsize=(5.5, 4), dpi=100); ax = fig.add_subplot(111)
        ax.plot(x_region, y_region, 'b.', label='Veri', markersize=5)
        ax.plot(x_region, fitted_curve_y, 'r-', label='Gaussian Fit')
        ax.set_xlabel("Dalgaboyu / Piksel"); ax.set_ylabel("Akı"); ax.legend(); ax.grid(True)
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=plot_frame); canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True); toolbar = NavigationToolbar2Tk(canvas, plot_frame)
        toolbar.update(); toolbar.pack(fill=tk.X)
        params_frame = ttk.LabelFrame(self, text="Fit Parametreleri", padding="10"); params_frame.pack(fill=tk.X, padx=10, pady=5)
        param_labels = [("Merkez (Mean):", f"{mean:.4f}"), ("Amplitüd:", f"{amplitude:.4f}"), ("Standart Sapma (StdDev):", f"{stddev:.4f}"), ("FWHM:", f"{fwhm:.4f}"), ("Offset (Süreklilik):", f"{offset:.4f}")]
        for i, (label, value) in enumerate(param_labels):
            ttk.Label(params_frame, text=label).grid(row=i, column=0, sticky=tk.W, padx=5, pady=2)
            ttk.Label(params_frame, text=value).grid(row=i, column=1, sticky=tk.W, padx=5, pady=2)
        ttk.Button(self, text="Kapat", command=self.destroy).pack(pady=10)
        self.protocol("WM_DELETE_WINDOW", self.destroy); self.grab_set(); self.wait_window(self)

# === Ana FITS Görüntüleyici Uygulaması ===
class FitsViewerApp:
    def __init__(self, root):
        self.root = root; self.root.title("ULANOM"); self.root.geometry("1100x750")
        self.hdulist = None; self.raw_image_data = None; self.image_data = None; self.raw_arc_image_data = None; self.arc_image_data = None
        self.science_header = None
        self.bias_data = None; self.dark_data = None; self.flat_data = None; self.bias_path = None; self.dark_path = None; self.flat_path = None
        self.science_exptime = None; self.dark_exptime = None
        self.current_vmin = None; self.current_vmax = None; self.current_stretch_type = 'linear'; self.current_limit_type = 'zscale'; self.current_percentile = 99.0; self.manual_vmin = 0.0; self.manual_vmax = 1.0; self.zscale_contrast = DEFAULT_ZSCALE_CONTRAST; self.zscale_nsamples = DEFAULT_ZSCALE_NSAMPLES
        self.stretch_map = {'linear': LinearStretch(), 'log': LogStretch(), 'sqrt': SqrtStretch(), 'squared': SquaredStretch(), 'histeq': HistEqStretch }
        self.rect_selector = None; self.selection_box = None; self.selection_patch = None; self.trace_peaks_x = None; self.trace_peaks_y = None; self.trace_coefficients = None; self.aperture_width = 10; self.polynomial_degree = 3; self.aperture_lower = None; self.aperture_upper = None; self.overlay_plots = []
        self.extracted_pixels = None; self.extracted_flux = None; self.arc_extracted_pixels = None; self.arc_extracted_flux = None;
        self.line_list_data = None; self.calib_points = []; self.dispersion_coefficients = None; self.wl_fit_degree = DEFAULT_WLCALIB_DEGREE;
        self.dispersion_covariance = None; self.science_wavelength_error = None
        self.background_regions = []; self.background_patches = []; self.subtract_bg_var = tk.BooleanVar(value=False)
        self.science_continuum_flux = None; self.science_normalized_flux = None; self.continuum_fit_degree = DEFAULT_CONTINUUM_DEGREE
        self.continuum_fit_type = tk.StringVar(value="Polynomial")
        self.science_plot_window = None; self.arc_plot_window = None; self.ew_selector = None; self.line_fit_selector = None;
        self.trace_coeffs_label_var = tk.StringVar(value="Trace Fit Katsayıları: -")
        self.poly_degree_label_var = tk.StringVar(value=f"Polinom Derecesi: {self.polynomial_degree}")
        self.ap_width_label_var = tk.StringVar(value=f"Açıklık Genişliği: {self.aperture_width} px")
        self.continuum_degree_label_var = tk.StringVar(value=f"Süreklilik Fit Derecesi: {self.continuum_fit_degree}")
        self.bias_status_var = tk.StringVar(value="Bias: Yok"); self.dark_status_var = tk.StringVar(value="Dark: Yok"); self.flat_status_var = tk.StringVar(value="Flat: Yok")
        self.show_wl_error_var = tk.BooleanVar(value=False) # DEĞİŞİKLİK: Dalgaboyu hata çubuklarını göster/gizle
        self.status_bar_var = tk.StringVar()
        self.line_fit_selector = None
        self._setup_menu(); self._setup_gui_layout(); self._update_menu_selection()
        self._update_calib_status_labels(); self.status_bar_var.set(" Durum: Hazır")
        self.spatial_profile_selector = None
        self.spatial_profile_line = None
        self.spatial_profile_window = None

    # --- _setup_menu metodu (DEĞİŞİKLİK) ---
    def _setup_menu(self):
        self.menu_bar = tk.Menu(self.root); self.root.config(menu=self.menu_bar)
        self.file_menu = tk.Menu(self.menu_bar, tearoff=0); self.menu_bar.add_cascade(label="Dosya", menu=self.file_menu)
        self.file_menu.add_command(label="Bilim Görüntüsü Aç...", command=self.open_fits_file)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Normalize Spektrumu Kaydet (.txt)...", command=self.save_normalized_spectrum)
        self.save_plot_menu = tk.Menu(self.file_menu, tearoff=0)
        self.file_menu.add_cascade(label="Grafiği Farklı Kaydet", menu=self.save_plot_menu)
        self.save_plot_menu.add_command(label="Ana Görüntüyü Kaydet (.png, .pdf)...", command=self.save_main_image_plot)
        self.save_plot_menu.add_command(label="Bilim Spektrumunu Kaydet (.png, .pdf)...", command=self.save_science_spectrum_plot)
        self.save_plot_menu.add_command(label="Ark Spektrumunu Kaydet (.png, .pdf)...", command=self.save_arc_spectrum_plot)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Çıkış", command=self.exit_app)
        self.calib_menu = tk.Menu(self.menu_bar, tearoff=0); self.menu_bar.add_cascade(label="Kalibrasyon", menu=self.calib_menu); self.calib_menu.add_command(label="Bias Yükle...", command=self.load_bias_frame); self.calib_menu.add_command(label="Dark Yükle...", command=self.load_dark_frame); self.calib_menu.add_command(label="Flat Yükle...", command=self.load_flat_frame); self.calib_menu.add_separator(); self.calib_menu.add_command(label="Kalibrasyonları Tekrar Uygula", command=self.reapply_calibrations); self.calib_menu.add_command(label="Yüklü Kalibrasyonları Temizle", command=self.clear_calibrations)
        self.scale_menu = tk.Menu(self.menu_bar, tearoff=0); self.menu_bar.add_cascade(label="Ölçek", menu=self.scale_menu); self.stretch_menu = tk.Menu(self.scale_menu, tearoff=0); self.scale_menu.add_cascade(label="Renk Dağılımı", menu=self.stretch_menu); self.stretch_var = tk.StringVar(value=self.current_stretch_type); stretches_for_menu = ['linear', 'log', 'sqrt', 'squared', 'histeq']
        for s_type in stretches_for_menu: self.stretch_menu.add_radiobutton(label=s_type.capitalize(), variable=self.stretch_var, value=s_type, command=self.apply_scaling_from_menu)
        self.limits_menu = tk.Menu(self.scale_menu, tearoff=0); self.scale_menu.add_cascade(label="Limitler", menu=self.limits_menu); self.limits_var = tk.StringVar(value=self.current_limit_type); self.limits_menu.add_radiobutton(label="Min/Max (100%)", variable=self.limits_var, value='minmax', command=self.apply_scaling_from_menu); self.percentile_menu = tk.Menu(self.limits_menu, tearoff=0); self.limits_menu.add_cascade(label="Yüzdelik", menu=self.percentile_menu); percentiles = [99.9, 99.5, 99.0, 98.0, 95.0, 90.0];
        for p in percentiles: self.percentile_menu.add_radiobutton(label=f"{p}%", variable=self.limits_var, value=f'percentile_{p}', command=partial(self.set_percentile_from_menu_and_apply, p))
        self.limits_menu.add_radiobutton(label="ZScale", variable=self.limits_var, value='zscale', command=self.apply_scaling_from_menu); self.limits_menu.add_radiobutton(label="ZMax", variable=self.limits_var, value='zmax', command=self.apply_scaling_from_menu)
        self.scale_menu.add_separator(); self.scale_menu.add_command(label="Ölçek Parametreleri...", command=self.show_scale_parameters_dialog); self.scale_menu.add_command(label="ZScale Parametreleri...", command=self.show_zscale_parameters_dialog); self.scale_menu.add_separator(); self.scale_menu.add_command(label="Piksel Histogram Grafiği...", command=self.show_histogram_graph)
        self.aperture_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.menu_bar.add_cascade(label="Açıklık & İndirgeme", menu=self.aperture_menu); self.aperture_menu.add_command(label="Ap: 1. Başlangıç Bölgesi Seç", command=self.activate_selection); self.aperture_menu.add_command(label="Ap: 2. İzi Bul ve Fit Et", command=self.find_and_fit_trace); self.aperture_menu.add_command(label="Ap: 3. Açıklığı Tanımla", command=self.define_aperture); self.aperture_menu.add_separator()
        self.aperture_menu.add_command(label="Arkaplan: Bölgeleri Tanımla...", command=self.start_background_selection); self.aperture_menu.add_checkbutton(label="Arkaplan Çıkar", variable=self.subtract_bg_var, onvalue=True, offvalue=False); self.aperture_menu.add_command(label="Temizle: Arkaplan Bölgeleri", command=self.clear_background_regions); self.aperture_menu.add_separator()
        self.aperture_menu.add_command(label="Ekstraksiyon: Spektrumu Çıkar ve Göster", command=self.run_extraction_and_show); self.aperture_menu.add_separator()
        self.aperture_menu.add_command(label="Normalizasyon: Süreklilik Düzeltmesi Yap", command=self.run_continuum_normalization);
        self.aperture_menu.add_command(label="Analiz: Çizgi Fiti Yap (Gaussian)...", command=self.start_line_fitting_mode)
        self.aperture_menu.add_separator()
        self.aperture_menu.add_command(label="DB Kal: Ark Lambası Aç...", command=self.load_arc_lamp); self.aperture_menu.add_command(label="DB Kal: Çizgi Listesi Yükle...", command=self.load_line_list); self.aperture_menu.add_command(label="DB Kal: Kalibrasyon Aracını Başlat...", command=self.start_wavelength_calibration);
        self.aperture_menu.add_command(label="DB Kal: Çözümü Kaydet...", command=self.save_wavelength_solution)
        self.aperture_menu.add_command(label="DB Kal: Çözümü Yükle...", command=self.load_wavelength_solution)
        self.aperture_menu.add_separator();
        self.aperture_menu.add_checkbutton(label="Grafikte Dalgaboyu Hatalarını Göster", variable=self.show_wl_error_var, command=self.toggle_errorbar_display) # DEĞİŞİKLİK
        self.aperture_menu.add_separator();
        self.aperture_menu.add_separator()
        self.aperture_menu.add_command(label="Görüntü: Uzaysal Profil Çiz...", command=self.start_spatial_profile_selection) # YENİ
        self.settings_menu = tk.Menu(self.aperture_menu, tearoff=0); self.aperture_menu.add_cascade(label="Ayarlar", menu=self.settings_menu)
        self.settings_menu.add_command(label="Açıklık Genişliği...", command=self.set_aperture_width_dialog);
        self.settings_menu.add_command(label="İz Fit Polinom Derecesi...", command=self.set_polynomial_degree_dialog);
        self.settings_menu.add_command(label="Süreklilik Fit Derecesi...", command=self.set_continuum_degree_dialog);
        self.continuum_fit_type_menu = tk.Menu(self.settings_menu, tearoff=0); self.settings_menu.add_cascade(label="Süreklilik Fit Tipi", menu=self.continuum_fit_type_menu)
        self.continuum_fit_type_menu.add_radiobutton(label="Polynomial", variable=self.continuum_fit_type, value="Polynomial", command=self.on_continuum_fit_type_change)
        spline_state = tk.NORMAL if SCIPY_AVAILABLE else tk.DISABLED; self.continuum_fit_type_menu.add_radiobutton(label="Spline", variable=self.continuum_fit_type, value="Spline", state=spline_state, command=self.on_continuum_fit_type_change)
        self.aperture_menu.add_separator(); self.aperture_menu.add_command(label="Temizle: Açıklık Çizimleri", command=lambda: self.clear_aperture_overlays(True, True, True, False))
        self.help_menu = tk.Menu(self.menu_bar, tearoff=0); self.menu_bar.add_cascade(label="Yardım", menu=self.help_menu); self.help_menu.add_command(label="Hakkında...", command=self.show_about_dialog)

    def _setup_gui_layout(self):
        self.aperture_menu.add_command(label="Normalizasyon: Süreklilik Düzeltmesi Yap", command=self.run_continuum_normalization)

        self.aperture_menu.add_separator()
        self.aperture_menu.add_command(label="İyileştirme: 1D Spektrumu Temizle (Kozmik Işın)...", command=self.run_cosmic_ray_rejection_1d)
        self.aperture_menu.add_separator()
        self.aperture_menu.add_command(label="Analiz: Eşdeğer Genişlik Hesapla...", command=self.start_ew_calculation_mode)
        self.aperture_menu.add_command(label="Analiz: Çizgi Fiti Yap (Gaussian)...", command=self.start_line_fitting_mode)
        self.aperture_menu.add_separator()
        self.right_panel = ttk.Frame(self.root); self.right_panel.pack(side=tk.RIGHT, fill=tk.Y, padx=5, pady=5);
        self.header_text = scrolledtext.ScrolledText(self.right_panel, wrap=tk.WORD, width=45, height=15)
        self.header_text.pack(side=tk.TOP, fill=tk.BOTH, expand=True); self.header_text.insert(tk.END, "FITS Header..."); self.header_text.config(state=tk.DISABLED)
        self.param_frame = ttk.LabelFrame(self.right_panel, text="Parametreler ve Durum");
        self.param_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=5, ipady=5)
        ap_subframe = ttk.Frame(self.param_frame); ap_subframe.pack(fill=tk.X, padx=5, pady=(0, 5))
        ttk.Label(ap_subframe, textvariable=self.poly_degree_label_var).pack(anchor=tk.W)
        ttk.Label(ap_subframe, textvariable=self.continuum_degree_label_var).pack(anchor=tk.W)
        ttk.Label(ap_subframe, textvariable=self.ap_width_label_var).pack(anchor=tk.W)
        ttk.Label(ap_subframe, textvariable=self.trace_coeffs_label_var).pack(anchor=tk.W)
        calib_subframe = ttk.LabelFrame(self.param_frame, text="Yüklü Kalibrasyonlar"); calib_subframe.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(calib_subframe, textvariable=self.bias_status_var).pack(anchor=tk.W)
        ttk.Label(calib_subframe, textvariable=self.dark_status_var).pack(anchor=tk.W)
        ttk.Label(calib_subframe, textvariable=self.flat_status_var).pack(anchor=tk.W)
        self.plot_frame = tk.Frame(self.root); self.plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5); self.figure = Figure(figsize=(7, 7), dpi=100); self.ax = self.figure.add_subplot(1, 1, 1); self.ax.set_facecolor('lightgrey'); self.ax.set_title("Görüntü alanı"); self.ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False); self.canvas = FigureCanvasTkAgg(self.figure, master=self.plot_frame); self.canvas_widget = self.canvas.get_tk_widget(); self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True); self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame); self.toolbar.update(); self.toolbar.pack(side=tk.BOTTOM, fill=tk.X); self.figure.tight_layout()
        self.status_bar = ttk.Label(self.root, textvariable=self.status_bar_var, relief=tk.SUNKEN, anchor=tk.W); self.status_bar.pack(side=tk.BOTTOM, fill=tk.X, pady=(2,0), padx=2)

    def start_spatial_profile_selection(self):
        if self.image_data is None:
            messagebox.showwarning("Veri Yok", "Uzaysal profil çizmek için önce bir bilim görüntüsü açın.", parent=self.root)
            return

        if self.rect_selector and self.rect_selector.active: self.rect_selector.set_active(False)
        if self.line_fit_selector and self.line_fit_selector.active: self.line_fit_selector.set_active(False)
        if self.ew_selector and self.ew_selector.active: self.ew_selector.set_active(False)
        if hasattr(self, 'spatial_profile_cid') and self.spatial_profile_cid: # Eğer olay bağlanmışsa, kopar
            self.canvas.mpl_disconnect(self.spatial_profile_cid)
            delattr(self, 'spatial_profile_cid')

        messagebox.showinfo("Profil Seçimi",
                             "Uzaysal profili çizilecek X-kolonunu (dikey çizgi) ana görüntü üzerinde tıklayarak seçin.",
                             parent=self.root)

        if self.spatial_profile_line:
            try:
                self.spatial_profile_line.remove()
            except ValueError: pass
            self.spatial_profile_line = None

        self.status_bar_var.set(" Durum: Uzaysal profil için X-kolonu seçin (tıklayın)...")
        self.spatial_profile_cid = self.canvas.mpl_connect('button_press_event', self._on_spatial_profile_column_selected)

    def _on_spatial_profile_column_selected(self, event):
        if event.inaxes != self.ax or event.xdata is None or self.image_data is None:
            return

        if hasattr(self, 'spatial_profile_cid') and self.spatial_profile_cid:
            self.canvas.mpl_disconnect(self.spatial_profile_cid)
            delattr(self, 'spatial_profile_cid')

        if self.spatial_profile_line:
            try: self.spatial_profile_line.remove()
            except ValueError: pass
            self.spatial_profile_line = None

        selected_x_column = int(round(event.xdata))
        if not (0 <= selected_x_column < self.image_data.shape[1]):
            messagebox.showwarning("Geçersiz Sütun", f"Seçilen X={selected_x_column} sütunu görüntü sınırları dışında.", parent=self.root)
            self.status_bar_var.set(" Durum: Hazır")
            return

        print(f"Uzaysal profil için X-kolonu seçildi: {selected_x_column}")

        min_y, max_y = self.ax.get_ylim()
        self.spatial_profile_line = self.ax.axvline(selected_x_column, color='cyan', linestyle='--', linewidth=1, label=f'Profil X={selected_x_column}')
        if self.spatial_profile_line not in self.overlay_plots:
            self.overlay_plots.append(self.spatial_profile_line)
        self._redraw_overlays()

        spatial_profile_data = self.image_data[:, selected_x_column]
        y_pixels = np.arange(self.image_data.shape[0])

        profile_title = f"Uzaysal Profil (X = {selected_x_column})"
        if self.science_file_path:
            profile_title = f"{os.path.basename(self.science_file_path)} - {profile_title}"

        if self.spatial_profile_window and self.spatial_profile_window.winfo_exists():
            self.spatial_profile_window.update_plot(y_pixels, spatial_profile_data,
                                                    x_axis_label="Y Piksel (Uzaysal Yön)",
                                                    y_axis_label="Akı / Counts",
                                                    title=profile_title)
            self.spatial_profile_window.lift()
        else:
            self.spatial_profile_window = SpectrumPlotWindow(self.root, self, y_pixels, spatial_profile_data,
                                                             title=profile_title,
                                                             x_axis_label="Y Piksel (Uzaysal Yön)")

        self.status_bar_var.set(f" Durum: Uzaysal profil X={selected_x_column} için çizildi.")

    def _save_figure_dialog(self, figure_to_save, default_filename_base):
        if figure_to_save is None: messagebox.showerror("Hata", "Kaydedilecek grafik bulunamadı.", parent=self.root); return
        filename = filedialog.asksaveasfilename(title="Grafiği Kaydet", initialfile=f"{default_filename_base}.png", defaultextension=".png", filetypes=[("PNG Dosyası", "*.png"), ("PDF Dosyası", "*.pdf"), ("JPEG Dosyası", "*.jpg"), ("SVG Dosyası", "*.svg"), ("Tüm Dosyalar", "*.*")])
        if filename:
            try: figure_to_save.savefig(filename, dpi=150); messagebox.showinfo("Başarılı", f"Grafik kaydedildi: {filename}", parent=self.root)
            except Exception as e: messagebox.showerror("Kaydetme Hatası", f"{e}", parent=self.root)

    def save_main_image_plot(self):
        if self.image_data is None: messagebox.showwarning("Veri Yok", "Kaydedilecek ana görüntü yok.", parent=self.root); return
        obj_name = "ana_goruntu"
        if self.science_header and 'OBJECT' in self.science_header: obj_name = self.science_header.get('OBJECT', obj_name).replace(' ', '_')
        elif self.hdulist and self.hdulist.filename(): obj_name = os.path.splitext(os.path.basename(self.hdulist.filename()))[0]
        self._save_figure_dialog(self.figure, obj_name)

    def save_science_spectrum_plot(self):
        if self.science_plot_window and self.science_plot_window.winfo_exists():
            obj_name = "bilim_spektrumu"
            if self.science_header and 'OBJECT' in self.science_header: obj_name = self.science_header.get('OBJECT', obj_name).replace(' ', '_')
            self._save_figure_dialog(self.science_plot_window.fig, f"{obj_name}_bilim_spek")
        else: messagebox.showwarning("Grafik Yok", "Kaydedilecek açık bilim spektrumu grafiği bulunamadı.", parent=self.root)

    def save_arc_spectrum_plot(self):
        if self.arc_plot_window and self.arc_plot_window.winfo_exists():
            obj_name = "ark_spektrumu"
            if self.science_header and 'OBJECT' in self.science_header: obj_name = self.science_header.get('OBJECT', obj_name).replace(' ', '_')
            elif self.arc_image_data is not None and hasattr(self.arc_image_data, 'meta') and 'OBJECT' in self.arc_image_data.meta : obj_name = self.arc_image_data.meta.get('OBJECT', obj_name).replace(' ', '_')
            self._save_figure_dialog(self.arc_plot_window.fig, f"{obj_name}_ark_spek")
        else: messagebox.showwarning("Grafik Yok", "Kaydedilecek açık ark spektrumu grafiği bulunamadı.", parent=self.root)

    def show_about_dialog(self):
        about_text = ("Ulanom\nSürüm: 1.5 (Geliştirme Sürecinde)\nYapımcı: Emre Bilgin\n\nTemel Özellikler:\n- FITS Görüntüleme\n- Temel Kalibrasyon (B/D/F) ve Dark Ölçekleme\n- Apertür Tanımlama ve İz Sürme\n- Arkaplan Çıkarma\n- Spektrum Ekstraksiyonu\n- Dalgaboyu Kalibrasyonu (Hata hesabı ile, Çözüm Kaydet/Yükle)\n- Süreklilik Normalizasyonu (Polynomial/Spline)\n- Gaussian Çizgi Fiti\n- Sonuçları Kaydetme (Metin/Resim)\n- FITS Header Bilgileri Kullanımı\n\n" f"SciPy (Spline/Gaussian Fit için): {'Kullanılabilir' if SCIPY_AVAILABLE else 'Yüklü Değil'}")
        messagebox.showinfo("Hakkında", about_text, parent=self.root)

    def _update_menu_selection(self):
        if hasattr(self, 'stretch_var'): self.stretch_var.set(self.current_stretch_type)
        if hasattr(self, 'limits_var'):
            if self.current_limit_type == 'percentile': menu_val = f'percentile_{self.current_percentile}'
            elif self.current_limit_type == 'manual': menu_val = ""
            else: menu_val = self.current_limit_type
            try: self.limits_var.set(menu_val)
            except tk.TclError: self.limits_var.set("")
    def open_fits_file(self):
        print("Bilim görüntüsü açma işlemi başlatıldı...")
        file_path = filedialog.askopenfilename(
            title="Bilim FITS Görüntüsü Seç",
            filetypes=(("FITS Dosyaları", "*.fits *.fit *.fts"), ("Tüm Dosyalar", "*.*"))
        )
        if not file_path:
            print("Dosya seçimi iptal edildi.")
            return

        try:
            self.cleanup_ui(clear_calib=False)

            self.science_file_path = file_path
            with fits.open(file_path) as hdul:
                self.hdulist = hdul
                science_hdu_index = -1
                for i, hdu_item in enumerate(hdul):
                    if isinstance(hdu_item, (fits.ImageHDU, fits.PrimaryHDU)) and \
                       hasattr(hdu_item, 'data') and hdu_item.data is not None and \
                       hdu_item.data.ndim == 2:
                        science_hdu_index = i
                        break

                if science_hdu_index == -1:
                    messagebox.showerror("Hata", "Dosyada uygun 2D görüntü verisi bulunamadı.", parent=self.root)
                    self.hdulist = None
                    return

                self.raw_image_data = hdul[science_hdu_index].data.astype(np.float32)
                self.science_header = hdul[science_hdu_index].header
                print(f"Ham bilim verisi '{os.path.basename(file_path)}' okundu (şekil: {self.raw_image_data.shape}).")

                try:
                    self.science_exptime = float(self.science_header.get('EXPTIME', None))
                    if self.science_exptime is not None:
                        print(f"Bilim çerçevesi poz süresi: {self.science_exptime} s")
                    else:
                        print("Uyarı: Bilim çerçevesi için 'EXPTIME' header'da bulunamadı veya geçersiz.")
                except (TypeError, ValueError):
                    self.science_exptime = None
                    print("Uyarı: Bilim çerçevesi için 'EXPTIME' değeri okunamadı/geçersiz.")

            self.header_text.config(state=tk.NORMAL)
            self.header_text.delete('1.0', tk.END)
            if self.science_header:
                try:
                    header_str = self.science_header.tostring(sep='\n', endcard=False, padding=False)
                    self.header_text.insert(tk.END, header_str)
                except Exception as e:
                    self.header_text.insert(tk.END, f"Header okunamadı: {e}")
            else:
                self.header_text.insert(tk.END, "Header bulunamadı.")
            self.header_text.config(state=tk.DISABLED)

            print("Yüklü kalibrasyonlar bilim verisine uygulanıyor (varsa)...")
            self.image_data = self._apply_calibrations(self.raw_image_data, is_science_frame=True)
            if self.image_data is None:
                self.image_data = self.raw_image_data.copy()
                print("Uyarı: Kalibrasyon uygulanamadı, ham veri kullanılıyor.")

            self.apply_scaling()

            self.root.title(f"FITS Görüntüleyici - {os.path.basename(file_path)}")
            messagebox.showinfo("Başarılı", f"Bilim görüntüsü yüklendi:\n{os.path.basename(file_path)}", parent=self.root)
            self.status_bar_var.set(f" Durum: {os.path.basename(file_path)} yüklendi.")

        except FileNotFoundError:
            messagebox.showerror("Hata", f"Dosya bulunamadı: {file_path}", parent=self.root)
            self.cleanup_ui(clear_calib=False)
        except Exception as e:
            messagebox.showerror("FITS Okuma Hatası", f"Dosya okunurken/işlenirken hata oluştu: {e}", parent=self.root)
            print(f"FITS Okuma Hatası: {e}")
            import traceback
            traceback.print_exc()
            self.cleanup_ui(clear_calib=False)

    def apply_scaling_from_menu(self):
        self.current_stretch_type = self.stretch_var.get(); selected_limit_menu_value = self.limits_var.get()
        if not selected_limit_menu_value.startswith('percentile_'):
            if selected_limit_menu_value != 'manual': self.current_limit_type = selected_limit_menu_value
            else: print("Uyarı: Manuel limitler menüden ayarlanamaz, diyalog kullanın."); self._update_menu_selection()
        print("Menüden ölçekleme uygulanıyor...")
        self.apply_scaling(redraw_overlays=True)

    def set_percentile_from_menu_and_apply(self, percentile_value):
        self.current_limit_type = 'percentile'; self.current_percentile = percentile_value; print(f"Percentile set: {self.current_percentile}%")
        self.apply_scaling(redraw_overlays=True)

    def _load_calibration_frame(self, frame_type):
        title_map = {'bias': "Bias Çerçevesi Seç", 'dark': "Dark Çerçevesi Seç", 'flat': "Flat Çerçevesi Seç"}
        file_path = filedialog.askopenfilename(title=title_map.get(frame_type, "Kalibrasyon Seç"), filetypes=(("FITS", "*.fits *.fit *.fts"), ("Hepsi", "*.*")))
        if not file_path: return False
        try:
            with fits.open(file_path) as hdul:
                calib_hdu_index = -1;
                for i, hdu in enumerate(hdul):
                    if isinstance(hdu, (fits.ImageHDU, fits.PrimaryHDU)) and hasattr(hdu, 'data') and hdu.data is not None and hdu.data.ndim == 2: calib_hdu_index = i; break
                if calib_hdu_index == -1: raise ValueError("Dosyada uygun 2D görüntü verisi bulunamadı.")
                calib_data = hdul[calib_hdu_index].data.astype(np.float32); calib_header = hdul[calib_hdu_index].header
                print(f"{frame_type.capitalize()} verisi '{os.path.basename(file_path)}' dosyasından okundu (şekil: {calib_data.shape}).")
                current_data_shape = None
                if self.raw_image_data is not None: current_data_shape = self.raw_image_data.shape
                elif self.raw_arc_image_data is not None: current_data_shape = self.raw_arc_image_data.shape
                if current_data_shape and calib_data.shape != current_data_shape: messagebox.showwarning("Boyut Uyuşmazlığı", f"{frame_type.capitalize()} ({calib_data.shape}) ve mevcut görüntü ({current_data_shape}) boyutları farklı!", parent=self.root)
                if frame_type == 'dark':
                    try: self.dark_exptime = float(calib_header.get('EXPTIME', None));
                    except (TypeError, ValueError): self.dark_exptime = None
                    if self.dark_exptime is not None: print(f"Dark çerçevesi poz süresi: {self.dark_exptime} s")
                    else: print("Uyarı: Dark çerçevesi için 'EXPTIME' header'da bulunamadı veya geçersiz.")
                elif frame_type == 'flat':
                    print("Flat normalizasyonu yapılıyor (medyan=1)...")
                    valid_flat_pixels = calib_data[np.isfinite(calib_data)]
                    if valid_flat_pixels.size > 0:
                        flat_median = np.nanmedian(valid_flat_pixels);
                        if np.isclose(flat_median, 0) or not np.isfinite(flat_median): messagebox.showerror("Flat Hatası", "Flat medyanı sıfır veya geçersiz."); calib_data = None
                        else: calib_data /= flat_median
                    else: messagebox.showerror("Flat Hatası", "Flat içinde geçerli piksel bulunamadı."); calib_data = None
                if calib_data is not None:
                    setattr(self, f"{frame_type}_data", calib_data); setattr(self, f"{frame_type}_path", file_path)
                    self._update_calib_status_labels(); messagebox.showinfo("Başarılı", f"{frame_type.capitalize()} yüklendi:\n{os.path.basename(file_path)}", parent=self.root)
                    if self.raw_image_data is not None or self.raw_arc_image_data is not None:
                        if messagebox.askyesno("Kalibrasyon Uygula", f"Yeni {frame_type} yüklendi. Mevcut verilere kalibrasyonları şimdi yeniden uygulamak ister misiniz?", parent=self.root): self.reapply_calibrations()
                    return True
                else:
                    setattr(self, f"{frame_type}_data", None); setattr(self, f"{frame_type}_path", None)
                    if frame_type == 'dark': self.dark_exptime = None
                    self._update_calib_status_labels(); return False
        except FileNotFoundError: messagebox.showerror("Hata", f"Dosya bulunamadı: {file_path}", parent=self.root); return False
        except Exception as e: messagebox.showerror(f"{frame_type.capitalize()} Yükleme Hatası", f"Hata: {e}", parent=self.root); return False

    def load_bias_frame(self): self._load_calibration_frame('bias')
    def load_dark_frame(self): self._load_calibration_frame('dark')
    def load_flat_frame(self): self._load_calibration_frame('flat')

    def _apply_calibrations(self, raw_data, is_science_frame=True):
        if raw_data is None: return None
        print("Kalibrasyonlar uygulanıyor...")
        corrected_data = raw_data.copy().astype(np.float64)
        if self.bias_data is not None:
            if self.bias_data.shape == corrected_data.shape: corrected_data -= self.bias_data; print("- Bias uygulandı.")
            else: print(f"! Uyarı: Bias boyutu uyuşmuyor. Atlanıyor.")
        if self.dark_data is not None:
            if self.dark_data.shape == corrected_data.shape:
                scaled_dark = self.dark_data.copy()
                current_frame_exptime = self.science_exptime if is_science_frame else None
                if self.dark_exptime and self.dark_exptime > 1e-6 and current_frame_exptime and current_frame_exptime > 1e-6:
                    if not np.isclose(self.dark_exptime, current_frame_exptime):
                        scale_factor = current_frame_exptime / self.dark_exptime
                        print(f"- Dark çerçevesi {scale_factor:.3f} ile ölçekleniyor (Hedef S.: {current_frame_exptime}s, Dark S.: {self.dark_exptime}s).")
                        scaled_dark *= scale_factor
                    else: print("- Dark ve Hedef poz süreleri aynı, ölçekleme yapılmadı.")
                elif self.dark_exptime and self.dark_exptime > 1e-6 and (current_frame_exptime is None or current_frame_exptime <= 1e-6) and is_science_frame: print(f"- Bilim poz süresi bilinmiyor/geçersiz, Dark ölçeklenemedi. Ham Dark çıkarılıyor.")
                elif (self.dark_exptime is None or self.dark_exptime <= 1e-6): print("- Dark poz süresi bilinmiyor/geçersiz, ölçekleme yapılamadı. Ham Dark çıkarılıyor.")
                corrected_data -= scaled_dark; print("- Dark uygulandı.")
            else: print(f"! Uyarı: Dark boyutu uyuşmuyor. Atlanıyor.")
        if self.flat_data is not None:
            if self.flat_data.shape == corrected_data.shape:
                safe_flat = self.flat_data.copy(); zero_mask = np.isclose(safe_flat, 0) | ~np.isfinite(safe_flat); num_zeros = np.sum(zero_mask)
                if num_zeros > 0: print(f"! Uyarı: Flat içinde {num_zeros} adet sıfır/geçersiz değer bulundu. Bu pikseller için bölme yapılmayacak."); safe_flat[zero_mask] = 1.0
                corrected_data /= safe_flat; print("- Flat uygulandı.")
            else: print(f"! Uyarı: Flat boyutu uyuşmuyor. Atlanıyor.")
        print("Kalibrasyon uygulaması tamamlandı."); return corrected_data.astype(np.float32)

    def reapply_calibrations(self):
        print("Kalibrasyonlar yeniden uygulanıyor..."); reapplied_science = False; reapplied_arc = False
        if self.raw_image_data is not None:
            print("Bilim verisine uygulanıyor..."); self.image_data = self._apply_calibrations(self.raw_image_data, is_science_frame=True)
            if self.image_data is not None: self.apply_scaling(); reapplied_science = True
            else: messagebox.showerror("Hata", "Bilim verisine kalibrasyon uygulanırken hata.", parent=self.root); self.image_data = self.raw_image_data; self.apply_scaling()
        if self.raw_arc_image_data is not None:
            print("Ark verisine uygulanıyor..."); self.arc_image_data = self._apply_calibrations(self.raw_arc_image_data, is_science_frame=False)
            if self.arc_image_data is not None:
                if self.aperture_lower is not None and self.aperture_upper is not None:
                    print("Ark spektrumu yeniden çıkarılıyor...")
                    arc_pix, arc_flux = self.extract_spectrum(image_data_to_extract=self.arc_image_data, aperture_lower=self.aperture_lower, aperture_upper=self.aperture_upper, subtract_background=False)
                    if arc_pix is not None:
                        self.arc_extracted_pixels = arc_pix; self.arc_extracted_flux = arc_flux
                        if self.arc_plot_window and self.arc_plot_window.winfo_exists(): self.show_1d_spectrum_window(self.arc_extracted_pixels, self.arc_extracted_flux, title="Ark Spektrumu (Yeniden Kalibre Edildi)", x_axis_label="Piksel (X)", y_axis_label="Akı / Counts")
                        print("Ark spektrumu güncellendi.")
                reapplied_arc = True
            else: messagebox.showerror("Hata", "Ark verisine kalibrasyon uygulanırken hata.", parent=self.root); self.arc_image_data = self.raw_arc_image_data
        if reapplied_science or reapplied_arc: messagebox.showinfo("Başarılı", "Kalibrasyonlar mevcut verilere yeniden uygulandı.", parent=self.root)

    def clear_calibrations(self):
        if (self.bias_data is None and self.dark_data is None and self.flat_data is None): messagebox.showinfo("Temiz", "Zaten yüklü kalibrasyon çerçevesi yok.", parent=self.root); return
        if messagebox.askyesno("Onay", "Yüklü tüm Bias, Dark ve Flat verilerini temizlemek istediğinizden emin misiniz?\n(Mevcut bilim/ark verisine kalibrasyonlar yeniden uygulanacak)", parent=self.root):
            print("Kalibrasyon verileri temizleniyor..."); self.bias_data = None; self.bias_path = None; self.dark_data = None; self.dark_path = None; self.dark_exptime = None; self.flat_data = None; self.flat_path = None
            self._update_calib_status_labels(); self.reapply_calibrations(); messagebox.showinfo("Temizlendi", "Tüm kalibrasyon verileri temizlendi.", parent=self.root)

    def _update_calib_status_labels(self):
        bias_txt = f"Bias: {os.path.basename(self.bias_path)}" if self.bias_path else "Bias: Yok"
        dark_txt = f"Dark: {os.path.basename(self.dark_path)}" if self.dark_path else "Dark: Yok"
        if self.dark_exptime is not None: dark_txt += f" ({self.dark_exptime}s)"
        flat_txt = f"Flat: {os.path.basename(self.flat_path)}" if self.flat_path else "Flat: Yok"
        if hasattr(self, 'bias_status_var'): self.bias_status_var.set(bias_txt)
        if hasattr(self, 'dark_status_var'): self.dark_status_var.set(dark_txt)
        if hasattr(self, 'flat_status_var'): self.flat_status_var.set(flat_txt)

    def activate_selection(self, selection_type='aperture'):
        if self.image_data is None: messagebox.showwarning("Veri Yok", "Önce bir bilim görüntüsü açın.", parent=self.root); return
        if self.rect_selector and self.rect_selector.active: messagebox.showinfo("Aktif", "Seçim aracı zaten aktif.", parent=self.root); return
        mode_text = "Başlangıç Açıklık Bölgesi" if selection_type == 'aperture' else "Arkaplan Bölgesi"; print(f"{mode_text} Seçimi Aktif.")
        messagebox.showinfo("Seçim Modu", f"{mode_text} seçmek için dikdörtgen çizin.\nBitirmek için fareyi bırakın.", parent=self.root)
        if self.rect_selector:
            try: self.rect_selector.set_active(False)
            except Exception as e: print(f"Eski seçici pasif yapılırken hata (önemsiz): {e}")
        callback = self.on_select_region if selection_type == 'aperture' else self.on_background_select
        props_dict = dict(facecolor='none', edgecolor='cyan' if selection_type == 'aperture' else 'red', linestyle='--')
        self.rect_selector = RectangleSelector(self.ax, callback, useblit=True, button=[1], minspanx=5, minspany=5, spancoords='data', interactive=True, props=props_dict)
        self.toolbar.set_message(f"{mode_text} seçimi aktif...")

    def on_select_region(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        if None in [x1, y1, x2, y2]:
            print("Seçim koordinatları alınamadı.")
            self.toolbar.set_message("")
            if self.rect_selector:
                self.rect_selector.set_active(False)
            return

        xmin, xmax = int(round(min(x1, x2))), int(round(max(x1, x2)))
        ymin, ymax = int(round(min(y1, y2))), int(round(max(y1, y2)))
        self.selection_box = (xmin, xmax, ymin, ymax)
        print(f"Açıklık Başlangıç Bölgesi seçildi: X={xmin}-{xmax}, Y={ymin}-{ymax}")

        self.clear_aperture_overlays(clear_selection=True, clear_trace=True, clear_aperture=True, clear_background=False)
        self._redraw_overlays()

        if self.rect_selector:
            self.rect_selector.set_active(False)

        self.toolbar.set_message("")

        current_mode = self.toolbar.mode
        if 'pan/zoom' in current_mode:
            print("Araç çubuğu PAN/ZOOM modundaydı, devre dışı bırakılıyor.")
            self.toolbar.pan()
        elif 'zoom rect' in current_mode:
            print("Araç çubuğu ZOOM modundaydı, devre dışı bırakılıyor.")
            self.toolbar.zoom()

        self.canvas.draw_idle()

        messagebox.showinfo("Bölge Seçildi", f"Seçilen Bölge:\nX: {xmin}-{xmax}\nY: {ymin}-{ymax}\n\nŞimdi 'İzi Bul ve Fit Et' adımına geçebilirsiniz.", parent=self.root)

    def _clean_cosmics_1d_iterative_median(self, flux_array_input, sigma_clip_threshold=5.0, box_width=5, iterations=2):
        """
        1D akı dizisindeki kozmik ışınları iteratif medyan filtresi ve sigma kırpma ile temizler.
        Tespit edilen kozmik ışınlar komşu piksellerin medyanı ile değiştirilir.
        """
        if flux_array_input is None or len(flux_array_input) < box_width:
            print("Kozmik ışın temizleme için yetersiz veri.")
            return flux_array_input

        flux_array = np.copy(flux_array_input) # Orijinal diziyi değiştirmemek için kopyala
        n_pixels = len(flux_array)
        half_box = box_width // 2

        for iteration in range(iterations):
            print(f"Kozmik ışın temizleme iterasyonu {iteration + 1}/{iterations}...")
            num_cleaned_in_iter = 0
            original_flux_for_iter = np.copy(flux_array) # Her iterasyonda o anki duruma göre işlem yap

            for i in range(n_pixels):
                start = max(0, i - half_box)
                end = min(n_pixels, i + half_box + 1)

                window_indices = list(range(start, i)) + list(range(i + 1, end)) # Ortadaki pikseli çıkarma
                if not window_indices:
                    continue

                local_median = np.nanmedian(original_flux_for_iter[window_indices])
                local_std = np.nanstd(original_flux_for_iter[window_indices])

                if local_std == 0 or not np.isfinite(local_std) or not np.isfinite(local_median):
                    continue

                if abs(original_flux_for_iter[i] - local_median) > sigma_clip_threshold * local_std:
                    flux_array[i] = local_median
                    num_cleaned_in_iter += 1

            print(f"İterasyon {iteration + 1}: {num_cleaned_in_iter} piksel temizlendi/düzeltildi.")
            if num_cleaned_in_iter == 0:
                break

        return flux_array

    def run_cosmic_ray_rejection_1d(self):
        print("1D Kozmik Işın Temizleme işlemi başlatılıyor...")
        if self.extracted_flux is None:
            messagebox.showwarning("Veri Yok", "Temizlenecek çıkarılmış 1D spektrum verisi bulunmuyor.\nLütfen önce spektrumu çıkarın.", parent=self.root)
            return

        if not messagebox.askyesno("Onay",
                                   "Mevcut çıkarılmış bilim spektrumu üzerindeki kozmik ışınlar temizlenecektir.\n"
                                   "Bu işlem, mevcut akı verilerini değiştirecektir.\nDevam etmek istiyor musunuz?",
                                   parent=self.root):
            print("Kozmik ışın temizleme işlemi iptal edildi.")
            return

        cleaned_flux = self._clean_cosmics_1d_iterative_median(self.extracted_flux, sigma_clip_threshold=5.0, box_width=5, iterations=2)

        num_diff_pixels = np.sum(self.extracted_flux != cleaned_flux)
        if num_diff_pixels > 0 :
            self.extracted_flux = cleaned_flux
            print(f"{num_diff_pixels} pikselde kozmik ışın temizlemesi uygulandı.")

            if self.science_normalized_flux is not None:
                self.science_normalized_flux = None
                self.science_continuum_flux = None
                print("Kozmik ışın temizlemesi yapıldı, mevcut normalizasyon kaldırıldı. Lütfen spektrumu yeniden normalize edin.")
                messagebox.showinfo("Normalizasyon Kaldırıldı",
                                     "Kozmik ışın temizlemesi sonrası spektrum değiştiği için mevcut normalizasyon kaldırıldı.\n"
                                     "Lütfen gerekirse spektrumu yeniden normalize edin.", parent=self.root)

            if self.science_plot_window and self.science_plot_window.winfo_exists():
                print("Bilim spektrumu grafiği temizlenmiş veri ile güncelleniyor...")
                current_x_data = self.science_plot_window.line.get_xdata() if self.science_plot_window.line else self.extracted_pixels
                current_x_label = self.science_plot_window.ax.get_xlabel()
                current_y_label = "Akı / Counts (Temizlenmiş)"
                current_title = self.science_plot_window.ax.get_title().replace(" (Normalize Edilmiş)", "")

                x_error_data_to_plot = None
                if "Dalga Boyu" in current_x_label and self.show_wl_error_var.get() and self.science_wavelength_error is not None and \
                   len(self.science_wavelength_error) == len(current_x_data):
                    x_error_data_to_plot = self.science_wavelength_error

                self.show_1d_spectrum_window(current_x_data, self.extracted_flux,
                                             title=current_title,
                                             x_axis_label=current_x_label,
                                             y_axis_label=current_y_label,
                                             x_error_data=x_error_data_to_plot)
            messagebox.showinfo("Başarılı", "1D spektrum kozmik ışınlardan temizlendi.", parent=self.root)
        else:
            messagebox.showinfo("Değişiklik Yok", "Kozmik ışın temizleme işlemi herhangi bir pikseli değiştirmedi.", parent=self.root)

    def draw_selection_box(self):
        if self.selection_box is None: return
        xmin, xmax, ymin, ymax = self.selection_box; width = xmax - xmin; height = ymax - ymin
        if self.selection_patch:
            try: self.selection_patch.remove()
            except ValueError: pass
        self.selection_patch = patches.Rectangle((xmin, ymin), width, height, linewidth=1, edgecolor='cyan', facecolor='none', linestyle='--', label='Seçim Bölgesi', zorder=11)
        self.ax.add_patch(self.selection_patch)
        if self.selection_patch not in self.overlay_plots: self.overlay_plots.append(self.selection_patch)

    def find_and_fit_trace(self):
        if self.image_data is None: messagebox.showerror("Hata", "İz bulmak için önce bir bilim görüntüsü açın.", parent=self.root); return
        if self.selection_box is None: messagebox.showerror("Hata", "İz bulmak için önce bir bölge seçin ('Ap: 1...').", parent=self.root); return
        xmin, xmax, ymin, ymax = self.selection_box; img_h, img_w = self.image_data.shape; xmin_c, xmax_c = max(0, xmin), min(img_w, xmax); ymin_c, ymax_c = max(0, ymin), min(img_h, ymax)
        if xmax_c <= xmin_c or ymax_c <= ymin_c: messagebox.showerror("Hata", "Seçilen bölge geçersiz veya görüntü dışında.", parent=self.root); return
        print(f"İz, X=[{xmin_c}, {xmax_c}), Y=[{ymin_c}, {ymax_c}) aralığında aranıyor..."); x_coords, y_peaks = [], []
        for x in range(xmin_c, xmax_c):
            column_data = self.image_data[ymin_c:ymax_c, x]; valid_indices = np.where(np.isfinite(column_data))[0]
            if valid_indices.size == 0: continue
            try: relative_peak_y_idx = np.argmax(column_data[valid_indices]); relative_peak_y = valid_indices[relative_peak_y_idx]; y_peaks.append(ymin_c + relative_peak_y); x_coords.append(x)
            except Exception as e: print(f"Sütun {x} için pik bulmada hata: {e}"); continue
        if len(x_coords) < self.polynomial_degree + 1: messagebox.showerror("Yetersiz Pik", f"Seçilen bölgede {self.polynomial_degree}. derece fit için yeterli pik bulunamadı ({len(x_coords)} pik bulundu, en az {self.polynomial_degree + 1} gerekli). Bölgeyi genişletmeyi veya polinom derecesini düşürmeyi deneyin.", parent=self.root); self.trace_peaks_x = None; self.trace_peaks_y = None; self.trace_coefficients = None; self.trace_coeffs_label_var.set("Trace Fit Katsayıları: - (Yetersiz Pik)"); self.clear_aperture_overlays(clear_selection=False, clear_trace=True, clear_aperture=True, clear_background=False); self._redraw_overlays(); return
        self.trace_peaks_x = np.array(x_coords); self.trace_peaks_y = np.array(y_peaks)
        try:
            self.trace_coefficients = np.polyfit(self.trace_peaks_x, self.trace_peaks_y, self.polynomial_degree); print(f"İz bulundu ve {self.polynomial_degree}. derece polinom fit edildi."); print(f"Fit Katsayıları: {self.trace_coefficients}")
            coeff_str = ", ".join([f"{c:.3e}" for c in self.trace_coefficients]); self.trace_coeffs_label_var.set(f"Trace Fit Katsayıları: [{coeff_str}]"); self.clear_aperture_overlays(clear_selection=False, clear_trace=True, clear_aperture=True, clear_background=False); self._redraw_overlays()
            messagebox.showinfo("Başarılı", f"{len(x_coords)} pik bulundu ve iz başarıyla fit edildi.", parent=self.root)
        except Exception as e: messagebox.showerror("Fit Hatası", f"İz fit edilirken hata oluştu: {e}", parent=self.root); self.trace_coefficients = None; self.trace_coeffs_label_var.set("Trace Fit Katsayıları: - (Fit Hatası)"); self.clear_aperture_overlays(clear_selection=False, clear_trace=True, clear_aperture=True, clear_background=False); self._redraw_overlays()

    def define_aperture(self):
        if self.image_data is None: messagebox.showerror("Hata", "Açıklık tanımlamak için önce bir bilim görüntüsü açın.", parent=self.root); return
        if self.trace_coefficients is None: messagebox.showerror("Hata", "Açıklık tanımlamak için önce izi fit edin ('Ap: 2...').", parent=self.root); return
        img_h, img_w = self.image_data.shape; fit_func = np.poly1d(self.trace_coefficients); x_range = np.arange(img_w); y_center = fit_func(x_range); half_width = self.aperture_width / 2.0;
        self.aperture_lower = y_center - half_width; self.aperture_upper = y_center + half_width; print(f"Açıklık sınırları {self.aperture_width} piksel genişlikle tanımlandı.")
        self.clear_aperture_overlays(clear_selection=False, clear_trace=False, clear_aperture=True, clear_background=False); self._redraw_overlays(); messagebox.showinfo("Başarılı", "Açıklık sınırları başarıyla tanımlandı.", parent=self.root)

    def set_aperture_width_dialog(self):
        new_width = simpledialog.askinteger("Açıklık Genişliği", "Yeni açıklık genişliğini piksel olarak girin:", parent=self.root, initialvalue=self.aperture_width, minvalue=1)
        if new_width is not None and new_width != self.aperture_width:
            self.aperture_width = new_width; self.ap_width_label_var.set(f"Açıklık Genişliği: {self.aperture_width} px"); print(f"Açıklık genişliği {self.aperture_width} piksel olarak ayarlandı.")
            if self.trace_coefficients is not None: print("Genişlik değişti, açıklık yeniden tanımlanıyor..."); self.define_aperture()

    def set_polynomial_degree_dialog(self):
        new_degree = simpledialog.askinteger("İz Fit Polinom Derecesi", "İz fit için yeni polinom derecesini girin:", parent=self.root, initialvalue=self.polynomial_degree, minvalue=0)
        if new_degree is not None and new_degree != self.polynomial_degree:
            self.polynomial_degree = new_degree; self.poly_degree_label_var.set(f"Polinom Derecesi: {self.polynomial_degree}"); print(f"İz fit polinom derecesi {self.polynomial_degree} olarak ayarlandı.")
            if self.selection_box is not None:
                print("Polinom derecesi değişti, iz yeniden fit ediliyor..."); self.trace_peaks_x = None; self.trace_peaks_y = None; self.trace_coefficients = None; self.aperture_lower = None; self.aperture_upper = None
                self.trace_coeffs_label_var.set("Trace Fit Katsayıları: - (Derece Değişti, Yeniden Fit Edin)"); self.clear_aperture_overlays(clear_selection=False, clear_trace=True, clear_aperture=True, clear_background=False); self._redraw_overlays()
                self.find_and_fit_trace()
            else: self.trace_coefficients = None; self.trace_coeffs_label_var.set("Trace Fit Katsayıları: - (Derece Değişti)")

    def clear_aperture_overlays(self, clear_selection=True, clear_trace=True, clear_aperture=True, clear_background=True):
        elements_to_remove = []; new_overlay_plots = []; needs_legend_update = False
        if clear_background:
            for patch in self.background_patches:
                try: patch.remove(); needs_legend_update = True
                except ValueError: pass
            self.background_patches.clear()
        for element in self.overlay_plots:
            remove_this = False; is_selection = element == self.selection_patch; is_scatter = hasattr(element, 'get_offsets'); is_line = isinstance(element, plt.Line2D); is_bg_patch = isinstance(element, patches.Rectangle) and element.get_edgecolor() == 'red'
            is_trace_line = False; is_aperture_line = False
            if is_line: style = element.get_linestyle(); color = element.get_color(); is_trace_line = style == '-' and color == 'lime'; is_aperture_line = style == '--' and color == 'yellow'
            if clear_selection and is_selection: remove_this = True
            if clear_trace and (is_scatter or is_trace_line): remove_this = True
            if clear_aperture and is_aperture_line: remove_this = True
            if clear_background and is_bg_patch: remove_this = True
            if remove_this:
                if not (clear_background and is_bg_patch):
                    try: element.remove(); elements_to_remove.append(element); needs_legend_update=True;
                    except (ValueError, AttributeError, TypeError): pass
            else:
                if not is_bg_patch or not clear_background: new_overlay_plots.append(element)
        self.overlay_plots = new_overlay_plots
        if clear_selection and self.selection_patch in elements_to_remove: self.selection_patch = None
        if needs_legend_update:
            current_legend = self.ax.get_legend()
            if current_legend: current_legend.remove()
            if self.overlay_plots:
                handles, labels = [], []
                for item in self.overlay_plots:
                    if hasattr(item, 'get_label') and item.get_label() and not item.get_label().startswith('_'):
                        handles.append(item); labels.append(item.get_label())
                if handles:
                    unique_handles_labels = {};
                    for h, l in zip(handles, labels):
                        if l not in unique_handles_labels and h is not None: unique_handles_labels[l] = h
                    if unique_handles_labels: self.ax.legend(unique_handles_labels.values(), unique_handles_labels.keys(), fontsize='x-small', loc='upper right')
            self.canvas.draw_idle()
        else: self.canvas.draw_idle()

    def start_background_selection(self):
        if self.image_data is None: messagebox.showerror("Hata", "Önce FITS açın."); return
        self.clear_background_regions()
        print("Arkaplan Bölgesi Seçimi Başladı."); messagebox.showinfo("Arkaplan Seçimi", "Arkaplanı temsil eden bir veya daha fazla dikdörtgen bölge çizin.\nHer bölge sonrası size sorulacak.", parent=self.root)
        self.activate_selection(selection_type='background')

    def on_background_select(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata; x2, y2 = erelease.xdata, erelease.ydata
        if None in [x1, y1, x2, y2]: print("Seçim başarısız."); self.toolbar.set_message("");

        xmin, xmax = int(round(min(x1, x2))), int(round(max(x1, x2))); ymin, ymax = int(round(min(y1, y2))), int(round(max(y1, y2)))
        if self.image_data is not None and self.aperture_lower is not None:
            ap_low_in_range = self.aperture_lower[max(0, xmin):min(self.image_data.shape[1], xmax)]
            ap_high_in_range = self.aperture_upper[max(0, xmin):min(self.image_data.shape[1], xmax)]
            if ap_low_in_range.size > 0 and ap_high_in_range.size > 0:
                ap_low_avg = np.nanmean(ap_low_in_range); ap_high_avg = np.nanmean(ap_high_in_range)
                if ymax > ap_low_avg and ymin < ap_high_avg:
                    if not messagebox.askyesno("Uyarı: Çakışma", "Seçilen arkaplan bölgesi açıklık ile çakışıyor.\n Yine de eklensin mi?", parent=self.root): self.toolbar.set_message("");
                else: print("Uyarı: Açıklık sınırları bu X aralığında hesaplanamadı.")
        new_region = (xmin, xmax, ymin, ymax)
        self.background_regions.append(new_region)
        print(f"Arkaplan bölgesi eklendi: {new_region}")
        self._redraw_overlays()
        self.toolbar.set_message("")

        if messagebox.askyesno("Devam?", "Başka bir arkaplan bölgesi eklemek ister misiniz?", parent=self.root):
            self.activate_selection(selection_type='background')
        else:
            print(f"Arkaplan bölgesi seçimi tamamlandı. Toplam {len(self.background_regions)} bölge.")
            if self.rect_selector:
                self.rect_selector.set_active(False)

    def draw_background_regions(self):
        for patch in self.background_patches:
            try:patch.remove()
            except ValueError: pass
        self.background_patches.clear()
        for region in self.background_regions:
            xmin, xmax, ymin, ymax = region; width = xmax - xmin; height = ymax - ymin
            patch = patches.Rectangle((xmin, ymin), width, height, linewidth=1, edgecolor='red', facecolor='none', linestyle=':', label='Arkaplan Bölgesi', zorder=7)
            self.ax.add_patch(patch); self.background_patches.append(patch)
            if patch not in self.overlay_plots: self.overlay_plots.append(patch)

    def clear_background_regions(self):
        self.background_regions.clear()
        self.clear_aperture_overlays(clear_selection=False, clear_trace=False, clear_aperture=False, clear_background=True)
        print("Arkaplan bölgeleri temizlendi.")

    def calculate_background(self, image_data, background_regions, aperture_lower, aperture_upper):
        """Verilen bölgelerden basit medyan arkaplan modeli hesaplar."""
        if not background_regions: return None;
        if image_data is None or aperture_lower is None or aperture_upper is None: return None
        img_h, img_w = image_data.shape; background_1d = np.full(img_w, np.nan)
        for x in range(img_w):
            bg_pixels_in_col = []; ap_y_start = int(round(aperture_lower[x])); ap_y_end = int(round(aperture_upper[x]))
            for r_xmin, r_xmax, r_ymin, r_ymax in background_regions:
                if x >= r_xmin and x < r_xmax:
                    reg_y_start = max(0, r_ymin); reg_y_end = min(img_h, r_ymax); y_indices = np.arange(reg_y_start, reg_y_end)
                    mask = ~((y_indices >= ap_y_start) & (y_indices < ap_y_end)); valid_y = y_indices[mask]
                    if valid_y.size > 0: bg_pixels_in_col.extend(image_data[valid_y, x])
            if bg_pixels_in_col:
                with warnings.catch_warnings(): warnings.simplefilter("ignore", category=RuntimeWarning); background_1d[x] = np.nanmedian(bg_pixels_in_col)
        print("Arkaplan modeli hesaplandı."); return background_1d

    def extract_spectrum(self, image_data_to_extract=None, aperture_lower=None, aperture_upper=None, subtract_background=False):
        """Belirtilen açıklığı kullanarak 1D spektrumu çıkarır. Opsiyonel olarak arkaplan çıkarır."""
        data = image_data_to_extract if image_data_to_extract is not None else self.image_data; lower = aperture_lower if aperture_lower is not None else self.aperture_lower; upper = aperture_upper if aperture_upper is not None else self.aperture_upper
        if data is None or lower is None or upper is None: messagebox.showerror("Hata", "Ekstraksiyon için FITS ve tanımlı açıklık gerekli."); return None, None
        print("Spektrum çıkarılıyor...")
        img_h, img_w = data.shape; pixels_x = np.arange(img_w); flux_1d = np.zeros(img_w) * np.nan
        background_model = None
        if subtract_background:
            if self.background_regions: background_model = self.calculate_background(data, self.background_regions, lower, upper)
            else: print("Uyarı: Arkaplan bölgeleri tanımlanmadığı için çıkarma yapılamıyor.")
            if background_model is None: print("Arkaplan hesaplanamadı, çıkarma yapılmayacak."); subtract_background = False
            else: print("Arkaplan modeli kullanılacak.")
        for x in pixels_x:
            y_low = int(round(lower[x])); y_high = int(round(upper[x])); y_start = max(0, y_low); y_end = min(img_h, y_high)
            if y_start >= y_end: continue
            aperture_slice = data[y_start:y_end, x].copy()
            if subtract_background and background_model is not None and np.isfinite(background_model[x]): aperture_slice -= background_model[x]
            flux_1d[x] = np.nansum(aperture_slice)
        bg_status = "(Arkaplan Çıkarıldı)" if subtract_background else "(Arkaplan Çıkarılmadı)"; print(f"Ekstraksiyon tamamlandı {bg_status}.")
        if data is self.image_data: self.extracted_pixels = pixels_x; self.extracted_flux = flux_1d
        return pixels_x, flux_1d

    def show_1d_spectrum_window(self, x_axis_data, y_axis_data, title="1D Spektrum", x_axis_label="Piksel/Dalga Boyu", y_axis_label="Akı / Counts", x_error_data=None):
        if x_axis_data is None or y_axis_data is None: messagebox.showerror("Hata", "Gösterilecek 1D spektrum verisi yok.", parent=self.root); return
        is_science_spectrum = "Bilim Spektrumu" in title; is_arc_spectrum = "Ark Spektrumu" in title
        target_window = None
        if is_science_spectrum: target_window = self.science_plot_window
        elif is_arc_spectrum: target_window = self.arc_plot_window
        if target_window and target_window.winfo_exists():
            print(f"Mevcut '{title}' penceresi güncelleniyor.")
            target_window.update_plot(x_axis_data, y_axis_data, x_axis_label=x_axis_label, y_axis_label=y_axis_label, title=title, x_error_data=x_error_data); target_window.lift()
        else:
            print(f"Yeni '{title}' penceresi oluşturuluyor.")
            new_window = SpectrumPlotWindow(self.root, self, x_axis_data, y_axis_data, title=title, x_axis_label=x_axis_label, x_error_data=x_error_data)
            if is_science_spectrum: self.science_plot_window = new_window
            elif is_arc_spectrum: self.arc_plot_window = new_window

    def run_extraction_and_show(self):
        if self.image_data is None or self.aperture_lower is None: messagebox.showerror("Hata", "Lütfen önce FITS açın ve açıklığı tanımlayın.", parent=self.root); return
        subtract_bg = self.subtract_bg_var.get()
        if subtract_bg and not self.background_regions: messagebox.showwarning("Arkaplan Eksik", "Arkaplan çıkarmak için önce bölgeleri tanımlayın."); subtract_bg = False
        pixels, flux = self.extract_spectrum(subtract_background=subtract_bg)
        if pixels is not None:
            title_suffix = " (BG Çıkarıldı)" if subtract_bg else ""
            object_name = self.science_header.get('OBJECT', 'Bilim Spektrumu') if self.science_header else 'Bilim Spektrumu'
            base_title = f"{object_name}{title_suffix}"
            if self.dispersion_coefficients is not None: self.apply_wavelength_solution(show_plot=True, title_extra=title_suffix)
            else: self.show_1d_spectrum_window(pixels, flux, title=base_title, x_axis_label="Piksel (X)", y_axis_label="Akı / Counts", x_error_data=None)

    def load_arc_lamp(self):
        if self.aperture_lower is None or self.aperture_upper is None: messagebox.showwarning("Apertür Eksik", "Ark lambası spektrumunu çıkarmak için önce ana bilim verisi üzerinde açıklığı tanımlayın.", parent=self.root); return
        file_path = filedialog.askopenfilename(title="Ark Lambası FITS Seç", filetypes=(("FITS", "*.fits *.fit *.fts"), ("Hepsi", "*.*")))
        if not file_path: return
        try:
            with fits.open(file_path) as hdul:
                arc_hdu_index = -1;
                for i, hdu in enumerate(hdul):
                    if isinstance(hdu, (fits.ImageHDU, fits.PrimaryHDU)) and hasattr(hdu, 'data') and hdu.data is not None and hdu.data.ndim == 2: arc_hdu_index = i; break
                if arc_hdu_index == -1: raise ValueError("Dosyada uygun 2D görüntü verisi bulunamadı.")
                self.raw_arc_image_data = hdul[arc_hdu_index].data.astype(np.float32); file_name = os.path.basename(file_path); print(f"Ham ark verisi '{file_name}' okundu.")
                if self.raw_image_data is not None and self.raw_arc_image_data.shape != self.raw_image_data.shape: messagebox.showwarning("Boyut Uyuşmazlığı", f"Ark ({self.raw_arc_image_data.shape}) ve bilim verisi ({self.raw_image_data.shape}) boyutları farklı!")
                print("Ark verisine kalibrasyonlar uygulanıyor (varsa)..."); self.arc_image_data = self._apply_calibrations(self.raw_arc_image_data, is_science_frame=False)
                if self.arc_image_data is None: messagebox.showerror("Ark Kalibrasyon Hatası", "Ark verisine kalibrasyon uygulanırken hata oluştu.", parent=self.root); self.raw_arc_image_data = None; return
                print("Kalibre edilmiş ark spektrumu çıkarılıyor...")
                arc_pix, arc_flux = self.extract_spectrum(image_data_to_extract=self.arc_image_data, aperture_lower=self.aperture_lower, aperture_upper=self.aperture_upper, subtract_background=False)
                if arc_pix is not None:
                    self.arc_extracted_pixels = arc_pix; self.arc_extracted_flux = arc_flux; messagebox.showinfo("Başarılı", f"Ark lambası yüklendi, kalibre edildi ve spektrumu çıkarıldı.\nDosya: '{file_name}'", parent=self.root)
                    self.show_1d_spectrum_window(self.arc_extracted_pixels, self.arc_extracted_flux, title="Ark Spektrumu (Kalibre Edilmiş)", x_axis_label="Piksel (X)", y_axis_label="Akı / Counts")
                else: self.arc_image_data = None; self.raw_arc_image_data = None; messagebox.showerror("Ark Ekstraksiyon Hatası", "Ark spektrumu çıkarılamadı.", parent=self.root)
        except FileNotFoundError: messagebox.showerror("Hata", f"Dosya bulunamadı: {file_path}", parent=self.root); self._reset_arc_data()
        except Exception as e: messagebox.showerror("Ark Lambası Hatası", f"Ark lambası yüklenirken/işlenirken hata: {e}", parent=self.root); self._reset_arc_data()

    def _reset_arc_data(self):
        self.raw_arc_image_data = None; self.arc_image_data = None; self.arc_extracted_pixels = None; self.arc_extracted_flux = None
        if self.arc_plot_window and self.arc_plot_window.winfo_exists(): self.arc_plot_window.destroy(); self.arc_plot_window = None

    def load_line_list(self):
        file_path = filedialog.askopenfilename(title="Çizgi Listesi Seç (.txt, .dat)",
                                               filetypes=(("Metin", "*.txt *.dat"), ("Hepsi", "*.*")))
        if not file_path:
            return
        try:
            wavelengths = np.loadtxt(file_path, comments='#', usecols=0)
            if wavelengths.ndim != 1:
                raise ValueError(f"Çizgi listesi tek boyutlu bir dizi olarak okunamadı. Okunan boyut: {wavelengths.ndim}")
            self.line_list_data = np.sort(wavelengths)
            messagebox.showinfo("Başarılı",
                                 f"{len(self.line_list_data)} çizgi '{os.path.basename(file_path)}' dosyasından yüklendi.",
                                 parent=self.root)
            print(f"Yüklenen Çizgi Listesi (ilk 10): {self.line_list_data[:10]}")
        except ValueError as ve:
            messagebox.showerror("Çizgi Listesi Okuma Hatası",
                                 f"Dosya okunurken değer hatası oluştu: {ve}\n\n"
                                 "Lütfen dosyanın formatını kontrol edin. Dosyanın ilk sütunu sayısal dalga boyu değerleri içermeli "
                                 "ve yorum satırları '#' ile başlamalıdır.",
                                 parent=self.root)
            print(f"Çizgi listesi yüklenirken ValueError: {ve}")
            self.line_list_data = None
        except Exception as e:
            messagebox.showerror("Çizgi Listesi Hatası",
                                 f"Çizgi listesi yüklenirken genel bir hata oluştu: {e}",
                                 parent=self.root)
            print(f"Çizgi listesi yüklenirken Exception: {e}")
            self.line_list_data = None

    def start_wavelength_calibration(self):
        if self.arc_extracted_flux is None or self.arc_extracted_pixels is None: messagebox.showerror("Hata", "Kalibrasyon için önce ark lambası spektrumunu yükleyip çıkarın.", parent=self.root); return
        if self.line_list_data is None: messagebox.showwarning("Uyarı", "Referans çizgi listesi yüklenmedi. Çizgi önerileri çalışmayacak.", parent=self.root)
        initial_covariance = self.dispersion_covariance
        calib_window = WavelengthCalibWindow(self.root, self, self.arc_extracted_pixels, self.arc_extracted_flux, self.line_list_data)
        calib_window.fitted_covariance = initial_covariance
        self.root.wait_window(calib_window)
        if hasattr(calib_window, 'applied') and calib_window.applied:
            if calib_window.fitted_coefficients is not None:
                self.dispersion_coefficients = calib_window.fitted_coefficients
                self.dispersion_covariance = calib_window.fitted_covariance
                self.calib_points = list(calib_window.calib_points_internal)
                self.wl_fit_degree = calib_window.fit_degree.get()
                print(f"Dalga boyu çözümü güncellendi (Derece {self.wl_fit_degree}): {self.dispersion_coefficients}")
                if self.dispersion_covariance is not None: print("Kovaryans matrisi alındı.")
                else: print("Uyarı: Kovaryans matrisi alınamadı.")
                if self.extracted_flux is not None: self.apply_wavelength_solution(show_plot=True)
                else: messagebox.showinfo("Çözüm Güncellendi", "Dalga boyu çözümü başarıyla güncellendi. Bilim spektrumu çıkarıldığında uygulanacak.", parent=self.root)
            else: print("Kalibrasyon penceresi 'Uygula' ile kapatıldı ama geçerli katsayı bulunamadı."); self.dispersion_coefficients = None; self.dispersion_covariance = None; self.science_wavelength_error = None
        else: print("Dalga boyu kalibrasyon penceresi 'İptal' ile veya kapatma tuşuyla kapatıldı.")

    # --- apply_wavelength_solution metodu (DEĞİŞİKLİK) ---
    def apply_wavelength_solution(self, show_plot=False, title_extra=""):
        if self.science_normalized_flux is not None:
            print("Dalga boyu çözümü uygulanıyor, önceki normalizasyon temizlendi.")
            self.science_normalized_flux = None
            self.science_continuum_flux = None

        if self.extracted_pixels is None or self.extracted_flux is None:
            print("Uygulanacak bilim spektrumu yok.")
            return

        self.science_wavelength_error = None # Hata dizisini sıfırla

        object_name = self.science_header.get('OBJECT', 'Bilim Spektrumu') if self.science_header else 'Bilim Spektrumu'
        base_title = f"{object_name}{title_extra}"

        if self.dispersion_coefficients is None:
            print("Uygulanacak dalga boyu çözümü yok.")
            if show_plot:
                self.show_1d_spectrum_window(self.extracted_pixels, self.extracted_flux,
                                            title=base_title, x_axis_label="Piksel (X)", y_axis_label="Akı / Counts")
            return

        try:
            disp_func = np.poly1d(self.dispersion_coefficients)
            wavelengths = disp_func(self.extracted_pixels)
            print("Dalga boyu çözümü bilim spektrumuna uygulandı.")

            final_x_label = "Dalga Boyu" # Varsayılan olarak "Dalga Boyu"
            unit_str = "" # Birim eklenecekse buraya eklenecek

            # Birim bilgisini header'dan almaya çalış
            if self.science_header:
                for key in ['CUNIT1', 'CUNIT2', 'CDELT1', 'CDELT2', 'CRVAL1', 'WAT1_001']: # Daha fazla anahtar kelime kontrol edilebilir
                    header_val = str(self.science_header.get(key, '')).strip().lower() # String'e çevir, boşlukları temizle, küçük harf yap

                    if 'angstrom' in header_val:
                        unit_str = "Angstrom"
                        break
                    elif 'nm' in header_val or 'nanometer' in header_val:
                        unit_str = "nm"
                        break
                    elif 'um' in header_val or 'micron' in header_val:
                        unit_str = r"$\mu$m" # LaTeX formatında mikron sembolü
                        break
                    elif 'wavelength' in header_val and 'unit' in header_val:
                        # Eğer 'WAT1_001' gibi bir keyword'de "units=..." varsa
                        parts = header_val.split('=')
                        if len(parts) > 1:
                            potential_unit = parts[-1].strip().replace('/', '') # "/A" gibi durumlarda "/" kaldır
                            if potential_unit:
                                unit_str = potential_unit.capitalize() # "Angstrom" veya "NM" olabilir
                                break

            if unit_str:
                final_x_label = f"Dalga Boyu ({unit_str})"
            else:
                final_x_label = "Dalga Boyu (Birim Bilinmiyor)" # Birim bulunamazsa

            # DEĞİŞİKLİK: Dalgaboyu hatalarını hesapla
            if self.dispersion_covariance is not None:
                print("Dalgaboyu hataları hesaplanıyor (kovaryans matrisinden)...")
                try:
                    degree = len(self.dispersion_coefficients) - 1
                    # Vandermonde matrisi, polinom fitindeki katsayıların nasıl değiştiğini açıklar
                    vandermonde_matrix = np.vander(self.extracted_pixels, N=degree + 1)
                    # Hata yayılımı formülü: Var(y) = X * Cov(beta) * X.T
                    # Burada y = f(x) = P(x) = sum(beta_i * x^i)
                    # X: [1, x, x^2, ..., x^degree]
                    # Cov(beta): katsayıların kovaryans matrisi
                    # Karekökünü alarak standart sapmayı buluruz
                    variance_wavelength = np.sum((vandermonde_matrix @ self.dispersion_covariance) * vandermonde_matrix, axis=1)
                    # Negatif varyans olmaması için sıfırdan küçükleri nan yap
                    variance_wavelength[variance_wavelength < 0] = np.nan
                    self.science_wavelength_error = np.sqrt(variance_wavelength)
                    print(f"Dalgaboyu hataları hesaplandı (örnek ilk 5: {self.science_wavelength_error[:5] if len(self.science_wavelength_error) > 5 else self.science_wavelength_error})")
                except Exception as e_err:
                    print(f"Hata: Dalgaboyu hatası hesaplanırken sorun oluştu: {e_err}")
                    import traceback
                    traceback.print_exc()
                    self.science_wavelength_error = None
            else:
                print("Uyarı: Kovaryans matrisi bulunamadığı için dalgaboyu hatası hesaplanamadı.")

            if show_plot:
                error_data_to_plot = None
                if self.show_wl_error_var.get() and self.science_wavelength_error is not None:
                    # Eğer hata verisinde çok fazla NaN veya inf varsa, çizim sorun olabilir
                    if np.sum(np.isfinite(self.science_wavelength_error)) / len(self.science_wavelength_error) > 0.5: # Yarıdan fazlası geçerliyse çiz
                         error_data_to_plot = self.science_wavelength_error
                    else:
                        print("Uyarı: Dalgaboyu hata verisi çoğunlukla geçersiz, hata çubukları çizilmiyor.")

                self.show_1d_spectrum_window(wavelengths, self.extracted_flux,
                                             title=f"{base_title} (Dalga Boyu Kalibreli)",
                                             x_axis_label=final_x_label,
                                             y_axis_label="Akı / Counts",
                                             x_error_data=error_data_to_plot)
        except Exception as e:
            messagebox.showerror("Dalga Boyu Uygulama Hatası", f"Hata: {e}")
            print(f"Dalga boyu uygulama hatası: {e}")
            import traceback
            traceback.print_exc()

    def _fit_and_normalize_spectrum(self, x_data, y_data):
        if x_data is None or y_data is None: return None, None
        if len(x_data) != len(y_data): return None, None
        finite_mask = np.isfinite(y_data)
        if not np.any(finite_mask): print("Hata: Akı verisinde sonlu değer yok."); return None, None
        x_finite = x_data[finite_mask]; y_finite = y_data[finite_mask]
        fit_type = self.continuum_fit_type.get(); degree = self.continuum_fit_degree
        if fit_type == "Polynomial" and len(x_finite) <= degree: print(f"Hata: Polinom Fit için yetersiz sonlu nokta ({len(x_finite)}), derece {degree} gerekli."); return None, None
        if fit_type == "Spline" and len(x_finite) <= 3: print(f"Hata: Spline Fit (k=3) için yetersiz sonlu nokta ({len(x_finite)}), en az 4 gerekli."); return None, None
        print(f"{fit_type} ile süreklilik fit ediliyor (sigma clipping ile)...")
        try:
            clipped_y = y_finite.copy(); current_mask = np.zeros_like(y_finite, dtype=bool); fit_params = None
            for i in range(5):
                if np.all(current_mask): print("Uyarı: Sigma clip sırasında tüm noktalar maskelendi."); break
                x_to_fit = x_finite[~current_mask]; y_to_fit = clipped_y[~current_mask]
                current_points = len(x_to_fit)
                if (fit_type == "Polynomial" and current_points <= degree) or \
                   (fit_type == "Spline" and current_points <= 3): print(f"Uyarı: Sigma clip iterasyon {i+1}'de fit için yetersiz nokta ({current_points})."); break
                continuum_est_iter = None
                try:
                    if fit_type == "Polynomial": fit_params = np.polyfit(x_to_fit, y_to_fit, degree); p = np.poly1d(fit_params); continuum_est_iter = p(x_finite)
                    elif fit_type == "Spline":
                        if not SCIPY_AVAILABLE: raise RuntimeError("Spline fit için SciPy kütüphanesi gerekli.")
                        fit_params = splrep(x_to_fit, y_to_fit, s=len(y_to_fit), k=3); continuum_est_iter = splev(x_finite, fit_params)
                except (np.linalg.LinAlgError, ValueError, RuntimeError, TypeError) as fit_err: print(f"Uyarı: Sigma clip iterasyon {i+1}'de {fit_type} fit hatası: {fit_err}. İterasyon durduruldu."); break
                if continuum_est_iter is None: break
                residuals = y_finite - continuum_est_iter; std_dev = np.nanstd(residuals[~current_mask])
                if not np.isfinite(std_dev) or std_dev == 0: print(f"Uyarı: Sigma clip iterasyon {i+1}'de std sapması geçersiz ({std_dev}). İterasyon durduruldu."); break
                new_outliers = np.abs(residuals) > (3.0 * std_dev); new_mask = current_mask | new_outliers
                if np.all(new_mask == current_mask): print(f"Sigma clip iterasyon {i+1}'de yakınsadı."); break
                current_mask = new_mask; print(f"Sigma clip iterasyon {i+1} tamamlandı, {np.sum(current_mask)} nokta maskelendi.")
            else: print("Sigma clip maksimum iterasyona ulaştı.")
            if fit_params is None: raise ValueError("Sigma clipping sonrası geçerli fit yapılamadı.")
            final_continuum = None
            if fit_type == "Polynomial": p = np.poly1d(fit_params); final_continuum = p(x_data)
            elif fit_type == "Spline": final_continuum = splev(x_data, fit_params)
            if final_continuum is None: raise ValueError("Final süreklilik hesaplanamadı.")
            safe_continuum = final_continuum.copy(); zero_mask = np.isclose(safe_continuum, 0) | ~np.isfinite(safe_continuum)
            if np.any(zero_mask): print(f"Uyarı: Fit edilen süreklilikte {np.sum(zero_mask)} sıfır/geçersiz değer bulundu. Bu noktalarda normalize akı NaN olacak."); safe_continuum[zero_mask] = np.nan
            normalized_flux = y_data / safe_continuum
            print("Süreklilik fit ve normalizasyon tamamlandı.")
            return normalized_flux.astype(np.float32), final_continuum.astype(np.float32)
        except Exception as e: messagebox.showerror("Normalizasyon Hatası", f"Süreklilik fit/normalizasyon sırasında hata: {e}", parent=self.root); print(f"Continuum fit/normalization error: {e}"); return None, None

    def run_continuum_normalization(self):
        print("Süreklilik düzeltmesi başlatılıyor...")
        target_x_for_norm = None
        target_y_for_norm = None
        final_x_label = None
        object_name = self.science_header.get('OBJECT', 'Bilim Spektrumu') if self.science_header else 'Bilim Spektrumu'
        bg_suffix = " (BG Çıkarıldı)" if self.subtract_bg_var.get() and self.background_regions else ""
        is_wavelength_calibrated = False
        error_data_to_plot = None

        if self.dispersion_coefficients is not None and self.extracted_pixels is not None and self.extracted_flux is not None:
            print("Dalga boyu kalibreli veri üzerinden normalizasyon yapılacak (orijinal akı kullanılıyor).")
            disp_func = np.poly1d(self.dispersion_coefficients)
            target_x_for_norm = disp_func(self.extracted_pixels)
            target_y_for_norm = self.extracted_flux

            unit_str = "Birim?"
            if self.science_header:
                for key in ['CUNIT1', 'CUNIT2', 'WAT1_001_UNITS', 'WAT2_001_UNITS']:
                    header_val_candidate = self.science_header.get(key, '').strip().lower()
                    if 'angstrom' in header_val_candidate: unit_str = "Angstrom"; break
                    elif 'nm' in header_val_candidate or 'nanometer' in header_val_candidate: unit_str = "nm"; break
                    elif 'um' in header_val_candidate or 'micron' in header_val_candidate: unit_str = r"$\mu$m"; break
            final_x_label = f"Dalga Boyu (Angstrom)"
            is_wavelength_calibrated = True
            if self.show_wl_error_var.get() and self.science_wavelength_error is not None and len(self.science_wavelength_error) == len(target_x_for_norm):
                error_data_to_plot = self.science_wavelength_error

        elif self.extracted_pixels is not None and self.extracted_flux is not None:
            print("Hafızadaki piksel bazlı orijinal bilim spektrumu normalizasyon için kullanılıyor.")
            target_x_for_norm = self.extracted_pixels
            target_y_for_norm = self.extracted_flux
            final_x_label = "Piksel (X)"
            is_wavelength_calibrated = False
        else:
            messagebox.showerror("Veri Yok", "Normalize edilecek bilim spektrumu bulunamadı.\nLütfen önce spektrumu çıkarın ('Ekstraksiyon' menüsü).", parent=self.root)
            return

        normalized_flux, continuum_flux = self._fit_and_normalize_spectrum(target_x_for_norm, target_y_for_norm)

        if normalized_flux is None:
            return

        self.science_normalized_flux = normalized_flux.astype(np.float32)
        self.science_continuum_flux = continuum_flux.astype(np.float32)

        final_plot_title = f"{object_name}{bg_suffix} (Normalize Edilmiş)"
        final_y_label = "Normalize Akı"

        self.show_1d_spectrum_window(target_x_for_norm,
                                     self.science_normalized_flux,
                                     title=final_plot_title,
                                     x_axis_label=final_x_label,
                                     y_axis_label=final_y_label,
                                     x_error_data=error_data_to_plot if is_wavelength_calibrated else None)

        messagebox.showinfo("Başarılı", "Süreklilik düzeltmesi uygulandı ve grafik güncellendi.", parent=self.root)

    def set_continuum_degree_dialog(self):
        new_degree = simpledialog.askinteger("Süreklilik Fit Derecesi", "Süreklilik (continuum) fit için polinom derecesini girin:", parent=self.root, initialvalue=self.continuum_fit_degree, minvalue=0)
        if new_degree is not None and new_degree != self.continuum_fit_degree:
            self.continuum_fit_degree = new_degree;
            self.continuum_degree_label_var.set(f"Süreklilik Fit Derecesi: {self.continuum_fit_degree}");
            print(f"Süreklilik fit derecesi {self.continuum_fit_degree} olarak ayarlandı.")
            self.science_normalized_flux = None; self.science_continuum_flux = None

    def on_continuum_fit_type_change(self):
        fit_type = self.continuum_fit_type.get(); print(f"Süreklilik fit tipi '{fit_type}' olarak değiştirildi.")
        self.science_normalized_flux = None; self.science_continuum_flux = None

    def save_normalized_spectrum(self):
        print("Normalize spektrumu kaydetme işlemi başlatıldı...")
        if self.science_normalized_flux is None: messagebox.showerror("Veri Yok", "Kaydedilecek normalize edilmiş spektrum verisi bulunamadı.\nLütfen önce 'Süreklilik Düzeltmesi Yap' işlemini çalıştırın.", parent=self.root); return
        x_data = None; x_label = "Bilinmeyen_Eksen"; wl_calibrated = False
        window_to_check = self.science_plot_window
        if window_to_check and window_to_check.winfo_exists() and "Dalga Boyu" in window_to_check.ax.get_xlabel():
            line = window_to_check.line
            if line is None and window_to_check.errorbar_container: line = window_to_check.errorbar_container.lines[0]
            if line: x_data = line.get_xdata(); x_label = "Dalgaboyu"; wl_calibrated = True; print(f"X ekseni olarak dalga boyu verisi kullanılacak ({len(x_data) if x_data is not None else 0} nokta).")
        elif self.extracted_pixels is not None:
            x_data = self.extracted_pixels; x_label = "Piksel"; print(f"X ekseni olarak piksel verisi kullanılacak ({len(x_data)} nokta).")
        if x_data is None: messagebox.showerror("Hata", "Spektrumun X ekseni verisi (piksel veya dalgaboyu) bulunamadı.", parent=self.root); return
        norm_flux = self.science_normalized_flux; orig_flux = self.extracted_flux; continuum = self.science_continuum_flux
        wl_error = self.science_wavelength_error if wl_calibrated else None
        data_arrays = {'x': x_data, 'norm': norm_flux, 'err': wl_error, 'orig': orig_flux, 'cont': continuum}
        valid_arrays = {k: v for k, v in data_arrays.items() if v is not None}
        if not all(len(arr) == len(x_data) for arr in valid_arrays.values()): messagebox.showerror("Veri Uyuşmazlığı", "Kaydedilecek veri sütunlarının uzunlukları eşleşmiyor.", parent=self.root); return
        default_base = 'spectrum'
        if self.science_header and 'OBJECT' in self.science_header: default_base = self.science_header['OBJECT'].replace(' ', '_')
        elif self.hdulist and self.hdulist.filename(): default_base = os.path.splitext(os.path.basename(self.hdulist.filename()))[0]
        initial_file = f"{default_base}_norm.txt"
        filename = filedialog.asksaveasfilename(title="Normalize Spektrumu Kaydet", initialfile=initial_file, defaultextension=".txt", filetypes=[("Metin Dosyası", "*.txt"), ("ASCII Dosyası", "*.asc"), ("Veri Dosyası", "*.dat"), ("Tüm Dosyalar", "*.*")])
        if not filename: print("Kaydetme işlemi iptal edildi."); return
        columns_to_save = [valid_arrays['x'], valid_arrays['norm']]; header_parts = [x_label, "Normalize_Aki"]
        if 'err' in valid_arrays: columns_to_save.insert(1, valid_arrays['err']); header_parts.insert(1, "Hata_" + x_label)
        if 'orig' in valid_arrays: columns_to_save.append(valid_arrays['orig']); header_parts.append("Orijinal_Aki")
        if 'cont' in valid_arrays: columns_to_save.append(valid_arrays['cont']); header_parts.append("Sureklilik")
        data_to_save = np.column_stack(columns_to_save); header_string = " ".join(header_parts)
        try:
            np.savetxt(filename, data_to_save, fmt="%.6f", header=header_string, comments='')
            messagebox.showinfo("Başarılı", f"Normalize spektrum başarıyla kaydedildi:\n{filename}", parent=self.root)
            print(f"Normalize spektrum şuraya kaydedildi: {filename}")
        except Exception as e: messagebox.showerror("Yazma Hatası", f"Dosya yazılırken hata oluştu:\n{e}", parent=self.root); print(f"Error saving spectrum to file: {e}")

    def save_wavelength_solution(self):
        if self.dispersion_coefficients is None: messagebox.showwarning("Çözüm Yok", "Kaydedilecek dalgaboyu kalibrasyon çözümü bulunmuyor.", parent=self.root); return
        filename = filedialog.asksaveasfilename(title="Dalgaboyu Çözümünü Kaydet", defaultextension=".json", filetypes=[("JSON Dosyası", "*.json"), ("Tüm Dosyalar", "*.*")])
        if not filename: print("DB Çözümü kaydetme işlemi iptal edildi."); return
        solution_data = {'coefficients': list(self.dispersion_coefficients), 'degree': self.wl_fit_degree, 'points': self.calib_points, 'covariance_matrix': self.dispersion_covariance.tolist() if self.dispersion_covariance is not None else None }
        try:
            with open(filename, 'w') as f: json.dump(solution_data, f, indent=4)
            messagebox.showinfo("Başarılı", f"Dalgaboyu çözümü kaydedildi:\n{filename}", parent=self.root); print(f"Dalgaboyu çözümü şuraya kaydedildi: {filename}")
        except Exception as e: messagebox.showerror("Yazma Hatası", f"Çözüm dosyası yazılırken hata oluştu:\n{e}", parent=self.root); print(f"Error saving wavelength solution: {e}")

    def load_wavelength_solution(self):
        filename = filedialog.askopenfilename(title="Dalgaboyu Çözümünü Yükle", filetypes=[("JSON Dosyası", "*.json"), ("Tüm Dosyalar", "*.*")])
        if not filename: print("DB Çözümü yükleme işlemi iptal edildi."); return
        try:
            with open(filename, 'r') as f: solution_data = json.load(f)
            if not all(k in solution_data for k in ['coefficients', 'degree', 'points']): raise ValueError("Dosya beklenen formatta değil (coefficients, degree, points eksik).")
            self.dispersion_coefficients = np.array(solution_data['coefficients']); self.wl_fit_degree = int(solution_data['degree']); self.calib_points = solution_data['points']
            if 'covariance_matrix' in solution_data and solution_data['covariance_matrix'] is not None: self.dispersion_covariance = np.array(solution_data['covariance_matrix'])
            else: self.dispersion_covariance = None
            self.science_wavelength_error = None
            print(f"Dalgaboyu çözümü '{filename}' dosyasından yüklendi."); print(f"  Katsayılar: {self.dispersion_coefficients}"); print(f"  Derece: {self.wl_fit_degree}"); print(f"  Kalibrasyon Noktaları: {len(self.calib_points)} adet")
            if self.dispersion_covariance is not None: print("  Kovaryans matrisi yüklendi.")
            if self.extracted_flux is not None:
                if messagebox.askyesno("Çözümü Uygula", "Yüklenen dalgaboyu çözümünü mevcut bilim spektrumuna uygulamak ister misiniz?", parent=self.root): self.apply_wavelength_solution(show_plot=True)
            else: messagebox.showinfo("Çözüm Yüklendi", "Dalgaboyu çözümü yüklendi. Bir bilim spektrumu çıkardığınızda uygulanacaktır.", parent=self.root)
        except FileNotFoundError: messagebox.showerror("Hata", f"Dosya bulunamadı: {filename}", parent=self.root)
        except json.JSONDecodeError: messagebox.showerror("Hata", f"'{filename}' geçerli bir JSON dosyası değil.", parent=self.root)
        except ValueError as ve: messagebox.showerror("Format Hatası", f"Dosya formatı hatalı: {ve}", parent=self.root)
        except Exception as e: messagebox.showerror("Yükleme Hatası", f"Çözüm dosyası okunurken hata oluştu:\n{e}", parent=self.root); print(f"Error loading wavelength solution: {e}")

    # --- toggle_errorbar_display metodu (DEĞİŞİKLİK) ---
    def toggle_errorbar_display(self):
        print(f"Dalgaboyu hata çubukları gösterilsin mi: {self.show_wl_error_var.get()}")
        if self.dispersion_coefficients is not None and self.science_plot_window and self.science_plot_window.winfo_exists():
            print("Grafik hata çubuğu durumu ile güncelleniyor...")
            bg_suffix = " (BG Çıkarıldı)" if self.subtract_bg_var.get() else ""
            self.apply_wavelength_solution(show_plot=True, title_extra=bg_suffix)
        elif self.dispersion_coefficients is not None and self.extracted_flux is not None:
            print("Bilim spektrum penceresi açık değil, hata çubuğu durumu değiştiğinde gösterilecek.")
        else:
            messagebox.showinfo("Bilgi", "Hata çubukları için yüklü bir dalgaboyu çözümü veya çıkarılmış spektrum yok.", parent=self.root)

    def start_ew_calculation_mode(self):
        print("Eşdeğer Genişlik hesaplama modu başlatılıyor...")
        if not (self.science_plot_window and self.science_plot_window.winfo_exists()):
            messagebox.showwarning("Grafik Yok", "Eşdeğer genişlik hesaplamak için önce bir bilim spektrumu grafiği açık olmalıdır.", parent=self.root)
            return

        if self.science_normalized_flux is None:
            if not messagebox.askyesno("Uyarı: Normalize Edilmemiş Spektrum",
                                       "Spektrum normalize edilmemiş görünüyor. Eşdeğer genişlik genellikle normalize edilmiş spektrumlar (süreklilik ~1) için daha anlamlıdır.\n"
                                       "Yine de devam etmek ve seçilen bölgenin uç noktalarını süreklilik olarak kullanmak ister misiniz?",
                                       parent=self.root):
                return
            print("Uyarı: Normalize edilmemiş spektrum üzerinden EW hesaplanacak.")
        else:
            print("Normalize edilmiş spektrum üzerinden EW hesaplanacak.")

        if self.line_fit_selector and hasattr(self.line_fit_selector, 'active') and self.line_fit_selector.active:
            try: self.line_fit_selector.set_active(False)
            except: pass
        if self.ew_selector and hasattr(self.ew_selector, 'active') and self.ew_selector.active:
            try: self.ew_selector.set_active(False)
            except: pass

        messagebox.showinfo("Çizgi Bölgesi Seçimi",
                             "Eşdeğer genişliği hesaplanacak çizginin bulunduğu bölgeyi grafikte fare ile seçin (başlangıçtan bitişe sürükleyin).",
                             parent=self.science_plot_window)

        self.science_plot_window.canvas_widget.focus_set()
        self.ew_selector = SpanSelector(self.science_plot_window.ax, self._on_ew_region_selected,
                                        'horizontal', useblit=True,
                                        props=dict(alpha=0.3, facecolor='lightcoral'),
                                        span_stays=True, minspan=3)
        print("Eşdeğer genişlik için bölge seçimi aktif.")

    def _on_ew_region_selected(self, x_min_selected, x_max_selected):
        if not (self.science_plot_window and self.science_plot_window.winfo_exists()):
            print("Hata: EW bölgesi seçildi ancak bilim spektrumu penceresi bulunamadı.")
            if self.ew_selector: self.ew_selector.set_active(False)
            return

        print(f"EW için bölge seçildi: X_min={x_min_selected:.4f}, X_max={x_max_selected:.4f}")

        if self.ew_selector:
            self.ew_selector.set_active(False)
            if hasattr(self.ew_selector, 'rect') and self.ew_selector.rect is not None:
                try:
                    self.ew_selector.rect.remove()
                    self.science_plot_window.canvas.draw_idle()
                except Exception as e_remove:
                    print(f"EW seçici dikdörtgeni silinirken hata: {e_remove}")

        line_obj = self.science_plot_window.line
        if line_obj is None and self.science_plot_window.errorbar_container:
            line_obj = self.science_plot_window.errorbar_container.lines[0]

        if line_obj is None:
            messagebox.showerror("Veri Hatası", "Grafikte çizgi verisi bulunamadı.", parent=self.science_plot_window)
            return

        line_data_x = line_obj.get_xdata()
        line_data_y = line_obj.get_ydata()

        selected_indices = np.where((line_data_x >= x_min_selected) & (line_data_x <= x_max_selected))[0]
        if len(selected_indices) < 2:
            messagebox.showerror("Yetersiz Nokta", "Eşdeğer genişlik hesaplamak için seçilen bölgede yeterli nokta (en az 2) bulunmuyor.", parent=self.science_plot_window)
            return

        x_region = line_data_x[selected_indices]
        y_region = line_data_y[selected_indices]

        continuum_level_region = None
        if self.science_normalized_flux is not None:
            continuum_level_region = np.ones_like(y_region)
            print("Normalize spektrumda süreklilik 1.0 kabul edildi.")
        else:
            if len(x_region) >= 2:
                m = (y_region[-1] - y_region[0]) / (x_region[-1] - x_region[0]) if (x_region[-1] - x_region[0]) != 0 else 0
                c = y_region[0] - m * x_region[0]
                continuum_level_region = m * x_region + c
                print(f"Normalize edilmemiş spektrum için yerel süreklilik (doğrusal fit): y = {m:.3f}x + {c:.3f}")
            else:
                messagebox.showerror("Hata", "Yerel süreklilik için yetersiz nokta.", parent=self.science_plot_window)
                return

        line_type_choice = simpledialog.askstring("Çizgi Tipi",
                                                 "Çizgi tipi nedir? ('absorpsiyon' veya 'emisyon' yazın):",
                                                 initialvalue='absorpsiyon', parent=self.science_plot_window)
        if line_type_choice:
            line_type_choice = line_type_choice.strip().lower()
        else:
            return

        if line_type_choice not in ['absorpsiyon', 'emisyon']:
            messagebox.showwarning("Geçersiz Tip", "Lütfen 'absorpsiyon' veya 'emisyon' girin.", parent=self.science_plot_window)
            return

        ew, unit = self._calculate_equivalent_width(x_region, y_region, continuum_level_region, line_type_choice)

        if ew is not None:
            messagebox.showinfo("Eşdeğer Genişlik Sonucu",
                                 f"Seçilen Bölge için Hesaplanan Eşdeğer Genişlik:\nEW = {ew:.4f} {unit}",
                                 parent=self.science_plot_window)
        else:
            messagebox.showerror("Hesaplama Hatası", "Eşdeğer genişlik hesaplanamadı.", parent=self.science_plot_window)

    def _calculate_equivalent_width(self, x_values, flux_values, continuum_values, line_type='absorpsiyon'):
        """
        Verilen x, akı ve süreklilik değerleri için eşdeğer genişliği hesaplar.
        x_values: Dalga boyu veya piksel değerleri dizisi.
        flux_values: Karşılık gelen akı değerleri dizisi.
        continuum_values: Her bir x_value için süreklilik seviyesi dizisi veya tek bir değer.
        line_type: 'absorpsiyon' veya 'emisyon'.
        """
        if len(x_values) < 2 or len(flux_values) < 2:
            print("EW hesaplamak için yetersiz veri noktası.")
            return None, ""

        if np.isscalar(continuum_values):
            continuum_values = np.full_like(flux_values, continuum_values)
        elif len(continuum_values) != len(flux_values):
            print("Hata: Akı ve süreklilik dizilerinin boyutları uyuşmuyor.")
            return None, ""

        integrand = None
        if line_type == 'absorpsiyon':
            safe_continuum = np.where(np.abs(continuum_values) < 1e-9, 1e-9, continuum_values)
            integrand = (safe_continuum - flux_values) / safe_continuum
        elif line_type == 'emisyon':
            safe_continuum = np.where(np.abs(continuum_values) < 1e-9, 1e-9, continuum_values)
            integrand = (flux_values - safe_continuum) / safe_continuum
        else:
            return None, ""

        ew = np.trapz(integrand, x_values)

        unit = ""
        if self.science_plot_window and self.science_plot_window.winfo_exists():
            xlabel_text = self.science_plot_window.ax.get_xlabel().lower()
            if "dalgaboyu" in xlabel_text or "wavelength" in xlabel_text:
                if "angstrom" in xlabel_text: unit = "Angstrom"
                elif "nm" in xlabel_text: unit = "nm"
                elif r"$\mu$m" in xlabel_text or "micron" in xlabel_text: unit = "µm"
                else: unit = "(Dalgaboyu Birimi)"
            elif "piksel" in xlabel_text or "pixel" in xlabel_text:
                unit = "Piksel"

        return ew, unit

    def start_line_fitting_mode(self):
        if not (self.science_plot_window and self.science_plot_window.winfo_exists()): messagebox.showwarning("Grafik Yok", "Çizgi fiti yapmak için önce bir bilim spektrumu grafiği açık olmalıdır.", parent=self.root); return
        if not SCIPY_AVAILABLE: messagebox.showerror("SciPy Eksik", "Gaussian çizgi fiti için SciPy kütüphanesi gereklidir.\nLütfen 'pip install scipy' ile yükleyin.", parent=self.root); return
        if self.line_fit_selector and hasattr(self.line_fit_selector, 'active') and self.line_fit_selector.active:
            try: self.line_fit_selector.set_active(False); print("Önceki çizgi fit seçici devre dışı bırakıldı.")
            except: pass
        messagebox.showinfo("Çizgi Seçimi", "Fit edilecek çizginin bulunduğu bölgeyi grafikte fare ile seçin (başlangıçtan bitişe sürükleyin).", parent=self.science_plot_window)
        self.science_plot_window.canvas_widget.focus_set()
        self.line_fit_selector = SpanSelector(self.science_plot_window.ax, self._on_line_region_selected, 'horizontal', useblit=False, props=dict(alpha=0.3, facecolor='skyblue'), span_stays=True, minspan=3)
        print("Çizgi fiti için bölge seçimi aktif.")

    def _on_line_region_selected(self, xmin, xmax):
        print(f"Çizgi bölgesi seçildi: X_min={xmin:.2f}, X_max={xmax:.2f}")
        if self.line_fit_selector:
            try: self.line_fit_selector.set_active(False)
            except Exception as e: print(f"SpanSelector devre dışı bırakılırken hata (önemsiz): {e}")
        if not (self.science_plot_window and self.science_plot_window.winfo_exists()): print("Hata: Bilim spektrumu penceresi bulunamadı."); return
        line_obj = self.science_plot_window.line
        if line_obj is None and self.science_plot_window.errorbar_container: line_obj = self.science_plot_window.errorbar_container.lines[0]
        if line_obj is None: messagebox.showerror("Hata", "Grafikte çizgi verisi bulunamadı."); return
        line_data_x = line_obj.get_xdata(); line_data_y = line_obj.get_ydata()
        try:
            indices = np.where((line_data_x >= xmin) & (line_data_x <= xmax))[0]
            if len(indices) < 4: messagebox.showerror("Yetersiz Nokta", "Gaussian fiti için seçilen bölgede yeterli nokta (en az 4) bulunmuyor.", parent=self.science_plot_window); return
            x_selected = line_data_x[indices]; y_selected = line_data_y[indices]
        except Exception as e: messagebox.showerror("Veri Hatası", f"Seçilen bölgeden veri alınırken hata: {e}", parent=self.science_plot_window); return
        fit_params, fitted_curve_y = self._perform_gaussian_fit(x_selected, y_selected)
        if fit_params: LineFitResultsDialog(self.root, fit_params, x_selected, y_selected, fitted_curve_y)
        else: messagebox.showerror("Fit Hatası", "Seçilen bölgeye Gaussian fiti yapılamadı.\nDaha dar veya daha uygun bir bölge seçmeyi deneyin.", parent=self.science_plot_window)

    def _perform_gaussian_fit(self, x, y):
        if not SCIPY_AVAILABLE: print("SciPy bulunamadığı için Gaussian fit yapılamıyor."); return None, None
        try:
            offset_guess = np.nanmin(y); y_shifted = y - offset_guess
            amplitude_guess = np.nanmax(y_shifted)
            if not np.isfinite(amplitude_guess) or amplitude_guess <= 0: amplitude_guess = 1.0
            mean_guess = x[np.nanargmax(y_shifted)]; stddev_guess = (x[-1] - x[0]) / 6.0
            if stddev_guess <= 1e-6: stddev_guess = 1.0
            p0 = [amplitude_guess, mean_guess, stddev_guess, offset_guess]
            print(f"[DEBUG] Gaussian initial guess (amp, mean, std, off): {p0}")
            bounds = ([-np.inf, np.min(x), 1e-4, -np.inf], [np.inf, np.max(x), (x[-1]-x[0]), np.inf])
            popt, pcov = curve_fit(gaussian_with_offset, x, y, p0=p0, bounds=bounds, maxfev=10000, ftol=1e-5, xtol=1e-5)
            amplitude, mean, stddev, offset = popt
            fwhm = 2 * np.sqrt(2 * np.log(2)) * np.abs(stddev)
            print(f"[DEBUG] Gaussian fit params: Amp={amplitude:.3f}, Mean={mean:.3f}, StdDev={stddev:.3f}, Offset={offset:.3f}, FWHM={fwhm:.3f}")
            fitted_curve_y = gaussian_with_offset(x, *popt)
            return (amplitude, mean, stddev, offset, fwhm), fitted_curve_y
        except RuntimeError: print("[DEBUG] Gaussian fit (curve_fit) yakınsamadı."); return None, None
        except Exception as e: print(f"[DEBUG] Gaussian fit sırasında genel hata: {e}"); return None, None

    def show_scale_parameters_dialog(self): ScaleDialog(self.root, self)
    def show_zscale_parameters_dialog(self): ZScaleDialog(self.root, self)
    def show_histogram_graph(self):
        if self.image_data is not None: HistogramDialog(self.root, self, self.image_data, self.current_vmin, self.current_vmax)
        else: messagebox.showwarning("Veri Yok", "Histogram için önce görüntü açın.")

    def cleanup_ui(self, clear_calib=True):
        print(f"UI temizleniyor (clear_calib={clear_calib})...")
        self.clear_aperture_overlays(True, True, True, True);
        self.selection_box = None; self.trace_peaks_x = None; self.trace_peaks_y = None; self.trace_coefficients = None; self.aperture_lower = None; self.aperture_upper = None
        self.current_vmin = None; self.current_vmax = None
        self.image_data = None; self.raw_image_data = None; self.extracted_pixels = None; self.extracted_flux = None
        self.arc_image_data = None; self.raw_arc_image_data = None; self.arc_extracted_pixels = None; self.arc_extracted_flux = None
        self.line_list_data = None; self.calib_points = []; self.dispersion_coefficients = None;
        self.dispersion_covariance = None; self.science_wavelength_error = None
        self.background_regions = []; self.background_patches = []; self.subtract_bg_var.set(False)
        self.science_continuum_flux = None; self.science_normalized_flux = None
        self.science_header = None; self.science_exptime = None;
        if clear_calib: self.dark_exptime = None
        if self.science_plot_window and self.science_plot_window.winfo_exists(): self.science_plot_window.destroy(); self.science_plot_window = None
        if self.arc_plot_window and self.arc_plot_window.winfo_exists(): self.arc_plot_window.destroy(); self.arc_plot_window = None
        if hasattr(self, 'hdulist') and self.hdulist:
            try: self.hdulist.close(); self.hdulist = None
            except Exception: pass
        if hasattr(self, 'header_text'): self.header_text.config(state=tk.NORMAL); self.header_text.delete('1.0', tk.END); self.header_text.insert(tk.END, "FITS Header..."); self.header_text.config(state=tk.DISABLED)
        if hasattr(self, 'poly_degree_label_var'): self.poly_degree_label_var.set(f"Polinom Derecesi: {self.polynomial_degree}")
        if hasattr(self, 'continuum_degree_label_var'): self.continuum_degree_label_var.set(f"Süreklilik Fit Derecesi: {self.continuum_fit_degree}")
        if hasattr(self, 'ap_width_label_var'): self.ap_width_label_var.set(f"Açıklık Genişliği: {self.aperture_width} px")
        if hasattr(self, 'trace_coeffs_label_var'): self.trace_coeffs_label_var.set("Trace Fit Katsayıları: -")
        if clear_calib:
            self.bias_data = None; self.dark_data = None; self.flat_data = None
            self.bias_path = None; self.dark_path = None; self.flat_path = None
            self.dark_exptime = None
        self._update_calib_status_labels()
        if hasattr(self, 'ax'): self.ax.cla(); self.ax.set_facecolor('lightgrey'); self.ax.set_title("Görüntü alanı"); self.ax.tick_params(axis='both',which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labelleft=False);
        if hasattr(self, 'canvas'): self.canvas.draw_idle()
        self.root.title("FITS Görüntüleyici")

    def exit_app(self):
        if self.hdulist:
            try: self.hdulist.close(); print("Açık FITS dosyası kapatıldı.")
            except Exception as e: print(f"FITS kapatılırken hata: {e}")
        self.root.quit(); self.root.destroy()
    def apply_scaling(self, redraw_overlays=True):
        if self.image_data is None:
            self.ax.cla()
            self.ax.set_facecolor('lightgrey')
            self.ax.set_title("Görüntü alanı")
            self.ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
            if hasattr(self, 'canvas'):
                self.canvas.draw_idle()
            return

        print(f"Görüntü ölçekleniyor: Stretch={self.current_stretch_type}, Limit={self.current_limit_type}")

        stretch_func = self.stretch_map.get(self.current_stretch_type, LinearStretch())
        if self.current_stretch_type == 'histeq':
            try:
                stretch_func = HistEqStretch(self.image_data)
            except Exception as hist_err:
                print(f"HistEqStretch hatası: {hist_err}, Linear stretch kullanılıyor.")
                stretch_func = LinearStretch()

        interval = None
        vmin, vmax = None, None
        limit_type_to_process = self.current_limit_type

        if limit_type_to_process == 'zmax':
            try:
                zscale_interval = ZScaleInterval(n_samples=self.zscale_nsamples, contrast=self.zscale_contrast)
                zmin_calc, _ = zscale_interval.get_limits(self.image_data)
                vmax_calc = np.nanmax(self.image_data)
                if np.isfinite(zmin_calc) and np.isfinite(vmax_calc) and vmax_calc > zmin_calc:
                    interval = ManualInterval(vmin=zmin_calc, vmax=vmax_calc)
                    print(f"ZMax kullanılıyor: Vmin={zmin_calc:.4g}, Vmax={vmax_calc:.4g}")
                else:
                    print(f"ZMax limitleri geçersiz ({zmin_calc=}, {vmax_calc=}), MinMax'a dönülüyor.")
                    interval = MinMaxInterval()
            except Exception as e:
                print(f"ZMax limitleri hesaplanırken hata ({e}), MinMax'a dönülüyor.")
                interval = MinMaxInterval()
        elif limit_type_to_process == 'manual':
            interval = ManualInterval(vmin=self.manual_vmin, vmax=self.manual_vmax)
        else:
            if limit_type_to_process == 'minmax': interval = MinMaxInterval()
            elif limit_type_to_process == 'percentile': interval = PercentileInterval(self.current_percentile)
            elif limit_type_to_process == 'zscale': interval = ZScaleInterval(n_samples=self.zscale_nsamples, contrast=self.zscale_contrast)
            else: print(f"Uyarı: Bilinmeyen limit tipi: {self.current_limit_type}. ZScale kullanılıyor."); interval = ZScaleInterval(n_samples=self.zscale_nsamples, contrast=self.zscale_contrast)

        try:
            if interval is not None:
                vmin, vmax = interval.get_limits(self.image_data)

            if not (vmin is not None and vmax is not None and np.isfinite(vmin) and np.isfinite(vmax) and vmax > vmin):
                print(f"Limit hesaplama başarısız veya geçersiz limitler ({vmin=}, {vmax=}). MinMax denemesi...")
                minmax_interval = MinMaxInterval()
                vmin, vmax = minmax_interval.get_limits(self.image_data)
                if not (vmin is not None and vmax is not None and np.isfinite(vmin) and np.isfinite(vmax) and vmax > vmin):
                    vmin = 0; vmax = 1; print("MinMax da başarısız veya geçersiz, [0, 1] kullanılıyor.")

            self.current_vmin = vmin; self.current_vmax = vmax;
            print(f"Hesaplanan Limitler: Vmin={vmin:.4g}, Vmax={vmax:.4g} ({interval.__class__.__name__ if interval else 'Fallback'})")

        except Exception as e:
            print(f"Limit hesaplama sırasında hata ({e}). MinMax denemesi...")
            try:
                minmax_interval = MinMaxInterval()
                vmin, vmax = minmax_interval.get_limits(self.image_data)
                if not (vmin is not None and vmax is not None and np.isfinite(vmin) and np.isfinite(vmax) and vmax > vmin): raise ValueError("MinMax failed too")
                self.current_vmin, self.current_vmax = vmin, vmax
                print(f"Hesaplanan Limitler (MinMax Fallback): Vmin={vmin:.4g}, Vmax={vmax:.4g}")
            except Exception as e_minmax:
                vmin, vmax = 0, 1; self.current_vmin, self.current_vmax = vmin, vmax
                print(f"MinMax hatası ({e_minmax}), [0, 1] kullanılıyor.")

        try:
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch_func, clip=True)

            if not self.ax.images:
                print("Adding new image to axes.")
                self.ax.imshow(self.image_data, cmap='gray', origin='lower', norm=norm, aspect='auto')
            else:
                print("Updating existing image on axes.")
                im = self.ax.images[0]
                im.set_data(self.image_data)
                im.set_norm(norm)

            img_title = self.root.title().replace("FITS Görüntüleyici - ", "")
            self.ax.set_title(img_title)

            if redraw_overlays:
                print("Calling _redraw_overlays from apply_scaling")
                self._redraw_overlays()
            else:
                self.canvas.draw_idle()

        except Exception as e:
            messagebox.showerror("Ölçekleme/Çizim Hatası", f"Görüntü çizilirken/ölçeklenirken hata: {e}", parent=self.root);
            print(f"Scaling/Drawing error details: {e}")
            self.ax.cla()
            self.ax.set_facecolor('lightgrey')
            self.ax.set_title("Görüntü alanı - Hata")
            self.ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
            self.canvas.draw_idle()


    def _redraw_overlays(self):
        """Mevcut tüm grafiksel kaplamaları (seçim kutusu, iz, açıklık, arka plan) yeniden çizer."""
        print("Kaplamalar (overlays) yeniden çiziliyor...")

        # Mevcut tüm overlayleri temizle (tekrar çizmek için)
        for artist in self.overlay_plots:
            try:
                artist.remove()
            except ValueError:
                pass
        self.overlay_plots = [] # Listeyi temizle

        # 1. Seçim Kutusunu Çiz (eğer varsa)
        if self.selection_box:
            xmin, xmax, ymin, ymax = self.selection_box
            width = xmax - xmin
            height = ymax - ymin
            # Rectangle patches'i her zaman yeni oluşturup listeye ekle, aksi halde eski objelere referans kalabilir.
            self.selection_patch = patches.Rectangle((xmin, ymin), width, height, linewidth=1, edgecolor='cyan', facecolor='none', linestyle='--', label='Seçim Bölgesi', zorder=11)
            self.ax.add_patch(self.selection_patch)
            self.overlay_plots.append(self.selection_patch)

        # 2. İz Noktalarını ve Fit Edilmiş İzi Çiz (eğer varsa)
        if self.trace_peaks_x is not None and self.trace_peaks_y is not None:
            trace_scatter = self.ax.scatter(self.trace_peaks_x, self.trace_peaks_y, c='red', marker='.', s=10, label='Bulunan İz Noktaları', zorder=12)
            self.overlay_plots.append(trace_scatter)
            if self.trace_coefficients is not None and self.image_data is not None:
                img_w = self.image_data.shape[1]
                x_fit_range = np.arange(img_w)
                fit_func = np.poly1d(self.trace_coefficients)
                y_fit_trace = fit_func(x_fit_range)
                trace_line, = self.ax.plot(x_fit_range, y_fit_trace, 'lime', linestyle='-', linewidth=1.5, label='Fit Edilmiş İz', zorder=11)
                self.overlay_plots.append(trace_line)

        # 3. Açıklık Sınırlarını Çiz (eğer varsa)
        if self.aperture_lower is not None and self.aperture_upper is not None and self.image_data is not None:
            img_w = self.image_data.shape[1]
            x_range = np.arange(img_w)
            # Etiketi sadece bir kere ekle
            ap_lower_line, = self.ax.plot(x_range, self.aperture_lower, 'yellow', linestyle='--', linewidth=1, label='Açıklık Sınırları', zorder=10)
            ap_upper_line, = self.ax.plot(x_range, self.aperture_upper, 'yellow', linestyle='--', linewidth=1, label='_nolegend_', zorder=10)
            self.overlay_plots.append(ap_lower_line)
            self.overlay_plots.append(ap_upper_line)

        # 4. Arka Plan Bölgelerini Çiz (eğer varsa)
        # draw_background_regions kendi içinde patch'leri yönetir ve self.background_patches'i doldurur.
        # Bu patch'ler de self.overlay_plots'a eklenmelidir.
        for region in self.background_regions: # Her zaman yeniden çizebilmek için
            xmin, xmax, ymin, ymax = region
            width = xmax - xmin
            height = ymax - ymin
            patch = patches.Rectangle((xmin, ymin), width, height, linewidth=1, edgecolor='red', facecolor='none', linestyle=':', label='Arkaplan Bölgesi', zorder=7)
            self.ax.add_patch(patch)
            self.overlay_plots.append(patch) # Her seferinde yeniden ekle

        # 5. Uzaysal Profil Çizgisini Çiz (eğer varsa)
        if self.spatial_profile_line:
            # Spatial profile line objesi zaten ax üzerinde olabilir, yeni bir tane çizme, sadece görünürlüğünü kontrol et
            # Ama eğer _redraw_overlays metodunu temizleme/yeniden oluşturma yaklaşımıyla kullanıyorsak:
            # O zaman burada da yeniden oluşturup eklemeliyiz.
            # Şu anki mantıkta `self.spatial_profile_line` objesi doğrudan `_on_spatial_profile_column_selected` içinde `ax.axvline` ile oluşturuluyor.
            # Dolayısıyla eğer sıfırdan çiziyorsak bu objeyi de yeniden yaratmalıyız.
            # Ancak, _redraw_overlays çağrısı genellikle mevcut çizimlerin üzerine yeni çizimler ekler veya eski objeleri günceller.
            # Eğer self.spatial_profile_line'ın remove() metodu çağrılmadıysa, tekrar oluşturmaya gerek yok, zaten eksendedir.
            # Bu durumdaki en iyi yaklaşım, bu objeyi de self.overlay_plots'a eklemek ve yukarıdaki döngüde kaldırılıp yeniden eklenmesini sağlamaktır.
            # _on_spatial_profile_column_selected'da zaten self.overlay_plots'a ekliyor.
            # Dolayısıyla, buradaki logic: eğer `self.spatial_profile_line` varsa ve hala eksendeyse, `overlay_plots` içinde kalır.
            # Eğer `_on_spatial_profile_column_selected` onu kaldırıp yeni bir tane oluşturursa, o da `overlay_plots`a eklenir.
            pass # Bu kısımda ekstra bir şey yapmaya gerek yok, yukarıdaki temizleme/yeniden oluşturma döngüsü onu zaten hallediyor.

        # 6. Efsaneyi (Legend) Güncelle
        unique_handles_labels = {}
        for artist in self.overlay_plots:
            if hasattr(artist, 'get_label'):
                label = artist.get_label()
                if label and not label.startswith('_') and label not in unique_handles_labels:
                    unique_handles_labels[label] = artist

        current_legend = self.ax.get_legend()
        if current_legend:
            current_legend.remove() # Önceki efsaneyi kaldır

        if unique_handles_labels:
            self.ax.legend(unique_handles_labels.values(), unique_handles_labels.keys(), fontsize='x-small', loc='upper right')

        self.canvas.draw_idle()

# === Ana Program Bloğu ===
if __name__ == "__main__":
    main_window = tk.Tk()
    style = ttk.Style(); available_themes = style.theme_names(); preferred_themes = ['vista', 'clam', 'alt', 'default']
    for theme in preferred_themes:
        if theme in available_themes:
            try: style.theme_use(theme); print(f"Using ttk theme: {theme}"); break
            except tk.TclError: continue
    else: print("Default ttk theme used.")

    # SciPy kontrolünü burada yap, hata mesajını Tkinter penceresinde göster
    if not SCIPY_AVAILABLE:
        messagebox.showwarning("SciPy Eksik",
                               "UYARI: SciPy kütüphanesi bulunamadı. Spline fit ve Gaussian çizgi merkezleme/fit özellikleri KULLANILAMAYACAKTIR.\n"
                               "Lütfen komut istemcisinde 'pip install scipy' yazarak yükleyin.",
                               parent=main_window)

    try:
        from astropy.stats import sigma_clip
    except ImportError:
        messagebox.showerror("Astropy Eksik",
                              "UYARI: Normalizasyon ve diğer astronomi işlemleri için 'astropy' gereklidir (sigma_clip bulunamadı).\n"
                              "Lütfen komut istemcisinde 'pip install astropy' yazarak yükleyin.",
                              parent=main_window)
        exit() # Astropy olmadan devam etmek anlamsız

    app = FitsViewerApp(main_window)
    main_window.protocol("WM_DELETE_WINDOW", app.exit_app)
    main_window.mainloop()