# K_Class Magnitude Plugin for SeisComP

This plugin provides a custom magnitude calculation for SeisComP, known as Class (or "Klass"), which is traditionally used in seismological practice in several countries.

---

## Disclaimer

This plugin is provided "as is" and is the result of an ongoing development and learning effort. It may contain bugs, errors, or inaccuracies.

No warranty of any kind is provided, express or implied. The user assumes all responsibility and risk for the use of this software. Thorough testing is strongly advised before deploying it in a production environment.

---

## Acknowledgments

The development of this plugin was heavily inspired by and adapted from the `MagnitudeProcessor_ML` (MLh) plugin, originally developed by the **Swiss Seismological Service (ETHZ/SED)**.

Significant portions of the class structure, configuration file handling, and processing logic are based on their work. The original authors are gratefully acknowledged for providing a robust and clear foundation that made this implementation possible.

---

## Developer Note: API Compatibility

This plugin has a known incompatibility between SeisComP API versions.

* **SeisComP 6.X:** The `setDefaults()` virtual function **does not exist** in the `Processing::MagnitudeProcessor` base class. To compile this plugin for SeisComP 6, you **must delete or comment out** the entire `setDefaults()` function in the `MagnitudeProcessor_K_Class` class.

* **SeisComP 7.X:** The `setDefaults()` function is part of the API and should be included.

The plugin will fail to compile on SeisComP 6 if the `setDefaults()` function is present.

---

## Amplitude

K_Class amplitude calculation uses a wavelet (Mexican hat) transform for the peak P wave amplitude search and the maximum value of an S wave at any horizontal component.
Since the wavelet is used only for P wave detection at the vertical component, its scale can be configured via the plugin configuration.

---

## Station Magnitude

The K_class plugin calculates the individual station magnitude by the formula obtained from (G.K. Aslanov et al., 2015):

$$
K_{\text{Class}} = 2.94 + 1.935(\log_{10}(a) + 1.734 \log_{10}(\text{hypdistkm}))
$$

Where $a$ is calculated as:

$$
a = (A_p + A_s) / V
$$

* $A_p$ and $A_s$ are the amplitudes of the P and S waves, respectively.
* $V$ is the gain of the instrument.
* `hypdistkm` is the distance from the sensor to the hypocenter in kilometers.

### Operational Ranges

* **Distance range:** 0 - 15 deg
* **Depth range:** 0 - 180 km

---

## Configuration

1.  Add the `K_Class` plugin to the existing plugins in the global configuration or your station profile.
2.  Set configurable parameters in the global bindings to compute K_class.