Amplitude
---------

K_Class amplitude calculation is using a wavelet (Mexican hat) transform for the peak P wave amplitude
search and maximum value of an S wave at any horizontal component.
Since the wavelet is used only for P wave detection at the vertical component its scale could be configurable via the plugin configuration. 

Station Magnitude
-----------------

The K_class plugin calculates the individual station magnitude by the formula obtained from (G.K. Aslanov et al., 2015):

.. math::

   K_Class = 2.94 + 1.935(\log10(a) + 1.734 \log10(hypdistkm))

Where a is calculated as

.. math::

   a = (Ap + As) / V

Where Ap and As are ampitudes of P and S wave respectivly and V is the gain of the instrument.
Hypdistkm is the distance from the sensor to the hypocenter in kilometers.

* Distance range: 0 - 15 deg
* Depth range: 0 - 180 km

Configuration
-------------

Add the *K_Class* plugin to the existing plugins in the global configuration.
Set configurable parameters in the global bindings to compute K_class.
