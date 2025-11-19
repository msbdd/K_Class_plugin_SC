Amplitude
---------

K_Class amplitude calculation is using a custom time window (between the P and S waves arrivals) for the maximum P wave amplitude
search on the vertical component and maximum value of an S wave at any horizontal component.

Station Magnitude
-----------------

The K_class plugin calculates the individual station magnitude by the formula obtained from (G.K. Aslanov et al., 2015):

TBD

Where Ap and As are ampitudes of P and S wave respectivly and V is the gain of the instrument.
Hypdistkm is the distance from the sensor to the hypocenter in kilometers.

* Distance range: 0 - 10 deg
* Depth range: 0 - 80 km

Configuration
-------------

Add the *K_Class* plugin to the existing plugins in the global configuration.
Set configurable parameters in the global bindings to compute K_class.
