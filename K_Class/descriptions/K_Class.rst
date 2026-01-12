Applicability
---------

* Depth: 0 - 80 km;
* Distance 0-1000 km;

Amplitude
---------

K_Class amplitude calculation is using a custom time window (between the P and S waves arrivals) for the maximum P wave amplitude
search on the vertical component and falls back for the Mlh amplitude on horizontals. 

Magnitude
---------

Maximum value of P wave on the vertical and S wave at any horizontal component is used. The calculation is done
by the formula used in the WSG software (A.P. Akimov et al., 2020):

.. math::

	K_{Class} = A \cdot \left( \log_{10}(Amp) + B(R) \right)

Where:

* :math:`A` is the main slope (default: 1.84);
* :math:`Amp` is the amplitude sum of P and S waves;
* :math:`B(R)` is a piecewise distance-correction term:

.. math::

     B(R) =
     \begin{cases}
     a_1 \cdot \log_{10}(R) + b_1, & R \le l_1 \\
     a_2 \cdot \log_{10}(R) + b_2, & l_1 < R \le l_2 \\
     a_3 \cdot \log_{10}(R) + b_3, & l_2 < R \le l_3 \\
     a_4 \cdot \log_{10}(R) + b_4, & R > l_3
     \end{cases}

where :math:`l_x` is hypocentral distance in km.

There is a python helper script provided for plotting the B(R)

.. figure:: B_R_plot.png
   :align: center
   :alt: B_R_plot image

Defaults
-------------

See SeisComP .xml file or refer to these values:

=========== =============
Coefficient Default Value
=========== =============
l1          75.0 (km)
l2          264.0 (km)
l3          800.0 (km)
A           1.84
a1          2.11
a2          1.1
a3          2.98
a4          0.0
b1          1.32
b2          3.21
b3          -1.34
b4          8.0
=========== =============

Configuration
-------------

Add the *K_Class* plugin to the existing plugins in the global configuration.
Set configurable coefficients depending on your region or start with the default set to compute K_class.
