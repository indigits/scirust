Interpolation Library
===========================

Interpolation is a two-phased operation in srinterp:

1. Create an interpolation scheme for the chosen algorithm and optimized for the given sample points. You get back a class that implements the _Interpolation_ trait.
2. Use this scheme to compute values at arbitrary points. Some interpolation algorithms also allow you to compute the derivative and the indefinite integral at that point.

The unit-like struct `Interpolate` provides simple factory methods to create the interpolation scheme in a simple method call:

<!-- * _RationalWithoutPoles_, creates a Floater-Hormann barycentric interpolation
* _RationalWithPoles_, creates a Bulirsch & Stoer rational interpolation -->
* _linear_, creates a linear spline interpolation

If unsure, we recommend using _RationalWithoutPoles_ for most cases.

Alternatively you can also use the algorithms directly. The following algorithms are available:


<!-- Interpolation on equidistant sample points
------------------------------------------

* *Polynomial*: Barycentric Algorithm -->


Interpolation on arbitrary sample points
----------------------------------------

<!-- * *Rational pole-free*: Barycentric Floater-Hormann Algorithm
* **Rational with poles**: Bulirsch & Stoer Algorithm
* *Neville Polynomial*: Neville Algorithm. Note that the Neville algorithm performs very badly on equidistant points. If you need to interpolate a polynomial on equidistant points, we recommend to use the barycentric algorithm instead. -->
* *Linear Spline*
<!-- * *Cubic Spline* with boundary conditions
* *Natural Cubic Spline*
* *Akima Cubic Spline* -->


Interpolation with additional data
----------------------------------

<!-- * *Generic Barycentric Interpolation*, requires barycentric weights -->
* *Generic Spline*, requires spline coefficients
<!-- * *Generic Cubic Hermite Spline*, requires the derivatives -->