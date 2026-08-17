# Distribution functions for M-spline baseline hazards

Density, distribution, quantile, hazard, cumulative hazard, and
restricted mean survival time functions for the M-spline baseline
hazards model.

## Usage

``` r
dmspline(x, basis, scoef, rate, log = FALSE)

pmspline(q, basis, scoef, rate, lower.tail = TRUE, log.p = FALSE)

qmspline(p, basis, scoef, rate, lower.tail = TRUE, log.p = FALSE)

hmspline(x, basis, scoef, rate, log = FALSE)

Hmspline(x, basis, scoef, rate, log = FALSE)

rmst_mspline(t, basis, scoef, rate, start = 0)
```

## Arguments

- x, q:

  Vector of quantiles

- basis:

  M-spline basis produced by
  [`splines2::mSpline()`](https://wwenjie.org/splines2/reference/mSpline.html)

- scoef:

  Vector (or matrix) of spline coefficients with length (or number of
  columns) equal to the dimension of `basis`

- rate:

  Vector of rate parameters

- log, log.p:

  Logical; if `TRUE`, probabilities `p` are given as \\\log(p)\\

- lower.tail:

  Logical; if `TRUE` (the default), probabilities are \\P(X \le x)\\,
  otherwise \\P(X \> x)\\

- p:

  Vector of probabilities

- t:

  Vector of times to which the restricted mean survival time is
  calculated

- start:

  Optional left-truncation time or times. The returned restricted mean
  survival will be conditioned on survival up to this time

## Value

`dmspline()` gives the density, `pmspline()` gives the distribution
function (CDF), `qmspline()` gives the quantile function (inverse-CDF),
`hmspline()` gives the hazard function, `Hmspline()` gives the
cumulative hazard function, and `rmst_mspline()` gives restricted mean
survival times.

## Details

Survival models with a flexible M-spline on the baseline hazard are
described by Phillippo et al. (2026) . Piecewise-exponential baseline
hazards are a special case where the degree of the M-spline polynomial
is 0.

The d/p/h/H functions are calculated from their definitions.
`qmspline()` uses numerical inversion via
[`flexsurv::qgeneric()`](http://chjackson.github.io/flexsurv-dev/reference/qgeneric.md).
`rmst_mspline()`uses numerical integration via
[`flexsurv::rmst_generic()`](http://chjackson.github.io/flexsurv-dev/reference/rmst_generic.md),
except for the special case of the piecewise-exponential hazard (i.e.
degree 0 M-splines) which uses the explicit formula from Royston and
Parmar (2013) .

Beyond the boundary knots, the hazard is assumed to be constant. (This
differs from the approach in
[`splines2::mSpline()`](https://wwenjie.org/splines2/reference/mSpline.html)
that extrapolates the polynomial basis functions, which is numerically
unstable and highly dependent on the data just before the boundary
knots.) As with all extrapolation, care should be taken when evaluating
the splines at times beyond the boundary knots (either directly through
the d/p/h/H/rmst functions, or indirectly by requesting quantiles with
`qmspline()` that correspond to times beyond the boundary knots). For
this reason evaluating the (unrestricted) mean survival time is not
generally recommended as this requires integrating over an infinite time
horizon (i.e. `rmst_mspline()` with `t = Inf`).

## References

Phillippo DM, Sadek A, Pedder H, Welton NJ (2026). “Network
Meta-Analysis of survival outcomes with non-proportional hazards using
flexible M-splines.” *Statistics in Medicine*, **45**(18-19), e70695.
ISSN 1097-0258.
[doi:10.1002/sim.70695](https://doi.org/10.1002/sim.70695) .  
  
Royston P, Parmar MKB (2013). “Restricted mean survival time: an
alternative to the hazard ratio for the design and analysis of
randomized trials with a time-to-event outcome.” *BMC Medical Research
Methodology*, **13**(1).
[doi:10.1186/1471-2288-13-152](https://doi.org/10.1186/1471-2288-13-152)
.
