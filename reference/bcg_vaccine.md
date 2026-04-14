# BCG vaccination

Data frame containing the results of 13 trials comparing BCG vaccination
to no vaccination for preventing tuberculosis (TB) (Dias et al. 2011;
Berkey et al. 1995) . The numbers of individuals diagnosed with TB in
each arm during the study follow-up period are recorded. The absolute
degrees latitude at which the study was conducted are also recorded.

## Usage

``` r
bcg_vaccine
```

## Format

A data frame with 26 rows and 6 variables:

- studyn:

  numeric study ID

- trtn:

  numeric treatment code

- trtc:

  treatment name

- latitude:

  absolute degrees latitude

- r:

  number diagnosed with TB

- n:

  sample size

## References

Berkey CS, Hoaglin DC, Mosteller F, Colditz GA (1995). “A random-effects
regression model for meta-analysis.” *Statistics in Medicine*,
**14**(4), 395–411.
[doi:10.1002/sim.4780140406](https://doi.org/10.1002/sim.4780140406) .  
  
Dias S, Sutton AJ, Welton NJ, Ades AE (2011). “NICE DSU Technical
Support Document 3: Heterogeneity: subgroups, meta-regression, bias and
bias-adjustment.” National Institute for Health and Care Excellence.
<https://sheffield.ac.uk/nice-dsu>.
