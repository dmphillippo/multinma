# Specify baseline connections

Helper function for the `connect_baseline` argument of
[`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md)
to specify how study baselines (intercepts) are linked, to connect a
disconnected network or single-arm studies.

## Usage

``` r
con(
  type = c("fixed", "random"),
  studies,
  prior_baseline = NULL,
  baseline_trt = NULL
)
```

## Arguments

- type:

  Type of connection, either "fixed" or "random".

- studies:

  Vector of study names.

- prior_baseline:

  Prior distribution for the baseline parameters of `studies` when
  `type = "random"`, as a
  [nma_prior](https://dmphillippo.github.io/multinma/dev/reference/nma_prior-class.md)
  object (see
  [priors](https://dmphillippo.github.io/multinma/dev/reference/priors.md)).

- baseline_trt:

  Which treatment the baseline prior distribution applies to, when
  `type = "random"`. By default, the random baseline prior will be
  placed on the network reference treatment. This should be changed if,
  for example, the prior was obtained from a baseline synthesis of a
  different treatment.

## Value

An object of class `nma_connect`.

## Details

A random baseline connection is specified with `type = "random"`. This
places an informative prior distribution on the baselines of studies
listed in `studies`, specified with the `prior_baseline` argument.

A fixed baseline connection is specified with `type = "fixed"`. This
shares the baseline parameter for studies listed in `studies`.
