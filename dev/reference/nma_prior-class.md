# The nma_prior class

The `nma_prior` class is used to specify prior distributions.

## Details

Objects of class `nma_prior` have the following components:

- `dist`:

  Distribution name

- `fun`:

  Name of constructor function, as string (e.g. `"normal"`)

- `...`:

  Parameters of the distribution

The distribution parameters, specified as named components in `...`,
match those in the constructor functions (see
[priors](https://dmphillippo.github.io/multinma/dev/reference/priors.md)).
