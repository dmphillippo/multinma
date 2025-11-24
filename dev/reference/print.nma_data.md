# Print `nma_data` objects

Print details of networks stored as
[nma_data](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
objects, as created by
[`set_ipd()`](https://dmphillippo.github.io/multinma/dev/reference/set_ipd.md),
[`set_agd_arm()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_arm.md),
[`set_agd_contrast()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_contrast.md),
[`set_agd_surv()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_surv.md),
or
[`combine_network()`](https://dmphillippo.github.io/multinma/dev/reference/combine_network.md).

## Usage

``` r
# S3 method for class 'nma_data'
print(x, ..., n = 10)

# S3 method for class 'mlnmr_data'
print(x, ..., n = 10)
```

## Arguments

- x:

  `nma_data` object

- ...:

  other options (not used)

- n:

  number of studies of each type to print

## Value

`x` is returned invisibly.
