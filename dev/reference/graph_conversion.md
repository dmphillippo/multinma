# Convert networks to graph objects

The method `as.igraph()` converts `nma_data` objects into the form used
by the [igraph](https://r.igraph.org/reference/aaa-igraph-package.html)
package. The method `as_tbl_graph()` converts `nma_data` objects into
the form used by the
[ggraph](https://ggraph.data-imaginist.com/reference/ggraph-package.html)
and
[tidygraph](https://tidygraph.data-imaginist.com/reference/tidygraph-package.html)
packages.

## Usage

``` r
# S3 method for class 'nma_data'
as.igraph(x, ..., collapse = TRUE)

# S3 method for class 'nma_data'
as_tbl_graph(x, ...)
```

## Arguments

- x:

  An
  [nma_data](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
  object to convert

- ...:

  Additional arguments

- collapse:

  Logical, collapse edges over studies? Default `TRUE`, only one edge is
  produced for each comparison (by IPD or AgD study type) with a
  `.nstudy` attribute giving the number of studies making that
  comparison. If `FALSE`, repeated edges are added for each study making
  the comparison.

## Value

An `igraph` object for `as.igraph()`, a `tbl_graph` object for
`as_tbl_graph()`.

## Examples

``` r
# Set up network of smoking cessation data
head(smoking)
#>   studyn trtn                   trtc  r   n
#> 1      1    1        No intervention  9 140
#> 2      1    3 Individual counselling 23 140
#> 3      1    4      Group counselling 10 138
#> 4      2    2              Self-help 11  78
#> 5      2    3 Individual counselling 12  85
#> 6      2    4      Group counselling 29 170

smk_net <- set_agd_arm(smoking,
                       study = studyn,
                       trt = trtc,
                       r = r,
                       n = n,
                       trt_ref = "No intervention")

# Print details
smk_net
#> A network with 24 AgD studies (arm-based).
#> 
#> ------------------------------------------------------- AgD studies (arm-based) ---- 
#>  Study Treatment arms                                                 
#>  1     3: No intervention | Group counselling | Individual counselling
#>  2     3: Group counselling | Individual counselling | Self-help      
#>  3     2: No intervention | Individual counselling                    
#>  4     2: No intervention | Individual counselling                    
#>  5     2: No intervention | Individual counselling                    
#>  6     2: No intervention | Individual counselling                    
#>  7     2: No intervention | Individual counselling                    
#>  8     2: No intervention | Individual counselling                    
#>  9     2: No intervention | Individual counselling                    
#>  10    2: No intervention | Self-help                                 
#>  ... plus 14 more studies
#> 
#>  Outcome type: count
#> ------------------------------------------------------------------------------------
#> Total number of treatments: 4
#> Total number of studies: 24
#> Reference treatment is: No intervention
#> Network is connected

# Convert to igraph object
igraph::as.igraph(smk_net)  # Edges combined by default
#> IGRAPH 8f4375f UN-- 4 6 -- 
#> + attr: name (v/c), .sample_size (v/n), .nstudy (e/n), .type (e/c)
#> + edges from 8f4375f (vertex names):
#> [1] No intervention       --Group counselling     
#> [2] No intervention       --Individual counselling
#> [3] Group counselling     --Individual counselling
#> [4] No intervention       --Self-help             
#> [5] Group counselling     --Self-help             
#> [6] Individual counselling--Self-help             
igraph::as.igraph(smk_net, collapse = FALSE)  # Without combining edges
#> IGRAPH 8133bd2 UN-- 4 28 -- 
#> + attr: name (v/c), .sample_size (v/n), .study (e/x), .type (e/c)
#> + edges from 8133bd2 (vertex names):
#>  [1] No intervention       --Group counselling     
#>  [2] No intervention       --Individual counselling
#>  [3] Group counselling     --Individual counselling
#>  [4] Group counselling     --Individual counselling
#>  [5] Group counselling     --Self-help             
#>  [6] Individual counselling--Self-help             
#>  [7] No intervention       --Individual counselling
#>  [8] No intervention       --Individual counselling
#> + ... omitted several edges

# Convert to tbl_graph object
tidygraph::as_tbl_graph(smk_net)  # Edges combined by default
#> # A tbl_graph: 4 nodes and 6 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 4 × 2 (active)
#>   name                   .sample_size
#>   <chr>                         <dbl>
#> 1 No intervention                7231
#> 2 Group counselling               555
#> 3 Individual counselling         7383
#> 4 Self-help                      1571
#> #
#> # Edge Data: 6 × 4
#>    from    to .nstudy .type
#>   <int> <int>   <int> <chr>
#> 1     1     2       2 AgD  
#> 2     1     3      15 AgD  
#> 3     2     3       4 AgD  
#> # ℹ 3 more rows
tidygraph::as_tbl_graph(smk_net, collapse = FALSE)  # Without combining edges
#> # A tbl_graph: 4 nodes and 28 edges
#> #
#> # An undirected multigraph with 1 component
#> #
#> # Node Data: 4 × 2 (active)
#>   name                   .sample_size
#>   <chr>                         <dbl>
#> 1 No intervention                7231
#> 2 Group counselling               555
#> 3 Individual counselling         7383
#> 4 Self-help                      1571
#> #
#> # Edge Data: 28 × 4
#>    from    to .study .type
#>   <int> <int> <fct>  <chr>
#> 1     1     2 1      AgD  
#> 2     1     3 1      AgD  
#> 3     2     3 1      AgD  
#> # ℹ 25 more rows
```
