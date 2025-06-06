% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc-examples.R
\name{example_smk_nodesplit}
\alias{example_smk_nodesplit}
\title{Example smoking node-splitting}
\description{
Calling \code{example("example_smk_nodesplit")} will run
node-splitting models with the smoking cessation data, using the code in
the Examples section below. The resulting \code{nma_nodesplit_df} object
\code{smk_fit_RE_nodesplit} will then be available in the global environment.
}
\details{
Smoking node-splitting for use in examples.
}
\examples{
# Set up network of smoking cessation data
head(smoking)

smk_net <- set_agd_arm(smoking,
                       study = studyn,
                       trt = trtc,
                       r = r,
                       n = n,
                       trt_ref = "No intervention")

# Print details
smk_net

\donttest{
# Fitting all possible node-splitting models
smk_fit_RE_nodesplit <- nma(smk_net, 
                            consistency = "nodesplit",
                            trt_effects = "random",
                            prior_intercept = normal(scale = 100),
                            prior_trt = normal(scale = 100),
                            prior_het = normal(scale = 5))
}

\dontshow{
if (requireNamespace("pkgdown", quietly = TRUE) && pkgdown::in_pkgdown()) {
  assign("smk_net", smk_net, .GlobalEnv)
  assign("smk_fit_RE_nodesplit", smk_fit_RE_nodesplit, .GlobalEnv)
}
}
}
\keyword{examples}
