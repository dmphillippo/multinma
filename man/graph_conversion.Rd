% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nma_data-class.R
\name{as.igraph.nma_data}
\alias{as.igraph.nma_data}
\alias{as_tbl_graph.nma_data}
\title{Convert networks to graph objects}
\usage{
\method{as.igraph}{nma_data}(x, ..., collapse = TRUE)

\method{as_tbl_graph}{nma_data}(x, ...)
}
\arguments{
\item{x}{An \link{nma_data} object to convert}

\item{...}{Additional arguments}

\item{collapse}{Logical, collapse edges over studies? Default \code{TRUE}, only
one edge is produced for each comparison (by IPD or AgD study type) with a
\code{.nstudy} attribute giving the number of studies making that comparison. If
\code{FALSE}, repeated edges are added for each study making the comparison.}
}
\value{
An \code{igraph} object for \code{as.igraph()}, a \code{tbl_graph} object for
\code{as_tbl_graph()}.
}
\description{
The method \code{as.igraph()} converts \code{nma_data} objects into the form used by
the \link[igraph:igraph-package]{igraph} package. The method \code{as_tbl_graph()}
converts \code{nma_data} objects into the form used by the
\link[ggraph:ggraph-package]{ggraph} and
\link[tidygraph:tidygraph-package]{tidygraph} packages.
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

# Convert to igraph object
igraph::as.igraph(smk_net)  # Edges combined by default
igraph::as.igraph(smk_net, collapse = FALSE)  # Without combining edges

# Convert to tbl_graph object
tidygraph::as_tbl_graph(smk_net)  # Edges combined by default
tidygraph::as_tbl_graph(smk_net, collapse = FALSE)  # Without combining edges
}
