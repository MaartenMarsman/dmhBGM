% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmhBGM-package.R
\docType{package}
\name{dmhBGM-package}
\alias{dmhBGM}
\alias{dmhBGM-package}
\title{dmhBGM: Bayesian Analysis of Graphical Models using DMH}
\description{
The \code{R} package \strong{dmhBGM} is a derivative of the \code{R} package
\strong{bgms} version \code{0.1.2}. The latter package provides tools for
Bayesian analysis of graphical models describing networks of variables. The
package uses Bayesian variable selection methods to model the underlying
network structure. The graphical model has an intractable normalizing
constant. \strong{bgms} circumvents this normalizing constant by employing a
joint pseudolikelihood. \strong{dmhBGM} uses the Double Metropolis Hastings
algorithm \insertCite{Liang_2010}{dmhBGM} instead.

The package uses Gibbs variable selection
\insertCite{GeorgeMcCulloch_1993}{dmhBGM} and uses a discrete spike and slab
prior on the pairwise interactions, which allows us to set the interactions
to exact zeroes. To account for the discontinuity at zero, we embed a
Metropolis approach for mixtures of mutually singular distributions
\insertCite{GottardoRaftery_2008}{dmhBGM} in a Gibbs sampler.
}
\references{
\insertAllCited{}
}
\author{
\strong{Maintainer}: Maarten Marsman \email{m.marsman@uva.nl} (\href{https://orcid.org/0000-0001-5309-7502}{ORCID})

Authors:
\itemize{
  \item Don van den Bergh (\href{https://orcid.org/0000-0002-9838-7308}{ORCID})
}

Other contributors:
\itemize{
  \item Nikola Sekulovski (\href{https://orcid.org/0000-0001-7032-1684}{ORCID}) [contributor]
}

}
\keyword{internal}
