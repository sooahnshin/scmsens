% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudo_inverse.R
\name{summary.lm_pseudo}
\alias{summary.lm_pseudo}
\title{Summary Method for lm_pseudo Objects}
\usage{
\method{summary}{lm_pseudo}(object, ...)
}
\arguments{
\item{object}{An object of class \code{lm_pseudo}.}

\item{...}{Additional arguments (ignored).}
}
\value{
A tibble containing the estimated coefficients with their corresponding terms.
}
\description{
Extracts the summary of an \code{lm_pseudo} object, returning a tibble with estimated coefficients.
}
\seealso{
\code{\link[=lm_pseudo]{lm_pseudo()}}, \code{\link[=predict.lm_pseudo]{predict.lm_pseudo()}}, \code{\link[=coef.lm_pseudo]{coef.lm_pseudo()}}, \code{\link[=tidy.lm_pseudo]{tidy.lm_pseudo()}}
}
