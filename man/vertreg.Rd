% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vertreg.R
\name{vertreg}
\alias{vertreg}
\title{Vertical Regression}
\usage{
vertreg(formula, data_pre, data_post, pseudo_inverse = FALSE)
}
\arguments{
\item{formula}{A string to specify the vertical regression model (e.g., "Y ~ X1 + X2" where "Y" is treated unit, "X1" and "X2" are control units).}

\item{data_pre}{Pre-period data where the rows are time points and the columns include treated/control units. See \code{synth_pre} for example.}

\item{data_post}{Post-period data. Should have identical columns to data_pre. See \code{synth_post} for example.}

\item{pseudo_inverse}{A logical to indicate whether to use the pseudo-inverse to fit the model.}
}
\value{
A vector of treatment effects.
}
\description{
Function to run a vertical regression.
}
\examples{
vertreg(
  formula = paste("Y ~ -1 + ", paste(paste0("X", 1:16), collapse = " + ")),
  data_pre = synth_pre,
  data_post = synth_post
)

}
\seealso{
\code{\link[=vertreg_stacked]{vertreg_stacked()}}
}
