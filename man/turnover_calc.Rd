% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port-perf.R
\name{turnover_calc}
\alias{turnover_calc}
\title{Turnover}
\usage{
turnover_calc(weights, rets_oos)
}
\arguments{
\item{weights}{a numerical nxp data matrix with portfolio weights.}

\item{rets_oos}{a numerical nxp data matrix with stock returns over the same observation period as the weights.}
}
\value{
a double, the average turnover rate over the specified observation period.
}
\description{
Calculates the average turnover rate of a portfolio strategy.
}
\examples{
set.seed(1234)
naive_port <- matrix(runif(200, 0, 1), 20, 10)
returns_oos <- matrix(rnorm(200, 0, 0.1), 20, 10)
turnover <- turnover_calc(naive_port, returns_oos)

}
