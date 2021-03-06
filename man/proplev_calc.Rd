% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port-perf.R
\name{proplev_calc}
\alias{proplev_calc}
\title{Proportional Leverage}
\usage{
proplev_calc(weights)
}
\arguments{
\item{weights}{a numerical nxp data matrix with portfolio weights.}
}
\value{
a double, the proportional leverage over the specified observation period.
}
\description{
Calculates the proportional leverage (% short sales) of a portfolio strategy.
}
\examples{
set.seed(1234)
naive_port <- matrix(runif(200, 0, 1), 20, 10)
turnover <- proplev_calc(naive_port)

}
