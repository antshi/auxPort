% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port-optim.R
\name{port_estim_gmv}
\alias{port_estim_gmv}
\title{Global Minimum-Variance Portfolio Optimization}
\usage{
port_estim_gmv(Sigma)
}
\arguments{
\item{Sigma}{a pxp covariance matrix of asset returns.}
}
\value{
a vector of length p with the weights of the portfolio.
}
\description{
Calculates the weights of a global minimum-variance portfolio strategy.
}
\examples{
data(sp500_rets)
example_rets <- sp500_rets[,2:11]
covMat <- var(example_rets)
port_estim_gmv(covMat)

}
