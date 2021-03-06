% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-clean.R
\name{find_nas}
\alias{find_nas}
\title{Find NAs}
\usage{
find_nas(data, all = TRUE)
}
\arguments{
\item{data}{a numerical vector or a numerical nxp data matrix.}

\item{all}{a logical, if TRUE (default), the function finds columns which consist only of NAs.
If FALSE, it finds the columns with any number of NAs.}
}
\value{
an index vector of the position of NAs.
}
\description{
Finds the NAs in the data and returns the respective column index.
}
\examples{
data(sp500)
index_nas <- find_nas(sp500, all=FALSE)

}
