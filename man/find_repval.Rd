% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-clean.R
\name{find_repval}
\alias{find_repval}
\title{Find Repeated Values}
\usage{
find_repval(data, maxrep = 5)
}
\arguments{
\item{data}{an nxp numerical data matrix.}

\item{maxrep}{an integer, indicating the maximal number of repetitions allowed. Default is 5 (one business week).}
}
\value{
an index vector of the columns with more than maxrep repetitions in the time series.
}
\description{
Finds repeated values in the data and returns the respective column index.
}
\examples{
data(sp500)
sp500[,1] <- as.Date(sp500[,1],format="\%d.\%m.\%Y",stringsAsFactors=FALSE)
sp500 <- sp500[,-which(substr(colnames(sp500),1,7)=="X.ERROR")]
NYSE_hol <- as.Date(timeDate::holidayNYSE(as.numeric(unique(format(sp500[,1],format="\%Y")))))
no_trades <- sort(c(NYSE_hol, as.Date(c("2001-09-11","2001-09-12",
"2001-09-13","2001-09-14"), format="\%Y-\%m-\%d")))
sp500 <- sp500[-match(no_trades,sp500[,1]),]
nonas <- which(apply(is.na(sp500[,-1]),2,sum)==0)
sp500 <- sp500[,c(1, nonas+1)]
repindex <- find_repval(sp500[,-1])

}
