#' Find NAs
#'
#' Finds the NAs in the data and returns the respective column index.
#'
#' @param data a numerical vector or a numerical nxp data matrix.
#' @param all a logical, if TRUE (default), the function finds columns which consist only of NAs.
#' If FALSE, it finds the columns with any number of NAs.
#'
#' @return an index vector of the position of NAs.
#'
#' @examples
#' data(sp500)
#' index_nas <- findNAs(sp500, all=FALSE)
#'
#' @export findNAs
#'
findNAs <- function(data, all = TRUE) {
  if (is.null(dim(data))) {
    indx <- which(is.na(data))
  } else{
    if (all) {
      indx <-
        which(sapply(apply(is.na(data), 2, function(x)
          which(all(x))), length) != 0)
    } else{
      indx <- which(sapply(apply(is.na(data), 2, which), length) != 0)
    }
  }
  return(indx)
}

#' Find Weekends
#'
#' Finds the weekends in a date vector and returns a respective index.
#'
#' @param dates a date vector.
#'
#' @return an index vector for the position of weekends.
#'
#' @examples
#' data(sp500)
#' dates <- as.Date(sp500[,1], format="%d.%m.%Y", stringsAsFactors=FALSE)
#' index_weekends <- findWeekends(dates)
#'
#' @export findWeekends
#'
findWeekends <- function(dates) {
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  indx <-
    which(weekdays(dates) == "Saturday" |
            weekdays(dates) == "Sunday")
  return(indx)
}


#' Find Repeated Values
#'
#' Finds repeated values in the data and returns the respective column index.
#'
#' @param data an nxp numerical data matrix.
#' @param maxrep an integer, indicating the maximal number of repetitions allowed. Default is 5 (one business week).
#'
#' @return an index vector of the columns with more than maxrep repetitions in the time series.
#'
#' @examples
#' data(sp500)
#' sp500[,1] <- as.Date(sp500[,1],format="%d.%m.%Y",stringsAsFactors=FALSE)
#' sp500 <- sp500[,-which(substr(colnames(sp500),1,7)=="X.ERROR")]
#' NYSE_hol <- as.Date(timeDate::holidayNYSE(as.numeric(unique(format(sp500[,1],format="%Y")))))
#' no_trades <- sort(c(NYSE_hol, as.Date(c("2001-09-11","2001-09-12","2001-09-13","2001-09-14"), format="%Y-%m-%d")))
#' sp500 <- sp500[-match(no_trades,sp500[,1]),]
#' nonas <- which(apply(is.na(sp500[,-1]),2,sum)==0)
#' sp500 <- sp500[,c(1, nonas+1)]
#' repindex <- findRepVal(sp500[,-1])
#'
#' @export findRepVal
#'
findRepVal <- function(data, maxrep = 5) {
  max_rep_number <- c()
  for (j in 1:ncol(data)) {
    values <- data[which(!is.na(data[, j])), j]
    counter <- 0
    max_rep_number[j] <- -1
    for (i in 2:length(values)) {
      diff_values <- values[i] - values[i - 1]
      if (diff_values == 0) {
        max_rep_number[j]  <- max(max_rep_number[j], counter)
        counter <- counter + 1
      } else{
        max_rep_number[j]  <- max_rep_number[j]
        counter <- 0
      }
    }
  }
  repval_indx <- which(max_rep_number >= maxrep)
  return(repval_indx)
}


#' Find Business Days
#'
#' Finds the first or last business days in a month and returns the respective date vector.
#'
#' @param dates a date vector.
#' @param type a character vector. type="first" (default) delivers the first business day in a month and
#' type="last" the last business day in a month.
#'
#' @return a date vector with the respective first or last business days for each month in dates.
#'
#' @examples
#' data(sp500)
#' dates <- as.Date(sp500[,1], format="%d.%m.%Y", stringsAsFactors=FALSE)
#' firstdays <- findBusinessDay(dates)
#' lastdays <- findBusinessDay(dates, type="last")
#'
#' @export findBusinessDay
#'
findBusinessDay <- function(dates, type = "first") {
  day <- format(dates, format = "%d")
  monthYr <- format(dates, format = "%Y-%m")
  if (type == "first") {
    y <- tapply(day, monthYr, min)
  } else if (type == "last") {
    y <- tapply(day, monthYr, max)
  } else{
    print("Not recognized type. Enter either first or last.")
  }
  businessDay <- as.Date(paste(row.names(y), y, sep = "-"))
  return(businessDay)
}

#' Returns Calculation
#'
#' Calculates the returns from a prices dataset.
#'
#' @param prices a numerical vector or an nxp data matrix with stock prices.
#' @param type a character vector, indicating which type of returns are to be calculated.
#' type="d" calculates the discrete returns (default) and type="c" calculates the continuous returns.
#'
#' @return a numerical vector or an (n-1)xp data matrix with returns.
#'
#' @examples
#' data(sp500)
#' sp500[,1] <- as.Date(sp500[,1],format="%d.%m.%Y",stringsAsFactors=FALSE)
#' sp500 <- sp500[,-which(substr(colnames(sp500),1,7)=="X.ERROR")]
#' NYSE_hol <- as.Date(timeDate::holidayNYSE(as.numeric(unique(format(sp500[,1],format="%Y")))))
#' no_trades <- sort(c(NYSE_hol, as.Date(c("2001-09-11","2001-09-12","2001-09-13","2001-09-14"),
#' format="%Y-%m-%d")))
#' sp500 <- sp500[-match(no_trades,sp500[,1]),]
#' nonas <- which(apply(is.na(sp500[,-1]),2,sum)==0)
#' sp500 <- sp500[,c(1, nonas+1)]
#' repindex <- findRepVal(sp500[,-1])
#' sp500_prices <- sp500[,-c(1, repindex+1)]
#' sp500_ret <- retsCalc(sp500_prices)
#'
#' @export retsCalc
#'
retsCalc <- function(prices, type = "d") {
  if (type == "d") {
    if (is.null(dim(prices))) {
      ret <- diff(prices) / prices[1:(length(prices) - 1)]
    } else{
      ret <- apply(prices, 2, diff) / prices[1:(nrow(prices) - 1), ]
    }
  } else if (type == "c") {
    if (is.null(dim(prices))) {
      ret <- diff(log(prices))
    } else{
      ret <- apply(log(prices), 2, diff)
    }
  } else{
    cat(
      "Invalid type for calculation of returns. Please enter either c for continuous or d for discrete returns."
    )
  }
  return(ret)
}
