#' Sample Expected Return Estimation
#'
#' Computes the sample estimator of the expected return.
#'
#' @param data an nxp data matrix
#' @return a vector of expected returns
#'
#' @details The sample estimator of the expected returns for an \eqn{n\times p} data matrix X are computed with the following formula:
#' \deqn{\hat{\mu}=\frac{1}{n}\sum_{i=1}^{n}X_{i}} (for \eqn{i=1,\ldots, n}).
#' @examples
#' data(sp500_rets)
#' sp_rets <- sp500_rets[,-1]
#' mu_sample <- mu_estim_sample(sp_rets)
#'
#' @export mu_estim_sample
#'

mu_estim_sample <- function(data) {
  data <- as.matrix(data)
  mu <- as.numeric(colMeans(data, na.rm = TRUE))

  return(mu)
}

#' Expected Return Estimation with Factors
#'
#' Computes the estimator of the expected return according to a factor model.
#'
#' @param data an nxp data matrix
#' @param rf a numeric vector or double for the risk-free return
#' @param factors an nxf data frame of factors, e.g. the Fama-French 5 Research Factors
#' @return a vector of expected returns
#'
#' @details The expected returns are calculated according to a factor model, estimated with an OLS regression.
#' For a set of factors F, the expected returns for an \eqn{n\times p} data matrix X are defined as
#' \deqn{\hat{\mu}= \hat{B}'E[F],} where \eqn{E[F]=\frac{1}{n}\sum_{i=1}^{n}F_{i}} and \eqn{\hat{B}} are the estimated
#' beta coefficients from the linear regression.
#'
#' @examples
#' data(sp500_rets)
#' data(ff_factors)
#' sp_rets <- sp500_rets[,-1]
#' ff_factors <- ff_factors[,-1]
#' rf <- ff_factors[,ncol(ff_factors)]
#' ff_factors <- ff_factors[,-ncol(ff_factors)]
#'
#' mu_sample <- mu_estim_fm(sp_rets, rf, ff_factors)
#'
#' @export mu_estim_fm
#'

mu_estim_fm <- function(data, rf, factors) {
  data <- as.matrix(data)
  factors <- as.data.frame(factors)
  rets_ex <- data - rf
  mu <- c()

  for (i in 1:ncol(rets_ex)) {
    fflm_data <- cbind(rets_ex[, i], factors)
    colnames(fflm_data) <- c("ExRet", colnames(factors))
    form_lm <-
      paste('ExRet ~', paste(colnames(factors), collapse = '+'))
    model_lm <- stats::lm(form_lm, data = fflm_data)
    betas <- model_lm$coefficients[-1]
    mu[i] <- as.numeric(betas %*% colMeans(factors))
  }

  return(mu)
}

#' Wrapper Function for Expected Returns Estimation
#'
#' Estimates the expected returns of a dataset
#' according to the user-defined function.
#'
#' @param data an nxp data matrix.
#' @param estim_func an estimation function.
#' @param ... additional parameters, parsed to estim_func.
#' @return a vector of expected returns
#'
#' @examples
#' data(sp500_rets)
#' sp_rets <- sp500_rets[,-1]
#' sigma_ml <- mu_estim_wrapper(sp_rets, estim_func=mu_estim_sample)
#'
#'
#' @export mu_estim_wrapper
#'
mu_estim_wrapper <- function(data, estim_func, ...) {
  mu <- estim_func(data, ...)
  return(mu)
}
