#' Naive Portfolio Diversification I
#'
#' Calculates the weights of a naive portfolio strategy.
#'
#' @param Sigma a pxp covariance matrix of asset returns.
#'
#' @return a vector of length p with the weights of the portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' port_estim_naive(covMat)
#'
#' @export port_estim_naive
#'
port_estim_naive <- function(Sigma) {
  n <- dim(Sigma)[1]
  weights <- as.numeric(rep.int(1 / n, n))

  return(weights)
}

#' Naive Portfolio Diversification II
#'
#' Calculates the weights of a naive portfolio strategy.
#'
#' @param p an integer, specifying the number of stocks.
#'
#' @return a vector of length p with the weights of the portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' port_estim_naive_slim(dim(example_rets)[2])
#'
#' @export port_estim_naive_slim
#'
port_estim_naive_slim <- function(p) {
  weights <- as.numeric(rep.int(1 / p, p))

  return(weights)
}

#' Global Minimum-Variance Portfolio Optimization
#'
#' Calculates the weights of a global minimum-variance portfolio strategy.
#'
#' @param Sigma a pxp covariance matrix of asset returns.
#'
#' @return a vector of length p with the weights of the portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' port_estim_gmv(covMat)
#'
#' @export port_estim_gmv
#'
port_estim_gmv <- function(Sigma) {
  Sigma <- as.matrix(Sigma)
  n <- dim(Sigma)[1]
  ones <- rep.int(1, n)
  SigmaInv <- solve(Sigma)

  weights <-
    as.numeric(SigmaInv %*% ones / as.numeric(ones %*% SigmaInv %*% ones))

  return(weights)
}

#' Tangency Portfolio Optimization
#'
#' Calculates the weights of a tangency portfolio strategy.
#'
#' @param Sigma a pxp matrix, the covariance matrix of asset returns.
#' @param mu a vector, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#'
#' @return a vector of length p with the weights of the portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' muVec <- colMeans(example_rets)
#' port_estim_tang(covMat, muVec)
#'
#' @export port_estim_tang
#'
port_estim_tang <- function(Sigma, mu, rf = 0) {
  Sigma <- as.matrix(Sigma)
  n <- dim(Sigma)[1]
  mueff <- mu - rf
  ones <- rep.int(1, n)
  SigmaInv <- solve(Sigma)

  weights <-
    as.numeric(SigmaInv %*% mueff / (as.numeric(ones %*% SigmaInv %*% mueff)))

  return(weights)
}

#' Markowitz Portfolio Diversification
#'
#' Calculates the weights of an optimal Markowitz portfolio strategy.
#'
#' @param Sigma a pxp covariance matrix of asset returns.
#' @param mu a vector, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#'
#' @return a vector of length p with the weights of the portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' muVec <- colMeans(example_rets)
#' port_estim_markowitz(covMat, muVec)
#'
#' @export port_estim_markowitz
#'
port_estim_markowitz <- function(Sigma,
                                 mu,
                                 rf = 0,
                                 gamma = 2) {
  Sigma <- as.matrix(Sigma)
  SigmaInv <- solve(Sigma)
  mueff <- mu - rf

  weights <- as.numeric((1 / gamma) * SigmaInv %*% mueff)

  return(weights)
}

#' Unbiased Markowitz Portfolio Diversification
#'
#' Calculates the weights of an optimal unbiased Markowitz portfolio strategy.
#'
#' @param Sigma a pxp covariance matrix of asset returns.
#' @param mu a vector, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#' @param size_sample an integer, the sample size (number of observations) of the returns data.
#'
#' @return a vector of length p with the weights of the portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' muVec <- colMeans(example_rets)
#' port_estim_markowitz_unb(covMat, muVec, size_sample=dim(example_rets)[1])
#'
#' @export port_estim_markowitz_unb
#'
port_estim_markowitz_unb <-
  function(Sigma,
           mu,
           rf = 0,
           gamma = 2,
           size_sample) {
    TT <- size_sample
    Sigma <- as.matrix(Sigma)
    n <- dim(Sigma)[1]
    unb <- (TT - n - 2) / TT

    weights <-
      as.numeric(unb * port_estim_markowitz(Sigma, mu, rf, gamma))

    return(weights)
  }

#' Eigenportfolio
#'
#' Calculates the weights of an eigenportfolio strategy.
#'
#' @param Sigma a pxp covariance matrix of asset returns.
#'
#' @return a vector of length p with the weights of the portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' port_estim_eigen(covMat)
#'
#' @export port_estim_eigen
#'
port_estim_eigen <- function(Sigma) {
  Sigma <- as.matrix(Sigma)
  p <- dim(Sigma)[2]
  volas <- as.numeric(sqrt(diag(Sigma)))
  Volamat_inv <- solve(diag(volas))
  Corrmat <- Volamat_inv %*% Sigma %*% Volamat_inv
  eigenvec <- eigen(Corrmat)$vectors
  eigenport_weights <- matrix(NA, p, p)

  for (i in 1:p) {
    eigenport_weights[, i] <- eigenvec[, i] / volas[i]
    eigenport_weights[, i] <-
      eigenport_weights[, i] / as.numeric(rep.int(1, p) %*% eigenport_weights[, i])
  }

  weights <- as.numeric(eigenport_weights[, 1])

  return(weights)
}


#' Nummerical Portfolio Optimization with constraints
#'
#' Calculates the weights of a GMV or tangency portfolio strategy with different constraints and/or penalties.
#'
#' @param Sigma a pxp matrix, the covariance matrix of asset returns.
#' @param mu a vector, the expected returns. Default value is NULL.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the assumed risk-averse parameter. Default value is 2.
#' @param Aeq a cxp equality constraint matrix, containing c constraints for p regressors. Default value is Aeq=NULL, no equality constraints.
#' @param beq a cx1 equality constraint vector. Default value is beq=NULL, no equality constraints.
#' @param A a cxp inequality constraint matrix, containing c constraints for p regressors. Default value is A=NULL, no inequality constraints.
#' @param b a cx1 inequality constraint vector. Default value is b=NULL, no inequality constraints.
#' @param lambda1 a tuning parameter value for the lasso penalty. Default value is lambda1=0.
#' @param penidx1 a logical px1 vector, indicating which coefficients are to be penalized with the lasso penalty lambda1. Default value is penidx1=NULL and imposes penalty on all p coefficients.
#' @param lambda2 a tuning parameter value for the ridge penalty. Default value is lambda2=0.
#' @param penidx2 a logical px1 vector, indicating which coefficients are to be penalized with the ridge penalty lambda2. Default value is penidx2=NULL and imposes penalty on all p coefficients.
#' @param gross_c a double, the requied gross exposure constraint. Default value is NULL (no constraint). Attention! Works only with GMV.
#' @param porttype a character string. Possible values are "Tang" (tangency portfolio) and "GMV" (global minimum variance portfolio). Default value is "GMV".
#' @param zero_tol a double, indicating the zero tolerance for the calculated weights. Default value is 1e-7.
#' @param res_all a logical. If TRUE, the result includes the calculated weights and the duals from the optimization. If FALSE, only the weights. Default value is FALSE.
#'
#' @details The portfolio optimization with constraints minimizes
#' \deqn{0.5||y - X \beta ||^2_2 + \lambda_1||\beta||_1 + \lambda_2||\beta||^2_2,}
#' subject to \eqn{Aeq \beta = beq} and \eqn{A \beta\le b}.
#'
#' @return a vector of length p with the weights of the portfolio.
#' @return a list with the weights vector and the duals from the optimization.
#'
#' @import CVXR
#' @import ROI
#' @import ROI.plugin.qpoases
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' muVec <- colMeans(example_rets)
#' port_estim_solver(covMat)
#' port_estim_solver(covMat, muVec, porttype="Tang")
#'
#' @export port_estim_solver
#'
port_estim_solver <-
  function(Sigma,
           mu = NULL,
           rf = 0,
           gamma = 2,
           Aeq = NULL,
           beq = NULL,
           A = NULL,
           b = NULL,
           lambda1 = 0,
           penidx1 = NULL,
           lambda2 = 0,
           penidx2 = NULL,
           gross_c = NULL,
           porttype = "GMV",
           zero_tol = 1e-7,
           res_all = FALSE) {
    Sigma <- as.matrix(Sigma)
    n <- dim(Sigma)[1]

    if (is.null(gross_c)) {
      Ident <- diag(1, n)

      #default values penidx1 and penidx2 (all entries are penalized)
      if (is.null(penidx1)) {
        penidx1 <- matrix(TRUE, n, 1)
      }
      dim(penidx1) <- c(n, 1)

      if (is.null(penidx2)) {
        penidx2 <- matrix(TRUE, n, 1)
      }
      dim(penidx2) <- c(n, 1)

      if (porttype == "GMV") {
        #default values Aeq and beq (sum constraint)
        if (is.null(Aeq)) {
          Aeq <- matrix(1, 1, n)
          beq <- rep(1, 1)
        }

        #number of equality constraints
        m1 <- dim(Aeq)[1]

        #default values A and b (no constraint)

        if (is.null(A)) {
          A <- matrix(NA, 0, n)
          b <- rep(0, 0)
        }

        #number of inequality constraints
        m2 <- dim(A)[1]

        #gamma value
        if (gamma != 2) {
          gamma <- 2
          cat("Gamma is different than 2 for GMV. Gamma set to 2.")
        }

        # quadratic coefficient
        H <-
          rbind(
            cbind((gamma / 2) * Sigma + lambda2 * Ident,
                  -(gamma / 2) * Sigma - lambda2 * Ident
            ),
            cbind(
              -(gamma / 2) * Sigma - lambda2 * Ident,
              (gamma / 2) * Sigma + lambda2 * Ident
            )
          )

        f <- lambda1 * rbind(penidx1, penidx1)

        #  all constraints (first m1 constraints are equality constraints)
        Amatrix <- rbind(cbind(Aeq, -Aeq), cbind(A, -A))
        bvector <- c(beq, b)

        # optimizer
        opt_problem <-
          ROI::OP(ROI::Q_objective(H, L = t(f)),
                  ROI::L_constraint(
                    L = Amatrix,
                    dir = c(rep("==", m1), rep("<=", m2)),
                    rhs = bvector
                  ))
        opt <- ROI::ROI_solve(opt_problem, solver = "qpoases")
        opt_sol <-  opt$message$primal_solution

        # duals
        duals <- as.numeric(opt$message$dual_solution[-(1:(2 * n))])

        # calculation of actual weights (wpos - wneg)
        betahat <-
          matrix(opt_sol[1:n] - opt_sol[(n + 1):length(opt_sol)], n, 1)

        # round the solutions
        betahat[which(abs(betahat) < zero_tol)] <- 0


      } else if (porttype == "Tang") {
        if (is.null(mu)) {
          stop(cat(
            "Please, enter a value for mu, the expected returns vector."
          ))
        }

        mueff <- as.matrix(mu - rf)

        #default values Aeq and beq (no constraint)
        if (is.null(Aeq)) {
          Aeq <- matrix(NA, 0, n)
          beq <- rep(0, 0)
        }

        #number of equality constraints
        m1 <- dim(Aeq)[1]

        #default values A and b (no constraint)

        if (is.null(A)) {
          A <- matrix(NA, 0, n)
          b <- rep(0, 0)
        }

        #number of inequality constraints
        m2 <- dim(A)[1]

        # quadratic coefficient
        H <-
          rbind(
            cbind((gamma / 2) * Sigma + lambda2 * Ident,
                  -(gamma / 2) * Sigma - lambda2 * Ident
            ),
            cbind(
              -(gamma / 2) * Sigma - lambda2 * Ident,
              (gamma / 2) * Sigma + lambda2 * Ident
            )
          )

        f <-
          rbind(-mueff, mueff) + lambda1 * rbind(penidx1, penidx1)

        #  all constraints (first m1 constraints are equality constraints)
        Amatrix <- rbind(cbind(Aeq, -Aeq), cbind(A, -A))
        bvector <- c(beq, b)

        # optimizer
        opt_problem <-
          ROI::OP(ROI::Q_objective(H, L = t(f)),
                  ROI::L_constraint(
                    L = Amatrix,
                    dir = c(rep("==", m1), rep("<=", m2)),
                    rhs = bvector
                  ))
        opt <- ROI::ROI_solve(opt_problem, solver = "qpoases")
        opt_sol <-  opt$message$primal_solution

        #duals
        duals <- as.numeric(opt$message$dual_solution[-(1:(2 * n))])

        # calculation of actual weights (wpos - wneg)
        betahat <-
          matrix(opt_sol[1:n] - opt_sol[(n + 1):length(opt_sol)], n, 1)

        # round the solutions
        betahat[which(abs(betahat) < zero_tol)] <- 0

        # normalization
        if (sum(betahat) != 0) {
          betahat <- betahat / sum(betahat)
        } else{
          betahat <- betahat
        }

      }

      if (res_all) {
        return(list(as.numeric(betahat), as.numeric(duals)))

      } else {
        return(as.numeric(betahat))
      }

    } else{
      w_var <- CVXR::Variable(n)
      obj_func <- CVXR::Minimize(CVXR::quad_form(w_var, Sigma))
      constr_gross <- sum(abs(w_var)) <= gross_c
      constr_sum <- sum(w_var) == 1
      problem_form <-
        CVXR::Problem(obj_func, constraints = list(constr_gross, constr_sum))
      result <- CVXR::solve(problem_form)
      weights <- result$getValue(w_var)
      return(as.numeric(weights))
    }
  }

#' Wrapper Function for Portfolio Optimization I
#'
#' Allows the execution of all included optimization functions.
#'
#' @param Sigma a pxp matrix, the covariance matrix of asset returns
#' @param estim_func a function for portfolio optimization
#' @param ... additional arguments to be passed to estim_func
#'
#' @return returns the value of the executed function estim_func
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' port_estim_wrapper(covMat, port_estim_gmv)
#'
#' muVec <- colMeans(example_rets)
#' port_estim_wrapper(covMat, port_estim_tang, muVec)
#'
#' @export port_estim_wrapper
#'
port_estim_wrapper <- function(Sigma, estim_func, ...) {
  return(estim_func(Sigma, ...))
}


#' Wrapper Function for Portfolio Optimization II
#'
#' Allows the execution of specific optimization functions.
#'
#' @param Sigma a pxp matrix, the covariance matrix of asset returns
#' @param est_type a function for portfolio optimization
#' @param ... additional arguments to be passed to the respective function, according to est_type.
#'
#' @return returns the weights, estimated according the est_type portfolio.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' port_estim(covMat, "GMV")
#'
#' muVec <- colMeans(example_rets)
#' port_estim(covMat, "Tang", muVec)
#'
#' @export port_estim
#'
port_estim <- function(Sigma, est_type = "Naive", ...) {
  if (est_type == "Naive") {
    weights <- port_estim_naive(Sigma)
  } else if (est_type == "GMV") {
    weights <- port_estim_gmv(Sigma)
  } else if (est_type == "Tang") {
    weights <- port_estim_tang(Sigma, ...)
  } else if (est_type == "Markowitz") {
    weights <- port_estim_markowitz(Sigma, ...)
  } else if (est_type == "Markowitz_unb") {
    weights <- port_estim_markowitz_unb(Sigma, ...)
  } else if (est_type == "Eigen") {
    weights <- port_estim_eigen(Sigma)
  } else if (est_type == "Solver") {
    weights <- port_estim_solver(Sigma, ...)
  }

  return(as.numeric(weights))
}

#' Estimation of the Efficient Frontier
#'
#' Calculates the weights and the standard deviations along the Efficient Frontier.
#'
#' @param Sigma a pxp covariance matrix of asset returns.
#' @param mu a vector, the expected returns.
#' @param mugrid a vector with length m as the grid of expected returns, along which the capital market line is to be estimated.
#' @param res_all a logical. If TRUE, the result includes the calculated weights and the standard deviations for the CML.
#' If FALSE, only the weights. Default value is FALSE.
#'
#' @return a pxm matrix with the weights of the CML portfolio along mugrid.
#' @return a vector of length m with the corresponding standard deviations.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' muVec <- colMeans(example_rets)
#' mugrid <- seq(0, 0.2, by=0.001)
#' results <- port_estim_efficient(covMat, muVec, mugrid, res_all=FALSE)
#'
#' @export port_estim_efficient
port_estim_efficient <-
  function(Sigma, mu, mugrid, res_all = FALSE) {
    n <- dim(Sigma)[1]
    ones <- rep.int(1, n)
    SigmaInv <- solve(Sigma)
    A <- as.numeric(ones %*% SigmaInv %*% mu)
    B <- as.numeric(mu %*% SigmaInv %*% mu)
    C <- as.numeric(ones %*% SigmaInv %*% ones)
    g <-
      as.numeric(1 / (C * B - A ^ 2) * (B * SigmaInv %*% ones  - A * SigmaInv %*% mu))
    h <-
      1 / (C * B - A ^ 2) * (C * SigmaInv %*% mu  - A * SigmaInv %*% ones)

    weights <- as.numeric(h %*% mugrid + g)
    sds <-
      as.numeric(sqrt(C / (C * B - A ^ 2) * (mugrid - A / C) ^ 2 + 1 /
                        C))
    if (res_all) {
      return(list(weights, sds))
    } else{
      return(weights)
    }
  }

#' Estimation of the Capital Market Line
#'
#' Calculates the weights and the standard deviations along the Capital Market Line (CML).
#'
#' @param Sigma a pxp covariance matrix of asset returns.
#' @param mu a vector, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param mugrid a vector with length m as the grid of expected returns, along which the capital market line is to be estimated.
#' @param res_all a logical. If TRUE, the result includes the calculated weights and the standard deviations for the CML.
#' If FALSE, only the weights. Default value is FALSE.
#'
#' @return a pxm matrix with the weights of the CML portfolio along mugrid.
#' @return a vector of length m with the corresponding standard deviations.
#'
#' @examples
#' data(sp500_rets)
#' example_rets <- sp500_rets[,2:11]
#' covMat <- var(example_rets)
#' muVec <- colMeans(example_rets)
#' mugrid <- seq(0, 0.2, by=0.001)
#' results <- port_estim_cml(covMat, muVec, rf=0, mugrid, res_all=FALSE)
#'
#' @export port_estim_cml
port_estim_cml <-
  function(Sigma,
           mu,
           rf = 0,
           mugrid,
           res_all = FALSE) {
    mueff <- mu - rf
    weights_tang <- port_estim_tang(Sigma, mu, rf)
    mu_tang <- as.numeric(weights_tang %*% mueff)
    sd_tang <-
      as.numeric(sqrt(t(weights_tang) %*% Sigma %*% weights_tang))
    cml_b <- sd_tang / (mu_tang - rf)
    cml_a <- -rf * cml_b
    sds <- cml_a + cml_b * mugrid
    weights <- as.matrix(mugrid / mu_tang) %*% t(weights_tang)

    if (res_all) {
      return(list(weights, sds))
    } else{
      return(weights)
    }
  }




# KanZhou-2Fund
port_estim_KZ2f <- function(Sigma,
                            mu,
                            rf = 0,
                            gamma = 2,
                            size_sample) {
  TT <- size_sample
  Sigma <- as.matrix(Sigma)
  n <- dim(Sigma)[1]
  SigmaInv <- solve(Sigma)
  mueff <- mu - rf
  sharpesquared <- as.numeric(mueff %*% SigmaInv %*% mueff)
  c3 <- (TT - n - 1) * (TT - n - 4) / (TT * (TT - 2))
  cstar <- c3 * sharpesquared / (sharpesquared + (n / TT))

  weights <- as.numeric((1 / gamma) * cstar * SigmaInv %*% mueff)

  return(weights)
}


# KanZhou-3Fund
port_estim_KZ3f <- function(Sigma,
                            mu,
                            rf = 0,
                            gamma = 2,
                            size_sample) {
  TT <- size_sample
  Sigma <- as.matrix(Sigma)
  n <- dim(Sigma)[1]
  ones <- rep.int(1, n)
  mueff <- mu - rf
  SigmaInv <- solve(Sigma)
  c3 <- (TT - n - 1) * (TT - n - 4) / (TT * (TT - 2))
  mug <-
    as.numeric(mueff %*% SigmaInv %*% ones) / as.numeric(ones %*% SigmaInv %*%
                                                           ones)
  psisquared <-
    as.numeric(mueff %*% SigmaInv %*% mueff) - (as.numeric((mueff %*% SigmaInv %*%
                                                              ones) ^ 2) / as.numeric(ones %*% SigmaInv %*% ones))

  weights <-
    as.numeric((c3 / gamma) * ((psisquared / (psisquared + n / TT)) * SigmaInv %*%
                                 mueff + ((n / TT) / (psisquared + n / TT)) * mug * SigmaInv %*% ones
    ))

  return(weights)
}

# TuZhou-Markowitz - unbiased
port_estim_TZMarkowitz <-
  function(Sigma,
           mu,
           rf = 0,
           gamma = 2,
           size_sample) {
    TT <- size_sample
    Sigma <- as.matrix(Sigma)
    weights_Markowitzunb <-
      port_estim_markowitz_unb(Sigma, mu, rf, gamma, TT)
    n <- dim(Sigma)[1]
    mueff <- mu - rf
    SigmaInv <- solve(Sigma)
    weights_naive <- port_estim_naive(Sigma)
    sharpesquared <- as.numeric(mueff %*% SigmaInv %*% mueff)
    pi1 <-
      as.numeric(
        as.numeric(weights_naive %*% Sigma %*% weights_naive) - ((2 / gamma) * as.numeric(weights_naive %*% mueff)) + ((1 /
                                                                                                                          (gamma ^ 2)) * sharpesquared)
      )
    c1 <- ((TT - 2) * (TT - n - 2)) / ((TT - n - 1) * (TT - n - 4))
    pi2 <-
      as.numeric(((1 / (gamma ^ 2)) * (c1 - 1) * sharpesquared) + (c1 / (gamma ^
                                                                           2)) * (n / TT))
    delta_TZMarkowitz <- pi1 / (pi1 + pi2)

    weights <-
      as.numeric((1 - delta_TZMarkowitz) * weights_naive + delta_TZMarkowitz * weights_Markowitzunb)

    return(weights)
  }


# TuZhou-KanZhou3Fund
port_estim_TZKZ3f <- function(Sigma,
                              mu,
                              rf = 0,
                              gamma = 2,
                              size_sample) {
  TT <- size_sample
  Sigma <- as.matrix(Sigma)
  weights_KZ3f <- port_estim_KZ3f(Sigma, mu, rf, gamma, TT)
  n <- dim(Sigma)[1]
  mueff <- mu - rf
  SigmaInv <- solve(Sigma)
  weights_naive <- port_estim_naive_slim(n)
  ones <- rep.int(1, n)
  mug <-
    as.numeric(mueff %*% SigmaInv %*% ones) / as.numeric(ones %*% SigmaInv %*%
                                                           ones)
  psisquared <-
    as.numeric(mueff %*% SigmaInv %*% mueff) - (as.numeric((mueff %*% SigmaInv %*%
                                                              ones) ^ 2) / as.numeric(ones %*% SigmaInv %*% ones))
  sharpesquared <- as.numeric(mueff %*% SigmaInv %*% mueff)
  c1 <- ((TT - 2) * (TT - n - 2)) / ((TT - n - 1) * (TT - n - 4))
  pi1 <-
    as.numeric(
      as.numeric(weights_naive %*% Sigma %*% weights_naive) - ((2 / gamma) * as.numeric(weights_naive %*% mueff)) + ((1 /
                                                                                                                        (gamma ^ 2)) * sharpesquared)
    )
  pi3 <-
    as.numeric(((1 / (gamma ^ 2)) * sharpesquared) - (1 / (c1 * (gamma ^ 2))) *
                 (sharpesquared - ((n / TT) * psisquared)))
  pi13 <- as.numeric(((1 / (gamma ^ 2)) * sharpesquared) -
                       ((1 / gamma) * as.numeric(weights_naive %*% mueff)) +
                       (1 / (gamma * c1)) * ((
                         psisquared * as.numeric(weights_naive %*% mueff) + (1 - psisquared) * mug * as.numeric(weights_naive %*% ones)
                       ) - (1 / gamma) * (
                         psisquared * as.numeric(mueff %*% SigmaInv %*% mueff) + (1 - psisquared) * mug * as.numeric(mueff %*% SigmaInv %*% ones)
                       )
                       ))
  delta_TZKZ3f <- (pi1 - pi13) / (pi1 + pi3 - 2 * pi13)

  weights <-
    as.numeric((1 - delta_TZKZ3f) * weights_naive + delta_TZKZ3f * weights_KZ3f)

  return(weights)
}
