
#' @title condGEE
#' @description Solves for the mean parameters (\eqn{theta}), the
#' variance parameter (\eqn{\sigma^2}), and their asymptotic variance 
#' in a conditional GEE for recurrent event gap times, as described by 
#' Clement, D. Y. and Strawderman, R. L. (2009) 
#' @author David Clement
#' @param data matrix of data with one row for each gap time; the first column 
#' should be a subject ID, the second column the gap time, the third column a 
#' completeness indicator equal to 1 if the gap time is complete and 0 if the 
#' gap time is censored, and the remaining columns the covariates for use in the 
#' mean and variance functions
#' @param start vector containing initial guesses for the unknown parameter vector
#' @param mu.fn the specification for the mean of the gap time; the default is 
#' a linear combination of the covariates; the function should take two arguments
#' (\code{theta}, and a matrix of covariates with each row corresponding to
#' one gap time) and it should return a vector of means
#' @param mu.d the derivative of \code{mu.fn} with respect to the parameter vector;
#' the default corresponds to a linear mean function
#' @param var.fn the specification for \eqn{V^2}, where the variance of the gap
#' time is \eqn{\sigma^2 V^2}; the default is a vector of ones; the function 
#' should take two arguments (\eqn{theta}, and a matrix of covariates with each
#' row corresponding to one gap time) and it should return a vector of variances
#' @param k1 the function to solve for the conditional mean length of the censored
#' gap times; its sole argument should be the vector of standardized (i.e.\
#'  \eqn{(Y-\mu)/(\sigma V)}) censored gap times; the default assumes the 
#'  standardized censored gap times follow a standard normal distribution, but 
#'  \code{K1.t3} and \code{K1.exp} are also provided in the package - they assume
#'   a standardized \emph{t} with 3 degrees of freedom and an exponential with 
#'   mean 0 and variance 1 respectively
#' @param k2 the function to solve for the conditional mean length of the square
#' of the censored gap times; its sole argument should be the vector of
#' standardized (i.e.\ \eqn{(Y-\mu)/(\sigma V)}) censored gap times; the default
#' assumes the standardized censored gap times follow a standard normal
#' distribution, but \code{K2.t3} and \code{K2.exp} are also provided in the
#' package - they assume a standardized \emph{t} with 3 degrees of freedom and 
#' an exponential with mean 0 and variance 1 respectively
#' @param robust logical, if \code{FALSE}, the mean and variance parameters are
#'  solved for simultaneously, increasing efficiency, but decreasing the leeway
#'   to misguess \code{start} and still find the root of the GEE
#' @param asymp.var logical, if \code{FALSE}, the function returns \code{NULL} 
#' for the asymptotic variance matrix
#' @param maxiter see \code{multiroot}; maximal number of iterations allowed
#' @param rtol see \code{multiroot}; relative error tolerance
#' @param atol see \code{multiroot}; absolute error tolerance 
#' @param ctol see \code{multiroot}; if between two iterations, the maximal 
#' change in the variable values is less than this amount, then it is assumed 
#' that the root is found
#' @param useFortran see \code{multiroot}; logical, if \code{FALSE}, then an
#'  \R implementation of Newton-Raphson is used
#' @importFrom rootSolve multiroot
#' @importFrom numDeriv jacobian
#' @return conditional expectation
#' @export
condGEE <- function(data, start, mu.fn = MU, mu.d = MU.d, var.fn = V, 
                    k1 = K1.norm, k2 = K2.norm, robust = TRUE, 
                    asymp.var = TRUE, maxiter = 100, rtol = 1e-6, 
                    atol = 1e-8, ctol = 1e-8, useFortran = TRUE) {
  
  p <- length(start)
  N <- length(data[, 1])
  n.covs <- length(data[1, ]) - 3
  
  uniq <- unique(data[, 1])
  n <- length(uniq)
  temp <- c(match(uniq, data[, 1]), length(data[, 1]) + 1)
  c <- temp[2:(n + 1)] - temp[1:n] 
  
  # default mean function is linear
  MU <- function(theta, covs)
   return(theta %*% t(cbind(rep(1,length(covs[, 1])), covs)))
  MU.d <- function(theta, covs)
   return(t(cbind(rep(1, length(covs[, 1])), covs)))
  
  # default V function is 1
  V <- function(theta, covs)
   return(rep(1, length(covs[, 1])))
  
  GEE.mean <- function(theta) {
    mu <- mu.fn(theta, matrix(data[, 4:(4 + n.covs - 1)], ncol = n.covs))
    var <- var.fn(theta, matrix(data[, 4:(4 + n.covs - 1)], ncol = n.covs))
    d.mu <- mu.d(theta, matrix(data[, 4:(4 + n.covs - 1)], ncol = n.covs))
    std.gaps <- as.numeric((data[, 2] - mu) / sqrt(sig2 * var))
    
    cens <- which(data[, 3] == 0)
    std.gaps[cens] <- k1(std.gaps[cens])
    
    tot <- as.vector(d.mu %*% diag(sqrt(1 / var), N) %*% std.gaps)
    
    return(tot / n)
  }
  
  GEE.var <- function(sig2) {
    mu <- mu.fn(theta, matrix(data[, 4:(4 + n.covs - 1)], ncol = n.covs))
    var <- var.fn(theta, matrix(data[, 4:(4 + n.covs - 1)], ncol = n.covs)) 
    std.gaps <- as.numeric((data[, 2] - mu) / sqrt(sig2 * var))
    std.gaps2 <- std.gaps ^ 2
    
    cens <- which(data[, 3] == 0)
    std.gaps2[cens] <- k2(std.gaps[cens])
    
    tot <- sum(std.gaps2) - N
    
    return(tot / n)
  }
  
  U <- function(eta) {
    theta <- eta[1:(p - 1)]
    sig2 <- eta[p]
    
    mu <- mu.fn(theta, matrix(data[,4:(4 + n.covs - 1)], ncol = n.covs))
    var <- var.fn(theta, matrix(data[, 4:(4 + n.covs - 1)], ncol = n.covs))
    d.mu <- mu.d(theta, matrix(data[, 4:(4 + n.covs - 1)], ncol = n.covs)) 
    std.gaps <- as.numeric((data[, 2] - mu) / sqrt(sig2 * var))
    std.gaps2 <- std.gaps ^ 2
    
    cens <- which(data[, 3] == 0)
    std.gaps2[cens] <- k2(std.gaps[cens])
    std.gaps[cens] <- k1(std.gaps[cens])
    
    tot.mean <- as.vector(d.mu %*% diag(sqrt(1 / var), N) %*% std.gaps)
    tot.var <- sum(std.gaps2) - N
    
    return(c(tot.mean, tot.var) / n)
  }
  
  var.S <- function(eta) {
    theta <- eta[1:(p - 1)]
    sig2 <- eta[p]
    tot <- matrix(rep(0, p ^ 2), p, p)
    
    for(i in 1:n) {
      stop.index <- sum(c[1:i])
      start.index <- stop.index - c[i] + 1
      ind <- (start.index:stop.index)
      len <- length(ind)
    
      mu <- mu.fn(theta, matrix(data[ind, 4:(4 + n.covs - 1)], ncol = n.covs))
      var <- var.fn(theta, matrix(data[ind, 4:(4 + n.covs - 1)], ncol = n.covs))        
      d.mu <- mu.d(theta, matrix(data[ind, 4:(4 + n.covs - 1)], ncol = n.covs))
     
      std.gaps <- as.numeric((data[ind, 2] - mu) / sqrt(sig2 * var))
      std.gaps2 <- std.gaps ^ 2
    
      cens <- which(data[ind, 3] == 0)
      if(length(cens) > 0) {
        std.gaps2[cens] <- k2(std.gaps[cens])
        std.gaps[cens] <- k1(std.gaps[cens])
      }
    
      tot.mean <- as.vector(d.mu %*% diag(sqrt(1 / var), len) %*% std.gaps)
      tot.var <- sum(std.gaps2) - len
      temp <- c(tot.mean, tot.var)
    
      tot <- tot + temp %*% t(temp)
    }
    
    return(tot / n)
  }
  
  # doing mean and var one by one is slower but more robust
  if(robust==TRUE) {
   conv <- 1
   theta <- start[1:(p-1)]
   sig2 <- start[p]
   old <- rep(0, p - 1)
   while(conv > 1e-4) {
     theta <- multiroot(GEE.mean, theta, maxiter, rtol, atol, ctol, useFortran)$root
     sig2 <- multiroot(GEE.var, sig2, maxiter, rtol, atol, ctol, useFortran)$root
     conv <- sum((theta-old) ^ 2)
     old <- theta
   }
   eta <- c(theta, sig2)
  } else { 
    eta <- multiroot(U, start, maxiter, rtol, atol, ctol, useFortran)$root 
  }
  
  a.var <- NULL
  if(asymp.var == TRUE) {
   bread.inv <- jacobian(U, eta)
   bread <- solve(bread.inv)
   meat <- var.S(eta)
  
   a.var <- bread %*% meat %*% t(bread) / n
  }
  
  return(list(eta = eta, a.var = a.var))
}
