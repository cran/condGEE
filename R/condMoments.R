

#' @title K1.norm
#' @description \code{E(Y|Y>w)} where Y is normal 
#' @author David Clement
#' @param w real value
#' @importFrom stats dnorm pnorm
#' @return conditional expectation
#' @export
K1.norm <- function(w) {
  return(dnorm(w) / (1 - pnorm(w)))
}


#' @title K2.norm
#' @description \code{E(Y^2|Y>w)} where Y is normal 
#' @author David Clement
#' @param w real value
#' @importFrom stats dnorm pnorm
#' @return conditional expectation
#' @export 
K2.norm <- function(w) {
  return(1 + w * dnorm(w) / (1 - pnorm(w))) 
}



#' @title K1.t3
#' @description \code{E(Y|Y>w)} where Y is t dist with 3 df 
#' @author David Clement
#' @param w real value
#' @importFrom stats pt
#' @return conditional expectation
#' @export
K1.t3 <- function(w) {
  df <- 3
  t.sd <- sqrt(df / (df - 2))
  temp <- 9 / (2 * ((t.sd * w) ^ 2 + 3)) * 0.367552597
  return(temp / (t.sd * (1 - pt(t.sd * w, df))))
}


#' @title K2.t3
#' @description \code{E(Y^2|Y>w)} where Y is t dist with 3 df 
#' @author David Clement
#' @param w real value
#' @importFrom stats pt
#' @return conditional expectation
#' @export
K2.t3 <- function(w) {
  df <- 3
  t.sd <- sqrt(df / (df - 2))
  yy <- t.sd * w
  sq3 <- sqrt(3)
  yy2 <- yy ^ 2
  temp <- atan(yy / sq3) * (2 * sq3 * yy2 + 6 * sq3) - 
    6 * yy - sq3 * pi * yy2 - 3 * sq3 * pi
  temp <- -3 * temp / (4 * (yy2 + 3)) * 0.367552597
  return(temp / (t.sd ^ 2 * (1 - pt(t.sd * w, df))))
}


#' @title K1.exp
#' @description \code{E(Y|Y>w)} where Y is exponential dist with mean 0
#' and variance 1
#' @author David Clement
#' @param w real value
#' @return conditional expectation
#' @export
K1.exp <- function(w) {
  return(pmax(w + 1, 0))
}
  
#' @title K2.exp
#' @description \code{E(Y^2|Y>w)} where Y is exponential dist with mean 0
#' and variance 1
#' @author David Clement
#' @param w real value
#' @return conditional expectation
#' @export
K2.exp <- function(w) {
  w <- pmax(w, -1)
  return(w ^ 2 + 2 * w + 2)
}
