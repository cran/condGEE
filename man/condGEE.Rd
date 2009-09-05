\name{condGEE}
\Rdversion{2.9.0}
\alias{condGEE}
\title{
  Parameter estimation in conditional GEE for recurrent event gap times
}
\description{
  Solves for the mean parameters (\eqn{\theta}), the variance parameter (\eqn{\sigma^2}), and their asymptotic variance in a conditional GEE for recurrent event gap times, as described by Clement, D. Y. and Strawderman, R. L. (2009) \emph{Biostatistics} \bold{10}, 451--467.  Makes a parametric assumption for the length of the censored gap time, and assumes gap times within subject are conditionally uncorrelated.
}
\usage{
condGEE(data, start, mu.fn=MU, mu.d=MU.d, var.fn=V, 
  k1=K1.norm, k2=K2.norm, robust=TRUE, asymp.var=TRUE, 
  maxiter=100, rtol=1e-6, atol=1e-8, ctol=1e-8, useFortran=TRUE)
}
\arguments{
  \item{data }{matrix of data with one row for each gap time; the first column should be a subject ID, the second column the gap time, the third column a completeness indicator equal to 1 if the gap time is complete and 0 if the gap time is censored, and the remaining columns the covariates for use in the mean and variance functions
  }
  \item{start }{vector containing initial guesses for the unknown parameter vector
  }
  \item{mu.fn }{the specification for the mean of the gap time; the default is a linear combination of the covariates; the function should take two arguments (\eqn{\theta}, and a matrix of covariates with each row corresponding to one gap time) and it should return a vector of means
  }
  \item{mu.d }{the derivative of \code{mu.fn} with respect to the parameter vector; the default corresponds to a linear mean function
  }
  \item{var.fn }{the specification for \eqn{V^2}, where the variance of the gap time is \eqn{\sigma^2 V^2}; the default is a vector of ones; the function should take two arguments (\eqn{\theta}, and a matrix of covariates with each row corresponding to one gap time) and it should return a vector of variances
  }
  \item{k1 }{the function to solve for the conditional mean length of the censored gap times; its sole argument should be the vector of standardized (i.e.\ \eqn{(Y-\mu)/(\sigma V)}) censored gap times; the default assumes the standardized censored gap times follow a standard normal distribution, but \code{K1.t3} and \code{K1.exp} are also provided in the package - they assume a standardized \emph{t} with 3 degrees of freedom and an exponential with mean 0 and variance 1 respectively
  }
  \item{k2 }{the function to solve for the conditional mean length of the square of the censored gap times; its sole argument should be the vector of standardized (i.e.\ \eqn{(Y-\mu)/(\sigma V)}) censored gap times; the default assumes the standardized censored gap times follow a standard normal distribution, but \code{K2.t3} and \code{K2.exp} are also provided in the package - they assume a standardized \emph{t} with 3 degrees of freedom and an exponential with mean 0 and variance 1 respectively

  }
  \item{robust}{logical, if \code{FALSE}, the mean and variance parameters are solved for simultaneously, increasing efficiency, but decreasing the leeway to misguess \code{start} and still find the root of the GEE
  }
  \item{asymp.var }{logical, if \code{FALSE}, the function returns \code{NULL} for the asymptotic variance matrix
  }
  \item{maxiter }{see \code{multiroot}; maximal number of iterations allowed
  }
  \item{rtol }{see \code{multiroot}; relative error tolerance 
  }
  \item{atol }{see \code{multiroot}; absolute error tolerance 
  }
  \item{ctol }{see \code{multiroot}; if between two iterations, the maximal change in the variable values is less than this amount, then it is assumed that the root is found
  }
  \item{useFortran }{see \code{multiroot}; logical, if \code{FALSE}, then an \R implementation of Newton-Raphson is used
  }
}
\details{
  Uses the function \code{multiroot} in the \code{rootSolve} package to solve the conditional GEE.  As in \code{multiroot}, there is no guarantee of finding the root.

  A monotone increasing transformation can be applied to the observed gap times before calling \code{condGEE}.

  When \code{robust=TRUE}, \eqn{\theta} and \eqn{\sigma^2} are solved for in an alternating fashion until convergence. Note that the estimating equation for the mean parameters depends on \eqn{\sigma^2} through the censored gap time.
}
\value{
  a list containing:
  \item{eta }{the parameter estimate \eqn{(\theta^T,\sigma^2)^T}
  }
  \item{a.var }{an estimate of the asymptotic variance matrix of the eta estimator
  }
}
\author{
  David Clement <dyc24@cornell.edu>
}
\references{
  Clement, D. Y. and Strawderman, R. L. 2009 \emph{Biostatistics} \bold{10}, 451--467.
}

