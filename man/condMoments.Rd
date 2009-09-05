\name{condMoments}
\Rdversion{2.9.0}
\alias{condMoments}
\alias{K1.norm}
\alias{K2.norm}
\alias{K1.t3}
\alias{K2.t3}
\alias{K1.exp}
\alias{K2.exp}
\title{
  First and second conditional moments for 3 distributions
}
\description{
  K1 provides \eqn{E(Y|Y>w)} and K2 provides \eqn{E(Y^2|Y>w)} for Y as standard normal, standardized \emph{t} with 3 degrees of freedom, or an exponential with mean 0 and variance 1.
}
\usage{
K1.norm(w)
K2.norm(w)
K1.t3(w)
K2.t3(w)
K1.exp(w)
K2.exp(w)
}
\arguments{
  \item{w }{a real-valued vector
  }
}
\value{a vector of conditional moments 
}
\author{
  David Clement <dyc24@cornell.edu>
}