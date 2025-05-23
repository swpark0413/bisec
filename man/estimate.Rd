% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{estimate}
\alias{estimate}
\alias{estimate.bpcr}
\alias{estimate.beig}
\title{Point and interval estimation of posterior distributions}
\usage{
estimate(object, ...)

\method{estimate}{bpcr}(object, prob = 0.95, orthogonal = FALSE, ...)

\method{estimate}{beig}(object, prob = 0.95, orthogonal = FALSE, ...)
}
\arguments{
\item{object}{An object from \strong{spikedEIG} and \strong{spikedBPCR}.}

\item{...}{Additional arguments passed to or from other methods.}

\item{prob}{A numeric scalar in the interval (0,1) giving the target probability content of the intervals. The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob.}

\item{orthogonal}{A logical value indicating whether to ensure the orthogonality of the posterior mean of the eigenvectors.}
}
\value{
An array (or matrix) containing the lower and upper bounds of the credible interval, as well as the posterior mean.
}
\description{
Compute the point estimate (mean) and equal-tailed credible interval to describe posterior distribution.
}
\examples{

\dontrun{

# generate a spiked covariance matrix:
n <- 50
p <- 100
K <- 4
leading <- c(1000,500,200,100)
remaining <- rep(0.1, p - K)
Sigma0 <- diag(c(leading, remaining), p)

# generate data
set.seed(413)
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)

# estimate eigenvalues and eigenvectors:
res <- bisec::spikedEIG(X = X, K = 4,prior = list(nu = (p+2), A = diag(2, p)), nsample = 100)
est <- bisec::estimate(res)
}

}
\seealso{
\code{spikedEIG} and \code{spikedBPCR}
}
\author{
Sewon Park
}
