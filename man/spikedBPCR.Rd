% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spikedBPCR.R
\name{spikedBPCR}
\alias{spikedBPCR}
\title{Bayesian principal component regression for a spiked covariance matrix}
\usage{
spikedBPCR(Y, X, K = NULL, prior = list(), nsample = 100)
}
\arguments{
\item{Y}{A response vector.}

\item{X}{A n \eqn{\times} p data matrix, where \code{n} is the number of observations and \code{p} is the number of predictors.}

\item{K}{The number of principal components.}

\item{prior}{A list containing prior information for the model. The list includes the following parameters (with default values in parentheses):
\describe{
  \item{\code{A}}{The positive definite scale matrix for the inverse-Wishart prior (default: \code{I}).}
  \item{\code{nu}}{The degrees of freedom for the inverse-Wishart prior (default: \code{p + K}).}
}}

\item{nsample}{A scalar indicating the number of posterior samples to generate.}
}
\value{
A list with the following elements:
\describe{
  \item{\code{post}}{A list of posterior samples for the parameters of interest.}
  \item{\code{p}}{The dimension of the covariance matrix.}
  \item{\code{K}}{The number of principal components retained.}
  \item{\code{prior}}{A list containing the prior information used in the model.}
  \item{\code{X}}{The original data matrix.}
  \item{\code{centering}}{A vector of means used to center both the response and the columns of the data matrix.}
}
}
\description{
Performs a Bayesian principal component regression on a given data matrix, assuming a spiked covariance structure.
}
\details{
Lee et al. (2024+) proposed a post-processed posterior method for Bayesian principal component regression (PCR).
Given that the data \eqn{\boldsymbol{Z}_1, \ldots, \boldsymbol{Z}_n} are generated from a multivariate normal distribution \eqn{N_{p+1}(\boldsymbol{0}_{p+1}, \boldsymbol{\Sigma})}, where the covariance matrix \eqn{\boldsymbol{\Sigma}} includes both the response variable \eqn{Y} and the predictors \eqn{\boldsymbol{X}}.
An inverse-Wishart prior \eqn{IW_{p+1}(\bm{A}_n, \nu_n)} is imposed on \eqn{\boldsymbol{\Sigma}}, leading to a posterior distribution from which initial samples are generated.

To estimate the top \eqn{K} eigenvectors \eqn{\bm{\Gamma}}, the coefficients of the principal scores \eqn{\bm{\eta}}, and the coefficients of the original variables \eqn{\bm{\beta}}, three functions are defined:
\itemize{
\item \eqn{f_{\bm{\Gamma}}(\bm{\Sigma})}: Extracts the top \eqn{K} eigenvectors from the submatrix \eqn{\bm{\Sigma}_{-1,-1}}, which represents the covariance of the predictors \eqn{\boldsymbol{X}} excluding the response \eqn{Y}.
\item \eqn{f_{\bm{\eta}}(\bm{\Sigma})}: Computes a vector based on the inverse of the top \eqn{K} eigenvalues and the first column of \eqn{\bm{\Sigma}} (excluding its first element), reflecting the relationship between the response and predictors.
\item \eqn{f_{\bm{\beta}}(\bm{\Sigma})}: Multiplies \eqn{f_{\bm{\Gamma}}(\bm{\Sigma})} and \eqn{f_{\bm{\eta}}(\bm{\Sigma})} to estimate the regression coefficients \eqn{\bm{\beta}}.
}
These functions transform the initial posterior samples of \eqn{\bm{\Sigma}} into samples for the posterior distributions of \eqn{\bm{\Gamma}}, \eqn{\bm{\eta}}, and \eqn{\bm{\beta}}, thereby providing estimates of these key components under the spiked covariance model.
For more details, see Lee et al. (2024+).
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
ed = eigen(Sigma0)
G = ed$vectors[,1:K]
eta = rnorm(K, 0, 2)
beta = G \%*\% matrix(eta,ncol = 1)
s2r = 100
BTSB = (t(beta) \%*\% Sigma0) \%*\% beta
sigma = sqrt(BTSB)/s2r
eps = rnorm(n, 0, sigma)
Y = X \%*\% beta + eps

# fit a Bayesian PCR model:
res <- spikedBPCR(Y, X, K, prior = list(nu = (p+2), A = diag(2, (p+1))))
}

}
\references{
Lee, K., Park, S., Kim, S. and Lee, J. (2024+), "Posterior asymptotics of high-dimensional spiked covariance model with inverse-Wishart prior",
	\emph{Arxiv}, URL: \url{https://arxiv.org/abs/2412.10753}.
}
\seealso{
\code{estimate} and \code{predict.bpcr}
}
\author{
Sewon Park
}
\keyword{Bayesian}
\keyword{PCR}
\keyword{and}
\keyword{covariance}
\keyword{spiked}
