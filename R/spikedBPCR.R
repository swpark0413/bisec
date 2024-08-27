#' Bayesian principal component regression for a spiked covariance matrix
#'
#' Performs a Bayesian principal component regression on a given data matrix, assuming a spiked covariance structure.
#'
#' Lee, Park, and Lee (2024+) proposed a post-processed posterior method for Bayesian principal component regression (PCR).
#' Given that the data \eqn{\boldsymbol{Z}_1, \ldots, \boldsymbol{Z}_n} are generated from a multivariate normal distribution \eqn{N_{p+1}(\boldsymbol{0}_{p+1}, \boldsymbol{\Sigma})}, where the covariance matrix \eqn{\boldsymbol{\Sigma}} includes both the response variable \eqn{Y} and the predictors \eqn{\boldsymbol{X}}.
#' An inverse-Wishart prior \eqn{IW_{p+1}(\bm{A}_n, \nu_n)} is imposed on \eqn{\boldsymbol{\Sigma}}, leading to a posterior distribution from which initial samples are generated.
#'
#' To estimate the top \eqn{K} eigenvectors \eqn{\bm{\Gamma}}, the coefficients of the principal scores \eqn{\bm{\eta}}, and the coefficients of the original variables \eqn{\bm{\beta}}, three functions are defined:
#' \itemize{
#' \item \eqn{f_{\bm{\Gamma}}(\bm{\Sigma})}: Extracts the top \eqn{K} eigenvectors from the submatrix \eqn{\bm{\Sigma}_{-1,-1}}, which represents the covariance of the predictors \eqn{\boldsymbol{X}} excluding the response \eqn{Y}.
#' \item \eqn{f_{\bm{\eta}}(\bm{\Sigma})}: Computes a vector based on the inverse of the top \eqn{K} eigenvalues and the first column of \eqn{\bm{\Sigma}} (excluding its first element), reflecting the relationship between the response and predictors.
#' \item \eqn{f_{\bm{\beta}}(\bm{\Sigma})}: Multiplies \eqn{f_{\bm{\Gamma}}(\bm{\Sigma})} and \eqn{f_{\bm{\eta}}(\bm{\Sigma})} to estimate the regression coefficients \eqn{\bm{\beta}}.
#' }
#' These functions transform the initial posterior samples of \eqn{\bm{\Sigma}} into samples for the posterior distributions of \eqn{\bm{\Gamma}}, \eqn{\bm{\eta}}, and \eqn{\bm{\beta}}, thereby providing estimates of these key components under the spiked covariance model.
#' For more details, see Lee, Park, and Lee (2024+).
#'
#' @param Y A response vector.
#' @param X A n \eqn{\times} p data matrix, where \code{n} is the number of observations and \code{p} is the number of predictors.
#' @param K The number of principal components.
#' @param prior A list containing prior information for the model. The list includes the following parameters (with default values in parentheses):
#' \describe{
#'   \item{\code{A}}{The positive definite scale matrix for the inverse-Wishart prior (default: \code{I}).}
#'   \item{\code{nu}}{The degrees of freedom for the inverse-Wishart prior (default: \code{p + K}).}
#' }
#' @param nsample A scalar indicating the number of posterior samples to generate.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{post}}{A list of posterior samples for the parameters of interest.}
#'   \item{\code{p}}{The dimension of the covariance matrix.}
#'   \item{\code{K}}{The number of principal components retained.}
#'   \item{\code{prior}}{A list containing the prior information used in the model.}
#'   \item{\code{X}}{The original data matrix.}
#'   \item{\code{centering}}{A vector of means used to center both the response and the columns of the data matrix.}
#' }
#'
#' @author Sewon Park
#' @seealso \code{estimate} and \code{predict.bpcr}
#' @keywords spiked covariance and Bayesian PCR
#'
#' @references Lee, K., Park, S., and Lee, J. (2024+), "Bayesian inference on spiked eigenstructure of high-dimensional covariances",
#' 	\emph{Arxiv}, URL: \url{https://arxiv.org/abs/2101.12179}.
#'
#' @importFrom hdbinseg get.factor.model
#' @importFrom CholWishart rInvWishart
#' @importFrom progress progress_bar
#' @importFrom purrr map
#' @importFrom magrittr `%>%`
#' @export
#' @examples
#'
#' \dontrun{
#'
#' # generate a spiked covariance matrix:
#' n <- 50
#' p <- 100
#' K <- 4
#' leading <- c(1000,500,200,100)
#' remaining <- rep(0.1, p - K)
#' Sigma0 <- diag(c(leading, remaining), p)
#'
#' # generate data
#' set.seed(413)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' ed = eigen(Sigma0)
#' G = ed$vectors[,1:K]
#' eta = rnorm(K, 0, 2)
#' beta = G %*% matrix(eta,ncol = 1)
#' s2r = 100
#' BTSB = (t(beta) %*% Sigma0) %*% beta
#' sigma = sqrt(BTSB)/s2r
#' eps = rnorm(n, 0, sigma)
#' Y = X %*% beta + eps
#'
#' # fit a Bayesian PCR model:
#' res <- spikedBPCR(Y, X, K, prior = list(nu = (p+2), A = diag(2, (p+1))))
#' }
#'


spikedBPCR <- function(Y, X, K = NULL, prior = list(), nsample = 100){

  eig_split = function(S, K, n, pb = NULL){

    pp1 = ncol(S)
    Sxx = S[2:pp1, 2:pp1]
    Sxy = S[2:pp1, 1, drop=F]
    Syy = S[1,1]

    p = pp1 - 1
    eigen_decomp = eig_decomp(Sxx, K = K)
    eig_vals = eigen_decomp$values[1:K]
    eig_vecs = eigen_decomp$vectors[,1:K, drop = F]



    D_lead = diag(1/(eig_vals), K, K)
    temp = arma_matmul(D_lead, t(eig_vecs))
    eta_hat = arma_matmul(temp, Sxy)
    beta_hat = arma_matmul(eig_vecs, eta_hat)

    if(!is.null(pb)){
      pb$tick()
    }

    return(list(beta = beta_hat, eta = eta_hat, evals = eig_vals, evecs = eig_vecs, Sxy = Sxy, Syy = Syy))
  }

  Z = cbind(Y, X)
  Z = scale(Z, scale = F)
  p <- ncol(Z)
  n <- nrow(Z)


  if(is.null(K)){
    K = hdbinseg::get.factor.model(X, ic="ah")$r.hat
  }

  ZZT = innerprod(Z, Z)
  Scov = ZZT / n
  privals = list(nu = p + K, A = diag(1,p))
  privals[names(prior)] <- prior
  A = privals$A
  nu = privals$nu

  # initial posterior step
  parameter_initposterior <- list(A = ZZT + A, nu = nu + n)


  # post-processing step

  pb <- progress::progress_bar$new(format = "MCMC iteration [:current / :total], Estimated time remaining: :eta",
                                   total = nsample,
                                   clear = FALSE,
                                   width= 80)

  ppp_save <- purrr::map(1:nsample,
                         ~(CholWishart::rInvWishart(1, parameter_initposterior$nu, parameter_initposterior$A)[,,1]) %>%
                           eig_split(S = ., K = K, n = n, pb = pb))

  ppp_save = splitList2(ppp_save)
  beta_samp = ppp_save$beta %>% do.call("cbind",.) %>% t()
  eta_samp = ppp_save$eta %>% do.call("cbind",.) %>% t()
  dimnames(beta_samp) = list(iterations = NULL, parameters = paste0('beta',1:ncol(beta_samp)))
  dimnames(eta_samp) = list(iterations = NULL, parameters = paste0('eta',1:ncol(eta_samp)))
  sxy_samp = ppp_save$sxy %>% do.call("cbind",.)
  syy_samp = ppp_save$syy %>% unlist()


  out <- list()
  out$p <- p
  out$K <- K
  out$post = list(beta = beta_samp, eta = eta_samp, syy = syy_samp, sxy = sxy_samp, evals =  ppp_save$evals, evecs = ppp_save$evecs)
  out$prior = privals
  out$X = X
  out$centering = attr(Z, "scaled:center")
  class(out) = 'bpcr'
  return(out)
}
