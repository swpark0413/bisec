#' Bayesian estimation of eigenstructure for a spiked covariance matrix
#'
#' Provides posterior samples for both eigenvalues and eigenvectors of a spiked covariance matrix.
#'
#' Lee et al. (2024+) proposed a Bayesian method for estimating the spiked eigenstructure of a covariance matrix,
#' focusing on the top \eqn{K} eigenvalues and eigenvectors. The method transforms the posterior distribution of the covariance matrix
#' into the posterior distributions of its top \eqn{K} eigenvalues \eqn{[\lambda_k(\boldsymbol{\Sigma}) \mid \mathbb{X}_n]}
#' and eigenvectors \eqn{[\boldsymbol{\xi}_k(\boldsymbol{\Sigma}) \mid \mathbb{X}_n]}, where \eqn{\mathbb{X}_n} denotes the \eqn{n} observations.
#'
#' Given an inverse-Wishart prior \eqn{\boldsymbol{\Sigma} \sim IW_p(\boldsymbol{A}_n, \nu_n)} on the covariance matrix,
#' the posterior is \eqn{\boldsymbol{\Sigma} \mid \mathbb{X}_n \sim IW_p(\boldsymbol{A}_n + \sum_{i=1}^n \boldsymbol{X}_i \boldsymbol{X}_i^T, \nu_n + n)}.
#' The algorithm for generating posterior samples from \eqn{[\lambda_k(\bm\Sigma) \mid \mathbb{X}_n]} and \eqn{[\bm\xi_k(\bm\Sigma) \mid \mathbb{X}_n]} is given below:
#' \itemize{
#' \item  Independently generate \eqn{\bm\Sigma_1,\ldots ,\bm\Sigma_N \sim IW(\bm{A}_n + \sum_{i=1}^n \bm{X}_i\bm{X}_i^T , \nu_n + n)}.
#' \item Compute the \eqn{k}-th eigenvalue and eigenvectors of \eqn{\bm\Sigma_1, \ldots, \bm\Sigma_N} and obtain
#' \deqn{
#' \lambda_k(\bm\Sigma_1),\ldots,\lambda_k(\bm\Sigma_N),\\
#' \bm\xi_k(\bm\Sigma_1),\ldots ,\bm\xi_k (\bm\Sigma_N),
#' } for \eqn{k=1,\ldots ,K}.
#' }
#' For more details, see Lee, Park, and Lee (2024+).
#'
#'
#' @param X An \eqn{n \times p} data matrix, where \code{n} is the number of observations and \code{p} is the number of predictors.
#' @param K The number of leading eigenvalues to be extracted.
#' @param prior A list containing prior information for the model. The list includes the following parameters (with default values in parentheses):
#' \describe{
#'   \item{\code{A}}{The positive definite scale matrix for the inverse-Wishart prior (default: \code{I}).}
#'   \item{\code{nu}}{The degrees of freedom for the inverse-Wishart prior (default: \code{p + k}).}
#' }
#' @param nsample A scalar value specifying the number of posterior samples to generate.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{post}}{A matrix of posterior samples for the parameters of interest.}
#'   \item{\code{p}}{The dimension of the covariance matrix.}
#'   \item{\code{K}}{The number of leading eigenvalues extracted.}
#'   \item{\code{prior}}{A list containing the prior information used in the model.}
#' }
#'
#' @author Sewon Park
#' @seealso \code{estimate}
#' @keywords spiked covariance
#'
#' @references Lee, K., Park, S., Kim, S. and Lee, J. (2024+), "Posterior asymptotics of high-dimensional spiked covariance model with inverse-Wishart prior",
#' 	\emph{Arxiv}, URL: \url{https://arxiv.org/abs/2412.10753}.
#'
#' @importFrom hdbinseg get.factor.model
#' @importFrom CholWishart rInvWishart
#' @importFrom progress progress_bar
#' @importFrom purrr map
#' @importFrom magrittr `%>%`
#' @export
#'
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
#'
#' # estimate eigenvalues and eigenvectors:
#' res <- bisec::spikedEIG(X = X, K = 4,prior = list(nu = (p+2), A = diag(2, p)), nsample = 100)
#' }
#'
spikedEIG = function(X, K = NULL, prior = list(), nsample = 100){

  eig_split = function(S, K, n,  pb = NULL){

    p = ncol(S)
    eigen_decomp = eig_decomp(S, K)
    eig_vecs = eigen_decomp$vectors
    eig_vals = eigen_decomp$values

    if(!is.null(pb)){
      pb$tick()
    }

    return(list(cov = get_upper(S), evals = eig_vals, evecs = eig_vecs))
  }

  p <- ncol(X)
  n <- nrow(X)

  if(is.null(K)){
    K = hdbinseg::get.factor.model(X, ic="ah")$r.hat
  }

  X = scale(X, scale = F)
  XXT = innerprod(X, X)
  Scov = XXT / n
  privals = list(nu = p + K, A = diag(1,p))
  privals[names(prior)] <- prior
  A = privals$A
  nu = privals$nu

  # initial posterior step
  parameter_initposterior <- list(A = (XXT) + A, nu = nu + n)


  # post-processing step
  pb <- progress::progress_bar$new(format = "MCMC iteration [:current / :total], Estimated time remaining: :eta",
                                   total = nsample,
                                   clear = FALSE,
                                   width= 80
  )


  Sigma_save <- purrr::map(1:nsample,
                           ~(CholWishart::rInvWishart(1, parameter_initposterior$nu, parameter_initposterior$A)[,,1]) %>%
                             eig_split(S = .,  K = K,n = n, pb = pb))

  Sigma_save = splitList1(Sigma_save)
  out <- list()
  out$post <- list(S = Sigma_save$cov, evals = Sigma_save$evals, evecs = Sigma_save$evecs)
  out$p <- p
  out$K <- K
  out$prior = privals
  class(out) = 'beig'
  return(out)
}




