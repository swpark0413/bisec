#' Point and interval estimation of posterior distributions
#'
#' Compute the point estimate (mean) and equal-tailed credible interval to describe posterior distribution.
#'
#' @param object An object from \strong{spikedEIG} and \strong{spikedBPCR}.
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the intervals. The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob.
#' @param orthogonal A logical value indicating whether to ensure the orthogonality of the posterior mean of the eigenvectors.
#'
#'
#' @return An array (or matrix) containing the lower and upper bounds of the credible interval, as well as the posterior mean.
#'
#' @author Sewon Park
#' @seealso \code{spikedEIG} and \code{spikedBPCR}
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
#' est <- bisec::estimate(res)
#' }
#'
estimate <- function(object, ...) {
  UseMethod("estimate")
}



#' @rdname estimate
#' @export
#'
estimate.bpcr <- function(object, prob = 0.95, orthogonal = FALSE, ...) {

  left_tail <- (1 - prob) / 2
  right_tail <- 1 - left_tail

  # posterior samples
  post_beta <- object$post$beta
  post_eta <- object$post$eta

  # posterior estimates for beta
  post_mean_beta <- colMeans(post_beta)
  credintv_beta = apply(post_beta, 2, function(x){quantile(x,c(left_tail, right_tail))})
  post_est_beta = rbind(credintv_beta, post_mean_beta)
  rownames(post_est_beta) <- c('lower','upper','mean')

  # posterior estimates for eta
  post_mean_eta <- colMeans(post_eta)
  credintv_eta = apply(post_eta, 2, function(x){quantile(x,c(left_tail, right_tail))})
  post_est_eta = rbind(credintv_eta, post_mean_eta)
  rownames(post_est_eta) <- c('lower','upper','mean')

  # posterior estimates for eigenvalues and eigenvectors
  post_est_eigen = get_est_eigen(object, prob, orthogonal)

  return(list(beta = post_est_beta, eta = post_est_eta, evals = post_est_eigen$evals, evecs = post_est_eigen$evecs))
}



#' @rdname estimate
#' @export
#'
estimate.beig <- function(object, prob = 0.95, orthogonal = FALSE, ...) {
  return(get_est_eigen(object, prob, orthogonal))
}



get_est_eigen <- function(object, prob = 0.95, orthogonal = FALSE){

  left_tail <- (1 - prob) / 2
  right_tail <- 1 - left_tail

  # posterior samples
  post_evecs <- object$post$evecs
  post_evals <- object$post$evals %>% do.call('rbind',.)

  # posterior estimates for eigenvalues
  post_mean_evals <- colMeans(post_evals)
  credintv_evals <- apply(post_evals, 2, function(x){quantile(x,c(left_tail, right_tail))})
  post_est_evals <- rbind(credintv_evals, post_mean_evals)
  rownames(post_est_evals) <- c('lower','upper','mean')

  # posterior estimates for eigenvectors
  temp <- lapply(splitList3(post_evecs), FUN = function(x){apply(x, 2, function(y){sign(sum(x[,1] * y)) * y})})
  post_mean_evecs = lapply(temp, rowMeans) %>% do.call(cbind,.) %>% apply(., 2, function(x){x/norm(x, type = '2')})
  if(orthogonal){
    post_mean_evecs <- post_mean_evecs %>% pracma::gramSchmidt() %>% .$Q
  }
  credintv_evecs <- apply(post_evecs %>% abind::abind(.,along = 0), c(2,3), function(x){quantile(x,c(left_tail, right_tail))})
  dim(post_mean_evecs) <- c(1, dim(post_mean_evecs))
  post_est_evecs <- abind::abind(credintv_evecs, post_mean_evecs, along=1)
  dimnames(post_est_evecs)[[1]] <- c('lower','upper','mean')

  return(list(evals = post_est_evals, evecs = post_est_evecs))
}
