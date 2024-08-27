#' Predict method for Bayesian PCR Model Fits
#'
#' Computes a predictive distribution for new data from a fitted Bayesian PCR model object.
#'
#'
#' @param object A fitted model object from \strong{spikedBPCR}.
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. If omitted, the model matrix is used. If newdata is provided and any variables were transformed (e.g. rescaled) in the data used to fit the model, then these variables must also be transformed in newdata.
#'
#' @return A \code{nsample} \eqn{\times} \code{nrow(newdata)} matrix of simulations from the posterior predictive distribution. Each row of the matrix is a vector of predictions generated using a single draw of the model parameters from the posterior distribution.
#' @author Sewon Park
#' @seealso \code{spikedBPCR}
#' @importFrom stats rnorm
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
#'
#' # estimate a predictive distribution for new data using the fitted Bayesian PCR model:
#' Xnew <- MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = Sigma0)
#' pred <- predict(res, Xnew)
#' }
#'
#'
predict.bpcr = function(object, newdata = NULL){
  y_scaled_mean = object$centering[1]
  y_scaled_std = 1
  scaled_mean = object$centering[-1]
  scaled_std = 1

  if(is.null(newdata)){
    newdata = t(object$X)
  } else {
    if(is.vector(newdata)){
      newdata = matrix(newdata, nrow = 1)
    }
    newdata <- apply(newdata, 1, function(x){(x - scaled_mean)/scaled_std}) %>% as.matrix()
  }

  beta = object$post$beta
  syy = object$post$syy
  sxy = object$post$sxy
  m = length(syy)
  n = dim(newdata)[2]

  results = matrix(NA, m, n)
  for(i in 1:m){
    b = beta[i, ,drop = F]
    mu = arma_matmul(b, newdata) %>% as.numeric()
    sigma2 = syy[i] - as.numeric(b %*% sxy[,i, drop = F])
    results[i,] = y_scaled_std * rnorm(n, mu, sqrt(sigma2)) + y_scaled_mean
  }
  return(results)

}
