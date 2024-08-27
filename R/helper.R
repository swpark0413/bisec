get_upper = function(X){
  dim_info = dim(X)
  if(dim_info[1] == dim_info[2]){
    return(X[upper.tri(X, diag = T)])
  } else {
    return(X)
  }
}


eig_decomp = function(S, K = NULL){
  eigen_decomp = Eigen_eig(S)
  p = ncol(S)
  if(is.null(K)){
    K = p
  }
  idx = seq(p, p-K+1, by = -1)
  eig_vecs = eigen_decomp$vectors[,idx ,drop = F]
  eig_vals = eigen_decomp$values[idx]
  return(list(values = eig_vals, vectors = eig_vecs))
}

