upper_sparse_cov <- function(x, y, id.var = "SampleID", lambda){
  
  nvar_x <- ncol(x)-length(id.var)
  nvar_y <- ncol(y)-length(id.var)
  blacklist <- 
    rbind(
      expand.grid(1:nvar_x,1:nvar_y),
      expand.grid((nvar_x + 1):(nvar_x + nvar_y),(nvar_x + 1):(nvar_x + nvar_y))
    ) %>% as.matrix
  
  xy_sparse <- sparsebnUtils::sparsebnData(
    select_at(inner_join(x,y, by = id.var),vars(-!!sym(id.var))), type ="continuous")
  
  xy_cov <- sparsebn::estimate.covariance(
    xy_sparse, 
    blacklist = blacklist, 
    lambdas = lambda
  )[[1]]
  
  return(xy_cov[1L:nvar_x,(nvar_x+1L):(nvar_x + nvar_y)])
}
