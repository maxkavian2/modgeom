'double_array2 <- function(x){
  y <- rep(0, times=length(x))
  .C("double_array", as.integer(y), as.integer(x), as.integer(length(x)))
}'


'test_dgels <- function(A,B){
  
  Av <- as.vector(A)
  Bv <- as.vector(B)
  m <- nrow(A)
  n <- ncol(A)
  nrhs <- ncol(B)
  
  result <- .C( "test_dgels", 
                as.integer(m),
                as.integer(n),
                as.integer(nrhs),
                as.double(Av),
                as.double(Bv)
              )
  #R <- matrix(result[[5]], ncol=nrhs, nrow=m)[1:n,]
  
  return(result[[5]]);
}'

'test_dsyrk <- function(M){
  
  A <- as.vector(M)
  n <- nrow(M)
  k <- ncol(M)
  R <- as.vector(matrix(rep(0,n*n), ncol=n, nrow=n))
  
  result <- .C( "test_dsyrk", 
                as.integer(n), 
                as.integer(k), 
                as.numeric(A), 
                as.numeric(R) )
  
  return( matrix(result[[4]], ncol=n,nrow=n) )
  
}'