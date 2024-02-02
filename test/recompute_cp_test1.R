

#
# tests the performance of bspline_find_interval_index
#


degree <- 3
knots <- seq(from=0, to=1, length.out=15)
support <- bspline_support(degree, knots=knots, native=TRUE)

n <- length(support)-degree-1
cps <- cbind( runif(n)*20, runif(n)*20, runif(n)*20 )
rm(n)

s1 <- Sys.time()
recompute_control_points(cps, support, degree=degree, native=FALSE)
s2 <- Sys.time()

# M <- matrix(rep(NA,times=nrow(cps)*ncol(cps)), ncol=ncol(cps), nrow=nrow(cps))
# for(i in 1:(nrow(cps)-1) )
#   M[i,] <- cps[i+1,] - cps[i,]
# M
# rm(M,i)

ss1 <- Sys.time()
recompute_control_points(cps, support, degree=degree, native=TRUE)
ss2 <- Sys.time()

print(s2 - s1)
print(ss2 - ss1)

