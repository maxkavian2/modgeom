

#
# test interpolate_bspline_function
#

u <- seq(from=0, to=1, length=10)
cp <- t(matrix(c(2,3,4,1,2,0,1,0,-1,-1,-1,-2,2,-3,0), ncol=5,nrow=3 ))
degree <- 4
knots = seq(from = 0, to = 1, length.out=nrow(cp) - degree + 1)
der.red <- 2

s1 <- Sys.time()
A1 <- interpolate_bspline(u, cp, knots=knots, degree=degree, derivate = der.red, native=TRUE)
s2 <- Sys.time()

ss1 <- Sys.time()
A2 <- interpolate_bspline(u, cp, knots=knots, degree=degree, derivate = der.red, native=FALSE)
ss2 <- Sys.time()

print(s2 - s1)
print(ss2 - ss1)

ratio <- as.numeric(ss2 - ss1) / as.numeric(s2 - s1)
print(paste("the ratio is:", ratio, "times faster"))




