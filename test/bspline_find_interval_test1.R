

#
# tests the performance of bspline_find_interval_index
#


degree <- 3
knots <- seq(from=0, to=1, length.out=10)
support <- bspline_support(degree, knots=knots, native=TRUE)

val <- .5

s1 <- Sys.time()
bspline(val, k=-1, U=support, native=FALSE)
s2 <- Sys.time()

ss1 <- Sys.time()
bspline(val, k=-1, U=support, native=TRUE)
ss2 <- Sys.time()

print(s2 - s1)
print(ss2 - s1)

