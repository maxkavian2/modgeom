

#
# tests the performance of bspline_find_interval_index
#


degree <- 8
knots <- seq(from=0, to=1, length.out=100)
support <- bspline_support(degree, knots=knots, native=TRUE)

vals <- seq(from=0, to=1, length.out = 10)

s1 <- Sys.time()
bspline_eval(vals, U=support, native=FALSE)
s2 <- Sys.time()

ss1 <- Sys.time()
bspline_eval(vals, U=support, native=TRUE)
ss2 <- Sys.time()

print(s2 - s1)
print(ss2 - ss1)

