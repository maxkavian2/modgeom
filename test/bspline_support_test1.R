
#
# tests the performance of bspline_support function
#

degree <- 50
knots <- seq(from=0, to=1, length.out=20000)
#knots[1001]  <- knots[1000]

s1 <- Sys.time()
bspline_support(degree, knots=knots)
s2 <- Sys.time()

ss1 <- Sys.time()
bspline_support(degree, knots=knots, native=TRUE)
ss2 <- Sys.time()

print(s2 - s1)
print(ss2 - ss1)
