

#
# tests the performance of bspline_find_interval_index
#


degree <- 3
knots <- seq(from=0, to=1, length.out=10)

support <- bspline_support(degree, knots=knots, native=TRUE)

val <- .78
bspline_find_interval_index(val , support, native=FALSE)
bspline_find_interval_index(val , support, native=TRUE)
