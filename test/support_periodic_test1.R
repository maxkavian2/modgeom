

# in this test we checked the periodic and non-periodic B-spline
# generation. 

library(modgeom)

# clamped B-spline evaluation (open) --------------------
knot_number <- 100

# periodic B-spline evaluation (closed) ---------------------------------
# the degree 
p <- 9

# the sequence of evaluation, 
#xseq <- seq(from = 0, to = .999999999999, length.out = 100);

# the knots
#user_knots <- c(0,sort(runif(knot_number)),1)
user_knots <- seq(from=0, to=1, length.out=knot_number+2)
sup <- bspline_support(degree = p, knots=user_knots, native=T)

s1 <- Sys.time()
sup_per <- transform_support_periodic(degree = p, U = sup, native=T)
s2 <- Sys.time()

ss1 <- Sys.time()
sup_per <- transform_support_periodic(degree = p, U = sup, native=F)
ss2 <- Sys.time()

print(s2 - s1)
print(ss2 - ss1)