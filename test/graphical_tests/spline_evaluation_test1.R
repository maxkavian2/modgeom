
"
  points evaluation of a b-spline.
  The result are the values of all the b-splines
"

library(modgeom)

# degree of the b-spline
degree <- 3

# used knots for this spline
user_knots <- c(0,.2,.4,.6,.8,1)

# the position which is interpolated
position <- c(0, .25, .75, 1)

# HERE STARTS THE FUN ...
# the current support
sup <- bspline_support(degree,knots=user_knots, native=TRUE)
print( sup ) 

# the values of all the b-splines at the current position
B <- bspline_eval(position , degree = degree, U = sup, native=TRUE)
print( B ) 
print( paste("sum of B =" ,sum(B)) )

print(bspline(position[length(position)], degree = degree, U = sup ))
