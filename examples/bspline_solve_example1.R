
# bspline_solve example

require(modgeom)

# the normal curve positions for drawing
position <- seq(from = 0, to = 1, length.out = 800)

# the control points
c_point <- t( rbind(c(0,1,2,3,3), c(0,3,3,3,3), c(2,1,1,3,4)) )

# the degree
deg <- 4

# the standard deviation of the randomization
normal.dev <- .3

# number of points of the dataset
rand_position <- runif(1050)

# COMPUTATION ----------------------------

# the correct model
R <- interpolate_bspline(position, c_point, degree=deg)

# the B-spline interpolation for the random
R_fuzzy <- interpolate_bspline(rand_position, c_point, degree=deg)
# randomization in the 3 coordinates
R_fuzzy <- rnorm(length(R_fuzzy), mean = 0, sd = normal.dev) + R_fuzzy

projection.angles <- projection_angles(rand_position, R_fuzzy, c_point, degree = deg )
print("projection angles ***")
print(projection.angles)

# TEST FITTING solution ---------------------
# N stores the new control points ...
N <- bspline_solve(R_fuzzy, rand_position, c_point, degree=deg, reg.factor = 1)

# the interpolation of the new control points
RN <- interpolate_bspline(position, N, degree=deg)


#DRAWING ----------------
require(rgl)



plot3d(R, type="l", 
       aspect=FALSE, xlab ="x", ylab="y", zlab ="z", 
       col=rgb(0,0,.3), lwd=2, main = "B-spline class 2 (clamped)" )

par3d(windowRect=c(184, 130, 813, 666))
#par3d(zoom = 1)

#draws the control points 
points3d(c_point, pch=2, col=rgb(1,0.5,0.5), size = 3)
lines3d(c_point, col=rgb(1,0,0), lty=2, lwd=.5)

#draws the data points
points3d(R_fuzzy, pch=2, col=rgb(0,1,0), size = 3)

#draws the curve and new control points
points3d(N, pch=2, col=rgb(1,0.5,1), size = 3)
lines3d(N, col=rgb(1,0,1), lty=2, lwd=.5)
lines3d(RN, col=rgb(0,0,1), lty=2, lwd=2)