
#
# sequential spline fitting using
# bspline_fit
#

# USER PARAMETERS -------------------------
require(modgeom)

# the normal curve positions for drawing
position <- seq(from = 0, to = 1, length.out = 800)

# the control points
c_point <- t( rbind(c(0,1,2,3,3), c(0,3,3,3,3), c(2,-1,7,3,4)) )

# the degree
deg <- 3

# the number of initial control points of the first guess
# this makes the result more complex
initial.control.points<- deg

# the standard deviation of the randomization
normal.dev <- .1

# number of points of the dataset
rand_position <- runif(300)
angle <- -20
zom <- 0.5


# COMPUTATION ----------------------------

# the correct model
R <- interpolate_bspline(position, c_point, degree=deg)

# the B-spline interpolation for the random
R_fuzzy <- interpolate_bspline(rand_position, c_point, degree=deg)
# randomization in the 3 coordinates
R_fuzzy <- rnorm(length(R_fuzzy), mean = 0, sd = normal.dev) + R_fuzzy

# TEST FITTING solution ---------------------
# finds the diagonal of the bounding box according to the current degree
nc_point <- bounding_diagonal_spline(R_fuzzy, initial.control.points, f=1.4)

w <- 3
nc_point <- bspline_fit(R_fuzzy, degree=deg, v=nc_point, 
                        tolerance=1e-7,
                        fp.tolerance=1e-3, step=1, #verbose=FALSE, 
                        algorithm="directional",
                        verb.progress=TRUE,
                        reg.factor=1, native=TRUE)
RN <- interpolate_bspline(position, nc_point, degree=deg,native=TRUE)



# DISPLAYS the final results, i.e. original curve and final curve.
require(rgl)



#draws the data points
plot3d(R_fuzzy, size=3,
       aspect=FALSE, xlab ="", ylab="", zlab ="", 
       col="gray50", lwd=2, main = "",
       xlim=c(-1,5),ylim=c(-1,5),zlim=c(-1,5), box=FALSE, axes=FALSE)
par3d(windowRect=c(184, 130, 813, 666))
#par3d(zoom = 1)
#rgl.viewpoint(25,angle, zoom =zom)
view3d(25,angle, zoom =zom)
lines3d(RN, col=rgb(0,0,1), lwd=w)
#lines3d(RNg, col=rgb(0,0,1), lwd=2)
#lines3d(RN0, col=rgb(0,0,1), lwd=1)



