
#
# Example of sequential spline fitting.
#
# This example does not show the performance of the function bspline_fit. Instead
# it is intended to graphically depict the sequential performance of the
# fitting algorithm. Check bspline_fit_example2.R in order to test a real 
# bspline_fit execution
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
RN <- interpolate_bspline(position, nc_point, degree=deg)
RN0 <- RN
ft_points <- bspline_footpoint(R_fuzzy, nc_point, degree = deg, 
                               tolerance=1e-1, verbose=F)
RN_ft <- interpolate_bspline(ft_points, nc_point, degree = deg)

#DRAWING ----------------
require(rgl)

#draws the data points
plot3d(R_fuzzy, size=3,
       aspect=FALSE, xlab ="", ylab="", zlab ="", 
       col="gray50", lwd=2, main = "",
       xlim=c(-1,5),ylim=c(-1,5),zlim=c(-1,5), box=FALSE, axes=FALSE)


#rgl.viewpoint(25,angle, zoom =zom)
view3d(25,angle, zoom =zom)

par3d(windowRect=c(184, 130, 813, 666))
#par3d(zoom = 1)

# draws stage 0 ab initio
points3d(nc_point, pch=2, col=rgb(0.7,0.7,0.7), size = 3)
lines3d(nc_point, col=rgb(0.5,0.5,1), lty=2, lwd=.5)

# solves several times the problem
lim<-50 # number of fixed iterations provided
RNg <- c()
for(j in 1:lim){
  #j<-1
  nc_point <- bspline_solve(R_fuzzy, ft_points, nc_point, 
                            degree=deg, reg.factor=1)
  # the interpolation of the new control points
  RN <- interpolate_bspline(position, nc_point, degree=deg, native=TRUE)
  ft_points <- bspline_footpoint(R_fuzzy, nc_point, u=ft_points,degree = deg, 
                                 tolerance=1e-6,verbose =F, step=1,native=TRUE)
  
  w<-1
  clr1 <- rgb(0.5,0.5,1);
  clr2 <- rgb(0.7,0.7, 0.7);
  if(j==lim){
    clr1 <- rgb(0,0,1);  
    clr2 <- rgb(0.7,0.7,0.7);
    w=3
  }
  
  # draws the diagonal spline (ab initio)
  points3d(nc_point, pch=2, col=clr2, size = 3)
  lines3d(nc_point, col=clr2, lty=2, lwd=w)
  # draws the current curve
  lines3d(RN, col=clr1, lwd=w)
  
  Sys.sleep(1)
  
  if(j== 1)
    RNg <- RN
}


# DISPLAYS the final results, i.e. original curve and final curve.
readline(prompt="Press [enter] to continue")
#draws the data points
plot3d(R_fuzzy, size=3,
       aspect=FALSE, xlab ="", ylab="", zlab ="", 
       col="gray50", lwd=2, main = "",
       xlim=c(-1,5),ylim=c(-1,5),zlim=c(-1,5), box=FALSE, axes=FALSE)

#rgl.viewpoint(25,angle, zoom =zom)
view3d(25,angle, zoom =zom)
lines3d(RN, col=rgb(0,0,1), lwd=w)
#lines3d(RNg, col=rgb(0,0,1), lwd=2)
lines3d(RN0, col=rgb(0,0,1), lwd=1)

# draws stage 0 ab initio
points3d(nc_point, pch=2, col=rgb(0.7,0.7,0.7), size = 3)
lines3d(nc_point, col=rgb(0.5,0.5,1), lty=2, lwd=.5)


adjs <- c(1,1)
text3d(x=RN[nrow(RN) %/% 2,][1],
       y=RN[nrow(RN) %/% 2,][2],
       z=RN[nrow(RN) %/% 2,][3]
       , expression(infinity), 
       adj = adjs, cex=2.5, 
       col="blue") 

#text3d(x=RNg[nrow(RN) %/% 2,][1],
#       y=RNg[nrow(RN) %/% 2,][2],
#       z=RNg[nrow(RN) %/% 2,][3]
#       , "1", adj = adjs, cex=2, 
#       col="blue") 

text3d(x=RN0[nrow(RN) %/% 2,][1],
       y=RN0[nrow(RN) %/% 2,][2],
       z=RN0[nrow(RN) %/% 2,][3]
       , "0", adj = adjs, cex=2, 
       col="blue") 

readline(prompt="Press [enter] to continue")
# it continue with the second example

#ANIMATES --------------------

#M <- par3d("userMatrix")
#if (!rgl.useNULL())
#  movie3d(  movie="/home/maxkavian/MAX_DOCUMENTS_debian/R_scripts/geometrical_modelling/PD_animate",
#            par3dinterp(time=(0:2)*0.75,userMatrix=list(M,
#                                                      rotate3d(M, pi/2, 1, 0, 0),
#                                                      rotate3d(M, pi/2, 0, 1, 0) ) ), 
#          duration=1.5,fps=24)



