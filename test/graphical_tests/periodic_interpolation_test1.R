
"
 periodic test
"

library(modgeom)

# the position which is interpolated
position <- seq(from = 0, to = 1, length.out=200)

# the degree
deg <- 4

# the matrix of control points here
#c_point <- rbind( c(0, 5, 0, 12, 10, 11, 10, 1, 1, 7),
#                  c(4, 12, -7, 3, -5, 10, 4, 8, 9, 9),
#                  c(1,  1,  1, 1,  -5,  1, 1, 1, 1, 10 ) )

c_point <- t( rbind( c(0, 5, 0, 0, 1), 
                     c(4, 12,-7,4, 1), 
                     c(0, 5, 0, 5, 0)  ) )

# adds the first control point to the last part
c_point1 <- rbind(c_point, c_point[1:1,])
#c_point2 <- rbind(c_point,c_point[1:(deg-1),])
#c_point2 <- rbind(c_point,c_point[1:deg,])


R1 <- interpolate_bspline(position, c_point, degree = deg, is.periodic = TRUE, native = TRUE)
R2 <- interpolate_bspline(position, c_point, degree = deg-1, is.periodic = TRUE, native = TRUE)
R3 <- interpolate_bspline(position, c_point, degree = deg-2, is.periodic = TRUE, native = TRUE)


#drawing -----
require(rgl)

par3d(windowRect=c(184, 130, 813, 666))
par3d(zoom = 1)

plot3d(R1, type="l", 
       aspect=T, xlab ="x", ylab="y", zlab ="z", 
       col=rgb(0,0,.3), lwd=2, main = "" )

#lines3d(t(R2), col = rgb(0,0,.5), lwd = 2 )
lines3d(R2, col = rgb(0,0,.5), lwd = 2 )
lines3d(R3, col = rgb(0,0,.9), lwd = 2 )

#draws the control points 
points3d( c_point1, pch=2, col=rgb(1,0.5,0.5,.5), size = 3)
lines3d(c_point1, col=rgb(1,0,0,1), lty=2, lwd=.5)

#animates ------
M <- par3d("userMatrix")
if (!rgl.useNULL())
  play3d( par3dinterp(time=(0:2)*1.5,userMatrix=list(M,
                                                     rotate3d(M, pi/2, 1, 0, 0),
                                                     rotate3d(M, pi/2, 0, 1, 0) ) ),          
          duration=6 )