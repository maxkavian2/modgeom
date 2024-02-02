
library(modgeom)
# the position which is interpolated
position <- seq(from = 0, to = 1-.005, length.out=200)

# the degrees
deg <- 4
c_point_number <- 7

par(mfrow=c(2,2))
for(i in 1:15){
  
  # the matrix of control points here
  a <- rnorm(c_point_number)*20-10;
  b <- rnorm(c_point_number)*20-10;
  c_point <- t (rbind( a, b ) )
  
  rm(a,b)
  
  R1 <- interpolate_bspline(position, c_point, degree = deg, native=TRUE)
  R1p <- interpolate_bspline(position, c_point, degree = deg, is.periodic=TRUE, native=TRUE)
  

  plot(c_point, pch="+", col="red", xlab = "x", ylab = "y",
       main="" )
  lines(c_point, col=rgb(1,0,0,.5), lty="dashed" )
  lines(R1, col=rgb(0,0,0), lwd =2 )

  plot(c_point, pch="+", col="red", xlab = "x", ylab = "y",
       main="" )
  lines(c_point, col=rgb(1,0,0,.5), lty="dashed" )
  lines(R1p, col=rgb(0,0,0), lwd =2 )

}

