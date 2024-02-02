
"
  test the bidimensional interpolation of points
"

library(modgeom)

# the position which is interpolated
position <- seq(from = 0, to = 0.99, length.out=200)

# the degrees
deg <- 5
c_point_number <- 6




#pdf("bsplines2_test_chart_test3.pdf", 9, 6.5)

par(mfrow=c(2,2))
for(i in 1:15){
  
  # the matrix of control points here
  a <- rnorm(c_point_number)*20-10; #a <- c(a,a[1:deg]);
  b <- rnorm(c_point_number)*20-10; #b <- c(b,b[1:deg]);
  c_point <- t (rbind( a, b ) )
  #print(c_point)
  rm(a,b)
  
  
  R1 <- interpolate_bspline(position, c_point, degree = deg, native=TRUE)
  R1p <- interpolate_bspline(position, c_point, degree = deg, is.periodic=TRUE, native=TRUE)
  #R2 <- interpolate_bspline(position, c_point, deg+1)
  #R3 <- interpolate_bspline(position, c_point, deg+2)
  
  
  
  plot(c_point, pch="+", col="red", xlab = "x", ylab = "y", 
       main="" )
  lines(c_point, col=rgb(1,0,0,.5), lty="dashed" )
  lines(R1, col=rgb(0,0,0), lwd =2 )
  
  plot(c_point, pch="+", col="red", xlab = "x", ylab = "y", 
       main="" )
  lines(c_point, col=rgb(1,0,0,.5), lty="dashed" )
  lines(R1p, col=rgb(0,0,0), lwd =2 )
  
  #lines(t(R2), col=rgb(0,0,.3), lwd = 2)
  #lines(t(R3), col=rgb(0,0,.5), lwd = 2)
  
}

#dev.off()