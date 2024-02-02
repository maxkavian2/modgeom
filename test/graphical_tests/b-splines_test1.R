
#
# Represents several B-splines along the knot sequence.
#

library(modgeom)

# choose the polinomial degree that you wish
p <- 11

# the sequence of evaluation
xseq <- seq(from = 0, to = 1, length.out = 100);

# the knots
user_knots <- c(0,.2,.4,.45,.50,1)
user_knots <- c(user_knots/2, 
                seq(from=.6, to=1.0, length.out=5))

# the graphical function
test_spline <- function(xseq, user_knots, p){
  
  for(j in 0:p){
    sup <- bspline_support(j, knots=user_knots, native=TRUE)   
    bbasis <- bspline_eval(xseq, degree=j, U=sup, native=TRUE)
    
    for(i in 1:ncol(bbasis)){
      if(i == 1){
        plot(xseq, bbasis[,1], type = "l", main=paste("p = ",j, sep=""),
             ylab = "B(x)", xlab="x", xaxt="n")      
      }else{
        points(xseq, bbasis[,i], type = "l")
      }
    }
    
    for(i in 1:length(sup))
      lines(c(sup[i],sup[i]),c(0,1), lty="dotted")
    
    axis(1, at=sup, labels=sup)
  }
  
}

# the chart generation -----------------
#pdf("bsplines2_test_chart.pdf", 8, 6)
par(mfrow=c(2,3))
test_spline(xseq, user_knots, p)
#dev.off()