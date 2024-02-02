
#
# In this test we checked the periodic and non-periodic B-spline
# generation. 
#

library(modgeom)

# clamped B-spline evaluation (open) --------------------
knot_number <- 8

# the graphical function
test_spline <- function(xseq, user_knots, p, 
                        sup = bspline_support(p, knots=user_knots)  ){
  
  bbasis <- bspline_eval(xseq, degree=p, U=sup)
  
  for(i in 1:ncol(bbasis)){
    if(i == 1){
      plot(xseq, bbasis[,1], type = "l", main=paste("p = ",p, sep=""),
           ylab = "B(x)", xlab="x", xaxt="n", ylim=c(0,1))      
    }else{
      points(xseq, bbasis[,i], type = "l")
    }
  }
  
  for(i in 1:length(sup))
    lines(c(sup[i],sup[i]),c(0,1), lty="dotted")
  
  axis(1, at=sup, labels=sup)
  
}

# choose the polinomial degree that you wish 
# splines will be evaluated from 0 to p
p <- 5

# the sequence of evaluation
xseq <- seq(from = 0, to = 0.99999999, length.out = 100);

# the knots
user_knots <- c(0,sort(runif(knot_number)),1)
#user_knots <- c(user_knots/2, 
#                seq(from=.6, to=1.0, length.out=5))



# the chart generation
#pdf("bsplines2_test_chart_test1b.pdf", 9.5, 6)
par(mfrow=c(2,3))
for(j in 0:p) test_spline(xseq, user_knots, j)
rm(j)
#dev.off()



# periodic B-spline evaluation (closed) ---------------------------------
# the degree 
p <- 11

# the sequence of evaluation
xseq <- seq(from = 0, to = .999999999999, length.out = 100);

# the knots
user_knots <- c(0,sort(runif(knot_number)),1)

#pdf("bsplines2_test_chart_periodic_test1b.pdf", 9.5, 6)
par(mfrow=c(2,3))
for(j in 0:p){ 
  sup <- bspline_support(j, knots=user_knots)
  sup_per <- transform_support_periodic(degree = j, U = sup, native=F)
  test_spline(xseq, NULL, j, sup=sup_per)
}
rm(j)
#dev.off()