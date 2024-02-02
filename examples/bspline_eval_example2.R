options(scipen = 5)
par(mgp=c(2.4,1,.4))
par(mar=c(4,5,2.5,2)+0.1) 
par(cex=.9)

# clamped B-spline evaluation (open) --------------------
knot_number <- 10

# the graphical function
test_spline <- function(xseq, user_knots, p, 
                        sup = bspline_support(p, knots=user_knots)  ){
  
  bbasis <- bspline_eval(xseq, degree=p, U=sup, native=TRUE)
  
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
  
  #axis(1, at=sup, labels=sup)
  sup.lab <- round(sup, digits=2)
  axis(1, at=sup, labels=FALSE)#sup.lab)
  text(x=sup, y=par()$usr[3]-0.09*(par()$usr[4]-par()$usr[3]),
       labels=sup.lab, srt=45, adj=1, xpd=TRUE)
  
}

# choose the polinomial degree that you wish 
# splines will be evaluated from 0 to p
p <- 11

# the sequence of evaluation
xseq <- seq(from = 0, to = 1, length.out = 200);

# the knots
user_knots <- c(0,sort(runif(knot_number)),1)

# the chart generation
par(mfrow=c(2,3))
for(j in 0:p) test_spline(xseq, user_knots, j)
rm(j)



# periodic B-spline evaluation (closed) ---------------------------------
# the degree 
p <- 11

# the sequence of evaluation
xseq <- seq(from = 0, to = 0.99999, length.out = 200);

# the knots
user_knots <- sort(c(0,sort(runif(knot_number)),1))

par(mfrow=c(2,3))
for(j in 0:p){ 
  sup <- bspline_support(j, knots=user_knots, native=TRUE)
  sup_per <- transform_support_periodic(degree = j, U = sup, native=TRUE)
  test_spline(xseq, NULL, j, sup=sup_per)
}
rm(j)
