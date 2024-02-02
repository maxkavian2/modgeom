

"
produces several curves of different kind, in order to check the differentiability 
depending on the degree of the spline.

first argument: the output folder (i.e. the folder where the images will be produced)
second argument: the degree

Note that for the current number of control points, the maximum degree allowed is 9
"


#args = commandArgs(trailingOnly=TRUE)

#setwd("/home/max/R_scripts/geometrical_modelling");
#source("bsplines2_functions.R");
library(modgeom)

# the position which is interpolated
position <- seq(from = 0, to = 1, length.out=2000)


# the matrix of control points here
c_point <- t(rbind( c(0, 5, 0, 12, 10, 11, 10, 1, 1, 7),
                    c(4, 12, -7, 3, -5, 10, 4, 8, 9, 10),
                    c(5,  10,  1, 1,  -5,  1, 1, 1, 1, 10 ) ))

#c_point <- t( rbind( c(0,5,0,0), c(4,12,-7,4), c(12,1,1,-12), c(2,5,0,0),c(12,5,0,0)  ) )

# the degree
deg <- 9 #as.numeric(args[2])

R <- list()
Rtext <- c()
for(dg in 0:(deg-1) ){
  R[[dg+1]] <- interpolate_bspline(position, c_point, degree = deg-dg, native=T)
  Rtext <- rbind(Rtext, interpolate_bspline(.16, c_point, degree = deg-dg, native=T))
}
Rtext[,1] = Rtext[,1]-0.6

#R2 <- interpolate_bspline(position, c_point, deg-1)
#R3 <- interpolate_bspline(position, c_point, degree = deg-2)
colscale <- c("red","gray")
colf <- colorRampPalette(colscale)
cols <- colf(length(R))

#drawing -----
require(rgl)

make.chart <- function(tag="", vpar=NULL){
  mar.red <- 0.1
  # graphical
  # changes the likelihood to switch to scientific notation
  options(scipen = 5)
  #set the distances of the axes labels to the plot
  par(mgp=c(2.4,1,.3))
  # sets the margins
  par(mar=c(4,5-mar.red,2.5-mar.red,2-mar.red)+0.1) 
  # sets the size of the character
  par(cex=1)
  
  # 3D pars
  par3d(windowRect=c(150, 150, 750, 750))
  rgl.viewpoint(vpar[1], vpar[2])
  #par3d(windowRect=c(184, 130, 813, 666))
  par3d(zoom = 0.7)
  
  plot3d(R[[1]], type="l", 
         aspect=T,  
         col=cols[1], lwd=1, main = "", xlab="x", ylab="y", zlab="z", box=F, 
         axes=F)
  
  axes3d(
    edges=c('x--', 'y--', 'z--'),
    labels=F,
    tick=T, at=c(-4:8)
    
  )
  
  
  
  #lines3d(t(R2), col = rgb(0,0,.5), lwd = 2 )
  if(length(R) > 1)
    for(i in 2:length(R) )
      lines3d(R[[i]], col = cols[i], lwd = 1 )
  
  #draws the control points 
  points3d( c_point, pch=2, col="gray", size = 6)
  #lines3d(c_point, col=rgb(0,1,0,0.5), lty=2, lwd=.5)
  #rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = c('pink'))
  
  # prints the texts
  labs <- c()
  for(i in deg:1)
    labs <- c(labs, paste("C",as.character(i-1), sep=""))
  
  rgl.texts(Rtext[,1], y =Rtext[,2], z=Rtext[,3], labs, col=cols  )
  
  #rgl.postscript( paste(dir,"/sample_",tag,"_.pdf", sep=""), fmt="pdf") 
  #rgl.snapshot( paste(dir,"/sample_",tag,"_.png", sep="") , fmt="png")
}

cc <- 0
sq <- c(0:90,89:1)
for(i in sq[sq %% length(sq) == 0] ){
  cc <- cc+1
  vpar1 <- c(i,i)
  cst <- as.character(cc)
  if(nchar(cst) == 1) cst <- paste("00",cst, sep="")
  if(nchar(cst) == 2) cst <- paste("0",cst, sep="")
  
  make.chart(tag=paste("angle_",cst,sep=""),vpar=vpar1)
}
rm(i,cc)


#animates ------
# M <- par3d("userMatrix")
# if (!rgl.useNULL())
#   play3d( par3dinterp(time=(0:2)*1.5,userMatrix=list(M,
#                                                       rotate3d(M, pi/2, 1, 0, 0),
#                                                       rotate3d(M, pi/2, 0, 1, 0) ) ),          
#           duration=6 )