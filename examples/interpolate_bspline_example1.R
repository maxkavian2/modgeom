library(rgl)
 
#the position which is interpolated
position <- seq(from = 0, to = 1, length.out=2000)

# the matrix of control points here
c_point <- t(rbind( c(0, 5, 0, 12, 10, 11, 10, 1, 1, 7),
                    c(4, 12, -7, 3, -5, 10, 4, 8, 9, 10),
                    c(5,  10,  1, 1,  -5,  1, 1, 1, 1, 10 ) ))
# the degree
deg <- 9 

R <- list()
Rtext <- c()
for(dg in 0:(deg-1) ){
  R[[dg+1]] <- interpolate_bspline(position, c_point, degree = deg-dg, native=TRUE)
  Rtext <- rbind(Rtext, interpolate_bspline(.16, c_point, degree = deg-dg, native=TRUE))
}
Rtext[,1] = Rtext[,1]-0.6

colscale <- c("red","gray")
colf <- colorRampPalette(colscale)
cols <- colf(length(R))

#drawing -----
make.chart <- function(tag="", vpar=NULL){
mar.red <- 0.1
options(scipen = 5)
par(mgp=c(2.4,1,.3))
par(mar=c(4,5-mar.red,2.5-mar.red,2-mar.red)+0.1) 
par(cex=1)
  
# 3D pars

  
plot3d(R[[1]], type="l", 
         aspect=TRUE,  
         col=cols[1], lwd=1, main = "", xlab="x", ylab="y", zlab="z", box=FALSE, 
         axes=FALSE)
  
par3d(windowRect=c(150, 150, 750, 750))
#rgl.viewpoint(vpar[1], vpar[2])
view3d(vpar[1], vpar[2])

#par3d(zoom = 0.7)

  axes3d(
    edges=c('x--', 'y--', 'z--'),
    labels=FALSE,
    tick=TRUE, at=c(-4:8)
    
  )
  

if(length(R) > 1)
    for(i in 2:length(R) )
      lines3d(R[[i]], col = cols[i], lwd = 1 )
  
#draws the control points 
points3d( c_point, pch=2, col="gray", size = 6)
  
# prints the texts
labs <- c()
  for(i in deg:1)
  labs <- c(labs, paste("C",as.character(i-1), sep=""))
  
rgl.texts(Rtext[,1], y =Rtext[,2], z=Rtext[,3], labs, col=cols  )


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