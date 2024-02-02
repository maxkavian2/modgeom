
# TEST 1 - resolution of optimiziation problems such as ||B - Ax|| -----------------------------
# Operative R version
A <- t(matrix(c(1,1,1,2,3,4,5,5,2,4,2,5,5,4,3),nrow=3,ncol=5))
B <- t(matrix(c(-10,-3,12,14,10,12,16,16,18,16), nrow=2,ncol=5))
R <- qr.solve(A,B, tol = 1e-12)
print(A)
print(B)
print(R)
# Native implementation here (currently deactivated)
# test_dgels(A,B)


# TEST 2 - computation of A %*% t(A) expressions ---------------------------
A <- t(matrix(c(1,1,1,2,3,4,5,5,2,4,2,5,5,4,3),nrow=3,ncol=5))
AA <- A %*% t(A)
print(AA)
# Native implementation here (currently deactivated)
#test_dsyrk(A)


# TEST 3 bspline solve with time assessment -------------------------
position <- seq(from = 0, to = 1, length.out = 100000)
c_point <- t( rbind(c(0,1,2,3,3,1,1,1,5,5,5,7,7,7,0,1,2,3,3,1,1,1,5,5,5,7,7,7), 
                    c(0,3,3,3,3,2,1,3,6,7,-1,2,3,4,0,3,3,3,3,2,1,3,6,7,-1,2,3,4), 
                    c(2,1,1,3,4,5,6,7,-2,-3,-4,-1,-2,-6,2,1,1,3,4,5,6,7,-2,-3,-4,-1,-2,-6)) )
deg <- 6
normal.dev <- .3
rand_position <- runif(100000)


R <- interpolate_bspline(position, c_point, degree=deg)
R_fuzzy <- interpolate_bspline(rand_position, c_point, degree=deg)
R_fuzzy <- rnorm(length(R_fuzzy), mean = 0, sd = normal.dev) + R_fuzzy


s1 <- Sys.time()
bspline_solve(R_fuzzy, rand_position, c_point, degree=deg, native = F)
s2 <- Sys.time()

ss1 <- Sys.time()
bspline_solve(R_fuzzy, rand_position, c_point, degree=deg, native=T)
ss2 <- Sys.time()

r1 <- s2 - s1
r2 <- ss2 - ss1

print(paste(r1," ** ",r2, sep=""))




# TEST 4 footpoint solve tests -------------------------------

# the positions which are interpolated
#position <- seq(from = 0, to = 1, length.out=800)

c_point <- t( rbind(c(0,1,2,3,3), c(0,3,3,3,3), c(2,1,1,3,4)) )

# the degree
deg <- 3

# the number of data points
n_data <- 5

# the data points
d_point <- t( rbind( runif(n_data)*(max(c_point[,1])-min(c_point[,1]))+min(c_point[,1]), 
                     runif(n_data)*(max(c_point[,2])-min(c_point[,2]))+min(c_point[,2]), 
                     runif(n_data)*(max(c_point[,3])-min(c_point[,3]))+min(c_point[,3])) )


# COMPUTATION 
# approximated parameters
s1 <- Sys.time()
R1 <- bspline_footpoint(d_point, c_point, 
                          degree=deg, step=1, tolerance=1e-9, verbose=F, prompt.mode="distance", prompt.begin = 0,
                          algorithm="directional", max.iterations=5000, fp.native=F)
s2 <- Sys.time()
ss1 <- Sys.time()
bspline_footpoint(d_point, c_point, 
                  degree=deg, step=1, tolerance=1e-9, verbose=F, prompt.mode="distance", prompt.begin =0,
                  algorithm="directional", max.iterations=5000, fp.native=T)
ss2 <- Sys.time()

r1 <- s2 - s1
r2 <- ss2 - ss1

print(paste(r1," ** ",r2, sep=""))


# param <- c(0.5,1) #seq(from=0, to=1, length.out=nrow(c_point)*8)
# R3 <- interpolate_bspline(param, c_point, degree= deg, derivate=0, knots = seq(from = 0, to = 1, length.out = nrow(c_point) - deg + 1), is.periodic=F,
#                     native=T)
# print(R3)




