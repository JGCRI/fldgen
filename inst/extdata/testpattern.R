## Run this in the package root directory to generate the test pattern input.

## grid parameters
nt <- 8                                 # number of time steps
nlat <- 4
nlon <- 4
ngrid <- nlat*nlon
fn <- 'inst/extdata/testpattern.nc'

## set up some basis vectors.  Some of these may need to be rethought if we
## change the grid size.
v1 <- rep(1,ngrid)
v2 <- c(rep(-1,ngrid/2), rep(1,ngrid/2))

vm <- matrix(0, nrow=nlat, ncol=nlon)
vm[1,] <- 1
vm[2,] <- -1
v3 <- as.vector(vm)

vm <- matrix(0, nrow=nlat, ncol=nlon)
vm[3,] <- 1
vm[4,] <- -1
v4 <- as.vector(vm)

## time vector
t <- seq(1,nt)
## add jitter to global mean
set.seed(8675309)
ts <- matrix(1.0e-2 * rnorm(nt), nrow=nt) %*% matrix(v1, nrow=1)

## for each of the other basis vectors, add periodic dependence that doesn't
## match up with the time grid.
x <- (t-1)/5 * 2*pi

ts <- ts + matrix(cos(2*x), nrow=nt) %*% matrix(v2, nrow=1)
ts <- ts + matrix(cos(3*x - pi/2), nrow=nt) %*% matrix(v3, nrow=1)
ts <- ts + matrix(cos(4*x - pi/4), nrow=nt) %*% matrix(v4, nrow=1)

## Create a griddata structure so we can write this to a netCDF file.
## We want flat geometry, so we will make all of the lat coordinates 0.  This
## will probably make the test pattern fail to display in a lot of netcdf
## visualization software.
griddata <- list(tas=ts, tgop=v1, lat=rep(0,nlat), lon=seq(1,nlon), time=t)

write.temperature(ts, fn, griddata, clobber=TRUE)

