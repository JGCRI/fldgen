# fieldgenr

Generate temperature fields with spatial and temporal correlation.


This package rovides functions to learn the spatial and temporal
correlation of the variability in an ESM and generate random
temperature fields with equivalent properties.

## Example

The data used in this example can be found in `tests/testthat/data`.

```R
library(fieldgenr)

griddata <- readdata('data/tann1.nc')
tgav <- readtgav('data/wgttann1.nc')
pscl <- pscl_analyze(griddata$tas, tgav)
reof <- eof_analyze(pscl$r, griddata$tgop)
Fx <- mvfft(reof$x)                     # Fourier transforms of the coordinates
                                        # of the basis functions
Fxmag <- abs(Fx)
Fxphase <- atan2(Im(Fx), Re(Fx))
tempgrids <- list()                     # Empty list to hold the temperature
                                        # realizations
length(tempgrids) <- 4

##  Run with the phases of the actual time series
meanfield <- pscl_apply(pscl, tgav)
tempgrids[[1]] <- reconst_fields(reof$rotation, mkcorrts(Fxmag, Fxphase) , meanfield)
## Run the rest with random phases
for(i in 2:4)
    ## If you wanted to use a smaller number of PCs, say 50, you could
    ## use reof$rotation[,1:50] and Fxmag[,1:50] in the next call.
    tempgrids[[i]] <- reconst_fields(reof$rotation, mkcorrts(Fxmag), meanfield)

```
