# fldgen: Generate temperature fields with spatial and temporal correlation.
[![Travis-CI Build Status](https://travis-ci.org/JGCRI/fldgen.svg?branch=master)](https://travis-ci.org/JGCRI/fldgen)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/JGCRI/fldgen?branch=master&svg=true)](https://ci.appveyor.com/project/JGCRI/fldgen)
[![Coverage Status](https://img.shields.io/codecov/c/github/JGCRI/fldgen/master.svg)](https://codecov.io/github/JGCRI/fldgen?branch=master)



The `fldgen` package rovides functions to learn the spatial and temporal
correlation of the variability in an ESM and generate random
temperature fields with equivalent properties.

## Installation

The easiest way to install `fldgen` is using `install_github` from the
`devtools` package.
```R
install_github('JGCRI/fldgen', ref='v1.0.0', build_vignettes=TRUE)
```

## Example

The data used in this example are installed with the package.  They
can be found in `system.file('extdata', package='fldgen')`.

```R
library(fldgen)

emulator <- train('extdata/tann1.nc')

tempgrids <- list()                     # Empty list to hold the temperature
                                        # realizations we are about to create.
length(tempgrids) <- 4

##  Run with the phases of the actual time series
meanfield <- pscl_apply(emulator$pscl, emulator$tgav)
tempgrids[[1]] <- reconst_fields(emulator$reof$rotation, mkcorrts(emulator) , meanfield)
## Run the rest with random phases
for(i in 2:4)
    tempgrids[[i]] <- reconst_fields(emulator$reof$rotation, mkcorrts(emulator), meanfield)

```

A more detailed example can be found in the tutorial vignette included
with the package.
```R
vignette('tutorial','fldgen')
```
