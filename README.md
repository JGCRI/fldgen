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
datadir <- system.file('extdata', package='fldgen')


infileT <- file.path(datadir, 'tas_annual_esm_rcp_r2i1p1_2006-2100.nc')
infileP <- file.path(datadir, 'pr_annual_esm_rcp_r2i1p1_2006-2100.nc')
emulator <- trainTP(c(infileT, infileP),
                    tvarname = "tas", tlatvar='lat', tlonvar='lon',
                    tvarconvert_fcn = NULL,
                    pvarname = "pr", platvar='lat', plonvar='lon',
                    pvarconvert_fcn = log)


residgrids <- generate.TP.resids(emulator, ngen = 3)

fullgrids <- generate.TP.fullgrids(emulator, residgrids,
                                   tgav = tgav,
                                   tvarunconvert_fcn = NULL,
                                   pvarunconvert_fcn = exp,
                                   reconstruction_function = pscl_apply)
```

A more detailed example can be found in the tutorial vignette included
with the package.
```R
vignette('tutorial2','fldgen')
```
