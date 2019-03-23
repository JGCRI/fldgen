#' Generate temperature fields with spatial and temporal correlation
#'
#' Provides functions to learn the spatial and temporal correlation of the
#' variability in an ESM and generate random temperature fields with equivalent
#' properties.
#'
#' @section Example:
#'
#' The data used in this example can be found in \code{system.file('extdata', package='fldgen')}
#'
#' \verb{
#' library(fldgen)
#' datadir <- system.file('extdata', package='fldgen')
#'
#'
#' infileT <- file.path(datadir, 'tas_annual_esm_rcp_r2i1p1_2006-2100.nc')
#' infileP <- file.path(datadir, 'pr_annual_esm_rcp_r2i1p1_2006-2100.nc')
#' emulator <- trainTP(c(infileT, infileP),
#'                     tvarname = "tas", tlatvar='lat', tlonvar='lon',
#'                     tvarconvert_fcn = NULL,
#'                     pvarname = "pr", platvar='lat', plonvar='lon',
#'                     pvarconvert_fcn = log)
#'
#'
#' residgrids <- generate.TP.resids(emulator, ngen = 3)
#'
#' fullgrids <- generate.TP.fullgrids(emulator, residgrids,
#'                                    tgav = tgav,
#'                                    tvarunconvert_fcn = NULL,
#'                                    pvarunconvert_fcn = exp,
#'                                    reconstruction_function = pscl_apply)
#' }
#'
#' For more examples, see the tutorial vignette (\code{vigentte('tutorial2', 'fldgen')}).
"_PACKAGE"
