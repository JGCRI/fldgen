#' Train the emulator on a set of input ESM data.
#'
#' This function runs all of the steps in the analysis that depend only on the
#' input data.  The results are saved in a structure of class \code{fldgen}.
#' This structure contains everything the emulator has learned about the model.
#' Therefore, it can be saved and reloaded in future sessions to bypass the
#' analysis step.
#'
#' The steps run by this procedure are \code{\link{read.temperatures}},
#' \code{\link{pscl_analyze}}, and \code{\link{eof_analyze}}.  Additionally, the
#' global mean temperature and the FFT of the EOF decomposition are calculated
#' and stored.  For backward compatibility, all of those functions remain part
#' of the public interface.
#'
#' At present the field generation functions still take as arguments the
#' individual components of the structure; however, in the next version they
#' will expect a \code{fldgen} object, which can be constructed manually using
#' \code{\link{fldgen_object}}.
#'
#' @param dat A single directory name, or a list of netCDF files.  If a
#' directory name is given, all netCDF files in the directory will be used.
#' @param latvar Name of the latitude coordinate variable in the netCDF files.
#' Normally this is \code{'lat'}, but occasionally it will be something
#' different, such as \code{'lat_2'}.
#' @param lonvar Name of the longitude coordinate variable in the netCDF files.
#' @return A \code{fldgen} object.
#' @export
train <- function(dat, latvar='lat', lonvar='lon')
{
    if(length(dat) == 1 && file.info(dat)$isdir) {
        ## This is a directory.  Replace with the list of netCDF files contained
        ## within.
        dat <- list.files(dat, '\\.nc$', full.names=TRUE)
    }

    readin <- function(fn) {read.temperatures(fn, latvar=latvar, lonvar=lonvar)}
    griddata <- concatGrids(lapply(dat, readin)) # information about the
                                        # temperature grid
    tgav <- griddata$tas %*% griddata$tgop # Global mean temperature

    pscl <- pscl_analyze(griddata$tas, tgav)

    reof <- eof_analyze(pscl$r, griddata$tgop)

    Fx <- mvfft(reof$x)

    phasecoef <- phase_eqn_coef(Fx)

    fldgen_object(griddata, tgav, pscl, reof, Fx, phasecoef)
}

#' Create a \code{fldgen} object from constituent parts
#'
#' Normally this code will be called internally, but it is available for
#' external use, so that data created using the old interface can be converted.
#'
#' @param griddata An object returned from \code{\link{read.temperatures}} or
#' \code{\link{concatGrids}}.
#' @param tgav Global mean temperature for the grids in griddata.
#' @param pscl Object returned from \code{\link{pscl_analyze}}.
#' @param reof Object returned from \code{\link{eof_analyze}}.
#' @param ftran The Fourier transform of the EOF projection coefficients.
#' @param phasecoef Object returned from \code{\link{phase_eqn_coef}}.
#' @export
#' @keywords internal
fldgen_object <- function(griddata, tgav, pscl, reof, ftran, phasecoef)
{
    Fx <- list(
        fx = ftran,
        mag = abs(ftran),
        phase = atan2(Im(ftran), Re(ftran)))
    fg <- list(griddata=griddata, tgav=tgav, pscl=pscl, reof=reof, fx=Fx, phasecoef=phasecoef)
    class(fg) <- 'fldgen'
    fg
}
