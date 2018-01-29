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
#' The names of the input files are recorded in the \code{infiles} field of the
#' returned object.  This is for the your information only; none of the code
#' here does anything with the recorded filenames.  You have the choice of
#' recording either absolute filenames (useful for determining exactly which
#' dataset was used) or relative filenames (useful if you might distribute a
#' saved model to other users who might not have the same directory structure as
#' you).
#'
#' @param dat A single directory name, or a list of netCDF files.  If a
#' directory name is given, all netCDF files in the directory will be used.
#' @param latvar Name of the latitude coordinate variable in the netCDF files.
#' Normally this is \code{'lat'}, but occasionally it will be something
#' different, such as \code{'lat_2'}.
#' @param lonvar Name of the longitude coordinate variable in the netCDF files.
#' @param meanfield Function to compute the mean temperature response field.
#' The default is a linear pattern scaling algorithm.
#' @param record_absolute If \code{TRUE}, record absolute paths for the input
#' files; otherwise, record relative paths.
#' @return A \code{fldgen} object.
#' @export
train <- function(dat, latvar='lat', lonvar='lon', meanfield=pscl_analyze,
                  record_absolute=FALSE)
{
    if(length(dat) == 1 && file.info(dat)$isdir) {
        ## This is a directory.  Replace with the list of netCDF files contained
        ## within.
        dat <- list.files(dat, '\\.nc$', full.names=TRUE)
    }

    if(record_absolute) {
        infiles <- normalizePath(dat)
    }
    else {
        infiles <- dat
    }

    readin <- function(fn) {read.temperatures(fn, latvar=latvar, lonvar=lonvar)}
    griddata <- concatGrids(lapply(dat, readin)) # information about the
                                        # temperature grid
    tgav <- griddata$tas %*% griddata$tgop # Global mean temperature

    pscl <- pscl_analyze(griddata$tas, tgav)

    reof <- eof_analyze(pscl$r, griddata$tgop)

    reof_l <- split_eof(reof, griddata)
    psd <- psdest(reof_l)
    Fx <- sqrt(psd)

    phasecoef <- phase_eqn_coef(Fx)

    fldgen_object(griddata, tgav, pscl, reof, Fx, phasecoef, infiles)
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
#' @param ftran The magnitude of the Fourier transform of the EOF projection coefficients.
#' @param phasecoef Object returned from \code{\link{phase_eqn_coef}}.
#' @param infiles Names of input files used to construct the data.
#' @export
#' @keywords internal
fldgen_object <- function(griddata, tgav, pscl, reof, ftran, phasecoef, infiles)
{
    Fx <- list(                         # We used to have other stuff in this
                                        # list, but it's obsolete now.
        mag = abs(ftran))
    fg <- list(griddata=griddata, tgav=tgav, pscl=pscl, reof=reof, fx=Fx,
               phasecoef=phasecoef, infiles=infiles)
    class(fg) <- 'fldgen'
    fg
}
