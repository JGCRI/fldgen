#' Train the emulator on a set of input ESM temperature and precipitation data.
#'
#' This function runs all of the steps in the analysis that depend only on the
#' input data.  The results are saved in a structure of class \code{fldgen}.
#' This structure contains everything the emulator has learned about the model.
#' Therefore, it can be saved and reloaded in future sessions to bypass the
#' analysis step.
#'
#' The steps run by this procedure are \code{\link{read.temperatures}},
#' \code{\link{read.precipitations}},
#' \code{\link{pscl_analyze}}, and \code{\link{eof_analyze}}.  Additionally,
#' the global mean temperature, global mean precipitation, and the FFT of the
#' EOF decomposition are calculated and stored.  For backward compatibility,
#' all of those functions remain part of the public interface.
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
#' dataset was used, if you have a lot of similarly named files) or relative
#' filenames (useful if you might distribute a saved model to other users who
#' might not have the same directory structure as you).
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
trainTP <- function(dat, varname = "tas", latvar='lat', lonvar='lon',
                    meanfield=pscl_analyze, record_absolute=FALSE)
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


    ### file pairing for T and P
    ### Will likely move to own function in the future, with more flexibility.

    # separate dat into list of precip files and temperature files. Relies on
    # CMIP5 naming conventions.
        pdat <- dat[grep("pr_", dat)]
    if(any(grepl("Aclim", pdat) == FALSE) & any(grepl("annual", pdat) == FALSE)){
        stop(paste("At least one precipitation file in", dat, "is not annual"))
    }

    tdat <- dat[grep("tas_", dat)]
    if(any(grepl("Aclim", tdat) == FALSE) & any(grepl("annual", tdat) == FALSE)){
        stop(paste("At least one temperature file in", dat, "is not annual"))
    }

    ### end file pairing


    readin <- function(fn) {read.temperatures(fn, varname = varname,
                                              latvar=latvar, lonvar=lonvar)}




    griddata <- concatGrids(lapply(tdat, readin)) # information about the
                                                  # temperature grid
    return(griddata)
    tgav <- griddata$tas %*% griddata$tgop # Global mean temperature

    meanfld <- meanfield(griddata$tas, tgav)

    reof <- eof_analyze(meanfld$r, griddata$tgop)

    reof_l <- split_eof(reof, griddata)
    psd <- psdest(reof_l)

    fldgen_object(griddata, tgav, pscl, reof, sqrt(psd$psd), psd$phases, infiles)
}

#' Create a \code{fldgen} object from constituent parts
#'
#' Normally this code will be called internally, but it is available for
#' external use, so that data created using the old interface can be converted.
#'
#' If there is only a single ESM run in the input data set, then \code{fxphase}
#' can be passed as a matrix of phases (instead of a list of matrices).  It will
#' be converted into a list automatically.
#'
#' @param griddata An object returned from \code{\link{read.temperatures}} or
#' \code{\link{concatGrids}}.
#' @param tgav Global mean temperature for the grids in griddata.
#' @param pscl Object returned from \code{\link{pscl_analyze}}.
#' @param reof Object returned from \code{\link{eof_analyze}}.
#' @param fxmag The magnitude of the Fourier transform of the EOF projection
#' coefficients.  This should be a matrix [Ntime x NEOF].  If using
#' \code{\link{psdest}}, note this is the square root of the power spectral
#' density returned by that function.
#' @param fxphase List of matrices [Ntime x NEOF] of phases of the Fourier
#' transform.  There should be one element in the list for each input ESM run.
#' @param infiles Names of input files used to construct the data.
#' @export
#' @keywords internal
fldgen_object <- function(griddata, tgav, pscl, reof, fxmag, fxphase, infiles)
{
    if(is.matrix(fxphase)) {
        phases <- list(fxphase)
    }
    else {
        phases <- fxphase
    }
    Fx <- list(mag = fxmag,
               phases = phases)
    fg <- list(griddata=griddata,
               tgav=tgav,
               pscl=pscl,
               reof=reof,
               fx=Fx,
               infiles=infiles)
    class(fg) <- 'fldgen'
    fg
}
