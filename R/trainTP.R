#' Train the emulator on a set of input ESM temperature and precipitation data.
#'
#' This function runs all of the steps in the analysis that depend only on the
#' input data.  The results are saved in a structure of class \code{fldgen}.
#' This structure contains everything the emulator has learned about the model.
#' Therefore, it can be saved and reloaded in future sessions to bypass the
#' analysis step.
#'
#' The steps run by this procedure are \code{\link{read.general}},
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
#' directory name is given, all netCDF files in the directory will be used.The
#' pairing of temperature and precipitaiton netCDF files in the directory
#' relies on the CMIP5 file naming conventions. Other naming conventions are
#' not currently supported.
#' @param tvarname Name of the temperature variable in the temperature netCDF
#' @param tlatvar Name of the latitude coordinate variable in the temperature
#' netCDF files.
#' Normally this is \code{'lat'}, but occasionally it will be something
#' different, such as \code{'lat_2'}. Can vary between paired T and P.
#' @param tlonvar Name of the longitude coordinate variable in the temperature
#' netCDF files.  Can vary between paired T and P.
#' @param pvarname Name of the precipitation variable in the precipitation
#' netCDF.
#' @param platvar Name of the latitude coordinate variable in the precipitation
#' netCDF files.
#' Normally this is \code{'lat'}, but occasionally it will be something
#' different, such as \code{'lat_2'}.  Can vary between paired T and P.
#' @param plonvar Name of the longitude coordinate variable in precipitation
#'  netCDF files.  Can vary between paired T and P.
#' @param meanfield Function to compute the mean temperature response field.
#' The default is a linear pattern scaling algorithm.
#' @param record_absolute If \code{TRUE}, record absolute paths for the input
#' files; otherwise, record relative paths.
#' @return A \code{fldgen} object.
#' @export
trainTP <- function(dat, tvarname = "tas", tlatvar='lat', tlonvar='lon',
                    pvarname = "pr", platvar='lat', plonvar='lon',
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
        pdat <- dat[grep(paste0(pvarname, "_"), dat)]
    if(any(grepl("Aclim", pdat) == FALSE) & any(grepl("annual", pdat) == FALSE)){
        stop(paste("At least one precipitation file in", dat, "is not annual"))
    }

    tdat <- dat[grep(paste0(tvarname, "_"), dat)]
    if(any(grepl("Aclim", tdat) == FALSE) & any(grepl("annual", tdat) == FALSE)){
        stop(paste("At least one temperature file in", dat, "is not annual"))
    }


    # get to scenario identifying info
    pdat %>%
        tibble::as_tibble() %>%
        dplyr::mutate(pfilename = value) %>%
        tidyr::separate(value, c("var", "timestep", "esm", "rcp", "run", "time"), sep = "_") %>%
        dplyr::select(-var, -timestep)  %>%
        tidyr::separate(time, c("startyr", "stopyr"), sep = "-") %>%
        dplyr::mutate(startyr = substr(startyr, 1, 4),
               stopyr = substr(stopyr, 1,4)) ->
        ptbl

    tdat %>%
        tibble::as_tibble() %>%
        dplyr::mutate(tfilename = value) %>%
        tidyr::separate(value, c("var", "timestep", "esm", "rcp", "run", "time"), sep = "_") %>%
        dplyr::select(-var, -timestep)  %>%
        tidyr::separate(time, c("startyr", "stopyr"), sep = "-") %>%
        dplyr::mutate(startyr = substr(startyr, 1, 4),
                      stopyr = substr(stopyr, 1,4)) ->
        ttbl

    ttbl %>%
        left_join(ptbl, by = c("esm", "rcp", "run", "startyr", "stopyr")) %>%
        dplyr::select(tfilename, pfilename) %>%
        na.omit ->
        paireddat

    if(nrow(paireddat) < 1){
        stop("You have no paired temperature and precipitation files in this
             directory. Use the function `train` to work on temperature files
             only. Or ensure that your temperature and precipitation netCDFs
             follow CMIP5 naming conventions.")
    }

    ### end file pairing


    ### Train

    # read the relevant Temperature and Precipitation inputs across all
    # scenarios and concatenate
    readinT <- function(fn) {read.general(fn, varname = tvarname,
                                          latvar=tlatvar, lonvar=tlonvar)}
    readinP <- function(fn) {read.general(fn, varname = pvarname,
                                          latvar=platvar, lonvar=plonvar)}



    griddataT <- concatGrids.general(lapply(paireddat$tfilename, readinT))
    griddataP <- concatGrids.general(lapply(paireddat$pfilename, readinP))


    # calculate global average T and global average P
    tgav <- griddataT$vardata %*% griddataT$globalop # Global mean temperature
    # pgav <- griddataP$vardata %*% griddataP$globalop # Global mean precip


    # Use this to calculate the mean field for T and the mean field for P
    # separately. We want the joint variability on the residuals, but we
    # want those residuals
    meanfldT <- meanfield(griddataT$vardata, tgav)
    meanfldP <- meanfield(griddataP$vardata, tgav)
    return(list(meanfldT, meanfldP))

    # normalize residuals for both because it won;t hurt probably



    # enforce pairing to bind columns of the matrices correctly to create
    # joint matrix of residuals

    joint_residuals <- as.matrix(c(meanfldT$r, meanfldP$r))
    ## need to figure out enforcing pairing based on tags. It should just happen
    ## automatically but it needs a test.
    return(joint_residuals)


    # update eof to include if's on length globalop and ncols meanfld$r so not
    # doing uneccesarily large block matrix calculations.
    reof <- eof_analyze(meanfld$r, # joint matrix of respective residuals
                        griddata$tgop)


    ## figure out needed updates
    reof_l <- split_eof(reof, griddata)
    psd <- psdest(reof_l)


    ## update this guy to general
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
