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
#' @section Supplying global mean temperatures:
#'
#' In some cases it may be desirable to supply global mean temperatures, rather
#' than calculating them from the grid.  For example, the grid may cover only
#' part of the grid, so that the average of the values would not be a true
#' global average.  In such cases the true global mean temperatures can be
#' supplied in a separate file.
#'
#' If this input file is used, it should be a plain text file with one
#' temperature value per line, each line corresponding to one time slice in the
#' corresponding netCDF file containing the grid data.  The name of the file
#' will be constructed by adding a suffix (specified in the
#' \code{globalAvg_file} argument) to the name of the netCDF file.  For example,
#' if the netCDF file is called \code{larry.nc}, and the suffix is
#' \code{curly.txt}, the name of the global mean temperature file will be
#' \code{larry.nccurly.txt}.  If multiple netCDF files are passed in, each one
#' will have its own file of global mean temperatures, constructed in the same
#' way.
#'
#' @param dat A single directory name, or a list of netCDF files.  If a
#' directory name is given, all netCDF files in the directory will be used.The
#' pairing of temperature and precipitation netCDF files in the directory
#' relies on the CMIP5 file naming conventions. Other naming conventions are
#' not currently supported.
#' @param tvarname Name of the temperature variable in the temperature netCDF
#' @param tlatvar Name of the latitude coordinate variable in the temperature
#' netCDF files.
#' Normally this is \code{'lat'}, but occasionally it will be something
#' different, such as \code{'lat_2'}. Can vary between paired T and P.
#' @param tlonvar Name of the longitude coordinate variable in the temperature
#' netCDF files.  Can vary between paired T and P.
#' @param tvarconvert_fcn The function used to transform the T variable prior
#' to training so that it has support on -infinity to infinity. Defaults to
#' NULL, as Temperature is effectively already supported on this range.
#' @param pvarname Name of the precipitation variable in the precipitation
#' netCDF.
#' @param platvar Name of the latitude coordinate variable in the precipitation
#' netCDF files.
#' Normally this is \code{'lat'}, but occasionally it will be something
#' different, such as \code{'lat_2'}.  Can vary between paired T and P.
#' @param plonvar Name of the longitude coordinate variable in precipitation
#'  netCDF files.  Can vary between paired T and P.
#' @param pvarconvert_fcn The function used to transform the P variable prior
#' to analyis so that it has support on -infinity to infinity. Defaults to
#' log() as precipitation values cannot be less than 0.
#' @param meanfield Function to compute the mean temperature response field.
#' The default is a linear pattern scaling algorithm.
#' @param globalAvg_file Optional string to use in constructing the name of the
#' file containing global mean temperatures, Tgav.  If omitted, calculate Tgav
#' as the mean of the grid cells in each frame. (See "Supplying global mean
#' temperatures" for more information.)
#' @param record_absolute If \code{TRUE}, record absolute paths for the input
#' files; otherwise, record relative paths.
#' @return A \code{fldgen} object.
#' @importFrom tibble as_tibble
#' @importFrom stats na.omit
#' @importFrom dplyr mutate select left_join %>%
#' @importFrom tidyr separate
#' @export
trainTP <- function(dat,
                    tvarname = "tas", tlatvar='lat', tlonvar='lon',
                    tvarconvert_fcn = NULL,
                    pvarname = "pr", platvar='lat', plonvar='lon',
                    pvarconvert_fcn = logPfloor,
                    meanfield=pscl_analyze,
                    globalAvg_file = NULL,
                    record_absolute=FALSE)
{

    ## silence package checks
    value <- var <- timestep <- time <- startyr <- stopyr <-
        tfilename <- pfilename <- NULL


    ## Make sure working with list of file names
    if(length(dat) == 1 && file.info(dat)$isdir) {
        ## This is a directory.  Replace with the list of netCDF files contained
        ## within.
        dat <- list.files(dat, '\\.nc$', full.names=TRUE)
    }


    ## pair the T and P files
    if(!is.null(dim(dat))){
        if(ncol(dat) == 2){
            ## This is a paired list of input files.
            ## It should follow the formalism
            ## column1 = temperature,
            ## column2 = precipitation file names.
            ## Entry should each be a netcdf name.
            paireddat <- tibble::tibble(tfilename = dat[,1],
                                        pfilename = dat[,2])
        }
        else{
            stop("Unrecognized input structure for training data.")
        }
    }
    else{
        ## Otherwise files should be name according to the CMIP5 convention so
        ## that they may be paired.
        paireddat <- file.pairer(dat, tvarname = tvarname, pvarname = pvarname)
    }


    ## Prepare the list of input files that contribute to the training for output
    infiles <- as.vector(as.matrix(paireddat))
    if(record_absolute) {
        infiles <- normalizePath(infiles)
    }

    ### Train

    # read the relevant Temperature and Precipitation inputs across all
    # scenarios and concatenate
    readinT <- function(fn) {read.general(fn, varname = tvarname,
                                          latvar=tlatvar, lonvar=tlonvar)}
    readinP <- function(fn) {read.general(fn, varname = pvarname,
                                          latvar=platvar, lonvar=plonvar)}



    griddataT <- concatGrids.general(lapply(paireddat$tfilename, readinT))
    griddataP <- concatGrids.general(lapply(paireddat$pfilename, readinP))


    # make sure supported on -infinity to infinity, add some output data to
    # facilitate the reversing later
    if(!is.null(tvarconvert_fcn)){
        griddataT$vardata_raw <- griddataT$vardata
        griddataT$vardata <- tvarconvert_fcn(griddataT$vardata)
    } else {
        griddataT$vardata_raw <- NULL
    }
    griddataT$tvarconvert_fcn <- tvarconvert_fcn


    if(!is.null(pvarconvert_fcn)){
        griddataP$vardata_raw <- griddataP$vardata
        griddataP$vardata <- pvarconvert_fcn(griddataP$vardata)
    } else {
        griddataP$vardata_raw <- NULL
    }
    griddataP$pvarconvert_fcn <- pvarconvert_fcn

    # Remove the gridcells that only contain NA values.
    griddataT <- drop_NAs(griddataT)
    griddataP <- drop_NAs(griddataP)

    # calculate global average T
    if(is.null(globalAvg_file)){

        # Calculate the global mean Temperature internally
        tgav <- griddataT$vardata %*% griddataT$globalop

    } else {

        tgav <- read_globalAvg(paireddat$tfilename,
                               globalAvg_file,
                               griddataT$vardata,
                               paireddat)

    }


    # Use this to calculate the mean field for T and the mean field for P
    # separately. We want the joint variability on the residuals, but we
    # want those residuals
    meanfldT <- meanfield(griddataT$vardata, tgav)
    meanfldP <- meanfield(griddataP$vardata, tgav)


    # characterize the distribution in each grid cell for each variable by
    # generating an empirical CDF function and an empirical quantile function
    # for each variable in each grid cell.
    tfuns <- characterize.emp.dist(meanfldT$r)
    pfuns <- characterize.emp.dist(meanfldP$r)


    # Use the empirical cdf to normalize residuals for both T and P
    normresidsT <- normalize.resids(inputresids = meanfldT$r,
                                    empiricalcdf = tfuns$cdf)
    normresidsP <- normalize.resids(inputresids = meanfldP$r,
                                    empiricalcdf = pfuns$cdf)


    # bind columns of the matrices to create joint matrix of residuals
    joint_residuals <- as.matrix(cbind(normresidsT$rn, normresidsP$rn))



    # EOF decomposition
    reof <- eof_analyze(joint_residuals, Ngrid = ncol(normresidsT$rn),
                        globop = griddataT$globalop)


    # split the EOFs by ESM run.
    # tags should be shared by T and P so doesn't matter which griddata you use
    reof_l <- split_eof(reof, griddataT)
    psd <- psdest(reof_l)

    # output a generalized fldgen object
    fldgen_object_TP(griddataT, griddataP,
                     tgav, meanfldT, meanfldP,
                     tfuns, pfuns,
                     reof, sqrt(psd$psd), psd$phases, infiles)
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
#' @param griddataT An object returned from \code{\link{read.temperatures}} or
#' \code{\link{concatGrids}}.
#' @param griddataP An object returned from \code{\link{read.temperatures}} or
#' \code{\link{concatGrids}}.
#' @param tgav Global mean temperature for the grids in griddata.
#' @param meanfldT Object returned from whatever mean field analysis is done,
#' currently \code{\link{pscl_analyze}}.
#' @param meanfldP Object returned from whatever mean field analysis is done,
#' currently \code{\link{pscl_analyze}}.
#' @param tfuns The empirical cdf and empirical quantile functions for each
#' grid cell, characterizing the distribution of residuals resulting from
#' meanfldT.
#' @param pfuns The empirical cdf and empirical quantile functions for each
#' grid cell, characterizing the distribution of residuals resulting from
#' meanfldP.
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
fldgen_object_TP <- function(griddataT, griddataP,
                             tgav, meanfldT, meanfldP,
                             tfuns, pfuns,
                             reof, fxmag, fxphase, infiles)
{
    if(is.matrix(fxphase)) {
        phases <- list(fxphase)
    }
    else {
        phases <- fxphase
    }
    Fx <- list(mag = fxmag,
               phases = phases)
    fg <- list(griddataT=griddataT,
               griddataP=griddataP,
               tgav=tgav,
               meanfldT=meanfldT,
               meanfldP=meanfldP,
               tfuns = tfuns,
               pfuns = pfuns,
               reof=reof,
               fx=Fx,
               infiles=infiles)
    class(fg) <- 'fldgen'
    fg
}

logPfloor <- function(x)
{
    ## Log transformation of precipitation, with a floor of 1e-9 (slightly less than
    ## 0.1mm for the year) on the pre-transform values
    x <- pmax(x, 1e-9)
    log(x)
}
