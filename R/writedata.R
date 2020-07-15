#### Output functions

#' Subset a trained emulator.
#'
#' A trained fldgen emulator features a large amount of data for both
#' using the emulator and rigorously validating an emulator.
#'
#' If one is just interested in the use of an emulator for generating
#' felds, this function can be called to reduce a trained emulator to
#' the bare essential list entries, which can then be saved and called
#' the same as an unreduced emulator by generate.TP.resids and
#' generate.TP.fullgrids
#'
#' Note that with this reduced emulator, there is NO way to reconstruct
#' the training data. A fully trained emulator contains a copy of the
#' training data, in addition to the training regressor values (tgav),
#' and the estimated linear model parameters and residuals
#'  (meanfieldT$b, w, r), which together can also reconstruct the data.
#'
#' Even though the coordinate information stored in an emulator$griddataT
#' is not needed directly to generate a new field of residuals or full data,
#' it is often needed in downstream use of the fields. Therfore an entry
#' reducedEmulator$griddataT$coord containg a matrix is saved in the
#' reducedEmulator. Each is a matrix of coordinates for each grid cell, with
#' cells in rows and latitude, longitude in the two columns. Keeping these
#' coordinate matrices for T and P is negligible.
#'
#' Finally, the reduced emulator produced by this function is specifically
#' meant for temperature and precipitation only, and is not robust to
#' extension to other variables.
#'
#' Finally finally, if a user is interested in a different subset of
#' list entries in a trained emulator, they are encouraged to subset and
#' save themself, as appropriate for their project.
#'
#' @param emulator A trained fldgen emulator, with all entries needed
#' for generating new residuals and for rigourously validating the
#' quality of the trained emulator
#'
#' @return reducedEmulator A trained fldgen emulator with only the list
#' entries needed by generate.TP.resids and generate.TP.fullgrids for
#' generating new fields:
#' \describe{
#' \item{griddataT}{Only the coordinate ids and set information.}
#' \item{griddataP}{Only the coordinate ids and set information, and the
#' function to convert from logP to P.}\
#' \item{tgav}{The Tgav data from training.}
#' \item{meanfldT}{the slope (w) and intercept (b) terms from the mean field
#' fit.}
#' \item{meanfldP}{the slope (w) and intercept (b) terms from the mean field
#' fit.}
#'\item{tfuns}{The empirical quantile functions for temperature, mapping
#'N(0,1) to the native distribution in each grid cell.}
#'\item{pfuns}{The empirical quantile functions for logP, mapping
#'N(0,1) to the native distribution in each grid cell.}
#'\item{reof}{The EOFs.}
#'\item{fx}{Time coefficients for each EOF from training data.}
#'\item{infiles}{The names of the files used for training the emulator.}
#' }
#'
#' @author ACS July 2020
#' @export
emulator_reducer <- function(emulator){

    if(length(names(emulator)) < 10){  # a full emulator has 10 list entries, check
        # to make sure that's showing up.
        stop('Your emulator is already reduced (missing at least one list entry)')
    }

    # This function reduces the size of the object while preserving the structure
    # expected by generate.TP.resids and generate.TP.fullgrids.
    list(griddataT = list(gridid_full = emulator$griddataT$gridid_full,
                          coord = emulator$griddataT$coord),
         griddataP = list(gridid_full = emulator$griddataP$gridid_full,
                          coord = emulator$griddataP$coord,
                          pvarconvert_fcn = emulator$griddataP$pvarconvert_fcn),
         tgav = emulator$tgav,
         # not reconstructing training data, don't need residuals in the
         # mean fields
         meanfldT = list(w = emulator$meanfldT$w,
                           b = emulator$meanfldT$b),
         meanfldP = list(w = emulator$meanfldP$w,
                           b = emulator$meanfldP$b),
         tfuns = list(quant = emulator$tfuns$quant),
         pfuns = list(quant = emulator$pfuns$quant),
         reof = emulator$reof,
         fx = emulator$fx,
         infiles = emulator$infiles) ->
        reducedEmulator

    return(reducedEmulator)
}


#' Write a temperature field as a netcdf file.
#'
#' Format a field as a netcdf file and write it to the specified file.  The lat,
#' lon, and time variables from the griddata input are used to recover the
#' original size of the grid, reversing the flattening procedure.
#'
#' @param fld A matrix [ntime x ngrid] containing the field to write out.
#' @param file Name of the output file.
#' @param griddata A griddata structure returned from the
#' \code{\link{read.temperatures}} function.
#' @param varname Name to use for the variable in the output file.
#' @param varunit Units of the variable in the field.
#' @param vardesc Description string to write into the \code{longname} attribute
#' in the output file.
#' @param tunit Units for the time variable.
#' @param clobber If \code{TRUE}, overwrite the file if it exists; otherwise
#' trying to write  a file that already exists is an error.
#' @export
write.temperature <- function(fld, file, griddata, varname='tas', varunit='K',
                              vardesc='2m air temperature', tunit = 'years since 2006',
                              clobber=FALSE)
{
    ## TODO: support writing multiple fields to a file and adding to (instead of
    ## overwriting) an existing file.
    if(clobber) {
        unlink(file)
    }
    else if(file.exists(file)) {
        stop('File ', file, ' exists and clobber=FALSE.')
    }

    ## Define the dimensions for the netcdf file
    lon <- ncdf4::ncdim_def('lon', 'degreesE', griddata$lon)
    lat <- ncdf4::ncdim_def('lat', 'degreesN', griddata$lat)
    nlon <- length(griddata$lon)
    nlat <- length(griddata$lat)
    t <- ncdf4::ncdim_def('time', tunit, griddata$time, unlim=TRUE)
    ntime <- length(griddata$time)

    ## There shouldn't be any missing data in fields we generate, so the fill
    ## value we use is arbitrary
    fillval <- -5555

    tas_def <- ncdf4::ncvar_def(varname, varunit, list(lon, lat, t), fillval,
                            vardesc, prec='float', compression=9)
    ncout <- ncdf4::nc_create(file, tas_def, force_v4=TRUE)

    ## Reorganize the data for output
    fld3d <- fld
    dim(fld3d) <- c(ntime, nlat, nlon)
    fld3d <- aperm(fld3d, c(3,2,1))

    ncdf4::ncvar_put(ncout, tas_def, fld3d)

    ncdf4::nc_close(ncout)
}
