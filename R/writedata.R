#### Output functions

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

#' Load and save emulator training data
#'
#' \code{savemodel} saves the results of training an emulator in a portable
#' format.  \code{loadmodel} loads a model from a file created this way and
#' returns it as a \code{fldgen} object.
#'
#' @param modeldata A \code{fldgen} object returned by either
#' \code{\link{train}} or \code{\link{fldgen_object}}.
#' @param file Name of the file to write the data to.
#' @param clobber Flag indicating whether it's ok to overwrite an existing file
#' @name saving_and_restoring
NULL

#' @rdname saving_and_restoring
#' @export
savemodel <- function(modeldata, file, clobber=FALSE)
{
    if(!inherits(modeldata, 'fldgen')) {
        stop('modeldata must be a fldgen object.')
    }

    if(!clobber && file.exists(file)) {
        stop('File ', file, ' exists, and noclobber is set.')
    }
    save(modeldata, file=file, compress='bzip2', compression_level=9)
}
