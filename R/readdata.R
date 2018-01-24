#' Read ESM temperature data
#'
#' The temperature data should be a netCDF file of annual grid values.  The data
#' is assumed to have its dimensions of time, lat, lon, where lon is the most
#' rapidly varying index, and time is the least.
#'
#' The output will be a list with six fields:
#' \describe{
#'   \item{\strong{tas}}{Matrix of temperature data.}
#'   \item{\strong{tgop}}{Vector operator for global mean temperature.}
#'   \item{\strong{lat}}{Vector of latitude values.  These are only replicated
#' once, so their length is nlat, not ngrid.  You have to use \code{rep(lat,
#' nlon)} to get latitudes for each grid cell.}
#'   \item{\strong{lon}}{Vector of longitude values.  These are only replicated
#' once, so their lingth is nlon, not ngrid.  Getting longitudes for each grid
#' cell is tricksy.  Try \code{as.vector(matrix(rep(lon, nlat), nrow=nlat,
#' byrow=TRUE))}.  Fortunately, you don't need them very often.}
#'   \item{\strong{time}}{Vector of time values, given as years since the base
#' year of the dataset.}
#'   \item{\strong{tags}}{A list of datasets that were concatenated to get this
#' structure.  The names of the list are the tags given to the dataset, and the
#' values are a vector of start-row, end-row pairs.}
#' }
#'
#' The data at each time is represented as a flattened vector of grid cells.
#' The flattening is performed by transposing to lat, lon ordering, so that lat
#' will be the most rapidly varying index (because R uses
#' Fortran-style indexing).  Then the spatial dimensions are discarded,
#' resulting in a 1D vector.  The time dimension is kept, resulting in a matrix
#' with years in rows and grid cells in columns.  The dimensions of the matrix
#' will be Nyear x Ngrid, where Nyear is the number of years in the data set,
#' and Ngrid is the number of grid cells.
#'
#' The \code{tgop} vector is a vector whose dot product with a flattened grid is
#' the area-weighted global mean for the grid.  It is stored as a column vector,
#' so \code{tas \%*\% tgop} is the time series of global mean temperatures.
#'
#' The lat and lon dimension variables from the input file  are also
#' stored in the structure.  These are primarily useful for writing out
#' generated grids as netCDF files.  The time dimension variable is converted to
#' integer years, starting at 0.
#'
#' Conventionally, we refer to the output list as \code{griddata}.  Notably, any
#' other function with a \code{griddata} argument is expecting one of these
#' structures.
#'
#' @param filename Name of the input netCDF file
#' @param len Maximum length of the time series to read.  If the data read is
#' longer, it will be trimmed.  (Default: read entire time series, regardless of
#' length.)
#' @param varname Name of the variable to read from the netcdf file.
#' Irrespective of the name in the netcdf file, it will be called 'tas' in the
#' return data.
#' @param tag A string identifying the name of the scenario.  If omitted, the
#' tag will default to the filename.
#' @param latvar Name of the latitude dimension variable in the input file.
#' @param lonvar Name of the longitude dimension variable in the input file.
#' @param timevar Name of the time dimension variable in the input file.
#' @return A \code{griddata} list (see details).
#' @importFrom assertthat assert_that
#' @export
read.temperatures <- function(filename, len=NULL, tag=basename(filename), varname='tas',
                              latvar='lat', lonvar='lon', timevar='time')
{
    tann <- ncdf4::nc_open(filename)

    ## tas3d should have dimensions of time x lat x lon in the netcdf file.
    ## Because R uses Fortran array layout, this gets reversed to lon x lat x time.
    tas3d <- ncdf4::ncvar_get(tann, var=varname)
    lat <- ncdf4::ncvar_get(tann, var=latvar)
    nlat <- length(lat)
    lon <- ncdf4::ncvar_get(tann, var=lonvar)
    nlon <- length(lon)
    time <- ncdf4::ncvar_get(tann, var=timevar)
    ntime <- length(time)
    timeout <- seq_along(time) - 1
    ncdf4::nc_close(tann)

    assert_that(all(dim(tas3d) == c(nlon, nlat, ntime)))

    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            time <- time[1:ntime]
            timeout <- timeout[1:ntime]
            tas3d <- tas3d[ , , 1:ntime]
        }
    }


    ## reorganize and flatten the 3-D array into a 2d array of ntime x ngrid
    ## As we do this, we will also make latitude the most rapidly varying index
    ## for the individual time slices.

    tas <- aperm(tas3d, c(3,2,1))
    dim(tas) <- c(ntime, nlat*nlon)

    ## create a grid cell area factor that has the same shape as the flattened
    ## grid.  Since latitude is the most rapidly varying dimension in the grid,
    ## we just need to replicate cos(latitude), nlon times
    areafac <- rep(cos(lat * pi/180.0), times=nlon)
    ## normalize areafac so that sum(area*grid) is the global weighted average
    areafac <- areafac/sum(areafac)

    ## The tag element of the output is a named list of 2-element vectors,
    ## giving the start and end of the data in the output.  When datasets are
    ## concatenated, this allows us to recover the originals
    tagout <- list(c(1,ntime))
    names(tagout) <- tag

    gd <- list(tas=as.matrix(tas), tgop=as.matrix(areafac), lat=lat, lon=lon,
               time=timeout, tags=tagout)
    class(gd) <- 'griddata'
    gd
}


#' @rdname saving_and_restoring
#' @export
loadmodel <- function(file)
{
    load(file)
    if(!exists('modeldata', inherits=FALSE)) {
        stop('No model data in file.')
    }
    if(!inherits(modeldata, 'fldgen')) {
        stop('Object loaded from file is not of type "fldgen".')
    }

    modeldata
}


#' Read and format global mean temperature
#'
#' Read global mean temperature from an input netCDF file and format for use
#' with the rest of the processing functions.
#'
#' Global mean temperature is a 1-D array, so storing in a netCDF file isn't
#' strictly necessary; however, some tools produce it in that format by
#' default.  A vector of tgav values from any source can be used with the other
#' functions in this package; all you have to do is convert it to an ntime x 1
#' matrix using \code{as.matrix}.
#'
#' @param tgavfilename netCDF file with global mean temperatures variable
#' @return A ntime x 1 matrix of global mean temperature
#' @export
readtgav <- function(tgavfilename)
{
    wgttann <- ncdf4::nc_open(tgavfilename)
    tgav <- ncdf4::ncvar_get(wgttann, var='tas')

    time <- ncdf4::ncvar_get(wgttann, var='time')
    ntime <- length(time)
    ncdf4::nc_close(wgttann)
    assert_that(dim(tgav) == ntime)

    as.matrix(tgav)
}


#' Concatenate a list of griddata objects
#'
#' Griddata objects are concatenated by appending their temperature observations
#' into a single grand matrix with all the time sices from the concatenated
#' objects.
#'
#' The time indices are likewise concatenated so that they will give accurate
#' time values.  The tags are concatenated and updated so that they point to the
#' new location of their original data.  Lat and lon values and the tgop
#' operator are unchanged, and they must be the same for all of the concatenated
#' grids, or an error is raised.
#'
#' Note: Time values are given as years since the start year for the data set.
#' If you concatenate data sets with different start years, the time values will
#' not be compatible.  For \emph{most} of what we do in this package, the
#' difference won't matter, but it's something to be aware of.
#'
#' @param gridlist List of grids to concatenate
#' @return A griddata (q.v. \code{link{read.temperatures}}) object with the concatenated
#' grid values.
#' @export
#' @keywords internal
concatGrids <- function(gridlist)
{
    tasl <- lapply(gridlist, function(x){x$tas})
    tgopl <- lapply(gridlist, function(x){x$tgop})
    latl <- lapply(gridlist, function(x){x$lat})
    lonl <- lapply(gridlist, function(x){x$lon})
    timel <- lapply(gridlist, function(x){x$time})
    tagsl <- lapply(gridlist, function(x){x$tags})

    ## verify that tgop, lat, and lon are compatible.
    checkfun <- function(cmpval) {
        function(x) {isTRUE(all.equal(x, cmpval))}
    }
    all(sapply(tgopl, checkfun(tgopl[[1]]))) ||
      stop('tgop values not compatible')
    all(sapply(latl, checkfun(latl[[1]]))) || stop('lat values not compatible')
    all(sapply(lonl, checkfun(lonl[[1]]))) || stop('lon values not compatible')

    ## Concatenate tas matrix by rows
    tas <- do.call(rbind, tasl)
    ## Concatenate time vectors
    time <- do.call(c, timel)

    ## the tags are a little trickier because we need to update the start and
    ## end rows.  On top of that, we have to deal with the case where some of
    ## the inputs were themselves concatenated.

    ## First, flatten the tags list
    tagsl <- flatten_tags(tagsl)
    ## find the cumulative number of rows to the start of each data set
    n <- cumsum(sapply(tagsl, function(x) {1+x[2]-x[1]}))
    length(n) <- length(n) -1
    n <- c(0,n)
    ## add the cumulative row count to each set
    tags <- Map(`+`, tagsl, n)

    ## Now reconstruct the output structure.
    gd <- list(tas=tas, tgop=tgopl[[1]], lat=latl[[1]], lon=lonl[[1]], time=time, tags=tags)
    class(gd) <- 'griddata'
    gd
}

#' Split a concatenated griddata object into its constituent grids
#'
#' This function undoes the effect of \code{\link{concatGrids}}.  It does not,
#' however, remember the path you took to get the concatenated grid, so
#' \code{splitGrids(c(g1, c(g2,g3)))} will give you \code{list(g1,g2,g3)}
#' instead of \code{list(g1, c(g2,g3))}, as you might have expected.
#'
#' @param griddata A griddata structure to split
#' @return A list of griddata structures, each containing a single data set
#' (i.e., a time series read in from a single file)
#' @export
#' @keywords internal
splitGrids <- function(griddata)
{
    lapply(names(griddata$tags),
           function(tname) {
               tag <- griddata$tags[[tname]] # tag is a vector of (startrow, stoprow)
               modtag <- tag - (tag[1]-1)
               out <-
                   list(tas = griddata$tas[tag[1]:tag[2], ],
                        tgop = griddata$tgop,
                        lat = griddata$lat,
                        lon = griddata$lon,
                        time = griddata$time[tag[1]:tag[2]],
                        tags = stats::setNames(list(modtag), tname))
               class(out) <- 'griddata'
               out
           })
}


#' Flatten a list of tags, restoring their start/end rows
#'
#' @param taglist List of tags, possibly nested
#' @keywords internal
flatten_tags <- function(taglist)
{
    taglist <- rlang::flatten(taglist)
    lapply(taglist, function(x) {
               diff <- x[1]-1
               x - diff
           })
}

#' c operator for griddata objects
#'
#' Writing \code{c(g1,g2)} is equivalent to
#' \code{\link{concatGrids}(list(g1,g2))}.
#'
#' @param ... One or more griddata objects
#' @export
#' @keywords internal
c.griddata <- function(...)
{
    concatGrids(list(...))
}
