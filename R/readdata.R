#' Read ESM data
#'
#' Read temperature data from the ESM that we are going to emulate.
#'
#' The temperature data should be a netCDF file of annual grid values.  The data
#' is assumed to have its dimensions of time, lat, lon, where lon is the most
#' rapidly varying index, and time is the least.
#'
#' The output will be a list with two fields:
#' \describe{
#'   \item{\strong{tas}}{Matrix of temperature data.}
#'   \item{\strong{tgop}}{Vector operator for global mean temperature.}
#' }
#'
#' The data at each time is represented as a flattened vector of grid cells.
#' The flattening is performed by transposing to lon, lat ordering, so that lon
#' will once again be the most rapidly varying index (because R uses
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
#' @return A \code{griddata} list (see details).
#' @importFrom assertthat assert_that
#' @export
read.ncdf <- function(filename)
{
    tann <- ncdf4::nc_open(filename)

    ## tas3d should have dimensions of time x lat x lon in the netcdf file.
    ## Because R uses Fortran array layout, this gets reversed to lon x lat x time.
    tas3d <- ncdf4::ncvar_get(tann, var='tas')
    lat <- ncdf4::ncvar_get(tann, var='lat')
    nlat <- length(lat)
    lon <- ncdf4::ncvar_get(tann, var='lon')
    nlon <- length(lon)
    time <- ncdf4::ncvar_get(tann, var='time')
    ntime <- length(time)
    timeout <- seq_along(time) - 1
    ncdf4::nc_close(tann)

    assert_that(all(dim(tas3d) == c(nlon, nlat, ntime)))


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

    list(tas=as.matrix(tas), tgop=as.matrix(areafac), lat=lat, lon=lon, time=timeout)
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


