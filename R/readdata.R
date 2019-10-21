#' Read ESM netCDF data for arbitrary variables
#'
#' The variable data should be a netCDF file of annual grid values.  The data
#' is assumed to have its dimensions of time, lat, lon, where lon is the most
#' rapidly varying index, and time is the least.
#'
#' The output will be a list with six fields:
#' \describe{
#'   \item{\strong{vardata}}{Matrix of variable data.}
#'   \item{\strong{globalop}}{Vector operator for global mean value of the variable.}
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
#' The \code{globalop} vector is a vector whose dot product with a flattened grid is
#' the area-weighted global mean for the grid.  It is stored as a column vector,
#' so \code{vardata \%*\% globalop} is the time series of global mean values.
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
#' @param tag A string identifying the name of the scenario.  If omitted, the
#' tag will default to the filename.
#' @param latvar Name of the latitude dimension variable in the input file.
#' @param lonvar Name of the longitude dimension variable in the input file.
#' @param timevar Name of the time dimension variable in the input file.
#' @return A \code{griddata} list (see details).
#' @importFrom assertthat assert_that
#' @export
read.general <- function(filename, len=NULL, tag=basename(filename), varname='tas',
                              latvar='lat', lonvar='lon', timevar='time')
{
    vann <- ncdf4::nc_open(filename)

    ## v3d should have dimensions of time x lat x lon in the netcdf file.
    ## Because R uses Fortran array layout, this gets reversed to lon x lat x time.
    v3d <- ncdf4::ncvar_get(vann, var=varname)
    lat <- ncdf4::ncvar_get(vann, var=latvar)
    nlat <- length(lat)
    lon <- ncdf4::ncvar_get(vann, var=lonvar)
    nlon <- length(lon)
    time <- ncdf4::ncvar_get(vann, var=timevar)
    ntime <- length(time)
    timeout <- seq_along(time) - 1
    ncdf4::nc_close(vann)

    assert_that(all(dim(v3d) == c(nlon, nlat, ntime)))

    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            time <- time[1:ntime]
            timeout <- timeout[1:ntime]
            v3d <- v3d[ , , 1:ntime]
        }
    }


    ## reorganize and flatten the 3-D array into a 2d array of ntime x ngrid
    ## As we do this, we will also make latitude the most rapidly varying index
    ## for the individual time slices.

    variable <- aperm(v3d, c(3,2,1))
    dim(variable) <- c(ntime, nlat*nlon)

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

    gd <- list(vardata=as.matrix(variable), globalop=as.matrix(areafac), lat=lat, lon=lon,
               time=timeout, tags=tagout)
    class(gd) <- 'griddata'
    gd
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
concatGrids.general <- function(gridlist)
{
    vardatal <- lapply(gridlist, function(x){x$vardata})
    globalopl <- lapply(gridlist, function(x){x$globalop})
    latl <- lapply(gridlist, function(x){x$lat})
    lonl <- lapply(gridlist, function(x){x$lon})
    timel <- lapply(gridlist, function(x){x$time})
    tagsl <- lapply(gridlist, function(x){x$tags})

    ## verify that globalop, lat, and lon are compatible.
    checkfun <- function(cmpval) {
        function(x) {isTRUE(all.equal(x, cmpval))}
    }
    all(sapply(globalopl, checkfun(globalopl[[1]]))) ||
        stop('globalop values not compatible')
    all(sapply(latl, checkfun(latl[[1]]))) || stop('lat values not compatible')
    all(sapply(lonl, checkfun(lonl[[1]]))) || stop('lon values not compatible')

    ## Concatenate tas matrix by rows
    vardata <- do.call(rbind, vardatal)
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
    gd <- list(vardata=vardata, globalop=globalopl[[1]], lat=latl[[1]], lon=lonl[[1]], time=time, tags=tags)
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
splitGrids.general <- function(griddata)
{
    lapply(names(griddata$tags),
           function(tname) {
               tag <- griddata$tags[[tname]] # tag is a vector of (startrow, stoprow)
               modtag <- tag - (tag[1]-1)
               out <-
                   list(vardata = griddata$vardata[tag[1]:tag[2], ],
                        globalop = griddata$globalop,
                        lat = griddata$lat,
                        lon = griddata$lon,
                        time = griddata$time[tag[1]:tag[2]],
                        tags = stats::setNames(list(modtag), tname))
               class(out) <- 'griddata'
               out
           })
}


#' Pair temperature and precipitation (or any other variable) netCDF files
#' based on file name, following CMIP5 convention.
#'
#' This function reads in a directory or list of netCDF files, with names
#' following CMIP5 scenario design conventions, and pairs the temperature
#' and precipitation files based on scenario. Technically any two variables
#' could have files paired using this function, just by changing the tvarname
#' and pvarname inputs.
#'
#' @param dat A single directory name, or a list of netCDF files.  If a
#' directory name is given, all netCDF files in the directory will be used.The
#' pairing of temperature and precipitation netCDF files in the directory
#' relies on the CMIP5 file naming conventions. Other naming conventions are
#' not currently supported.
#' @param tvarname The name of the temperature variable in CMIP5 ('tas' default).
#' @param varname The name of the precipitation variable in CMIP5 ('pr' default).
#' @return A tibble of paired file names in two columns: tfilename, pfilename.
#' Each row of the tibble is a single scenario.
#' @export
#' @keywords internal
file.pairer <- function(dat, tvarname = 'tas', pvarname = 'pr')
{
    # silence package checks
    value <- var <- timestep <- time <- startyr <- stopyr <-
        tfilename <- pfilename <- NULL

    # separate dat into list of precip files and temperature files. Relies on
    # CMIP5 naming conventions.
    pdat <- dat[grep(paste0(pvarname, "_"), dat)]
    if(any(grepl("aclim", tolower(pdat) == FALSE) &
           grepl("annual", tolower(pdat)) == FALSE)) {
      stop(paste("At least one precipitation file in", dat, "is not annual"))
    }
    pdirs <- dirname(pdat)
    pdat <- basename(pdat)

    tdat <- dat[grep(paste0(tvarname, "_"), dat)]
    if(any(grepl("aclim", tolower(tdat)) == FALSE &
           grepl("annual", tolower(tdat)) == FALSE)) {
      stop(paste("At least one temperature file in", dat, "is not annual"))
    }
    tdirs <- dirname(tdat)
    tdat <- basename(tdat)

    ## This pattern extracts the metadata from the filenames of the input files.  In order,
    ## the fields are variable, realm/resolution, model, scenario, runid, time range.
    ## The model name can have '-' characters, and the time range has the pattern NNNN-NNNN
    ## (where the number of digits is arbitrary).  All the others have to be alphanumeric.
    ptrn <- "([[:alnum:]]+)_([[:alnum:]]+)_([-[:alnum:]]+)_([[:alnum:]]+)_?(.*)_([0-9]+-[0-9]+)\\.nc"
    pmeta <- stringr::str_match(pdat, ptrn)
    tmeta <- stringr::str_match(tdat, ptrn)
    if(any(is.na(pmeta[,1])) || any(is.na(tmeta[,1]))) {
      ## One or more of the files doesn't conform to the naming convention, throw an error
      bad <- function(m,d) {d[is.na(m[,1])]}
      badp <- bad(pmeta, pdat)
      badt <- bad(tmeta, tdat)
      flist <- paste(paste(badp, collapse='\n'), paste(badt, collapse='\n'), sep='\n')
      stop('The following files had malformed filenames:\n', flist)
    }
    fields <- c('pfilename','var','realmres', 'model','scenario','runid','time')
    ptbl <- structure(as.data.frame(pmeta, stringsAsFactors=FALSE), names=fields)
    ptbl$pfilename <- file.path(pdirs, ptbl$pfilename)            # restore full path
    fields[1] <- 'tfilename'
    ttbl <- structure(as.data.frame(tmeta, stringsAsFactors=FALSE), names=fields)
    ttbl$tfilename <- file.path(tdirs, ttbl$tfilename)
    paireddat <- dplyr::full_join(ttbl, ptbl, by=fields[c(-1,-2)])

    ## Check for failed matches
    fail <- is.na(paireddat$tfilename) | is.na(paireddat$pfilename)
    if(all(fail)) {
      stop('None of the input files had matches.\ntfilenames:\n',
           paste(tmeta, collapse='\n'), '\npfilenames:\n', paste(pmeta, collapse='\n'))
    }
    if(any(fail)) {
      paireddat <- paireddat[!fail,]
      ngood <- nrow(paireddat)
      warning('Some temperature files had no match in precipitation, or vice versa.  ',
              'Proceeding with matched data (', ngood, ' pairs.)')
    }

    return(paireddat[c('tfilename','pfilename')])

}


#' Read in global average temperature time series
#'
#' Use the string defined by the globalAvg_file \code{trainTP} argument to
#' import a global average temperature time series from a corresponding txt
#' file.
#'
#' @param nc_files the path and name to the nc file the global average is
#' paired with.
#' @param globalAvg_file The string end added onto the nc_file name for the txt
#' file to read in.
#' @param vardata The data frame to use to check the dimensions of the data
#' imported by the function.
#' @param paireddat The paireddat data frame, this is used to make sure that
#' the number of rows returned matches what is expected.
#' @return A matrix of global averages, the number of rows should match the
#' number of rows in the vardata data frame.
#' @importFrom utils read.table
#' @export
#' @keywords internal
read_globalAvg <- function(nc_files, globalAvg_file, vardata, paireddat){

    # silence package checks
    X <- NULL


    # TODO currently, because of the paireddat element, this is not back-
    # compatible with the T-only training. When creating a single, streamlined
    # train function that works on T or TP, update. Based on plan for that
    # (only support the input of a list of file names, one var or paired),
    # will likely just come down to renaming paireddat here and in the new
    # master train function.


    # TODO I think that theres needs to be some changes to this function so
    # that it takes in a global average data frame and checks to make sure
    # that the time is compatible with the the training files! other wise
    # there could be some serious problems.
    tgav_list <- lapply(nc_files, function(nc_file = X){

        # Check to make sure that only one nc file is read in at a time
        if(length(nc_file) > 1) {
            stop('more than one file read into read_globalAvg')
            }


        # Read in the global average from the txt file.
        global_file <- paste0(nc_file, globalAvg_file)


        # Throw an error if the file is missing
        if(!file.exists(global_file)) {
            stop('could not find: ', basename(global_file))
            }


        # Import the txt file
        tgav <- as.matrix(read.table(file = global_file))


        # Check the imported tgav dimensions
        expected_rows <- nrow(vardata) / nrow(paireddat)
        if(nrow(tgav) != expected_rows) {
            stop(nrow(tgav), ' rows in ', basename(global_file), ' when ', expected_rows, ' rows expected.')
            }

        tgav
        })


    # Format into a matrix
    matrix(data = unlist(tgav_list), ncol = 1)
}



###########################################################################################


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


#' Read ESM precipitation data
#'
#' The precipitation data should be a netCDF file of average annual grid values.
#' The data is assumed to have its dimensions of time, lat, lon, where lon is
#' the most rapidly varying index, and time is the least.
#'
#' The output will be a list with six fields:
#' \describe{
#'   \item{\strong{pr}}{Matrix of precipitation data.}
#'   \item{\strong{pgop}}{Vector operator for global mean precipitation.}
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
#' The \code{pgop} vector is a vector whose dot product with a flattened grid is
#' the area-weighted global mean for the grid.  It is stored as a column vector,
#' so \code{pr \%*\% pgop} is the time series of global mean precipitations.
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
read.precipitations <- function(filename, len=NULL, tag=basename(filename), varname='pr',
                              latvar='lat', lonvar='lon', timevar='time')
{
    pann <- ncdf4::nc_open(filename)

    ## pr3d should have dimensions of time x lat x lon in the netcdf file.
    ## Because R uses Fortran array layout, this gets reversed to lon x lat x time.
    pr3d <- ncdf4::ncvar_get(pann, var=varname)
    lat <- ncdf4::ncvar_get(pann, var=latvar)
    nlat <- length(lat)
    lon <- ncdf4::ncvar_get(pann, var=lonvar)
    nlon <- length(lon)
    time <- ncdf4::ncvar_get(pann, var=timevar)
    ntime <- length(time)
    timeout <- seq_along(time) - 1
    ncdf4::nc_close(pann)

    assert_that(all(dim(pr3d) == c(nlon, nlat, ntime)))

    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            time <- time[1:ntime]
            timeout <- timeout[1:ntime]
            pr3d <- pr3d[ , , 1:ntime]
        }
    }


    ## reorganize and flatten the 3-D array into a 2d array of ntime x ngrid
    ## As we do this, we will also make latitude the most rapidly varying index
    ## for the individual time slices.

    pr <- aperm(pr3d, c(3,2,1))
    dim(pr) <- c(ntime, nlat*nlon)

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

    gd <- list(pr=as.matrix(pr), pgop=as.matrix(areafac), lat=lat, lon=lon,
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
        modeldata <- NULL               # silence check notes.
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


