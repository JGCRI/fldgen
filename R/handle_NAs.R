#' Remove NA grid cells from processing in the emulator training process
#'
#' This function removes grid cells with missing data, NA values for all
#' timesteps, from being processed by \code{trainTP}. This enables \code{fldgen}
#' to work with output that has NAs in it.
#'
#' The parts of the griddata structure that are affected are:
#' \describe{
#' \item{vardata}{Missing grid cells are removed}
#' \item{globalop}{Entries corresponding to missing grid cells are removed.}
#' }
#'
#' The following items are added, \emph{only if} there were missing grid cells
#' to drop.
#' \describe{
#' \item{ncol_full}{Number of columns in the original matrix, before missing
#' values were dropped.}
#' \item{gridid_full}{Column ids in the original matrix, for the grid cells that
#' were not dropped.}
#' \item{coord}{Matrix of latitude and longitude for each grid cell. This is distinct from
#' the \emph{lat} and \emph{lon} entries, which gives the latitude coordinate
#' extracted from the original netCDF file.} #'
#' }
#'
#' @param griddata output from \code{concatGrids.general}
#' @return A new griddata structure, with missing grid cells removed
#' @export
#' @keywords internal
drop_NAs <- function(griddata){

    nmiss <- apply(griddata$vardata, 2, function(col) {sum(is.na(col))})
    if(all(nmiss == 0)) {
        ## IF there are no NA values then return the griddata structure
        ## unmodified.
        return(griddata)

    }

    ## Check for partially missing grid cells.  We require that cells are either
    ## always missing, or never.
    if(any(nmiss > 0 & nmiss != nrow(griddata$vardata))) {
        stop('Inconsistent grid cell entries, both NAs and numeric values for ', names(griddata$tags))
    }

    coord_full <- coord_array(griddata$lat, griddata$lon)

    ## Record the original number of columns
    griddata$ncol_full <- ncol(griddata$vardata)

    ## Find the grid cells with valid data
    valid_cells <- which(nmiss == 0)

    ## Drop missing cells from vardata and globalop
    griddata$vardata <- griddata$vardata[ , valid_cells]
    griddata$globalop <- griddata$globalop[valid_cells, ]
    ## Need to drop missing cells from the coordinate arrays too.
    coord <- coord_full[valid_cells , ]

    ## Record the other data we need
    griddata$gridid_full <- valid_cells # valid_cells has the original
                                        # column ids.
    griddata$coord <- coord

    griddata
}


#' Add the NA grid cells that may have been removed during the the emulator training process
#'
#' This function adds NA gridcells that were removed during the training process in order to convert
#' the results to same dimensions as the input grid.
#'
#' @param data Matrix of data to which to restore missing grid cells
#' @param griddata Griddata structure from the training procedure
#' @return Matrix of data with grid cells restored.
#' @export
#' @keywords internal
add_NAs <- function(data, griddata)
{
    if(is.null(griddata$ncol_full)) {
        ## Nothing was dropped from this dataset, so nothing to restore.
        return(data)
    }

    nrow <- nrow(data)
    ncol <- griddata$ncol_full

    output <- matrix(nrow=nrow, ncol=ncol)
    output[ , griddata$gridid_full] <- data

    output
}

#' Convert lat/lon coordinates to a coordinate array
#'
#' NetCDF assumes the grid is rectangular, so it stores latitude by row and
#' longitude by column.  We are potentially handling irregular coordinate grids
#' here, so we need to be able to covert to lat/lon for each grid cell.
#'
#' @param lat Latitude coordinate array for the full rectangular grid
#' @param lon Longitude coordinate array for the full rectangular grid
#' @return Matrix [ngrid x 2], latitude in the first column, longitude in the
#' second.
#' @export
coord_array <- function(lat, lon)
{
    nlat <- length(lat)
    nlon <- length(lon)

    ## Latitude is the most rapidly varying coordinate
    latgrid <- rep(lat, times=nlon)
    ## Longitude is trickier.  We'll build it as a matrix with the lon values
    ## running in rows.
    longrid <- as.vector(
        matrix(rep(lon, times=nlat), nrow=nlat, byrow=TRUE))

    coord <- matrix(c(latgrid, longrid), ncol=2)
    colnames(coord) <- c('lat','lon')
    coord
}
