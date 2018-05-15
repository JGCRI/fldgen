#### Functions for working with input and output data

#' Create an operator to extract grid cells from a field by latitude and
#' longitude.
#'
#' This function creates an operator for extracting a selection of grid cells
#' identified by a latitude and longitude bounding box.  The result is a logical
#' vector with \code{TRUE} values in the indices corresponding to grid cells
#' inside the bounding box.
#'
#' This operator can be used as an index to extract the desired grid cells, for
#' example:
#'
#' \code{mask <- extract_box; fld_subset <- fld[ ,mask]}.
#'
#' It can also
#' be used to calculate averages over a specified area:
#'
#' \code{meantemp <- fld \%*\% as.numeric(mask)}.
#'
#' Or, you could use it to mask out grid cells that are not in the region of
#' interest:
#'
#' \code{fld[ ,!mask] <- NA}
#'
#' @param lat A vector of length 2 giving the minimum and maximum latitude for
#' the bounding box.
#' @param lon A vector of length 2 giving the minimum and maximum longitude for
#' the bounding box.
#' @param griddata A \code{griddata} structure, which contains the latitude and
#' longitude coordinates.
#' @return Logical vector indicating which grid cells are in the bounding box.
#' @export
extract_box <- function(lat, lon, griddata)
{
    nlat <- length(griddata$lat)
    nlon <- length(griddata$lon)

    ## Latitude is the most rapidly varying coordinate.
    latvec <- rep(griddata$lat, times=nlon)
    ## Longitude is the least rapidly varying.  Use the byrow argument to the
    ## matrix constructor to get values in the right order.
    lonvec <- as.vector(matrix(rep(griddata$lon, times=nlat),
                               ncol=nlon, byrow=TRUE))

    ## the operator is just the elementwise & of all the conditions
    latvec >= lat[1] & latvec <= lat[2] & lonvec >= lon[1] & lonvec <= lon[2]
}
