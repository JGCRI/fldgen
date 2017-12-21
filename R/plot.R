#### functions for plotting fields

#' Reorganize a spatial field into a data frame
#'
#' The easiest way to plot this sort of data in ggplot is to organize it into a
#' frame with lat, lon, and value columns.  To do this you have to replicate the
#' latitude and longitude values as required.  These functions convert the
#' matrix representation of fields into a data frame with variables suitable for
#' plotting.  \code{fld2df} operates on a single time slice, while
#' \code{fldts2df} works on an entire time sequence.
#'
#' In addition to reformatting the data as a data frame, these functions remap
#' 0-360 longitude to -180, 180.
#'
#' @param fld Vector of ngrid=nlat*nlon values: a single time slice of the
#' field.
#' @param griddata The griddata structure returned from
#' \code{\link{read.temperatures}}.
#' @param ti The time index.  The time variable in the output will be set to
#' \code{griddata$time[ti]}.  If the time index is omitted, then the time
#' variable will be omitted from the output.
#' @return Data frame with variables lat, lon, value, and optionally t.
#' @export
fld2df <- function(fld, griddata, ti=NULL)
{
    nlat <- length(griddata$lat)
    nlon <- length(griddata$lon)
    assert_that(is.vector(fld))
    assert_that(length(fld) == nlat*nlon)

    ## Odds are
    if(all(griddata$lon >= 0) && any(griddata$lon > 180)) {
        loncorrect <- dplyr::if_else(griddata$lon > 180, griddata$lon - 360,
                                  griddata$lon)
    }
    else {
        loncorrect <- griddata$lon
    }

    ## The fields are in lat, lon order; i.e., lat varies most rapidly.
    lat <- rep(griddata$lat, times=nlon)
    lon <- as.vector(matrix(rep(loncorrect, times=nlat), nrow=nlat,
                            byrow=TRUE))
    if(is.null(ti)) {
        tibble::tibble(lat=lat, lon=lon, value=as.vector(fld))
    }
    else {
        tibble::tibble(lat=lat, lon=lon, t=griddata$time[ti], value=as.vector(fld))
    }
}


#' Reorganize a time sequence of spatial fields into a data frame
#'
#'
#' @rdname fld2df
#' @param fldts Matrix [ntime, ngrid] of spatial fields
#' @export
fldts2df <- function(fldts, griddata)
{
    ntime <- dim(fldts)[1]
    assert_that(ntime == length(griddata$time))

    tidxs <- 1:ntime
    dplyr::bind_rows(
        lapply(tidxs, function(i) {fld2df(fldts[i,], griddata, i)}) )
}



#' Plot a single field in matrix form
#'
#' Transform the field into a data frame using \code{\link{fld2df}} and plot
#' using the \code{gcammaptools} package.  If \code{gcammaptools} isn't
#' available, return \code{NULL}.  (TODO: make an ersatz plot if gcammaptools
#' isn't there)
#'
#' @param fld Vector of ngrid=nlat*nlon values: a single time slice of the
#' field.
#' @param griddata The griddata structure returned from
#' \code{\link{read.temperatures}}.
#' @param nb Number of breaks in the color scale.  If nb < 2, use a smooth
#' gradient.
#' @param minval Lower limit of the color scale.  The default value was chosen
#' to work well for fields of residuals from the mean temperature response.
#' @param maxval Upper limit of the color scale.  The default value was chosen
#' to work well for fields of residuals from the mean temperature response.
#' @export
plot_field <- function(fld, griddata, nb=6, minval=-3.5, maxval=3.5)
{
    if(requireNamespace('gcammaptools')) {
        tdf <- fld2df(fld, griddata)
        ## TODO: make these options a little more customizable
        if(nb < 2) {
            gcammaptools::plot_GCAM_grid(tdf, col='value', extent=gcammaptools::EXTENT_WORLD,
                                         legend=TRUE) +
              ggplot2::scale_fill_distiller(palette='RdYlBu', direction=-1,
                                            limits=c(minval, maxval), oob=scales::squish)
        }
        else {
            ## Discretize the output values.
            tdf$value <- minval +
              findInterval(tdf$value, seq(minval, maxval, length.out=nb))/nb * (maxval-minval)
            gcammaptools::plot_GCAM_grid(tdf, col='value', extent=gcammaptools::EXTENT_WORLD,
                                         legend=TRUE) +
              ggplot2::scale_fill_distiller(palette='RdYlBu', direction=-1, limits=c(minval,maxval))
        }
    }
    else {
        NULL
    }
}
