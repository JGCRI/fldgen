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
#' @param legendstr String to use for the title of the graph legend
#' @param palettestr String with the name of the ggplot2 color palette to use.
#' Defaults to 'RdYlBu'. Options are Diverging: 'BrBG', 'PiYG', 'PRGn', 'PuOr',
#' 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral'; Qualitative: 'Accent',
#' 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3';
#' Sequential: 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
#' 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn',
#' 'YlGnBu', 'YlOrBr', 'YlOrRd'.
#' @param palettedir Numerical value, only accepts 1 or -1, determining the
#' order of the color palette used. Defaults to -1.
#' @importFrom dplyr mutate %>% if_else
#' @importFrom ggplot2 aes element_blank
#' @export
plot_field <- function(fld, griddata, nb=6, minval=-3.5, maxval=3.5, legendstr="Temperature (K)", palettestr = 'RdYlBu', palettedir = -1 )
{
        # silence package checks
        lat <- lon <- long <- group <- value <- countries <- NULL


        # reshape data to be easier to plot
        tdf <- fld2df(fld, griddata)

        # # Get outlines of countries to add to map to make them easier to interpret.
        #
        # in order to plot correctly, the data frame of country outline
        # information needs a `value` column. Add one.
        fldgen::pckg_country_outlines_ggplot %>%
            mutate(value = 0) ->
            countries

        # Tweak the easter tip of Russia to get the map to break evenly across
        # longitudes.
        maxgroup <- max(countries$group)

        countries %>%
            mutate(group = if_else(long > 180, maxgroup + 1, group),
                   long = if_else(long > 180, long - 360, long)) ->
            countries


        ## TODO: make these options a little more customizable
        if(nb < 2) {
            ggplot2::ggplot(tdf, aes(x=lon, y = lat, fill = value)) +
                ggplot2::geom_raster() +
                ggplot2::scale_fill_distiller(palette=palettestr, direction=palettedir,
                                              limits=c(minval, maxval), oob=scales::squish,
                                              guide=ggplot2::guide_colorbar(title=legendstr, title.position='top')) +
                ggplot2::scale_x_continuous(expand=c(0,0)) +
                ggplot2::scale_y_continuous(expand=c(0,0)) +
                ggplot2::theme(axis.title=element_blank(),
                              axis.text=element_blank(),
                              axis.ticks=element_blank(),
                              legend.position = 'bottom') +
                # black outlines of countries
                ggplot2::geom_path(data = countries,
                                   aes(x = long, y = lat, group = group),
                                   size = 0.25, color = 'black')

            # # Legacy code using gcammaptools
            # gcammaptools::plot_GCAM_grid(tdf, col='value', extent=gcammaptools::EXTENT_WORLD,
            #                              legend=TRUE) +
            #     ggplot2::scale_fill_distiller(palette=palettestr, direction=palettedir,
            #                                   limits=c(minval, maxval), oob=scales::squish,
            #                                   guide=ggplot2::guide_colorbar(title=legendstr, title.position='top'))

        }
        else {
            ## Discretize the output values.
            tdf$value <- minval +
              findInterval(tdf$value, seq(minval, maxval, length.out=nb))/nb * (maxval-minval)

            ggplot2::ggplot(tdf, aes(x=lon, y = lat, fill = value)) +
                ggplot2::geom_raster() +
                ggplot2::scale_fill_distiller(palette=palettestr, direction=palettedir,
                                              limits=c(minval, maxval), oob=scales::squish,
                                              guide=ggplot2::guide_colorbar(title=legendstr, title.position='top')) +
                ggplot2::scale_x_continuous(expand=c(0,0)) +
                ggplot2::scale_y_continuous(expand=c(0,0)) +
                ggplot2::theme(axis.title=element_blank(),
                               axis.text=element_blank(),
                               axis.ticks=element_blank(),
                               legend.position = 'bottom') +
                # black outlines of countries
                ggplot2::geom_path(data = countries,
                                   aes(x = long, y = lat, group = group),
                                   size = 0.25, color = 'black')
        }

}


#' Plot a single year's temperature for every generated gridded time series
#' in a list.
#'
#' Takes in a list of temperature fields (residual or full data), a year
#' index, and the underlying temperature gridded data from a trained
#' emulator, as well as a few arguments relevant to plotting
#'
#' @param fieldlist List of residual
#' @param yearindex Index of the year to be plotted. For example, in data
#' ranging from 2006 to 2100, the yearindex for 2006 is 1.
#' @param emulatorgriddata The griddataT structure returned in a trained
#' emulator.
#' @param minval Lower limit of the color scale.  The default value was chosen
#' to work well for fields of residuals from the mean temperature response.
#' @param maxval Upper limit of the color scale.  The default value was chosen
#' to work well for fields of residuals from the mean temperature response.
#' @param legendstr String to use for the title of the graph legend
#' @param palettestr String with the name of the ggplot2 color palette to use.
#' Defaults to 'RdYlBu'. Options are Diverging: 'BrBG', 'PiYG', 'PRGn', 'PuOr',
#' 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral'; Qualitative: 'Accent',
#' 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3';
#' Sequential: 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
#' 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn',
#' 'YlGnBu', 'YlOrBr', 'YlOrRd'.
#' @param palettedir Numerical value, only accepts 1 or -1, determining the
#' order of the color palette used. Defaults to -1.
#'
#' @author ACS May 2020
#' @export
plotglobalfieldsT <- function(fieldlist, yearindex, emulatorgriddata,
                              minval=-3.5, maxval=3.5, legendstr,
                              palettestr = 'RdYlBu', palettedir = -1) {


    Ngrid <- ncol(emulatorgriddata[[1]])

    ## Extract a single example field from each series (yearindex determines which one) and create a plot
    lapply(fieldlist, function(g) {
        suppressWarnings(
            plot_field(g[yearindex, 1:Ngrid], emulatorgriddata, 14,  # 14 is the number of color levels in the plot
                       minval, maxval,
                       legendstr = legendstr, palettestr = palettestr, palettedir = palettedir) +
                ggplot2::guides(title=legendstr, title.position='top', title.hjust=0.5)
        )
    })
}



#' Plot a single year's precipitation for every generated gridded time series
#' in a list.
#'
#' Takes in a list of precipitation fields (residual or full data), a year
#' index, and the underlying precipitation gridded data from a trained
#' emulator, as well as a few arguments relevant to plotting
#'
#' @param fieldlist List of residual
#' @param yearindex Index of the year to be plotted. For example, in data
#' ranging from 2006 to 2100, the yearindex for 2006 is 1.
#' @param emulatorgriddata The griddataP structure returned in a trained
#' emulator.
#' @param minval Lower limit of the color scale.  The default value was chosen
#' to work well for fields of residuals from the mean precipitation response.
#' @param maxval Upper limit of the color scale.  The default value was chosen
#' to work well for fields of residuals from the mean precipitation response.
#' @param legendstr String to use for the title of the graph legend
#' @param palettestr String with the name of the ggplot2 color palette to use.
#' Defaults to 'RdYlBu'. Options are Diverging: 'BrBG', 'PiYG', 'PRGn', 'PuOr',
#' 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral'; Qualitative: 'Accent',
#' 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3';
#' Sequential: 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
#' 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn',
#' 'YlGnBu', 'YlOrBr', 'YlOrRd'.
#' @param palettedir Numerical value, only accepts 1 or -1, determining the
#' order of the color palette used. Defaults to -1.
#'
#' @author ACS May 2020
#' @export
plotglobalfieldsP <- function(fieldlist, yearindex, emulatorgriddata,
                              minval=-1.5, maxval=1.5, legendstr,
                              palettestr = 'BrBG', palettedir = 1) {

    Ngrid <- ncol(emulatorgriddata[[1]])

    ## Extract a single example field from each series (yearindex determines which one) and create a plot
    lapply(fieldlist, function(g) {
        suppressWarnings(
            plot_field(g[yearindex, (Ngrid + 1):(2*Ngrid)], emulatorgriddata, 14,  # 14 is the number of color levels in the plot
                       minval, maxval,
                       legendstr = legendstr, palettestr = palettestr, palettedir = palettedir) +
                ggplot2::guides(title=legendstr, title.position='top', title.hjust=0.5)
        )
    })
}



