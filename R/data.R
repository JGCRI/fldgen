#' pckg_country_outlines_ggplot
#'
#' A data frame containing the data from ggplot2::map_data('world'),
#' outlines of every country for use in provided convenience functions.
#'
#' @format A data frame with 99338 observations on the following 6 variables.
#' \describe{
#' \item{\code{long}}{a numeric vector}
#' \item{\code{lat}}{a numeric vector}
#' \item{\code{group}}{a numeric vector}
#' \item{\code{order}}{a numeric vector}
#' \item{\code{region}}{a character vector}
#' \item{\code{subregion}}{a character vector}
#' }
#' @details{
#' Including a call to map_data, however, threw test errors for github
#' actions in the PR.
#' Instead, we ran  the following offline to save this data:
#' pckg_country_outlines_ggplot <- ggplot2::map_data('world')
#' usethis::use_data(pckg_country_outlines_ggplot, overwrite = TRUE, compress = 'xz')
#' This results in pckg_country_outlines_ggplot being available
#' in the package.
#' }
'pckg_country_outlines_ggplot'
