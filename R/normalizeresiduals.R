#' Normalize residuals
#'
#' A vector of residuals must be normalized for the fldgen algorithm to work
#' properly. This function takes a vector of residuals (T or P), calculates the
#' quantiles of each residual, and maps the residual to the corresponding
#' quantile value in a normal distribution.
#'
#' The output will be a list with Three fields:
#' \describe{
#'   \item{\strong{r}}{Column vector of the original, input residuals for a grid cell}
#'   \item{\strong{quants}}{Column vector of the quantiles of each input residual.}
#'   \item{\strong{rn}}{Column vectorof the new, normally distributed residuals.}
#' }
#'
#' UPDATE: The data at each time is represented as a flattened vector of grid cells.
#' The flattening is performed by transposing to lat, lon ordering, so that lat
#' will be the most rapidly varying index (because R uses
#' Fortran-style indexing).  Then the spatial dimensions are discarded,
#' resulting in a 1D vector.  The time dimension is kept, resulting in a matrix
#' with years in rows and grid cells in columns.  The dimensions of the matrix
#' will be Nyear x Ngrid, where Nyear is the number of years in the data set,
#' and Ngrid is the number of grid cells.
#'
#'
#' Conventionally, we refer to the output list as \code{griddata}.  Notably, any
#' other function with a \code{griddata} argument is expecting one of these
#' structures.
#'
#' @param r Column vector of the original, input residuals for a grid cell
#' @param len Maximum length of the time series to read.  If the data read is
#' longer, it will be trimmed.  (Default: read entire time series, regardless of
#' length.)
#' @export
normalize.resids <- function(r, len=NULL){

    ntime <- length(r)
    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            r <- r[1:ntime]
        }
    }

    # Get the quantiles of every residual for each grid cell
    # ecdf(r)(r) -> quantiles
    # quanfun <- apply(r, 2, ecdf)
    get.quans <- function(x){stats::ecdf(x)(x)}
    quantiles <- apply(r,2, get.quans)


    q2 <-   quantiles


    # get the values of the normal distribution corresponding to these quantiles.
    normresids <- apply(q2, 2, stats::qnorm)
    # normresids[is.infinite(normresids)] <- 10 # correct way to address Inf

    # return
    output <- list(r = r, quants = q2, rn = normresids)
    class(output) <- 'quantilemapping'

    return(output)

}


#' unNormalize residuals
#'
#' A vector of residuals must be normalized for the fldgen algorithm to work
#' properly. This function takes a vector of residuals (T or P), calculates the
#' quantiles of each residual, and maps the residual to the corresponding
#' quantile value in a normal distribution.
#'
#' The output will be a list with Three fields:
#' \describe{
#'   \item{\strong{r}}{Column vector of the original, input residuals for a grid cell - need because CDF of grid cell unknown}
#'   \item{\strong{rn}}{Column vector of the input, randomized, normally distributed residuals.}
#'   \item{\strong{rnew}}{Column vector of the quantiles of each input residual.}
#' }
#'
#' UPDATE: The data at each time is represented as a flattened vector of grid cells.
#' The flattening is performed by transposing to lat, lon ordering, so that lat
#' will be the most rapidly varying index (because R uses
#' Fortran-style indexing).  Then the spatial dimensions are discarded,
#' resulting in a 1D vector.  The time dimension is kept, resulting in a matrix
#' with years in rows and grid cells in columns.  The dimensions of the matrix
#' will be Nyear x Ngrid, where Nyear is the number of years in the data set,
#' and Ngrid is the number of grid cells.
#'
#'
#' Conventionally, we refer to the output list as \code{griddata}.  Notably, any
#' other function with a \code{griddata} argument is expecting one of these
#' structures.
#'
#' @param r Column vector of the original, input residuals for a grid cell
#' @param rn Column vector of the input, randomized, normally distributed residuals
#' @param len Maximum length of the time series to read.  If the data read is
#' longer, it will be trimmed.  (Default: read entire time series, regardless of
#' length.)
#' @export
unnormalize.resids <- function(r, rn, len=NULL){

    ntime <- length(r)
    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            r <- r[1:ntime]
        }
    }

    # Get the quantiles of every random, normal residual for each grid cell
    newquantiles <- apply(rn, 2, pnorm)


    # get the new quantile values of the original distribution
    vec.quantile <- function(x){as.vector(quantile(x, newquantiles, type = 4))}
    rnew <- apply(r, 2, vec.quantile)


    # return
    output <- list(r = r, rn = normresids, rnew = rnew)
    class(output) <- 'invquantilemapping'

    return(output)

}
