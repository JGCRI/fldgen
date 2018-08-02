#' Normalize residuals
#'
#' A vector of residuals in a single grid cell must be normalized for the
#' fldgen algorithm to work properly. This function takes a matrix of
#' residuals (each grid cell is a column, T or P), calculates the quantiles
#' of each residual, and maps to the corresponding value in a normal
#' distribution.
#'
#' The output will be a list with Three fields:
#' \describe{
#'   \item{\strong{r}}{Matrix of the original, input residuals for each grid
#'   cell.}
#'   \item{\strong{quants}}{Matrix of the quantiles of each input residual.}
#'   \item{\strong{rn}}{Matrix of the new, normally distributed residuals.}
#' }
#'
#'
#'
#' Conventionally, we refer to the output list as \code{quantilemapping}.
#' Notably, any other function with a \code{quantilemapping} argument is
#' expecting one of these structures.
#'
#' @param r Matrix of the original, input residuals for each grid cell (columns)
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

    # deal with the fact that the max point gets set to the 1.00 quantile
    offset <- function(x){
        x[which.max(x)] <- max(x) - 1/nrow(r)
        return(x)}

    q2 <- apply(quantiles, 2, offset)

    # get the values of the normal distribution corresponding to these quantiles.
    normresids <- apply(q2, 2, stats::qnorm)
    #normresids[is.infinite(normresids)] <- 10 # address Inf

    # return
    output <- list(r = r, quants = q2, rn = normresids)
    class(output) <- 'quantilemapping'

    return(output)

}


#' unNormalize residuals
#'
#' A vector of normally distributed residuals in a single grid cell, likely
#' generated during the field generation step of this package, must be mapped
#' back to the original distribution of ESM residuals for the generated fields
#' to have any meaning. The original distribution of ESM residuals (T or P) is
#' empirically characterized by the vector of residuals in the grid cell.
#' fldgen algorithm to work properly. This function takes a matrix of normally
#' distributed residuals (each grid cell is a column, T or P), calculates the
#' quantiles of the normal distribution of each residual, and maps to the
#' corresponding value in the original, empirically characterized distribution
#' of residuals.
#'
#' The output will be a list with Three fields:
#' \describe{
#'   \item{\strong{r}}{Matrix of the original, input residuals for a grid cell,
#'   empirically characterizing the distribution.}
#'   \item{\strong{rn}}{Matrix of the input, randomized, normally distributed
#'   residuals.}
#'   \item{\strong{rnew}}{Matrix of the new residuals, sampled from the
#'   original distribution of residuals according to the quantiles of the
#'   input normally distributed random residuals.}
#' }
#'
#' Conventionally, we refer to the output list as \code{invquantilemapping}.
#' Notably, any other function with a \code{invquantilemapping} argument is
#' expecting one of these structures.
#'
#' @param r Matrix of the original, input residuals for a grid cell
#' @param rn Matrix of the input, randomized, normally distributed residuals
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
