#' characterize empirical distributions
#'
#' Given a vector of residuals in a single grid cell that empirically
#' characterizes the distribution of residuals, return a CDF and
#' quantile function for that grid cell.
#'
#' The output will be a list with two fields:
#' \describe{
#'   \item{\strong{cdf}}{List of empirical CDF functions for each grid cell}
#'   \item{\strong{quant}}{List of empirical quantile functions for each grid
#'   cell.}
#' }
#'
#'
#'
#' Conventionally, we refer to the output list as \code{empiricaldistribution}.
#' Notably, any other function with a \code{empiricaldistribution} argument is
#' expecting one of these structures.
#'
#' @param inputresids Matrix of the original, input residuals for each grid
#'  cell (columns)
#' @param len Maximum length of the time series to read.  If the data read is
#' longer, it will be trimmed.  (Default: read entire time series, regardless of
#' length.)
#' @export
characterize.emp.dist <- function(inputresids, len=NULL){

    ntime <- length(inputresids)
    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            inputresids <- inputresids[1:ntime]
        }
    }

    # add a dummy point to the input residuals so that the largest residual
    # observed in each gridcell is now the penultimate, and won't get assigned
    # p = 1 from ecdf
    maxresid <- 2 * max(apply(inputresids, 2, max))
    newrow <- t(matrix(rep(maxresid, ncol(inputresids))))

    inputresids2 <- rbind(inputresids, newrow)

    # get the empirical cdf function for each grid cell (column of inputresids)
    empCDFs <- apply(inputresids2, 2, ecdf)


    # define a function that can create the inverse of an empirical cdf
    # function. Such an inverse is the quantile function.
    # This exploits the fact that an object of type ecdf is a step function.
    # This inverse is created with a linear interpolation call via `approxfun`.
    # This corresponds to an argument of `type = 4` in the typical call to
    # `quantile`.
    #
    # The stepfunction would technically be the true inverse of the empirical
    # cdfs returned by `ecdf`. However, poor performance at edge cases leads
    # to unacceptably large errors. The interpolated quantile function is a
    # more accurate inverse to the stepfunction empirical cdfs.
    #
    # Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages,
    # American Statistician 50, 361–365. doi: 10.2307/2684934.
    inv.ecdf <- function(f) {
        # need to address value < min(value) has probability 0.
        # smallest p empirical cdf returned from ecdf can handle is typically
        # 0.0025ish
        # calculate the last line segment and use the slope to extrapolate
        # back to p = 0, then add that point to be fitted
        minval <- min(knots(f))
        minp <- f(minval)

        minval2 <- sort(knots(f),partial=2)[2]
        minp2 <- f(minval2)

        slope <- 0.01 * (minval - minval2) / (minp - minp2)
        newminval <- minval - slope * minp


        # define the new points and interpolate
        newy <- c(newminval, knots(f)) # values from the distributions
        newx <- c(0, f(knots(f))) # probabilities


        approxfun(x = newx, y = newy)


    }

    ##### QUESTION ######
    # if we're going ahead and interpolating to create our quantile,
    # why not just interpolate for our empirical cdf too?
    # harder to justify why not using built in R functions where we can
    # but also more consistent. Could also be a cleaner place to address the
    # 0 and 1 quantile issue
    #
    #
    # define a function that can create an interpolated empirical cdf rather than
    # a step function
    interp.ecdf <- function(f) {
        minval <- min(knots(f))
        minp <- f(minval)

        minval2 <- sort(knots(f),partial=2)[2]
        minp2 <- f(minval2)

        slope <- 0.01 * (minval - minval2) / (minp - minp2)
        newminval <- minval - slope * minp


        # define the new points and interpolate
        newx <- c(newminval, knots(f)) # values from the distributions
        newy <- c(0, f(knots(f))) # probabilities

        approxfun(x = newx, y = newy)
    }

    #######################


    # apply this function to the empirical cdfs to get empirical quantile
    # functions for every grid cell.
    empQuants <- lapply(empCDFs, inv.ecdf)

    empCDFs2 <- lapply(empCDFs, interp.ecdf)


    # return
    output <- list(cdf = empCDFs2, quant = empQuants)
    class(output) <- 'empiricaldistribution'

    return(output)
}


#' Normalize residuals
#'
#' A vector of residuals in a single grid cell must be normalized for the
#' fldgen algorithm to work properly. This function takes a matrix of
#' residuals (each grid cell is a column, T or P), calculates the quantiles
#' of each residual, and maps to the corresponding value in a normal
#' distribution.
#'
#' The output will be a list with four fields:
#' \describe{
#'   \item{\strong{inputresids}}{Matrix of the original, input residuals for
#'   each grid cell.}
#'   \item{\strong{empiricalcdf}}{List of the empirical cdf functions for each
#'   grid cell.}
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
#' @param inputresids Matrix of the original, input residuals for each grid
#'  cell (columns)
#' @param empiricalcdf List of the empirical cdf functions for each grid cell.
#' @param len Maximum length of the time series to read.  If the data read is
#' longer, it will be trimmed.  (Default: read entire time series, regardless of
#' length.)
#' @export
normalize.resids <- function(inputresids, empiricalcdf, len=NULL){

    ntime <- length(inputresids)
    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            inputresids <- inputresids[1:ntime]
        }
    }


    ## assert that inputresids has same number of columns as entries in
    ##empiricalcdf



    # Get the quantiles of every residual for each grid cell.
    get.quans <- function(index){
        empiricalcdf[[index]](inputresids[,index])
        }


    q2 <- sapply(1:ncol(inputresids), get.quans)



    # get the values of the normal distribution corresponding to these quantiles.
    normresids <- apply(q2, 2, stats::qnorm)

    # return
    output <- list(inputresids = inputresids, empiricalcdf = empiricalcdf,
                   quants = q2, rn = normresids)
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
#' The output will be a list with four fields:
#' \describe{
#'   \item{\strong{empiricalquant}}{List of the empirical quantile functions
#'   for each grid cell.}
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
#' @param empiricalquant List of the empirical quantile functions for each grid
#' cell.
#' @param rn Matrix of the input, randomized, normally distributed residuals
#' @param len Maximum length of the time series to read.  If the data read is
#' longer, it will be trimmed.  (Default: read entire time series, regardless of
#' length.)
#' @export
unnormalize.resids <- function(empiricalquant, rn, len=NULL){

    ntime <- length(rn)
    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            rn <- rn[1:ntime]
        }
    }

    # Get the quantiles of every random, normal residual for each grid cell
    newquantiles <- apply(rn, 2, stats::pnorm)


    # get the values of the original distribution corresponding to the new
    # quantile values, using the input empiricalquant list of functions.
    vec.quantile <- function(index){
        empiricalquant[[index]](newquantiles[,index])

    }

    rnew <- sapply(1:length(empiricalquant), vec.quantile)


    # return
    output <- list(empiricalquant = empiricalquant, rn = rn, rnew = rnew)
    class(output) <- 'invquantilemapping'

    return(output)

}
