#' Characterize empirical distributions
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
#' @importFrom stats approxfun
#' @export
characterize.emp.dist <- function(inputresids, len=NULL){


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



    # we're going to be adding two dummy points at p = 0, p = 1
    # we have ntime observations between.
    # So we're going from 0 to 1 inclusive in N+1 steps (for N+2 points).
    ntime <- nrow(inputresids)
    offset <- 1/(ntime + 1)


    # helper functions
    maxN <- function(x, N=2){
        len <- length(x)
        if(N>len){
            warning('N greater than length(x).  Setting N=length(x)')
            N <- length(x)
        }
        sort(x,partial=len-N+1)[len-N+1]
    }

    minN <- function(x, N=2){
        len <- length(x)
        if(N>len){
            warning('N greater than length(x).  Setting N=length(x)')
            N <- length(x)
        }
        sort(x,partial=N)[N]
    }


    # function to create the new max for a column of residuals according
    # to 8-20-18 notes/documentation
    newmax <- function(residual.column){
        x.max1 <- max(residual.column)         # largest
        x.max2 <- maxN(residual.column, N = 2) # 2nd largest
        # x.max3 <- maxN(residual.column, N = 3) # 3rd largest

        # corresponding probabilities
        p.max0 <- 1
        p.max1 <- (offset-1) * offset
        p.max2 <- (offset-2) * offset
        # p.max3 <- (offset-3) * offset

        # slope between two segments
        slope1 = (p.max1 - p.max2) / (x.max1 - x.max2)
        # slope2 = (p.max2 - p.max3) / (x.max2 - x.max3)

        # want to use k*slope1 to extrapolate to P=1 and our new x.max value
        # want to capture the saturating nature of the tails of the CDF
        # between segment 1 and 2, the slope has decreased by k = slope1/slope2
        # use this k for our new slope.
        k <- 100 #slope1 / slope2
        slope0 <- k * slope1 # less steep than slope1

        # calculate the new point
        x.max0 <- (1 / slope0) * (p.max0 - p.max1) + x.max1
        return(x.max0)

    }


    newmin <- function(residual.column){
        x.min1 <- min(residual.column)         # smallest
        x.min2 <- minN(residual.column, N = 2) # 2nd smallest
        #x.min3 <- minN(residual.column, N = 3) # 3rd smallest

        # corresponding probabilities
        p.min0 <- 0
        p.min1 <- 1 * offset
        p.min2 <- 2 * offset
       # p.min3 <- 3 * offset

        # slope between two segments
        slope1 = (p.min1 - p.min2) / (x.min1 - x.min2)
        #slope2 = (p.min2 - p.min3) / (x.min2 - x.min3)

        # want to use k*slope1 to extrapolate to P=1 and our new x.min value.
        # want to capture the saturating nature of the tails of the CDF
        # between segment 1 and 2, the slope has decreased by k = slope1/slope2
        # use this k for our new slope.
        k <- 100 # slope1 / slope2
        slope0 <- k * slope1 # less steep than slope1

        # calculate the new point
        x.min0 <- (1 / slope0) * (p.min0 - p.min1) + x.min1
        return(x.min0)

    }


    # add a dummy point to the input residuals so that the largest residual
    # observed in each gridcell is now the penultimate, and won't get assigned
    # p = 1 from ecdf
    maxresid <-  (1 + 0.0001) * apply(inputresids, 2, max)
    newrow1 <- t(matrix(maxresid))

    # dummy min for p=0
    minresid <-  apply(inputresids, 2, newmin)
    newrow2 <- t(matrix(minresid))


    # add the new rows
    inputresids2 <- rbind(inputresids, newrow1, newrow2)



    # define a function that can create an interpolated empirical cdf rather than
    # a step function
    interp.ecdf <- function(column.resids) {

        # for one column of the residuals, get the order of points.
        # runs 1...Ntime+2 because we have added two dummy rows.
        # needs to run 0...Ntime+1 for the probabilities
        order.vector <- rank(column.resids)

        # cumulative probabilities are order-1 / Ntime+1
        probabilities <- (order.vector - 1) * offset

        # return the ecdf for the function
        approxfun(x = column.resids, y = probabilities, rule=2)
    }

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
    inv.ecdf <- function(column.resids) {

        # for one column of the residuals, get the order of points.
        order.vector <- rank(column.resids)

        # cumulative probabilities are order-1 / Ntime+1
        probabilities <- (order.vector - 1) * offset

        # return the ecdf for the function
        approxfun(x = probabilities, y = column.resids, rule=2)

    }



    #######################


    # apply this function to each column of residuals to get empirical quantile
    # functions for every grid cell.
    empQuants <- apply(inputresids2, 2, inv.ecdf)

    empCDFs <-apply(inputresids2, 2, interp.ecdf)


    # return
    output <- list(cdf = empCDFs, quant = empQuants)
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
