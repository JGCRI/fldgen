#' Using a trained emulator, generate new full fields.
#'
#' This function takes in a trained emulator - a structure of class
#' \code{fldgen}.
#' This structure contains everything the emulator has learned about the model,
#' and is used to generate new fields of residuals. Also taking in a list of
#' generated residual fields and a global average yield, a global gridded
#' mean field is constructed accoring to the input reconstruction_function.
#' The mean field and residual fields are added to return a list of different
#' realizations of full fields.
#'
#'
#' @param emulator A trained \code{fldgen} temperature precipitation joint
#' emulator.
#' @param residgrids A list of new residual fields, each entry in the list is a
#' new realization, a matrix that is [Nyears x 2 * Ngrid]; the first 1:Ngrid
#' cols are the temperature residuals and columns (Ngrid + 1):(2*Ngrid) are the
#' precipitation residuals.
#' @param tgav  A vector (or N x 1 matrix) of global mean temperatures, used to
#' calculate the mean warming response field.
#' @param tvarunconvert_fcn The function to undo any transformation done to the
#' input training data in \code{trainTP} to correct the support. This should be
#' the inverse function of the tvarconvert_fcn argument to \code{trainTP}. This
#' is stored in a \code{trainTP} \code{emulator} under
#' \code{emulator$griddattaT$tvarconvert_fcn}.
#' Defaults to NULL as temperature typically doesn't need to be transformed to
#' correct the support. WARNING: currently rely on the user to define the
#' correct inverse function here, though we do some checks and throw warnings
#' if it looks like there may be an issue.
#' @param pvarunconvert_fcn The function to undo any transformation done to the
#' input training data in \code{trainTP} to correct the support. This should be
#' the inverse function of the pvarconvert_fcn argument to \code{trainTP}. This
#' is stored in a \code{trainTP} \code{emulator} under
#' \code{emulator$griddattaP$pvarconvert_fcn}.
#' Defaults to exp() as precipitation is usually log-transformed in ordere to
#' correct the support. WARNING: currently rely on the user to define the
#' correct inverse function here, though we do some checks and throw warnings
#' if it looks like there may be an issue.
#' @param reconstruction_function A function for constructing a mean field from
#' trained pattern scaling result + a vector of global annual mean temperatures.
#' @return A list of:
#' 1) fullgrids = A list of new full fields, each entry in
#' the list is a new realization, a matrix that is [Nyears x 2 * Ngrid]; the
#' first 1:Ngrid cols are the temperature field and columns
#' (Ngrid + 1):(2*Ngrid) are the precipitation field.
#' 2) meanfieldT = the reconstructed, pattern scaled temperature mean field.
#' 3) meanfieldP = the reconstructed, pattern scaled precipitation mean field.
#' @export

generate.TP.fullgrids <- function(emulator, residgrids, tgav,
                                  tvarunconvert_fcn = NULL, pvarunconvert_fcn = exp,
                                  reconstruction_function = pscl_apply){

    # dbl check that convert and unconvert functions are actual inverses on some test data.
    # throw a warning if this is not the case but it is up to the user to track down.

    if(!is.null(emulator$griddataT$tvarconvert_fcn)) {
        if(is.null(tvarunconvert_fcn)) {
            stop('Temperature transform function was applied, but no inverse function supplied.')
        }

        if(!chkinverse(emulator$griddataT$tvarconvert_fcn, tvarunconvert_fcn,
                       range=c(250,310))) {
            warning('Your functions for transforming and inverse transforming the support of T may not be inverses of each other 1.')
        }
    }

    if(!is.null(emulator$griddataP$pvarconvert_fcn)) {
        if(is.null(pvarunconvert_fcn)) {
            stop('Precipitation transform function was applied, but no inverse function supplied.')
        }

        if(!chkinverse(emulator$griddataP$pvarconvert_fcn, pvarunconvert_fcn,
                       range=c(1e-9, 1e-5), logspace=TRUE)) {
            warning('Your functions for transforming and inverse transforming the support of P may not be inverses of each other 1.')
        }
    }


    ## Convert the tgav input to a column vector if necessary
    assertthat::assert_that(is.numeric(tgav))
    if(!is.matrix(tgav) || ncol(tgav) != 1) {
        tgav <- matrix(tgav, ncol=1)
    }

    # save a copy of the number of grids cells for each variable
    Ngrid <- ncol(residgrids[[1]])/2

    meanfieldT <- reconstruction_function(emulator$meanfldT, tgav)
    meanfieldP <- reconstruction_function(emulator$meanfldP, tgav)


    # Add the meanfield values to the residual grids
    lapply(residgrids, function(matrix, gridcells = Ngrid){

        # Separate the tas and pr data from one another.
        tas <- matrix[ , 1:Ngrid]
        pr  <- matrix[ , (Ngrid + 1):(2 * Ngrid)]

        # Add the meanfield to the data
        tas[ , 1:Ngrid] <- tas[ , 1:Ngrid] + meanfieldT
        pr[ , 1:Ngrid]  <- pr[ , 1:Ngrid] + meanfieldP



        # convert from (-inf, inf) support to natural support.
        if( !is.null(emulator$griddataT$tvarconvert_fcn)){

            tas <- tvarunconvert_fcn(tas)

        } else{
            message('Generated T full fields not being transformed to a different support. Up to user to know if this is desirable.')
        }


        if(!is.null(emulator$griddataP$pvarconvert_fcn)){

            pr <- pvarunconvert_fcn(pr)

        } else{
            message('Generated P full fields not being transformed to a different support. Up to user to know if this is desirable.')
        }



        # Return output
        return(list(tas = tas, pr = pr))


    }) ->
        fullgrids


    return(list(fullgrids = fullgrids, meanfieldT = meanfieldT, meanfieldP =
                  meanfieldP))

}


#' Check if two functions are inverses of one another over a given range of inputs.
#'
#' This is a quick and dirty check; it will definitely fail if the input range
#' is not in the domain of both functions.
#'
#' @param f1 First function to test
#' @param f2 Second function to test
#' @param range Range of test values
#' @param ntest Number of values to test in the range
#' @param logspace Flag indicating whether to space test values logarithmically.
#' @return logical
#' @keywords internal
chkinverse <- function(f1, f2, range=c(1.0, 5.0), ntest=10, logspace=FALSE)
{
    if(logspace)
        range <- log(range)

    xvals <- seq(range[1], range[2], length.out = ntest)

    if(logspace)
        xvals <- exp(xvals)

    xt1 <- f1(f2(xvals))
    xt2 <- f2(f1(xvals))

    isTRUE(all.equal(xt1, xvals)) && isTRUE(all.equal(xt2, xvals))
}
