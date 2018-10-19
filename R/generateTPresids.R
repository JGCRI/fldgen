#' Using a trained emulator, generate new fields of residuals.
#'
#' This function takes in a trained emulator - a structure of class
#' \code{fldgen}.
#' This structure contains everything the emulator has learned about the model,
#' and is used to generate new fields of residuals.
#'
#' First, new residuals are generated. The distribution of residuals over time
#' for each variable in each grid cell is N(0,1). These generated residuals are
#' then mapped back to the native distribution for each variable's residuals
#' in each grid cell.
#'
#'
#' @param emulator A trained \code{fldgen} temperature precipitation joint
#' emulator.
#' @param ngen The number of new fields to generate
#' @param method The algorithm used to generate new EOF coefficients, defaults
#' to 1.
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
#' @return A list containing:
#' 1)residgrids = A list of new residual fields, each entry in the list is a new
#' realization, a matrix that is [Nyears x 2 * Ngrid]; the first 1:Ngrid cols
#' are the temperature residuals and columns (Ngrid + 1):(2*Ngrid) are the
#' precipitation residuals.
#' 2) tvarunconvert_fcn = The function to inverse transform T values from
#' (-infinity, infinity) to the original support, if necessary.
#' 3) pvarunconvert_fcn = The function to inverse transform P values from
#' (-infinity, infinity) to the original support, if necessary.
#' @export
generate.TP.resids <- function(emulator, ngen, method = 1, tvarunconvert_fcn = NULL, pvarunconvert_fcn = exp)
{
    # dbl check that convert and unconvert functions are actual inverses on some test data.
    # throw a warning if this is not the case but it is up to the user to track down.
    testx <- 1:5

    if(!is.null(emulator$griddataT$tvarconvert_fcn) &
       !(is.null(tvarunconvert_fcn))){
        testx %>%
            emulator$griddataT$tvarconvert_fcn() %>%
            tvarunconvert_fcn() ->
            x1

        testx %>%
            tvarunconvert_fcn() %>%
            emulator$griddataT$tvarconvert_fcn() ->
            x2

        if((max(abs(x1-testx))>1e-14) | (max(abs(x2-testx))>1e-14) ){
            warning('Your functions for transforming and inverse transforming the support of T may not be inverses of each other 1.')
        }

    } else if((!is.null(emulator$griddataT$tvarconvert_fcn) & is.null(tvarunconvert_fcn)) |
              (is.null(emulator$griddataT$tvarconvert_fcn) & !is.null(tvarunconvert_fcn))){
        warning('Your functions for transforming and inverse transforming the support of T may not be inverses of each other 2.')
    }

    if(!is.null(emulator$griddataP$pvarconvert_fcn) &
       !(is.null(pvarunconvert_fcn))){
        testx %>%
            emulator$griddataP$pvarconvert_fcn() %>%
            pvarunconvert_fcn() ->
            x1

        testx %>%
            pvarunconvert_fcn() %>%
            emulator$griddataP$pvarconvert_fcn() ->
            x2

        if((max(abs(x1-testx))>1e-14) | (max(abs(x2-testx))>1e-14) ){
            warning('Your functions for transforming and inverse transforming the support of P may not be inverses of each other 1.')
        }

    } else if((!is.null(emulator$griddataP$pvarconvert_fcn) & is.null(pvarunconvert_fcn)) |
              (is.null(emulator$griddataP$pvarconvert_fcn) & !is.null(pvarunconvert_fcn))){
        warning('Your functions for transforming and inverse transforming the support of P may not be inverses of each other 2.')
    }

    ####
    Ngrid <- ncol(emulator$meanfldT$r)

    # Generate the new residual fields in the normal space.
    newgrids <- lapply(1:ngen,
                       function(x) {
                           ## It takes the full fldgen object, but doesn't appear
                           ## to use anything other than phase and magnitude,
                           ## which by design contains both T and P. May be able
                           ## to improve implementation in future.
                           # This uses a RNG, it will be different every call:
                           bcoord <- mkcorrts(emulator, method=method)
                           ## Because the generated residuals are generated to
                           ## follow N(0, 1), we can't add back to the meanfield
                           ## as is. They have to be transformed back to the
                           ## distribution of T or P residuals that are
                           ## empirically characterized for each grid cell in
                           ## emulator$tfuns$quant, emulator$pfuns$quant.
                           reconst_fields(emulator$reof$rotation, bcoord)
                       })

    residgrids <- lapply(newgrids, function(g) {
        g[, 1:Ngrid] <- unnormalize.resids(empiricalquant = emulator$tfuns$quant,
                                           rn = g[ ,1:Ngrid])$rnew

        g[, (Ngrid+1):(2*Ngrid)] <- unnormalize.resids(empiricalquant = emulator$pfuns$quant,
                                                       rn = g[ , (Ngrid+1):(2*Ngrid)])$rnew

        return(g)}
    )


    # check the conversion
    if(!is.null(emulator$griddataT$tvarconvert_fcn)){

         residgrids <- lapply(residgrids, function(g){
            g[, 1:Ngrid] <- tvarunconvert_fcn(g[, 1:Ngrid])
            return(g)
        })
    } else{
        print('Generated T residuals are not being transformed to a different support. Up to user to know if this is desirable.')
    }


    if(!is.null(emulator$griddataP$pvarconvert_fcn)){

        residgrids <- lapply(residgrids, function(g){
            g[, (Ngrid+1):(2*Ngrid)] <- pvarunconvert_fcn(g[, (Ngrid+1):(2*Ngrid)])
            return(g)
        })
    } else{
        print('Generated P residuals are not being transformed to a different support. Up to user to know if this is desirable.')
    }



    ### ADD THE NA COLUMNS BACK FOR ISIMIP DATA

    ### return the generated residuals
    return(list(residgrids = residgrids, tvarunconvert_fcn = tvarunconvert_fcn, pvarunconvert_fcn = pvarunconvert_fcn))
}
