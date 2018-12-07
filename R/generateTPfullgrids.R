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
#' @param tgavdf  A data frame with two columns. Column \code{tgav} = the vector
#' of global annual mean temperatures for constructing a mean field with
#' \code{reconstruction_function}.To be added to each list entry of residual
#' fields in \code{residgrids}. And a column \code{time}, used to subset the
#' trained emulator so that pattern scaling information from the appropriate
#' years is applied to tgav time series. Note: grids cannot be generated for
#' tgav years missing from the trained emulator.
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
#' @param addNAs A logical TRUE/FALSE, whether NA values need to be added back
#' to the full grids of data (that is, if ISIMIP data, restricting to land cells
#' only, was used for training rather than data with full global coverage.)
#' Defaults to FALSE.
#' @return A list of:
#' 1) fullgrids = A list of new full fields, each entry in
#' the list is a new realization, a matrix that is [Nyears x 2 * Ngrid]; the
#' first 1:Ngrid cols are the temperature field and columns
#' (Ngrid + 1):(2*Ngrid) are the precipitation field.
#' 2) coordinates = something KD needs to fill in, relevant to using ISIMIP data
#' in the grand experiment.
#' 3) time = the input time column used to subset the trained emulator so that
#' pattern scaling information from the appropriate years is applied to tgav
#' time series. Note: grids cannot be generated for tgav years missing from
#' the trained emulator.
#' 4) meanfieldT = the reconstructed, pattern scaled temperature mean field.
#' 5) meanfieldP = the reconstructed, pattern scaled precipitation mean field.
#' 6) tvarunconvert_fcn = the tvarunconvert_fcn input.
#' 7) pvarunconvert_fcn = the pvarunconvert_fcn input.
#' @importFrom stringr str_split
#' @export

generate.TP.fullgrids <- function(emulator, residgrids, tgavdf,
                                  tvarunconvert_fcn = NULL, pvarunconvert_fcn = exp,
                                  reconstruction_function = pscl_apply, addNAs = FALSE){

    # silence package checks
    column_index <- lat <- lon <- NULL

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


    # Check inputs
    if(!is.data.frame(tgavdf)){ stop('tgavdf must be a data frame') }

    req_cols <- c('tgav', 'time')
    missing  <- !req_cols %in% names(tgavdf)
    if(any(missing)){ stop('tgav missing the following columns: ',
                           paste(req_cols[missing], collapse = ', ')) }


    # TODO eventually we are going to want to check time that is actually pulled from the traing
    # data instead of parsing it out from the tag names.
    startYr_endYr <- gsub(".nc", "", stringr::str_split(basename(emulator$infiles), "_")[[1]][6])
    startYr <- as.integer(stringr::str_split(startYr_endYr, '-')[[1]][1])
    endYr   <- as.integer(stringr::str_split(startYr_endYr, '-')[[1]][2])

    nc_time <- startYr:endYr
    # Compare the emulator time with the tgavdf time
    missing_time <- !tgavdf$time %in% nc_time

    if(any(missing_time)){
        stop('emulator was not trained with data for the following years: ',
             paste(tgavdf$time[missing_time], collapse = ', '))
        }


    if(any(!nc_time %in% tgavdf$time)){

        # TODO
        stop('Need to add code that will subset the resid and emulator training data for years in common')

    }


    # save a copy of the number of grids cells for each variable
    Ngrid <- ncol(residgrids[[1]])/2


    # Define the time and tgav
    time <- tgavdf$time
    tgav <- matrix(tgavdf$tgav, ncol = 1)

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
            print('Generated T full fields not being transformed to a different support. Up to user to know if this is desirable.')
        }


        if(!is.null(emulator$griddataP$pvarconvert_fcn)){

            pr <- pvarunconvert_fcn(pr)

        } else{
            print('Generated P full fields not being transformed to a different support. Up to user to know if this is desirable.')
        }



        # Return output
        return(list(tas = tas, pr = pr))


    }) ->
        fullgrids




    # Create the croodinates mapping file.
    # TODO there is probably a better way to do this
    # Create a coordinate mapping data frame. Since
    # the lat is the most rapidly varying index left
    # join the lat to a lon tibble.
    tibble::tibble(lat = emulator$griddataT$lat) %>%
        dplyr::mutate(join = 1) ->
        lat_tibble

    tibble::tibble(lon = emulator$griddataT$lon) %>%
        dplyr::mutate(join = 1) %>%
        dplyr::left_join(lat_tibble, by = "join") %>%
        dplyr::mutate(column_index = 1:ncol(emulator$griddataT$vardata)) %>%
        dplyr::select(column_index, lat, lon) ->
        coordinates

    # If griddata mapping files exist then then add NA grid cells will have to be added back into
    # the grids in order to return output that has the same grid dimensions.
    if(!is.null(emulator$griddataT$mapping) & addNAs){

        # Since the T and P fields are concatenated together in a single data frame the mapping data must be updated to
        # reflect the single data frame.
        emulator$griddataT$mapping %>%
            dplyr::bind_rows(tibble::tibble(original_index = emulator$griddataP$mapping$original_index + max(emulator$griddataT$mapping$original_index),
                                            new_index = emulator$griddataP$mapping$new_index + max(emulator$griddataT$mapping$new_index, na.rm = TRUE))) ->
            mapping

        # Use the new mapping file to insert the NAs into the fullgrids
        fullgrids <- lapply(fullgrids, FUN = add_NAs, mapping = mapping)

    } else if(!is.null(emulator$griddataT$mapping) & !addNAs){

        # If the NAs were removed but the user has selected NOT to add the NAs back
        # to the grid then discard the NA cells from the coordinates mapping file.
        emulator$griddataT$mapping %>%
            na.omit %>%
            dplyr::select(column_index = new_index, lat, lon) ->
            coordinates

    }




    return(list(fullgrids = fullgrids, coordinates = coordinates, time = time,
                meanfieldT = meanfieldT, meanfieldP = meanfieldP,
                tvarunconvert_fcn = tvarunconvert_fcn, pvarunconvert_fcn = pvarunconvert_fcn))

}
