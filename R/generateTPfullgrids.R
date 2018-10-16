#' Using a trained emulator, generate new full fields.
#'
#' This function takes in a trained emulator - a structure of class
#' \code{fldgen}.
#' This structure contains everything the emulator has learned about the model,
#' and is used to generate new fields of residuals.
#'
#' First, ...
#'
#' @param emulator A trained \code{fldgen} temperature precipitation joint
#' emulator.
#' @param residgrids Descript
#' @param tgav  Descript
#' @param reconstruction_function Descript
#' @param addNAs descript
#' @param tvarunconvert_fcn The function to undo any transformation done to the
#' input training data in \code{trainTP} to correct the support. This should be
#' the inverse function of the tvarconvert_fcn argument to \code{trainTP}.
#' Defaults to NULL as temperature typically doesn't need to be transformed to
#' correct the support. WARNING: currently rely on the user to define the
#' correct inverse function here.
#' @param pvarunconvert_fcn The function to undo any transformation done to the
#' input training data in \code{trainTP} to correct the support. This should be
#' the inverse function of the pvarconvert_fcn argument to \code{trainTP}.
#' Defaults to exp() as precipitation is usually log-transformed in ordere to
#' correct the support. WARNING: currently rely on the user to define the
#' correct inverse function here.
#' @return Descript - can crib some of the language from generateTPresids
#' @export

generate.TP.fullgrids <- function(emulator, residgrids, tgav, reconstruction_function = pscl_apply, addNAs = FALSE,
                                  tvarunconvert_fcn = NULL, pvarunconvert_fcn = exp){
    # TODO kalyn you need to add documentation and testing for this function!!


    # Check inputs
    if(!is.data.frame(tgav)){ stop('tgav must be a data frame') }

    req_cols <- c('tgav', 'time')
    missing  <- !req_cols %in% names(tgav)
    if(any(missing)){ stop('tgav missing the following columns: ', paste(req_cols[missing], collapse = ', ')) }


    # TODO eventually wer are going to want to check time that is actually pulled from the traing
    # data instead of parsing it out from the tag names.
    startYr_endYr <- gsub(".nc", "", stringr::str_split(basename(emulator$infiles), "_")[[1]][6])
    startYr <- as.integer(stringr::str_split(startYr_endYr, '-')[[1]][1])
    endYr   <- as.integer(stringr::str_split(startYr_endYr, '-')[[1]][2])

    nc_time <- startYr:endYr
    # Compare the emulator time with the tgav time
    missing_time <- !tgav$time %in% nc_time

    if(any(missing_time)){stop('emulator was not trained with data for the following years: ', paste(tgav$time[missing_time], collapse = ', '))}


    if(any(!nc_time %in% tgav$time)){

        # TODO
        stop('Need to add code that will subset the resid and emulator training data for years in common')

    }

    # Define the time and tgav
    time <- tgav$time
    tgav <- matrix(tgav$tgav, ncol = 1)

    meanfieldT <- reconstruction_function(emulator$meanfldT, tgav)
    meanfieldP <- reconstruction_function(emulator$meanfldP, tgav)

    # save a copy of the number of grids cells for each variable
    Ngrid <- ncol(residgrids[[1]])/2


    # Add the meanfield values to the residual grids
    lapply(residgrids, function(matrix, gridcells = Ngrid){

        # Separate the tas and pr data from one another.
        tas <- matrix[ , 1:Ngrid]
        pr  <- matrix[ , (Ngrid + 1):(2 * Ngrid)]

        # Add the meanfield to the data
        tas[ , 1:Ngrid] <- tas[ , 1:Ngrid] + meanfieldT
        pr[ , 1:Ngrid]  <- pr[ , 1:Ngrid] + meanfieldP

        # Return output
        list(tas = tas, pr = pr)

    }) ->
        fullgrids


    ## deal with the support conversion if necessary

    if(!is.null(emulator$griddataT$tvarconvert_fcn)){

        fullgrids <- lapply(fullgrids, function(g){
            g[, 1:Ngrid] <- tvarunconvert_fcn(g[, 1:Ngrid])
            return(g)
        })
    } else{
        print('Generated T full fields are not being touched')
    }


    if(!is.null(emulator$griddataP$pvarconvert_fcn)){

        fullgrids <- lapply(fullgrids, function(g){
            g[, (Ngrid+1):(2*Ngrid)] <- pvarunconvert_fcn(g[, (Ngrid+1):(2*Ngrid)])
            return(g)
        })
    } else{
        print('Generated P full fields are not being touched')
    }



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

    return(list(fullgrids = fullgrids, coordinates = coordinates, time = time))

}
