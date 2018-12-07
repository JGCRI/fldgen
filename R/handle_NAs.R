#' Remove NA grid cells from processing in the emulator training process
#'
#' This function removes grid cells with missing data, NA values for all timesteps,
#' from being processed by \code{trainTP}. This enables \code{fldgen} to work with output that has NAs in it.
#'
#'
#' @param griddata_input output from \code{concatGrids.general}
#' @return A list similar in structure and content to \code{concatGrids.general} output but without the
#' NA gridcells if applicable and a mapping data frame of the original and new grid cell indices.
#' @importFrom dplyr filter select mutate group_by left_join %>% n summarise
#' @importFrom stats na.omit
#' @export
#' @keywords internal
drop_NAs <- function(griddata_input){

    # Silence package checks
    remove_index <- . <- original_index <- new_index <- lat <- lon <- NULL

    # IF there are no NA values then return a NULL mapping file and the input griddata.
    if(!TRUE %in% is.na(griddata_input$vardata)){

        griddata_output <- griddata_input
        griddata_output[['mapping']] <- NULL

    } else {

        # Pull out the array indices for the vardata values that are NAs.
        NA_ind <- as.data.frame(which(is.na(griddata_input$vardata), arr.ind = TRUE))

        # Count the number of NAs in each column. We expect the entire column
        # to equal NA, if there is a mix of NAs and numeric values then return an error.

        NA_ind %>%
            dplyr::group_by(col) %>%
            dplyr::summarise(NA_count = n()) ->
            NA_count

        NA_count$remove_index <- TRUE

        if(any(NA_count$NA_count != nrow(griddata_input$vardata))){

            stop('Inconsistent grid cell entries, both NAs and numeric values for ', names(griddata_input$tags))

        } else {

            tibble::tibble(lat = griddata_input$lat) %>%
                dplyr::mutate(join = 1) ->
                lat_tibble

            # Create a mapping file that contains the a Boolean to remove the original_index and the
            # new index number after the NA columns have been removed.
            tibble::tibble(lon = griddata_input$lon) %>%
                dplyr::mutate(join = 1) %>%
                dplyr::left_join(lat_tibble, by = "join") %>%
                dplyr::mutate(original_index = 1:ncol(griddata_input$vardata)) %>%
                dplyr::select(original_index, lat, lon)  %>%
                dplyr::left_join(NA_count, by = c("original_index" = "col")) ->
                original_NA_info

            # Create a df of the original index and the new index for the indices that
            original_NA_info %>%
                dplyr::filter(is.na(remove_index)) %>%
                dplyr::mutate(new_index = 1:nrow(.)) %>%
                dplyr::select(original_index, new_index, lat, lon) ->
                new_index_mapping

            # Combine the new_index_mapping file with the original and NA information df to
            # create a mapping file all of the original index values and the new index values.
            # If the new_index is NA then it means the column should be removed from the griddata_list
            original_NA_info %>%
                dplyr::left_join(new_index_mapping, by = c("original_index", "lat", "lon")) %>%
                dplyr::select(original_index, new_index, lat, lon) ->
                mapping

            # Create output list
            griddata_output <- griddata_input

            # Now remove the NA values from the list.
            griddata_output[['vardata']]  <- griddata_input$vardata[ , which(!is.na(mapping$new_index))]
            griddata_output[['globalop']] <- griddata_input$globalop[which(!is.na(mapping$new_index)), ]
            griddata_output[['mapping']]  <- data.frame(mapping)

        }

    }

    griddata_output

}


#' Add the NA grid cells that may have been removed during the the emulator training process
#'
#' This function adds NA gridcells that were removed during the training process in order to convert
#' the results to same dimensions as the input grid.
#'
#' @param data a matrix of the time by grid cell data that is missing NA gridcells
#' @param mapping a mapping data frame created by \code{drop_NAs} that contains the original and new index information for each grid cell.
#' @return a complete data frame of NA grid cells and the generated fields ordered by grid cell
#' @importFrom dplyr filter select mutate group_by left_join %>%
#' @export
#' @keywords internal
add_NAs <- function(data, mapping){

    # Silence package
    new_index <- NULL

    # Using information from the mapping file determine which indices from the original gird
    # exist in the data frame and which ones need to be fille with NAs.
    currently_exist   <- stats::na.omit(mapping)
    currently_missing <- dplyr::filter(mapping, is.na(new_index))

    # Use the original index to name the columns that exist in the data.
    data <- as.data.frame(data)
    names(data) <- currently_exist$original_index

    # Use the missing original indices as the column names for the data frame of NA values to add
    # to the data
    data_NA <- as.data.frame(matrix(nrow = nrow(data), ncol = nrow(currently_missing)))
    names(data_NA) <- currently_missing$original_index

    # Combine the data and the NA data frames
    output <- cbind(data, data_NA)

    # Reorder the columns to reflect the original index
    output <- output[as.character(mapping$original_index)]

    # Remove the column names from the output
    names(output) <- NULL

    # Return output as matrix -- may need to be a matrix!
    #as.matrix(output)
    output
}

