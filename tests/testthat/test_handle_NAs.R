

context('Handle NAs')


# Read in the nc file
inputfile     <- system.file('extdata/tas_annual_esm_rcp_r2i1p1_2006-2100.nc', package='fldgen')
original_data <- read.general(inputfile, len = 2, tag = 'g1', varname='tas',
                                                  latvar='lat_2', lonvar='lon_2', timevar='time')

test_that('When there are no NAs there are no changes.', {

    out1 <- drop_NAs(original_data)

    expect_identical(out1, original_data)

})


test_that('drop_NAs throws errors.', {

    # Change some of the entries to NAs but not all
    data_bad <- original_data
    data_bad$vardata[ 1, 7:10] <- NA

    expect_error(drop_NAs(data_bad), 'Inconsistent grid cell entries, both NAs and numeric values for g1')

    })

test_that('drop_NAs & add_NAs are inverses.', {

    # Change some of the columns to NAs
    data_w_NAs <- original_data
    data_w_NAs$vardata[ , 7:10] <- NA

    # Use the function to remove the NAs, we expect this to drop columns 7-10 (4 columns)
    data_wo_NAs <- drop_NAs(data_w_NAs)

    expected_cols <- ncol(data_w_NAs$vardata) - 4
    expect_equal(ncol(data_wo_NAs$vardata), expected_cols)

    expect_false(is.null(data_wo_NAs$ncol_full))
    expect_false(is.null(data_wo_NAs$gridid_full))
    expect_false(is.null(data_wo_NAs$coord))

    # Check to see that when the NAs are added back the output matches the original data with the NAs.
    data_NAs_added <- add_NAs(data_wo_NAs$vardata, data_wo_NAs)

    expect_equal(data_NAs_added, data_w_NAs$vardata)
})












