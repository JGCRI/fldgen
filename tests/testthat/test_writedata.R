context('Writing data')

griddata <- read.temperatures(system.file('extdata/tann1.nc', package='fldgen'))

test_that('data written to netCDF is identical to original.',
{
    d <- tempdir()
    file <- file.path(d,'test.nc')
    if(file.exists(file)) {
        unlink(file)
    }
    write.ncdf(griddata$tas, file, griddata)
    expect_true(file.exists(file))

    griddata2 <- read.temperatures(file)
    expect_equal(griddata2$tas, griddata$tas)

    unlink(file)
})


test_that('clobber argument to write.ncdf works',
{
    ## create a dummy file
    d <- tempdir()
    file <- file.path(d, 'test.nc')
    write('NOT to leave the room, even if you come and get him.',
          file)

    expect_error(write.ncdf(griddata$tas, file, griddata), 'exists and clobber')

    expect_silent(write.ncdf(griddata$tas, file, griddata, clobber=TRUE))

    ## verify that we can read the result
    expect_silent(read.temperatures(file))

    unlink(file)
})
