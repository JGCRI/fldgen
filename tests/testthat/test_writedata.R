context('Writing data')

griddata <- read.temperatures(system.file('extdata/tann1.nc', package='fldgen'))

test_that('data written to netCDF is identical to original.',
{
    d <- tempdir()
    file <- file.path(d,'test.nc')
    if(file.exists(file)) {
        unlink(file)
    }
    write.temperature(griddata$tas, file, griddata)
    expect_true(file.exists(file))

    griddata2 <- read.temperatures(file)
    expect_equal(griddata2$tas, griddata$tas)

    unlink(file)
})


test_that('clobber argument to write.temperature works',
{
    ## create a dummy file
    d <- tempdir()
    file <- file.path(d, 'test.nc')
    write('NOT to leave the room, even if you come and get him.',
          file)

    expect_error(write.temperature(griddata$tas, file, griddata), 'exists and clobber')

    expect_silent(write.temperature(griddata$tas, file, griddata, clobber=TRUE))

    ## verify that we can read the result
    expect_silent(read.temperatures(file))

    unlink(file)
})


## Test loading and saving models.  Technically the former is reading rather
## than writing data, but it's a lot more convenient to test these things
## together.

# test_that('reading and writing models works',
# {
#     ## These tests are kind of long, since the fldgen objects are rather large.
#     ## Therefore, skip them on the CI environments, which are somewhat
#     ## time-sensitive.  And since we're skipping them on the CI, I also haven't
#     ## included the data file in the repository (it is quite large for the time
#     ## being).  We can revisit this once we find a way to shrink the memory
#     ## footprint of the fldgen objects.
#     skip_on_travis()
#     skip_on_appveyor()
#     skip_on_covr()
#
#     ## read data in legacy format
#     oldfile <- 'data/testmodel-oldfmt.rda'
#     expect_silent(testmodel1 <- loadmodel(oldfile, oldfmt=TRUE))
#     expect_is(testmodel1, 'fldgen')
#
#     ## write in new format
#     newfile <- tempfile(fileext='.rds')
#     expect_silent(savemodel(testmodel1, newfile, compress=FALSE))
#     expect_true(file.exists(newfile))
#
#     ## read in new format
#     expect_silent(testmodel2 <- loadmodel(newfile))
#     expect_is(testmodel2, 'fldgen')
#     expect_equal(testmodel1, testmodel2)
#
#     ## trying to read with mismatched formats won't work
#     expect_error(suppressWarnings(tm1 <- loadmodel(newfile, oldfmt=TRUE), 'magic number'))
#     expect_error(tm2 <- loadmodel(oldfile, oldfmt=FALSE), 'unknown input format')
#     unlink(newfile)
# })
