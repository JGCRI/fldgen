context('Reading data')

inputfile <- system.file('extdata/tann1.nc', package='fldgen')
ntime <- 95
nlat <- 192
nlon <- 288
ngrid <- nlat*nlon

## This will be used in two tests below
griddata <- read.ncdf(inputfile)

test_that('Full data read works.',
{
    expect_equal(dim(griddata$tas), c(ntime, ngrid))
    expect_equal(dim(griddata$tgop), c(ngrid, 1))
    expect_equal(length(griddata$lat), nlat)
    expect_equal(length(griddata$lon), nlon)
    expect_equal(length(griddata$time), ntime)
})


test_that('Trim time series length works.',
{
    ntime.trim <- 10
    gdtrim <- read.ncdf(inputfile, ntime.trim)

    expect_equal(dim(gdtrim$tas), c(ntime.trim, ngrid))
    expect_equal(dim(gdtrim$tgop), c(ngrid, 1))
    expect_equal(length(gdtrim$lat), nlat)
    expect_equal(length(gdtrim$lon), nlon)
    expect_equal(length(gdtrim$time), ntime.trim)

    expect_equal(gdtrim$tas, griddata$tas[1:ntime.trim, ])

})
