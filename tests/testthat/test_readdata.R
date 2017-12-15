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

    expect_equal(griddata$tags, list(tann1.nc=c(1,ntime)))

    expect_equal(class(griddata), 'griddata')
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

    expect_equal(gdtrim$tags, list(tann1.nc=c(1,ntime.trim)))

    expect_equal(gdtrim$tas, griddata$tas[1:ntime.trim, ])

})


test_that('Grids can be concatenated',
{
    g1 <- read.ncdf(inputfile, 2, tag='g1')
    g2 <- read.ncdf(inputfile, 3, tag='g2')
    g3 <- read.ncdf(inputfile, 4, tag='g3')

    g23 <- c(g2, g3)
    expect_equal(g23$tas, rbind(g2$tas, g3$tas))
    expect_equal(g23$tgop, g2$tgop)
    expect_equal(g23$lat, g3$lat)
    expect_equal(g23$lon, g3$lon)
    expect_equal(g23$time, c(g2$time, g3$time))
    expect_equal(g23$tags,
                 list(g2=c(1,3), g3=c(4,7)))

    g123 <- c(g1, g23)
    expect_equal(g123$tas, rbind(g1$tas, g2$tas, g3$tas))
    expect_equal(g123$tgop, g1$tgop)
    expect_equal(g123$lat, g2$lat)
    expect_equal(g123$lon, g3$lon)
    expect_equal(g123$time, c(g1$time, g2$time, g3$time))
    expect_equal(g123$tags,
                 list(g1=c(1,2), g2=c(3,5), g3=c(6,9)))
})
