context('Reading data')


## Legacy tests of temperature specific files
inputfile <- system.file('extdata/tann1.nc', package='fldgen')
ntime <- 95
nlat <- 192
nlon <- 288
ngrid <- nlat*nlon

## This will be used in two tests below
griddata <- read.temperatures(inputfile)

test_that('legacy T - Full data read works.',
{
    expect_equal(dim(griddata$tas), c(ntime, ngrid))
    expect_equal(dim(griddata$tgop), c(ngrid, 1))
    expect_equal(length(griddata$lat), nlat)
    expect_equal(length(griddata$lon), nlon)
    expect_equal(length(griddata$time), ntime)

    expect_equal(griddata$tags, list(tann1.nc=c(1,ntime)))

    expect_equal(class(griddata), 'griddata')
})


test_that('legacy T - Trim time series length works.',
{
    ntime.trim <- 10
    gdtrim <- read.temperatures(inputfile, ntime.trim)

    expect_equal(dim(gdtrim$tas), c(ntime.trim, ngrid))
    expect_equal(dim(gdtrim$tgop), c(ngrid, 1))
    expect_equal(length(gdtrim$lat), nlat)
    expect_equal(length(gdtrim$lon), nlon)
    expect_equal(length(gdtrim$time), ntime.trim)

    expect_equal(gdtrim$tags, list(tann1.nc=c(1,ntime.trim)))

    expect_equal(gdtrim$tas, griddata$tas[1:ntime.trim, ])

})


test_that('legacy T - Grids can be concatenated and split',
{
    g1 <- read.temperatures(inputfile, 2, tag='g1')
    g2 <- read.temperatures(inputfile, 3, tag='g2')
    g3 <- read.temperatures(inputfile, 4, tag='g3')

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

    sg <- splitGrids(g123)
    expect_equal(length(sg), 3)
    expect_equal(sg[[1]], g1)
    expect_equal(sg[[2]], g2)
    expect_equal(sg[[3]], g3)
})


test_that('legacy T - Splitting a grid with a single element is ok.',
{
    gdl <- splitGrids(griddata)
    expect_equal(length(gdl), 1)
    expect_equal(gdl[[1]], griddata)
})


## Updated tests for temperature and precipitation
inputTfile <- system.file('extdata/tas_annual_esm_rcp_r2i1p1_startyr-endyr.nc', package='fldgen')
inputPfile <- system.file('extdata/pr_annual_esm_rcp_r2i1p1_startyr-endyr.nc', package='fldgen')
ntime <- 95
nlat <- 192
nlon <- 288
ngrid <- nlat*nlon


## T first
## This will be used in two tests below
griddata <- read.general(filename = inputTfile,  varname='tas',
                         latvar='lat_2', lonvar='lon_2',
                         timevar='time')

test_that('Full T data read works.',
          {
              expect_equal(dim(griddata$vardata), c(ntime, ngrid))
              expect_equal(dim(griddata$globalop), c(ngrid, 1))
              expect_equal(length(griddata$lat), nlat)
              expect_equal(length(griddata$lon), nlon)
              expect_equal(length(griddata$time), ntime)

              expect_equal(griddata$tags, list(`tas_annual_esm_rcp_r2i1p1_startyr-endyr.nc`=c(1,ntime)))

              expect_equal(class(griddata), 'griddata')
          })


test_that('Trim T time series length works.',
          {
              ntime.trim <- 10
              gdtrim <- read.general(filename = inputTfile, len = ntime.trim,
                                     varname='tas',
                                     latvar='lat_2', lonvar='lon_2',
                                     timevar='time')

              expect_equal(dim(gdtrim$vardata), c(ntime.trim, ngrid))
              expect_equal(dim(gdtrim$globalop), c(ngrid, 1))
              expect_equal(length(gdtrim$lat), nlat)
              expect_equal(length(gdtrim$lon), nlon)
              expect_equal(length(gdtrim$time), ntime.trim)

              expect_equal(gdtrim$tags, list(`tas_annual_esm_rcp_r2i1p1_startyr-endyr.nc`=c(1,ntime.trim)))

              expect_equal(gdtrim$vardata, griddata$vardata[1:ntime.trim, ])

          })


test_that('T Grids can be concatenated and split',
          {
              g1 <- read.general(filename = inputTfile, len = 2,
                                 tag = 'g1', varname='tas',
                                 latvar='lat_2', lonvar='lon_2',
                                 timevar='time')
              g2 <- read.general(filename = inputTfile, len = 3,
                                 tag = 'g2', varname='tas',
                                 latvar='lat_2', lonvar='lon_2',
                                 timevar='time')
              g3 <- read.general(filename = inputTfile, len = 4,
                                 tag = 'g3', varname='tas',
                                 latvar='lat_2', lonvar='lon_2',
                                 timevar='time')

              g23 <- concatGrids.general(list(g2, g3))
              expect_equal(g23$vardata, rbind(g2$vardata, g3$vardata))
              expect_equal(g23$globalop, g2$globalop)
              expect_equal(g23$lat, g3$lat)
              expect_equal(g23$lon, g3$lon)
              expect_equal(g23$time, c(g2$time, g3$time))
              expect_equal(g23$tags,
                           list(g2=c(1,3), g3=c(4,7)))

              g123 <- concatGrids.general(list(g1, g23))
              expect_equal(g123$vardata, rbind(g1$vardata, g2$vardata, g3$vardat))
              expect_equal(g123$globalop, g1$globalop)
              expect_equal(g123$lat, g2$lat)
              expect_equal(g123$lon, g3$lon)
              expect_equal(g123$time, c(g1$time, g2$time, g3$time))
              expect_equal(g123$tags,
                           list(g1=c(1,2), g2=c(3,5), g3=c(6,9)))

              sg <- splitGrids.general(g123)
              expect_equal(length(sg), 3)
              expect_equal(sg[[1]], g1)
              expect_equal(sg[[2]], g2)
              expect_equal(sg[[3]], g3)
          })


test_that('Splitting a T grid with a single element is ok.',
          {
              gdl <- splitGrids.general(griddata)
              expect_equal(length(gdl), 1)
              expect_equal(gdl[[1]], griddata)
          })



## P second
## This will be used in two tests below
griddata <- read.general(filename = inputPfile,  varname='pr',
                         latvar='lat', lonvar='lon',
                         timevar='time')

test_that('Full P data read works.',
          {
              expect_equal(dim(griddata$vardata), c(ntime, ngrid))
              expect_equal(dim(griddata$globalop), c(ngrid, 1))
              expect_equal(length(griddata$lat), nlat)
              expect_equal(length(griddata$lon), nlon)
              expect_equal(length(griddata$time), ntime)

              expect_equal(griddata$tags, list(`pr_annual_esm_rcp_r2i1p1_startyr-endyr.nc`=c(1,ntime)))

              expect_equal(class(griddata), 'griddata')
          })


test_that('Trim P time series length works.',
          {
              ntime.trim <- 10
              gdtrim <- read.general(filename = inputPfile, len = ntime.trim,
                                     varname='pr',
                                     latvar='lat', lonvar='lon',
                                     timevar='time')

              expect_equal(dim(gdtrim$vardata), c(ntime.trim, ngrid))
              expect_equal(dim(gdtrim$globalop), c(ngrid, 1))
              expect_equal(length(gdtrim$lat), nlat)
              expect_equal(length(gdtrim$lon), nlon)
              expect_equal(length(gdtrim$time), ntime.trim)

              expect_equal(gdtrim$tags, list(`pr_annual_esm_rcp_r2i1p1_startyr-endyr.nc`=c(1,ntime.trim)))

              expect_equal(gdtrim$vardata, griddata$vardata[1:ntime.trim, ])

          })


test_that('P Grids can be concatenated and split',
          {
              g1 <- read.general(filename = inputPfile, len = 2,
                                 tag = 'g1', varname='pr',
                                 latvar='lat', lonvar='lon',
                                 timevar='time')
              g2 <- read.general(filename = inputPfile, len = 3,
                                 tag = 'g2', varname='pr',
                                 latvar='lat', lonvar='lon',
                                 timevar='time')
              g3 <- read.general(filename = inputPfile, len = 4,
                                 tag = 'g3', varname='pr',
                                 latvar='lat', lonvar='lon',
                                 timevar='time')

              g23 <- concatGrids.general(list(g2, g3))
              expect_equal(g23$vardata, rbind(g2$vardata, g3$vardata))
              expect_equal(g23$globalop, g2$globalop)
              expect_equal(g23$lat, g3$lat)
              expect_equal(g23$lon, g3$lon)
              expect_equal(g23$time, c(g2$time, g3$time))
              expect_equal(g23$tags,
                           list(g2=c(1,3), g3=c(4,7)))

              g123 <- concatGrids.general(list(g1, g23))
              expect_equal(g123$vardata, rbind(g1$vardata, g2$vardata, g3$vardat))
              expect_equal(g123$globalop, g1$globalop)
              expect_equal(g123$lat, g2$lat)
              expect_equal(g123$lon, g3$lon)
              expect_equal(g123$time, c(g1$time, g2$time, g3$time))
              expect_equal(g123$tags,
                           list(g1=c(1,2), g2=c(3,5), g3=c(6,9)))

              sg <- splitGrids.general(g123)
              expect_equal(length(sg), 3)
              expect_equal(sg[[1]], g1)
              expect_equal(sg[[2]], g2)
              expect_equal(sg[[3]], g3)
          })


test_that('Splitting a P grid with a single element is ok.',
          {
              gdl <- splitGrids.general(griddata)
              expect_equal(length(gdl), 1)
              expect_equal(gdl[[1]], griddata)
          })
