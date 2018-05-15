context('Data handling and processing')

inputfile <- system.file('extdata/tann1.nc', package='fldgen')
griddata <- read.temperatures(inputfile)

test_that('Data can be extracted by latitude/longitude', {
    nlat <- length(griddata$lat)
    nlon <- length(griddata$lon)
    ntime <- length(griddata$time)
    ## Whole grid
    op1 <- extract_box(c(-90, 90), c(0, 360), griddata)

    expect_equal(length(op1), nlat*nlon)
    expect_true(all(op1))

    ## Nonexistent box
    op2 <- extract_box(c(-90, 90), c(-180,-170), griddata)
    expect_equal(length(op2), nlat*nlon)
    expect_false(any(op2))

    ## 5deg^2 box.  This winds up to be 5x5 grid cells
    op3 <- extract_box(c(0,5), c(0,5), griddata)
    tas_subset <- griddata$tas[,op3]
    expect_equal(dim(tas_subset), c(nrow(griddata$tas), 25))

    ## box-average (unweighted)
    tmean <- griddata$tas %*% (op3/sum(op3))
    expect_equal(dim(tmean), c(ntime, 1))

    tmean_manual <-
        sapply(1:ntime,
               function(t) {
                   ## turn into vector so we can address by row and column
                   tf <- matrix(griddata$tas[t,], nrow=nlat)
                   mean(tf[97:101, 1:5])
               })
    expect_equal(tmean, matrix(tmean_manual, ncol=1))

    ## weighted version
    tgop_box <- griddata$tgop * op3
    wttmean <- griddata$tas %*% (tgop_box/sum(tgop_box))
    expect_equal(dim(tmean), c(ntime,1))

    tgop_box_manual <- matrix(griddata$tgop, nrow=nlat)[97:101, 1:5]
    norm <- sum(tgop_box_manual)
    wttmean_manual <-
        sapply(1:nrow(griddata$tas),
               function(r) {
                   tf <- matrix(griddata$tas[r,], nrow=nlat)
                   sum(tf[97:101, 1:5]  * tgop_box_manual) / norm
               })
    expect_equal(wttmean, matrix(wttmean_manual, ncol=1))
})


