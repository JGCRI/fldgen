context('Field generation on partial input grids')

library(fldgen)

## Input file locations
datadir <- file.path(getwd(), 'data')
tasfile <- file.path(datadir,'tas_annual_ipsl-cm5a-lr_rcp8p5_xxx_2006-2099.nc')
prfile <- file.path(datadir,'pr_annual_ipsl-cm5a-lr_rcp8p5_xxx_2006-2099.nc')
flist <- c(tasfile, prfile)
coordfile <- file.path(datadir, 'coord-ref.csv')
tgav_suffix <- 'GlobalAvg.txt'
flist <- c(tasfile, prfile)

ngrid <- 67420                          # number of valid grid cells in these datasets
ntime <- 94                             # number of time steps in these datasets

## Train
emu <- trainTP(flist, tvarname = "tas", pvarname = "pr", globalAvg_file =
                 tgav_suffix)

## Generate residuals
ngen <- 2
resids <- generate.TP.resids(emu, ngen=ngen)

## Generate full fields
tgav <- emu$tgav
pvarunconvertfn <- exp
fullgrids <- generate.TP.fullgrids(emu, resids, tgav, NULL, pvarunconvertfn)

## read the reference coordinate file
coord_ref <- readr::read_csv(coordfile)

test_that('Emulator was trained correctly.', {
    ## Check that the global mean operator is properly normalized
    expect_equal(sum(emu$griddataT$globalop), 1)
    expect_equal(sum(emu$griddataP$globalop), 1)

    ## Check tgav values are a ntime x 1 matrix
    expect_true(is.matrix(emu$tgav))
    expect_equal(dim(emu$tgav), c(ntime,1))

    ## Check the dimension of the matrix of EOFs.
    expect_equal(dim(emu$reof$rotation), c(2*ngrid, ntime))

    ## Check the dimensions of the components of the griddata vars
    for(gd in list(emu$griddataT, emu$griddataP)) {
        expect_equal(dim(gd$vardata), c(ntime, ngrid))
        expect_equal(length(gd$globalop), ngrid)
        expect_equal(gd$ncol_full, length(gd$lat)*length(gd$lon))
        expect_equal(length(gd$gridid_full), ngrid)

        ## Check the generated coordinates against the reference (note that
        ## passing this test perforce means that we dropped the correct grid
        ## cells)
        expect_equal(gd$coord[ , 'lat'], coord_ref$lat)
        expect_equal(gd$coord[ , 'lon'], coord_ref$lon)
    }
})


test_that('Residuals have the right size.', {
    expect_equal(length(resids), ngen)
    for(resid_grid in resids) {
        expect_true(is.matrix(resid_grid))
        expect_equal(dim(resid_grid), c(ntime, 2*ngrid)) # 2*ngrid because the
                                        # grid has both variables.
    }

    ## Check that the generated residuals are different from one another.
    for(i in seq(1, ngen-1)) {
        for(j in seq(2, ngen)) {
            expect_false(isTRUE(all.equal(resids[[i]], resids[[j]])))
        }
    }
})


test_that('Fullgrids structure is constructed correctly.', {

    for(meanfld in list(fullgrids$meanfieldT, fullgrids$meanfieldP)) {
        expect_true(is.matrix(meanfld))
        expect_equal(dim(meanfld), c(ntime, ngrid))
    }

    for(r in fullgrids$fullgrids) {
        for(grd in c('tas','pr')) {
            m <- r[[grd]]
            expect_true(is.matrix(m))
            expect_equal(dim(m), c(ntime, ngrid))
        }
    }
})
