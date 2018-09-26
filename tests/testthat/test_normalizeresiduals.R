context('Native-normal mappings inverse check')

###############################
### Forward: Native to Normal to Native

### Get the residuals for the pair of example TP files
emulator <- suppressWarnings(trainTP(c(system.file('extdata/tas_annual_esm_rcp_r2i1p1_startyr-endyr.nc', package='fldgen'),
                      system.file('extdata/pr_annual_esm_rcp_r2i1p1_startyr-endyr.nc', package='fldgen')),
                    tvarname = "tas", tlatvar = "lat_2", tlonvar = "lon_2",
                    pvarname = "pr", platvar = "lat", plonvar = "lon"))

residT <- emulator$meanfldT$r
tfuns <- emulator$tfuns
tfuns2 <- characterize.emp.dist(residT)
residP <- emulator$meanfldP$r
pfuns <- emulator$pfuns
pfuns2 <- characterize.emp.dist(residP)

rm(emulator)

test_that('Emulator is calling characterize.emp.dist correctly for T residuals',
          {
              expect_equal(length(tfuns$cdf)  , length(tfuns2$cdf))
              expect_equal(length(tfuns$quant), length(tfuns2$quant))

              expect_equal(length(tfuns$cdf), length(tfuns$quant))
              expect_equal(length(tfuns2$cdf), length(tfuns2$quant))

              expect_equal(length(tfuns$cdf), ncol(residT))


              testdata <- t(as.matrix(runif(length(tfuns$cdf), 0.1, 6)))


              get.quans <- function(index){
                  tfuns$cdf[[index]](testdata[,index])
              }

              q <- sapply(1:ncol(testdata), get.quans)


              get.quans2 <- function(index){
                  tfuns2$cdf[[index]](testdata[,index])
              }

              q2 <- sapply(1:ncol(testdata), get.quans2)

              expect_equal(q, q2)

          })

test_that('Emulator is calling characterize.emp.dist correctly for P residuals',
          {
              expect_equal(length(pfuns$cdf)  , length(pfuns2$cdf))
              expect_equal(length(pfuns$quant), length(pfuns2$quant))

              expect_equal(length(pfuns$cdf),  length(pfuns$quant))
              expect_equal(length(pfuns2$cdf), length(pfuns2$quant))

              expect_equal(length(pfuns$cdf), ncol(residP))


              testdata <- t(as.matrix(runif(length(pfuns$cdf), 0.1, 6)))


              get.quans <- function(index){
                  pfuns$cdf[[index]](testdata[,index])
              }

              q <- sapply(1:ncol(testdata), get.quans)


              get.quans2 <- function(index){
                  pfuns2$cdf[[index]](testdata[,index])
              }

              q2 <- sapply(1:ncol(testdata), get.quans2)

              expect_equal(q, q2)

          })




test_that('T residuals are preserved through normalize and unnormalize',
          {
              # convert the T residuals from their native distribution in each
              # grid cell to N(0,1)
              t1 <- normalize.resids(inputresids = residT, empiricalcdf = tfuns$cdf)$rn

              # Convert T residuals in N(0,1) back to the native distribution
              # in each grid cell
              t2 <- unnormalize.resids(empiricalquant = tfuns$quant, rn = t1)$rnew

              # compare t2 with the original residuals
              maxdiff <- max(abs(t2-residT))
              ftol <- 1.0e-6 * min(abs(residT))
              expect_lte(maxdiff, ftol)
          })


test_that('P residuals are preserved through normalize and unnormalize',
          {
              # convert the T residuals from their native distribution in each
              # grid cell to N(0,1)
              p1 <- normalize.resids(inputresids = residP, empiricalcdf = pfuns$cdf)$rn

              # Convert T residuals in N(0,1) back to the native distribution
              # in each grid cell
              p2 <- unnormalize.resids(empiricalquant = pfuns$quant, rn = p1)$rnew

              # compare p2 with the original residuals
              maxdiff <- max(abs(p2-residP))
              ftol <- 1.0e-6 * min(abs(residP))
              expect_lte(maxdiff, ftol)
          })


test_that('Normal residuals are preserved through unnormalize to T distributions and normalize',
          {
              set.seed(1)
              # create some normally distributed residuals
              residN <- as.matrix(rnorm(100, 0, 1))
              for(j in 1:((ncol(residT)/4)-1)){
                  residN <- cbind(residN, as.matrix(rnorm(100, 0, 1)))
              }

              residN <- cbind(residN, residN, residN, residN)
              expect_equal(ncol(residN), ncol(residT))


              # convert to T distributions
              n1 <- unnormalize.resids(empiricalquant = tfuns$quant, rn = residN)$rnew

              # Convert T residuals in N(0,1) back to the native distribution
              # in each grid cell
              n2 <- normalize.resids(inputresids = n1,
                                     empiricalcdf = tfuns$cdf)$rn

              # compare n2 with the original residuals
              maxdiff <- max(abs(n2-residN))
              ftol <- 1.0e-6
              expect_lte(maxdiff, ftol)
          })


test_that('Normal residuals are preserved through unnormalize to P distributions and normalize',
          {
              set.seed(1)
              # create some normally distributed residuals
              residN <- as.matrix(rnorm(100, 0, 1))
              for(j in 1:((ncol(residP)/4)-1)){
                  residN <- cbind(residN, as.matrix(rnorm(100, 0, 1)))
              }

              residN <- cbind(residN, residN, residN, residN)
              expect_equal(ncol(residN), ncol(residP))


              # convert to T distributions
              n1 <- unnormalize.resids(empiricalquant = pfuns$quant, rn = residN)$rnew

              # Convert T residuals in N(0,1) back to the native distribution
              # in each grid cell
              n2 <- normalize.resids(inputresids = n1,
                                     empiricalcdf = pfuns$cdf)$rn

              # compare n2 with the original residuals
              maxdiff <- max(abs(n2-residN))
              ftol <- 1.0e-6
              expect_lte(maxdiff, ftol)
          })

