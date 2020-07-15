library('fldgen')

train_models <- function(models, tasvar='tasAdjust', prvar='prAdjust',
                         datadir='./training-data') {

    ## The following would give you the complete set of models:
    ## models <- c('GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC5')

    for (model in models) {
        datafiles <- list.files(path=datadir, pattern=model, full.names=TRUE)
        cat('Processing model ', model, '  datafiles:\n', paste(datafiles, collapse='\n'),'\n')
        emu <- trainTP(datafiles, tvarname=tasvar, pvarname=prvar)
        emu$griddataP$vardata_raw <- NULL
        outfilename <- paste0('fldgen-',model, '.rds')
        coord <- emu$griddataT$coord
        coord[67382, ] <- c(-49.75, 178.75)
        emu$griddataT$coord <- coord
        emu$griddataP$coord <- coord

        saveRDS(emu, outfilename)


        emulator <- emulator_reducer(emu)
        outfilename <- paste0('fldgen-',model, '_reducedEmulator.rds')
        saveRDS(reducedEmulator, outfilename)
    }
}


