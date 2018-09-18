#' Using a trained emulator, generate new fields of residuals.
#'
#' This function takes in a trained emulator - a structure of class
#' \code{fldgen}.
#' This structure contains everything the emulator has learned about the model,
#' and is used to generate new fields of residuals.
#'
#' First, new residuals are generated. The distribution of residuals over time
#' for each variable in each grid cell is N(0,1). These generated residuals are
#' then mapped back to the native distribution for each variable's residuals
#' in each grid cell.
#'
#'
#' @param emulator A trained \code{fldgen} temperature precipitation joint
#' emulator.
#' @param ngen The number of new fields to generate
#' @param method The algorithm used to generate new EOF coefficients, defaults
#' to 1.
#' @return A list of new residual fields, each entry in the list is a new
#' realization, a matrix that is [Nyears x 2 * Ngrid]; the first 1:Ngrid cols
#' are the temperature residuals and columns (Ngrid + 1):(2*Ngrid) are the
#' precipitation residuals.
#' @export
generate.TP.resids <- function(emulator, ngen, method = 1)
{
    Ngrid <- ncol(emulator$meanfldT$r)

    # Generate the new residual fields in the normal space.
    newgrids <- lapply(1:ngen,
                       function(x) {
                           ## It takes the full fldgen object, but doesn't appear
                           ## to use anything other than phase and magnitude,
                           ## which by design contains both T and P. May be able
                           ## to improve implementation in future.
                           # This uses a RNG, it will be different every call:
                           bcoord <- mkcorrts(emulator, method=method)
                           ## Because the generated residuals are generated to
                           ## follow N(0, 1), we can't add back to the meanfield
                           ## as is. They have to be transformed back to the
                           ## distribution of T or P residuals that are
                           ## empirically characterized for each grid cell in
                           ## emulator$tfuns$quant, emulator$pfuns$quant.
                           reconst_fields(emulator$reof$rotation, bcoord)
                       })

    residgrids <- lapply(newgrids, function(g) {
        g[, 1:Ngrid] <- unnormalize.resids(empiricalquant = emulator$tfuns$quant,
                                           rn = g[ ,1:Ngrid])$rnew

        g[, (Ngrid+1):(2*Ngrid)] <- unnormalize.resids(empiricalquant = emulator$pfuns$quant,
                                                       rn = g[ , (Ngrid+1):(2*Ngrid)])$rnew

        return(g)}
    )



    ### ADD THE NA COLUMNS BACK FOR ISIMIP DATA

    ### return the generated residuals
    return(residgrids)
}
