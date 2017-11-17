#### Functions for generating field variability

#' Perform EOF analysis on residuals
#'
#' Decompose residuals from the mean field into orthogonal components.  One of
#' these components (conventionally the first component) is always the
#' normalized global mean temperature operator.  This convention isolates the
#' time variability of the global mean in a single components and ensures that
#' the remaining components describe variability that is purely spatial.
#'
#' @section Value:
#'
#' The return value of this function is in the same format as the structure
#' returned by \code{\link[stats]{princomp}}.  The most important fields are:
#'
#' \describe{
#'   \item{rotation}{A matrix [ngrid x N] containing orthogonal basis vectors
#'   for the residuals.  Each column is a basis vector.  The number of basis
#'   vectors, N, will be less than or equal to the number of time slices in the
#'   input.}
#'   \item{x}{A matrix [ntime x N] containing the coordinates of the residuals
#'   in the coordinate system defined by the basis vectors.  Each row is a time
#'   slice.  Thus, \eqn{r(t) = \sum_{i=1}^{N} x[t,i] * rotation[:,i]}.}
#' }
#'
#' @param resids A matrix [ntime x ngrid] of residuals from the mean field
#' @param tgop A column vector [ngrid x 1] containing the global mean operator.
#' This is produced as the \code{tgop} field of the structure returned by
#' \code{\link{readdata}}.
#' @importFrom stats prcomp sd
#' @export
eof_analyze <- function(resids, tgop)
{
    ## pull the residuals from the pattern scaling.  Also, set xh0 to be the
    ## tgop, divided by a normalization factor
    xh0 <- tgop / sqrt(sum(tgop*tgop))

    ## Project resids onto xh0, and subtract that projection from the original
    ## residuals.  This is the same as calculating:
    ##   resids - resids %*% xh0 %*% t(xh0)
    ##
    ## First calculate projection coefficients.  We will need them again later.
    proj0 <- resids %*% xh0

    ## For each residual vector, subtract its projection onto xh0 from the
    ## residual itself
    resids <- resids - proj0 %*% t(xh0)

    ## Now all of the remaining residuals are in the subspace where the global
    ## mean is zero.  We can compute principal components on those, and the PCs
    ## will all have zero global mean.

    ## I'm not sure whether we want to scale to unit variance or not.  The
    ## prototype version says we want to do this because there is higher
    ## variance at high latitudes, and the man page for this function says it is
    ## "advisable", but to me it seems like it adds a lot of effort without much
    ## payoff.  We'll try it sans scaling for now.
    pcs <- prcomp(resids, retx=TRUE, center=FALSE, scale=FALSE)

    ## It's pretty common to get some degenerate basis vectors.  The signature
    ## of these is having very low singular values (stored in sdev).  We check
    ## for this and eliminate them from the basis.
    sdmin <- 1.0e-8 * max(pcs$sdev)
    pckeep <- pcs$sdev > sdmin
    ## drop the components that don't qualify
    pcs$sdev <- pcs$sdev[pckeep]
    pcs$rotation <- pcs$rotation[,pckeep]
    pcs$x <- pcs$x[,pckeep]

    ## next we want to graft the zeroth basis vector onto our EOF results.
    pcnames <- c('PC0', colnames(pcs$rotation))
    pcs$rotation <- cbind(xh0, pcs$rotation)
    colnames(pcs$rotation) <- pcnames
    ## the first column of the 'x' matrix will just be the projection
    ## coefficients that we computed above.
    pcs$x <- cbind(proj0, pcs$x)
    colnames(pcs$x) <- pcnames

    ## We also need to add the sdev of the zeroth components to the sdev vector
    pcs$sdev <- c(sd(proj0), pcs$sdev)

    ## If we decide to use scaling in the future, we will need to graft onto
    ## that one as well.  For now it just contains the scalar FALSE, so nothing
    ## needs to be done.

    pcs
}

#' Reconstitute temperature fields
#'
#' Reconstitute temperature residual fields from the basis vectors and a
#' matrix of coefficients.  If meanfield is supplied, add the residual field
#' to the mean field to get a reconstituted temperature field.
#'
#' @param basis A matrix [ngrid x N] of basis vectors, one vector per column.
#' @param bcoord A matrix [ntime x N] of basis coordinates, one time step in
#' each row, with coordinates across the columns.
#' @param meanfield An optional matrix [ntime x ngrid] of mean field values, one
#' time slice per row. If supplied, it will be added to the residuals calculated
#' from \code{basis} and \code{bcoord}.
#' @return A matrix [ntime x ngrid] of residual values (or temperature values,
#' if \code{meanfield} was supplied).
#' @export
reconst_fields <- function(basis, bcoord, meanfield=NULL)
{
    resid <- bcoord %*% t(basis)
    if(is.null(meanfield))
        resid
    else
        meanfield + resid
}


#' Make time series with specified autocorrelation properties
#'
#' The time series produced will have the power spectrum given in the columns of
#' \code{Fxmag}.
#'
#' It usually will not be necesary to provide the \code{phase} input, since the
#' intent is usually to generate random phases; however, it is occasionally
#' useful to the input ESM data, or some other data set for which the phases are
#' known.  To make this sort of exercise easier, the function will accept a
#' matrix containing phases for all frequencies, including the negative ones,
#' despite the fact that only the phases for positive frequencies are needed.
#'
#' @section ToDo:
#'
#' It should be possible to supply a partial matrix of phases, with some values
#' specified and the rest set to NA.  This would allow us, for example, to
#' specify phases for the global mean component while the other components are
#' randomized.
#'
#' @param Fxmag A matrix [ntime x N] in which each column is the absolute value
#' of the fourier transform of a (real-valued) time series.  This should include
#' the entire Fourier transform, both positive and negative frequencies.
#' @param phase An optional matrix of phases.  These should be uniform random
#' deviates in [0,2pi).  We only need phases for the positive frequencies.  For
#' an even number of samples, that's N/2.  For an odd number, that's (N-1)/2.
#' Excess phases will be ignored, so it's ok to just dimension phase the same as
#' Fxmag.  If phase is omitted, then it will be generated using runif.
#' @param complexout The inverse FFT produces complex-valued results; however,
#' the imaginary parts should all be zero.  By default we return Re(rslt), but
#' setting this flag causes the result to be left in complex form.
#' @importFrom stats mvfft runif
#' @export
mkcorrts <- function(Fxmag, phase=NULL, complexout=FALSE)
{
    N <- nrow(Fxmag)
    M <- ncol(Fxmag)
    if(N %% 2 == 0) {
        ## even number of samples
        Nplus <- N/2 - 1                # number of +ve frequencies
        plusrows <- seq(2, length.out=Nplus) # rows containing +ve frequencies
        Nc <- N/2 + 1                   # row containing the Nyquist frequency
        Nphase <- Nc                    # number of phases
        minusrows <- seq(Nc+1, length.out=Nplus) # rows containing -ve frequencies
        ## total rows:  N/2-1 (+ve) + N/2-1 (-ve) + 1 (zero) + 1 (Nyquist) == N
    }
    else {
        ## odd number of rows
        Nplus <- (N-1)/2                # number of +ve frequencies
        Nphase <- Nplus + 1
        plusrows <- seq(2, length.out=Nplus) # rows containing +ve freqs
        minusrows <- seq(Nplus+2, length.out=Nplus) # rows containing -ve freqs
        ## total rows: (N-1)/2 (+ve) + (N-1)/2 (-ve) + 1 (zero) == N
    }

    if(is.null(phase)) {
        phase <- matrix(2*pi*runif(Nphase*M), ncol=M)
    }
    else if(nrow(phase) > Nphase) {
        ## We need Nplus phases for each time series.  We also need to look at
        ## the phases for the zero and (if present) Nyquist components.
        ## However, one thing we might do occasionally is to copy the phases
        ## from another time series.  If the number of rows is greater than
        ## Nplus, assume that's what we are doing and copy just the relevant
        ## phases.
        phase <- phase[1:Nphase,]
    }

    ## Assign the magnitudes to the output array.  Force them to be complex.
    Fxout <- Fxmag * 1+0i
    ## The zero and Nyquist components don't get a phase as such, but they can
    ## be either positive or negative.  If the phase is between pi/2 and 3pi/2,
    ## we make them negative; otherwise they are positive.
    Fxout[1,] = Fxout[1,]*dplyr::if_else(abs(phase[1,]-pi) < pi/2, -1.0, 1.0)
    if(N %% 2 == 0) {
        ## The Nyquist frequency is only present if N is even.
        Fxout[Nc,] = Fxout[Nc,]*dplyr::if_else(abs(phase[Nc,]-pi) < pi/2, -1.0,
             1.0)
    }

    ## plus rows get new phases
    Fxout[plusrows,] <- Fxout[plusrows,] * (cos(phase[plusrows,]) +
                                                  1i*sin(phase[plusrows,]))
    ## minus rows have to be the conjugates of the plus rows
    Fxout[minusrows,] <- Conj(Fxout[rev(plusrows),])

    xout <- mvfft(Fxout, inverse=TRUE)
    if(complexout)
        xout / N
    else
        Re(xout) / N
}
