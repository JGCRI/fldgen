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
#' \code{\link{read.temperatures}}.
#' @importFrom stats prcomp sd
#' @export
#' @keywords internal
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
#' For a mean field generated using the default linear pattern scaling, one can
#' get the mean field by running \code{\link{pscl_apply}} on the \code{pscl}
#' member of a \code{\link[=train]{fldgen}} object.
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
#'
#' The \code{phase} argument allows a user to provide fixed phases to be used in
#' the field construction.  This \emph{usually} is not desirable because an
#' arbitrary set of phases will not generally produce properly uncorrelated
#' outputs.  However, it may occasionally be useful to reproduce the input ESM
#' data, or some other data set for which the phases are known.  To make this
#' sort of exercise easier, the function will accept a matrix containing phases
#' for all frequencies, including the negative ones, despite the fact that only
#' the phases for positive frequencies are needed.
#'
#' @section Methods:
#'
#' Method I assignes uniform random values in [0,2pi) to all phases.  This
#' method is simple and fast, but it produces a slight bias in the covariance of
#' the projection coefficients for the EOFs.  (The covariance of the projection
#' coefficients for two distinct EOFs should be zero.)  Although this produces
#' biases in the grid cell statistics, in many cases they are small enough to be
#' ignorable.
#'
#' Method II is more conservative in its phase randomization.  It uses the
#' Fourier transform of one of the input data sets as a prototype and for each
#' frequency component it generates a random phase shift and adds it to that
#' frequency component for all of the EOFs.  This guarantees that the covariance
#' will be zero, but at the cost of making the space of possible realizations
#' much lower dimension, roughly 1/2 N instead of 1/2 M*N (where M is the number
#' of EOFs and N is the number of basis functions).
#'
#' @section ToDo:
#'
#' It should be possible to supply a partial matrix of phases, with some values
#' specified and the rest set to NA.  This would allow us, for example, to
#' specify phases for the global mean component (EOF-0) while the other
#' components are randomized.
#'
#' @param fldgen A \code{\link[=train]{fldgen}} object.
#' @param phase An optional matrix of phases.  See notes in the Details section.
#' @param complexout The inverse FFT produces complex-valued results; however,
#' the imaginary parts should all be zero.  By default we return Re(rslt), but
#' setting this flag causes the result to be left in complex form.  This is
#' mostly useful for testing.
#' @param method Integer specifying method 1 or method 2 for generating the
#' phases.  Ignored if \code{phase} is set.
#' @importFrom stats mvfft runif
#' @export
mkcorrts <- function(fldgen, phase=NULL, method=1, complexout=FALSE)
{
    Fxmag <- fldgen$fx$mag

    Nt <- nrow(Fxmag)
    Nf <- nphase(Nt)
    M <- ncol(Fxmag)

    plusrows <- 1:Nf
    minusrows <- find_minusf_coord(plusrows,Nt)

    if(is.null(phase)) {
        if(method==2) {
            ## generate phase shifts for each frequency component
            dtheta <- runif(Nf, 0, 2*pi)
            dtheta[1] <- phsym(dtheta[1])
            if(Nt %% 2 == 0) {
                dtheta[Nf] <- phsym(dtheta[Nf])
            }

            ## pick one of the phase templates at random and apply the phase
            ## shifts to the plusrows.
            n <- length(fldgen$fx$phases)
            phase <- fldgen$fx$phases[[sample.int(n, 1)]][plusrows,]
            ## add thetai to each column to get the final
            phase <- broadcast_apply_col(phase, dtheta, `+`)
        }
        else if(method==1) {
            phase <- matrix(runif(Nf*M, 0, 2*pi), ncol=M)
            phase[1,] <- phsym(phase[1,]) # symmetrize f==0 mode.
            if(Nt %% 2 == 0) {
                ## For even N the Nyquist mode is present and must be symmetrized.
                phase[Nf,] <- phsym(phase[Nf,])
            }
        }
        else {
            stop('Unknown phase method: ', method)
        }
    }
    else {
        phase <- phase[1:Nf,]
    }

    ## Assign the magnitudes to the output array.  Force them to be complex.
    Fxout <- Fxmag * 1+0i
    ## The zero and Nyquist components don't get a phase as such, but they can
    ## be either positive or negative.  If the phase is between pi/2 and 3pi/2,
    ## we make them negative; otherwise they are positive.

    ## plus rows get new phases
    Fxout[plusrows,] <- Fxout[plusrows,] * (cos(phase[plusrows,]) +
                                                  1i*sin(phase[plusrows,]))
    ## minus rows have to be the conjugates of the plus rows
    Fxout[minusrows,] <- Conj(Fxout[plusrows,])

    xout <- mvfft(Fxout, inverse=TRUE)
    if(complexout)
        xout / Nt
    else
        Re(xout) / Nt
}

#' Estimate the power spectral density from a list of ESM runs
#'
#' This is a simplistic estimate of the PSD; it's just the mean of the
#' squared-magnitude of the DFT, averaged across the runs in the dataset.  There
#' are more sophisticated things you could do, but a crude estimate is more than
#' good enough for what we are trying to do here.
#'
#' The matrix returned contains the PSD estimates for \emph{all}
#' frequencies (including negative frequencies).  Each column corresponds to one of the
#' principal components (PC0, PC1, ... PCN), and each row corresponds to a
#' frequency bin ($0, 1/N_t, 2/N_t, \ldots, 1/(2N_t), \ldots, -2/N_t, -1/N_t$,
#' in units of yr$^{-1}$).  The \emph{square root} of this matrix should be passed
#' as the first argument to \code{\link{mkcorrts}}.
#'
#' When doing the averaging we assume that all of the input time series have the
#' same lengths, so that all of the frequency bins correspond.  Passing series
#' of unequal length is an error.
#'
#' @param prcomp_l A list of principal components structures returned by
#' eof_analyze, or a single such structure.
#' @return A matrix [ntime, numPC] of PSD estimates (see details).
#' @export
#' @keywords internal
psdest <- function(prcomp_l)
{
    if(inherits(prcomp_l, 'prcomp')) {
        ## User passed a single griddata instead of a list.  Wrap it in a list
        ## and continue.
        prcomp_l <- list(prcomp_l)
    }

    ## Check to see that all time series are the same length.
    lens <- sapply(prcomp_l, function(pc) {nrow(pc$x)})
    if(any(lens != lens[1])) {
        stop('All input prcomp structures must have equal length')
    }

    ffts <- lapply(prcomp_l,
                   function(pc) {
                       fx <- stats::mvfft(pc$x)
                       list(Mod(fx)^2,
                            atan2(Im(fx), Re(fx)))
                   })

    psds <- lapply(ffts, function(x) {x[[1]]})
    phases <- lapply(ffts, function(x) {x[[2]]})

    list(psd=Reduce(`+`, psds) / length(psds),
         phases=phases)
}


#' Split an EOF structure made of multiple time series into a list
#'
#' Each element of the output list will be an EOF structure for a single time
#' series (i.e., a single ESM run).  This will \emph{not} be the same as the
#' structure you would have gotten had you run the single ESM run through the
#' process by itself.  It will have more, and possibly different principal
#' components because the EOF analysis was done in conjunction with the other
#' ESM runs.
#'
#' @param reof Structure of residual EOFs, as returned from
#' \code{\link{eof_analyze}}.
#' @param griddata Merged griddata structure used to perform the mean field and
#' EOF analyses
#' @return List of residual EOF structures, with one element for each input ESM
#' run.
#' @export
#' @keywords internal
split_eof <- function(reof, griddata)
{
    lapply(griddata$tags,
           function(tag) {
               rr <- reof
               strt <- tag[1]
               end <- tag[2]
               rr$x <- rr$x[strt:end,]
               rr
           })
}



#' Generate the coefficients in the equations that determine the phase constraints.
#'
#' The result will be a structure containing:
#' \describe{
#'   \item{A}{A matrix [Nf x Nb], where Nf is the number of
#' nonnegative frequencies in the Fourier transform and Nb is the number of
#' basis functions (i.e., EOFs).  A[k,j] = (a_k^i a_k^j), where i
#' is the index of the basis function used for the reference (conventionally 2,
#' which is the first spatial basis function.  Don't use the time-like basis
#' function, i=1, since it will generally have small coefficients).}
#'   \item{i}{The index of the reference basis function.}
#' }
#'
#' These equations are given as equation XX in the paper.
#'
#' @param Fx The Fourier transform of the EOF projection coefficients.
#' @param i  The index of the basis function to use for the reference.  Don't use
#' EOF-0 (i=1).
#' @export
#' @keywords internal
phase_eqn_coef <- function(Fx, i=2)
{
    Fxmag <- abs(Fx)
    Nt <- nrow(Fxmag)
    Nf <- nphase(Nt)
    krows <- 1:Nf

    ## The j,k dependence is a_k^j.  The a's are in Fxmag.
    A <- Fxmag[krows,]

    ## Now, each column of A needs to be multiplied by a_k^i.
    A <- broadcast_apply_col(A, Fxmag[krows,i], `*`)

    ## The rows that are neither f=0 nor f=fc need to be multiplied by 2.  There
    ## may or may not be a row for fc, depending on the parity of Nt
    if(Nt%%2 == 0) {
        ## even N => there is an fc
        plusrows <- seq(2, Nf-1)
    }
    else {
        ## odd N => no row for fc
        plusrows <- seq(2, Nf)
    }
    A[plusrows,] <- 2.0 * A[plusrows,]

    list(A=A/Nt^2, i=i)
}
