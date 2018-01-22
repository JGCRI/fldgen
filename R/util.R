################################################################
#### Utility Functions
################################################################

## In a 1D FFT grid, given the i coordinate of a positive frequency f,
## find the i' coordinate corresponding to -f.  As with all R arrays, the
## coordinates are unit-offset.
##
## i == the input coordinates.  A vector of values may be given
## N == the length of the FFT dimension.
find_minusf_coord <- function(i, N)
{
    dplyr::if_else(i==1, 1, N-i+2)
}

## Find the number of independent phases for an FFT of length N
nphase <- function(N)
{
    N <- as.numeric(N)   # N-1 will be numeric, even if N is integer, so make sure types will agree.
    dplyr::if_else(N %% 2 == 0, N, N-1) / 2 + 1
}


## Apply a function (element-wise) to each column of a matrix A and a vector v.
## E.g., broadcast_apply_col(A, v, `*`) would multiply each column of A
## element-wise with v and return the resulting matrix.
broadcast_apply_col<- function(A, v, func)
{
    B <- A
    for(j in seq(1,ncol(A))) {
        B[,j] <- func(B[,j], v)
    }
    B
}


## Generate the coefficients in the equations that determine the phase
## condition.  These are given as equation XX in the paper.
##
## The result will be a matrix A[Nf x Nb], where Nf is the number of
## nonnegative frequencies in the Fourier transform and Nb is the number of
## basis functions (i.e., EOFs).  A[k,j] = (a_k^i a_k^j / f_k) * x_k^j, where i
## is the index of the basis function used for the reference (conventionally 2,
## which is the first spatial basis function.  Don't use the time-like basis
## function, i=1, since it will generally have small coefficients).
##
## Fx: The Fourier transform object created by train_fldgen
## i:  The index of the basis function to use for the reference
phase_eqn_coef <- function(Fx, i=2)
{
    Nt <- nrow(Fx$mag)
    Nc <- nphase(Nt)
    krows <- 1:Nc

    ## TODO: We only want the coefficients here, not the x values.
    ## Find the x values, x_k^j = cos(phase(k,j)-phase(k,i))
    x <- cos(broadcast_apply_col(Fx$phase[krows,], Fx$phase[krows,i], `-`))
    ## The j,k dependence is a_k^j * x_k^j.  The a are in Fx$mag.
    A <- Fx$mag[krows,] * x

    ## The other terms have k-dependence only.  f_k is 2pi *k/Nt; however, there
    ## is a factor of pi that cancels, so we're going to use f_k = 2*k/Nt.
    ## NB: These k values are numbered starting at 0.
    k <- seq(0,Nc-1)
    fk <- 2*k/Nt
    ## Also, k=0 is special.  Its coefficient is a fixed value of 2.  Since fk
    ## will appear in the denominator, fk(0) = 1/2.  (k=0 is in array slot 1)
    fk[1] <- 0.5

    ## There is a factor of 1/Nt^3 in the equations that doesn't affect anything
    ## (because these equations will all be set equal to zero and solved), but
    ## including it makes the magnitudes of the coefficients more reasonable.
    ## Since fk goes in the denominator, we multiply by Nt^3
    fk <- fk*(Nt^3)

    ## Now, each column of A needs to be multiplied by a_k^i/fk.  This will be our
    ## return value.
    broadcast_apply_col(A, Fx$mag[krows,i]/fk, `*`)
}
