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

## Generate random phases that satisfy the zero cross-correlation condition.
##
## For each j > 0, j != A$i, select a random set of deltaij such that
##    sum( A$A[j] * cos(deltaij) ) == 0.
## We do this by selecting N-1 of the cos(deltaij) randomly and solving for the
## remaining one.  If the resulting cos(deltaij) is invalid (i.e., because it is
## not in [-1,1]), we make another random choice and try again.
##
## Astruct: the structure returned from phase_eqn_coef
## Nt:  the total number of time steps in the data
gen_random_delta <- function(Astruct, Nt)
{
    ## Find the x values, x_k^j = cos(phase(k,j)-phase(k,i))
    ## x <- cos(broadcast_apply_col(Fx$phase[krows,], Fx$phase[krows,i], `-`))

    A <- Astruct$A
    i <- Astruct$i

    jvals <- seq(1,ncol(A))
    Nf <- nrow(A)

    ## perform the calculation for a single j
    delta1j <- function(j) {
        if(j==1) {
            ## The first component (EOF-0) is special and shouldn't be used with
            ## this procedure.  We'll fill the column with zeros so that the
            ## other column lines up right.
            return(rep(0,Nf))
        }
        if(j==i) {
            ## For j==i, the deltaij are zero, by definition
            message('j= ', j, ':  returning deltaij = 0 because i==j.')
            return(rep(0,Nf))
        }

        ## Pick an x to solve for.  It can't be the first one (f==0), nor can it
        ## be the last (f==fc) because of symmetry conditions.  (Techincally, if
        ## Nt is odd, it *could* be the last one because that isn't fc, but it's
        ## simpler just to exclude the last frequency, and it doesn't really
        ## cost us anything.)
        ##
        ## It's arbitrary which one we solve for, so pick the one with the
        ## largest coefficient.
        aij <- A[,j]
        kpivot <- which.max(aij[-c(1,Nf)]) + 1 # +1 b/c we cut off the first element
        apivot <- aij[kpivot]

        message('j= ', j,': solving for k= ', kpivot)

        iter <- 0
        ITMAX <- 1000
        while (iter < ITMAX) {
            iter <- iter+1
            ## Choose a random vector of x values (where x==cos(deltaij))
            x <- runif(Nf, -1, 1)
            ## f=0 must have x=1 or x=-1
            if(x[1] < 0)
                x[1] <- -1
            else
                x[1] <- 1

            ## If Nt is even, then x[Nf] must also be +/- 1
            if(Nt%%2 == 0) {
                if(x[Nf] < 0)
                    x[Nf] <- -1
                else
                    x[Nf] <- 1
            }

            ## zero out the x that we are solving for
            x[kpivot] <- 0

            ## Now, solve aij[p] * x[p] = -sum(aij*x)
            ## We can include the x that we are solving for in the sum because
            ## we just set it to zero.
            x[kpivot] <- -sum(aij*x) / aij[kpivot]

            ## Check to see if x is a valid cosine
            if(x[kpivot] >= -1 && x[kpivot] <= 1) {
                ## success
                message('    j= ', j,': Found result after ', iter, ' iterations.')
                break
            }
        }

        if(iter > ITMAX) {
            ## Didn't find an acceptable solution in the allowed number of
            ## iterations.
            warning('    j= ', j,
                    ': ITMAX reached with no result.  Setting result to best guess.')
            ## return the best approximation we can get.
            if(x[kpivot] < 0)
                x[kpivot] <- -1
            else
                x[kpivot] <- 1
        }
        ## convert back to phase angle differences
        acos(x)
    }

    ## do it for all j
    sapply(jvals, delta1j)
}
