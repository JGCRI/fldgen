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


## Convert arbitrary phases into the nearest self-symmetric (i.e., 0 or pi)
## phases.
## NB: phases must be in [0,2pi] for this to work properly
phsym <- function(phase)
{
    dplyr::if_else(abs(phase-pi) < pi/2, pi, 0.0)
}
