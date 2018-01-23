################################################################
#### Functions to hand test the steps in the derivation of the formulae used in
#### the method.
####
#### These functions are neither needed by nor included in the package.  They
#### will stick around in extdata for as long as they seem useful.
################################################################

## compute phi(t) from the vector of Fourier coefficients (this is essentially
## an inverse Fourier transform).  x is not used; it's just there for comparison
## during the calc.
mkphi1 <- function(a,x)
{
    Nt <- length(a)
    t <- seq(0, Nt-1)

    fourterm <- function(t) {
        k <- seq(0,Nt-1)
        expterm <- 2*pi*k/Nt*t * 1i
        summand <- a * exp(expterm)
        sum(summand)
    }

    phit <- sapply(t, fourterm)
    phit/Nt
}


get_plusf_idx <- function(Nt, Nf)
{
    if(Nt%%2 == 0) {
        ## There is a kc term
        plusf <- seq(2,Nf-1)
    }
    else {
        ## No kc term
        plusf <- seq(2,Nf)
    }
}


## compute phi(t) using the equivalent cosine series.
mkphi2 <- function(acplx, x)
{
    Nt <- length(acplx)
    Nf <- fldgen:::nphase(Nt)
    a <- abs(acplx[1:Nf])
    theta <- atan2(Im(acplx[1:Nf]), Re(acplx[1:Nf]))
    t <- seq(0, Nt-1)

    tterm <- function(t) {
        k <- seq(0, Nf-1)
        cosarg <- 2*pi*k/Nt*t + theta[k+1]
        costerm <- a * cos(cosarg)

        ## apply a factor of two to all rows that are not k==0 or k==kc
        plusf <- get_plusf_idx(Nt, Nf)
        costerm[plusf] <- 2.0*costerm[plusf]

        ## Return the sum of these.
        sum(costerm)
    }

    phit <- sapply(t, tterm)
    phit/Nt
}

## compute sum(phi_i(t) * phi_j(t)) over t using orthogonal properties of the
## cos series
phiprod1 <- function(ai, aj, thetai, thetaj)
{
    Nt <- length(ai)
    Nf <- fldgen:::nphase(Nt)
    t <- seq(0, Nt-1)

    ai <- ai[1:Nf]
    aj <- aj[1:Nf]
    thetai <- thetai[1:Nf]
    thetaj <- thetaj[1:Nf]

    tterm <- function(t) {
        k <- seq(0, Nf-1)
        cosft <- 2*pi*k/Nt * t
        cosargi <- cosft + thetai
        cosargj <- cosft + thetaj

        costerms <- ai * aj * cos(cosargi) * cos(cosargj)

        ## apply factor of 4 to all rows that are not k==0 or k==kc
        plusf <- get_plusf_idx(Nt, Nf)
        costerms[plusf] <- 4.0 * costerms[plusf]

        ## Return the sum of these
        sum(costerms)
    }

    tterms <- sapply(t, tterm)
    sum(tterms) / Nt^3
}


## compute sum(phi_i(t) * phi_j(t)) over t using integral properties
phiprod2 <- function(fxmag, fxphase, i, j)
{
    Nt <- nrow(fxmag)
    Nf <- fldgen:::nphase(Nt)
    k <- seq(0,Nf-1)

    ai <- fxmag[k+1,i]
    aj <- fxmag[k+1,j]

    thetai <- fxphase[k+1,i]
    thetaj <- fxphase[k+1,j]
    deltaij <- thetaj - thetai

    ## coefficient for rows that aren't k==0 or k==kc is 4
    kterms <- ai * aj * cos(deltaij)
    plusf <- get_plusf_idx(Nt, Nf)
    kterms[plusf] <- 2.0 * kterms[plusf]

    list(Aij=kterms, detaij=deltaij, ai=ai, aj=aj)
}
