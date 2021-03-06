% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train.R
\name{fldgen_object}
\alias{fldgen_object}
\title{Create a \code{fldgen} object from constituent parts}
\usage{
fldgen_object(griddata, tgav, pscl, reof, fxmag, fxphase, infiles)
}
\arguments{
\item{griddata}{An object returned from \code{\link{read.temperatures}} or
\code{\link{concatGrids}}.}

\item{tgav}{Global mean temperature for the grids in griddata.}

\item{pscl}{Object returned from \code{\link{pscl_analyze}}.}

\item{reof}{Object returned from \code{\link{eof_analyze}}.}

\item{fxmag}{The magnitude of the Fourier transform of the EOF projection
coefficients.  This should be a matrix \link{Ntime x NEOF}.  If using
\code{\link{psdest}}, note this is the square root of the power spectral
density returned by that function.}

\item{fxphase}{List of matrices \link{Ntime x NEOF} of phases of the Fourier
transform.  There should be one element in the list for each input ESM run.}

\item{infiles}{Names of input files used to construct the data.}
}
\description{
Normally this code will be called internally, but it is available for
external use, so that data created using the old interface can be converted.
}
\details{
If there is only a single ESM run in the input data set, then \code{fxphase}
can be passed as a matrix of phases (instead of a list of matrices).  It will
be converted into a list automatically.
}
\keyword{internal}
