% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LIK.density.R, R/LSCV.density.R
\name{LIK.density}
\alias{LIK.density}
\alias{LSCV.density}
\title{Cross-validation bandwidths for spatial kernel density estimates}
\usage{
LIK.density(
  pp,
  hlim = NULL,
  hseq = NULL,
  resolution = 64,
  edge = TRUE,
  auto.optim = TRUE,
  type = c("fixed", "adaptive"),
  seqres = 30,
  parallelise = NULL,
  zero.action = 0,
  verbose = TRUE,
  ...
)

LSCV.density(
  pp,
  hlim = NULL,
  hseq = NULL,
  resolution = 64,
  edge = TRUE,
  auto.optim = TRUE,
  type = c("fixed", "adaptive"),
  seqres = 30,
  parallelise = NULL,
  zero.action = 0,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{pp}{An object of class \code{\link[spatstat.geom]{ppp}} giving the observed
2D data to be smoothed.}

\item{hlim}{An optional vector of length 2 giving the limits of the
optimisation routine with respect to the bandwidth. If unspecified, the
function attempts to choose this automatically.}

\item{hseq}{An optional increasing sequence of bandwidth values at which to
manually evaluate the optimisation criterion. Used only in the case
\code{(!auto.optim && is.null(hlim))}.}

\item{resolution}{Spatial grid size; the optimisation will be based on a
[\code{resolution} \eqn{\times}{x} \code{resolution}] density estimate.}

\item{edge}{Logical value indicating whether to edge-correct the density
estimates used.}

\item{auto.optim}{Logical value indicating whether to automate the numerical
optimisation using \code{\link{optimise}}. If \code{FALSE}, the optimisation
criterion is evaluated over \code{hseq} (if supplied), or over a seqence of
values controlled by \code{hlim} and \code{seqres}.}

\item{type}{A character string; \code{"fixed"} (default) performs classical leave-one-out
cross-validation for the fixed-bandwidth estimator. Alternatively, \code{"adaptive"} utilises
multiscale adaptive kernel estimation (Davies & Baddeley, 2018) to run the cross-validation
in an effort to find a suitable global bandwidth for the adaptive estimator. Note that data points are not `left out' of
the pilot density estimate when using this option (this capability is currently in development). See also the entry for \code{...}.}

\item{seqres}{Optional resolution of an increasing sequence of bandwidth
values. Only used if \code{(!auto.optim && is.null(hseq))}.}

\item{parallelise}{Numeric argument to invoke parallel processing, giving
the number of CPU cores to use when \code{!auto.optim} \bold{and} \code{type = "fixed"}. Experimental. Test
your system first using \code{parallel::detectCores()} to identify the
number of cores available to you.}

\item{zero.action}{A numeric integer, either \code{-1}, \code{0} (default), \code{1} or \code{2} controlling how the function should behave in response to numerical errors at very small bandwidths, when such a bandwidth results in one or more zero or negative density values during the leave-one-out computations. See `Details'.}

\item{verbose}{Logical value indicating whether to provide function progress
commentary.}

\item{...}{Additional arguments controlling pilot density estimation and multi-scale bandwidth-axis
resolution when \code{type = "adaptive"}. Relevant arguments are \code{hp}, \code{pilot.density},
\code{gamma.scale}, and \code{trim} from \code{\link{bivariate.density}}; and \code{dimz} from 
\code{\link{multiscale.density}}. If \code{hp} is missing and required, the function makes a (possibly recursive)
call to \code{LSCV.density} to set this using fixed-bandwidth LSCV. The remaining defaults are \code{pilot.density = pp},
\code{gamma.scale = "geometric"}, \code{trim = 5}, and \code{dimz = resolution}.}
}
\value{
A single numeric value of the estimated bandwidth (if
  \code{auto.optim = TRUE}). Otherwise, a \eqn{[}\code{seqres} \eqn{x} 2\eqn{]} matrix 
  giving the bandwidth sequence and corresponding CV
  function value.
}
\description{
Isotropic fixed or global (for adaptive) bandwidth selection for standalone 2D density/intensity
based on either unbiased least squares cross-validation (LSCV) or likelihood (LIK) cross-validation.
}
\details{
This function implements the bivariate, edge-corrected versions of fixed-bandwidth least squares cross-validation and likelihood cross-validation
as outlined in Sections 3.4.3 and 3.4.4 of Silverman (1986) in order to select an optimal fixed smoothing bandwidth. With \code{type = "adaptive"} it may also be used to select the global bandwidth
for adaptive kernel density estimates, making use of multi-scale estimation (Davies and Baddeley, 2018) via \code{\link{multiscale.density}}.
Note that for computational reasons, the leave-one-out procedure is not performed on the pilot density in the adaptive setting; it
is only performed on the final stage estimate. Current development efforts include extending this functionality, see \code{\link{SLIK.adapt}}. See also `Warning' below.

Where \code{LSCV.density} is based on minimisation of an unbiased estimate of the mean integrated squared error (MISE) of the density, \code{LIK.density} is based on
maximisation of the cross-validated leave-one-out average of the log-likelihood of the density estimate with respect to \eqn{h}.

In both functions, the argument \code{zero.action} can be used to control the level of severity in response to small bandwidths that result (due to numerical error) in at least one density value being zero or less.
When \code{zero.action = -1}, the function strictly forbids bandwidths that would result in one or more \emph{pixel} values of a kernel estimate of the original data (i.e. anything over the whole region) being zero or less---this is the most restrictive truncation. With \code{zero.action = 0} (default), the function
automatically forbids bandwidths that yield erroneous values at the leave-one-out data point locations only. With \code{zero.action = 1}, the minimum machine value (see \code{.Machine$double.xmin} at the prompt) is
used to replace these individual leave-one-out values. When \code{zero.action = 2}, the minimum value of the valid (greater than zero) leave-one-out values is used to replace any erroneous leave-one-out values.
}
\section{Warning}{
 Leave-one-out CV for bandwidth selection in kernel
density estimation is notoriously unstable in practice and has a tendency to
produce rather small bandwidths, particularly for spatial data. Satisfactory bandwidths are not guaranteed
for every application; \code{zero.action} can curb adverse numeric effects for very small bandwidths during the optimisation procedures. This method can also be computationally expensive for
large data sets and fine evaluation grid resolutions. The user may also need to
experiment with adjusting \code{hlim} to find a suitable minimum.
}

\examples{

data(pbc)
pbccas <- split(pbc)$case

LIK.density(pbccas)
LSCV.density(pbccas)

\donttest{
#* FIXED 

# custom limits
LIK.density(pbccas,hlim=c(0.01,4))
LSCV.density(pbccas,hlim=c(0.01,4))

# disable edge correction
LIK.density(pbccas,hlim=c(0.01,4),edge=FALSE)
LSCV.density(pbccas,hlim=c(0.01,4),edge=FALSE)

# obtain objective function
hcv <- LIK.density(pbccas,hlim=c(0.01,4),auto.optim=FALSE)
plot(hcv);abline(v=hcv[which.max(hcv[,2]),1],lty=2,col=2)

#* ADAPTIVE
LIK.density(pbccas,type="adaptive")
LSCV.density(pbccas,type="adaptive")
 
# change pilot bandwidth used
LIK.density(pbccas,type="adaptive",hp=2)
LSCV.density(pbccas,type="adaptive",hp=2)
} 

}
\references{
Davies, T.M. and Baddeley A. (2018), Fast computation of
spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.

Silverman, B.W. (1986), \emph{Density Estimation for Statistics
and Data Analysis}, Chapman & Hall, New York.

Wand, M.P. and Jones,
C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
}
\seealso{
\code{\link{SLIK.adapt}} and functions for bandwidth selection in package
  \code{spatstat}: \code{\link[spatstat.explore]{bw.diggle}};
  \code{\link[spatstat.explore]{bw.ppl}}; \code{\link[spatstat.explore]{bw.scott}};
  \code{\link[spatstat.explore]{bw.frac}}.
}
\author{
T. M. Davies
}
