\name{LSCV.risk}
\Rdversion{1.1}
\alias{LSCV.risk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leave-one-out least-squares cross-validation (LSCV) bandwidths for the relative risk function
}
\description{
Attempts to estimate a jointly optimal, common case-control fixed bandwidth for use in the kernel-smoothed relative risk function via leave-one-out least-squares cross-validation (LSCV). The user can choose between two methods described in Kelsall and Diggle (1995a;b) and Hazelton (2008).
}
\usage{
LSCV.risk(cases, controls, hlim = NULL,
	 method = c("kelsall-diggle", "hazelton"), res = 128, 
	 WIN = NULL, edge = TRUE, comment = TRUE)
}
\arguments{
  \item{cases}{
An object of type \code{\link{data.frame}}, \code{\link{list}}, \code{\link{matrix}}, or \code{\link[spatstat]{ppp}} describing the observed case data from which we wish to calculate the LSCV bandwidth. See `Details' for further information.
}
\item{controls}{
As for \code{cases}, but for the control observations. Both \code{cases} and \code{controls} must be of the same object class.
}
  \item{hlim}{
A numeric vector of length 2 giving the interval over which to search for the common bandwidth which minimises the selection criterion. If \code{NULL} (default), the function attempts to automatically select an appropriate range based on multiples of Stoyan and Stoyan's (1994) rule-of-thumb. The user is strongly recommended to supply their own \code{hlim}.
}
  \item{method}{
A character vector giving the specific selection criterion to minimise; see either Kelsall and Diggle (1995b) or Hazelton (2008). See `Details'. Defaults to \code{"kelsall-diggle"}.
}
  \item{res}{
Single integer giving the square grid resolution over which evaluation of the selection criterion takes place. Defaults to a 128 by 128 grid.
}
  \item{WIN}{
A polygonal \code{\link[spatstat]{owin}} object giving the study region. Ignored if \code{data} is already a \code{\link[spatstat]{ppp.object}}.
}
  \item{edge}{
Boolean. Whether or not to employ edge-correction in the calculations. Defaults to \code{TRUE}.
}
  \item{comment}{
Boolean. Whether or not to print function progress during execution. Defaults to \code{TRUE}.
}
}
  
\details{
This function calculates a `jointly optimal', common isotropic LSCV bandwidth for the (Gaussian) kernel-smoothed relative risk function (case-control density-ratio). If the \code{cases}, \code{controls} arguments are \code{\link{data.frame}} or \code{\link{matrix}} objects, these must each have exactly two columns containing the x (\code{[,1]}) and y (\code{[,2]}) data values. Should they be \code{\link{list}}s, these must have two vector components of equal length named \code{x} and \code{y}. Alternatively, \code{cases} and \code{controls} may be objects of class \code{\link[spatstat]{ppp}} (see \code{\link[spatstat]{ppp.object}}), and the argument \code{WIN} can be ignored.\cr

It can be shown that choosing a bandwidth that is equal for both case and control density estimates is preferable to computing `separately optimal' bandwidths (Kelsall and Diggle, 1995a). Setting \code{method = "kelsall-diggle"}, \code{LSCV.risk} computes the common bandwidth which minimises the approximate mean integrated squared error of the log-transformed risk surface (see specifically Kelsall and Diggle, 1995b).\cr

Alternatively, the user has the option of computing the common case-control bandwidth which minimises a \emph{weighted} mean integrated squared error of the (raw) relative risk function (see Hazelton, 2008). Generally, this author has found the Kelsall-Diggle method to provide more stable performance.
}

\value{
A single numeric value of the estimated bandwidth. The user may need to experiment with adjusting \code{hlim} to find a suitable minimum.
}

\references{
Kelsall, J.E. and Diggle, P.J. (1995a), Kernel estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.\cr\cr
Kelsall, J.E. and Diggle, P.J. (1995b), Non-parametric estimation of spatial variation in relative risk, \emph{Statistics in Medicine}, \bold{14}, 2335-2342.\cr\cr
Hazelton, M. L. (2008), Letter to the editor: Kernel estimation of risk surfaces without the need for edge correction, \emph{Statistics in Medicine}, \bold{27}, 2269-2272.\cr\cr
Stoyan, D. and Stoyan, H. (1994), \emph{Fractals, Random Shapes and Point Fields.} Wiley, Great Britain. ISBN 0-471-93757-6.
}
\author{
T.M. Davies
}

\section{Warning}{
Leave-one-out LSCV for jointly optimal, common bandwidth selection in the kernel-smoothed risk function is even more unstable (in terms of high variability) than the standalone density version. Caution is advised; not all applications will yield a successful result (this is termed ``a breakdown of the methodology'' by Kelsall and Diggle, 1995a). Undersmoothing has been noted in this author's personal experience. This method can also be computationally expensive for large data sets and fine evaluation grid resolutions.
}

\seealso{
\code{spatstat}'s function \code{\link[spatstat]{bw.relrisk}}
}

\examples{
data(chorley)

LSCV.risk(cases = split(chorley)[[1]], controls = split(chorley)[[2]],
 hlim = c(0.1,2))
}

