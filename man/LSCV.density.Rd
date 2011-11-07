\name{LSCV.density}
\Rdversion{1.1}
\alias{LSCV.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leave-one-out least-squares cross-validation (LSCV) for bivariate KDE bandwidths
}
\description{
Provides an isotropic LSCV bandwidth estimate for use in 2-dimensional kernel density estimation (see e.g. Bowman and Azzalini, 1997).
}
\usage{
LSCV.density(data, hlim = NULL, res = 128, edge = TRUE,
	 WIN = NULL, quick = TRUE, comment = TRUE)
}
\arguments{
  \item{data}{
An object of type \code{\link{data.frame}}, \code{\link{list}}, \code{\link{matrix}}, or \code{\link[spatstat]{ppp}} describing the observed data from which we wish to calculate the LSCV bandwidth. See `Details' for further information.
}
  \item{hlim}{
A numeric vector of length 2 giving the interval over which to search for the bandwidth that minimises the selection criterion. If \code{NULL} (default), the function attempts to automatically select an appropriate range based on multiples of Stoyan and Stoyan's (1994) rule-of-thumb. The user is strongly recommended to supply their own \code{hlim}.
}
  \item{res}{
Single integer giving the square grid resolution over which evaluation of the selection criterion takes place. Defaults to a 128 by 128 grid.
}
  \item{edge}{
Boolean. Whether or not to employ edge-correction in the calculations. Defaults to \code{TRUE}.
}
  \item{WIN}{
A polygonal \code{\link[spatstat]{owin}} object giving the study region. Ignored if \code{data} is already a \code{\link[spatstat]{ppp.object}}.
}
  \item{quick}{
Intended for advanced use; users are recommended not to change the default \code{TRUE}. Setting \code{quick = FALSE} forces the function to individually evaluate the CV objective function at each of \code{seq(hlim[1], hlim[2], length = 50)} bandwidths, returning the corresponding values. Can be useful for diagnostic purposes.
}
  \item{comment}{
Boolean. Whether or not to print function progress during execution. Defaults to \code{TRUE}.
}


}
  
\details{
This function calculates a LSCV smoothing bandwidth for kernel density estimates of 2-dimensional (bivariate) data. If the \code{data} argument is a \code{\link{data.frame}} or a \code{\link{matrix}}, this must have exactly two columns containing the x (\code{[,1]}) and y (\code{[,2]}) data values. Should \code{data} be a \code{\link{list}}, this must have two vector components of equal length named \code{x} and \code{y}. Alternatively, \code{data} may be an object of class \code{\link[spatstat]{ppp}} (see \code{\link[spatstat]{ppp.object}}).
}

\value{
A single numeric value of the estimated bandwidth (if \code{quick = FALSE}, this value is named \code{hopt}; additionally returned are the objective function values (\code{lscv}) and the index of the minimum value (\code{ind})). The user may need to experiment with adjusting \code{hlim} to find a suitable minimum.
}

\references{
Bowman, A.W. and Azzalini, A. (1997), \emph{Applied Smoothing Techniques for Data Analysis: The Kernel Approach with S-Plus Illustrations.} Oxford University Press Inc., New York. ISBN 0-19-852396-3.\cr\cr
Stoyan, D. and Stoyan, H. (1994), \emph{Fractals, Random Shapes and Point Fields.} Wiley, Great Britain. ISBN 0-471-93757-6.\cr\cr
}

\author{
T.M. Davies
}

\section{Warning}{
Leave-one-out LSCV for bandwidth selection in kernel density estimation is notoriously unstable in practice and has a tendency to produce rather small bandwidths. Satisfactory bandwidths are not guaranteed for every application. This method can also be computationally expensive for large data sets and fine evaluation grid resolutions.
}

\seealso{
\code{spatstat}'s function \code{\link[spatstat]{bw.relrisk}}
}

\examples{
data(PBC)

##PBC cases
LSCV.density(split(PBC)[[1]],hlim=c(10,400))

##PBC controls
LSCV.density(split(PBC)[[2]],hlim=c(10,400))

}

