\name{CV.sm}
\Rdversion{1.1}
\alias{CV.sm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-squares cross-validation (LSCV) for bivariate KDE bandwidths
}
\description{
Provides an isotropic LSCV bandwidth estimate for use in bivariate kernel density estimation, taken from the function \code{\link[sm]{h.select}} in the package \code{\link[sm]{sm}} (see Bowman and Azzalini, 1997; 2010).
}
\usage{
CV.sm(data)
}
\arguments{
  \item{data}{
An object of type \code{\link{data.frame}}, \code{\link{list}}, \code{\link{matrix}}, or \code{\link[spatstat]{ppp}} describing the observed data from which we wish to calculate the LSCV bandwidth. See `Details' for further information.
}
}
  
\details{
This function calculates a LSCV smoothing bandwidth for kernel density estimates of bivariate data. If the \code{data} argument is a \code{\link{data.frame}} or a \code{\link{matrix}}, this must have exactly two columns containing the x (\code{[,1]}) and y (\code{[,2]}) data values. Should \code{data} be a \code{\link{list}}, this must have two vector components of equal length named \code{x} and \code{y}. Alternatively, \code{data} may be an object of class \code{\link[spatstat]{ppp}} (see \code{\link[spatstat]{ppp.object}}).
}

\value{
A single numeric value of the estimated bandwidth.
}

\references{
Bowman, A.W. and Azzalini, A. (1997), \emph{Applied Smoothing Techniques for Data Analysis: The Kernel Approach with S-Plus Illustrations.} Oxford University Press Inc., New York. ISBN 0-19-852396-3.\cr\cr
Bowman, A.W. and Azzalini, A. (2010), R package `sm': nonparametric smoothing methods (version 2.2-4), URL: \url{http://www.stats.gla.ac.uk/~adrian/sm}; \url{http://azzalini.stat.unipd.it/Book_sm}.
}
\author{
T.M. Davies
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
data(PBC)

##PBC cases
CV.sm(PBC$data[PBC$data[,3]==1,1:2])

##PBC controls
CV.sm(PBC$data[PBC$data[,3]==0,1:2])

}

