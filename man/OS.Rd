\name{OS}
\Rdversion{1.1}
\alias{OS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximal smoothing principle (oversmoothing) for bivariate KDE bandwidths
}
\description{
Provides an (isotropic) bandwidth estimate for use in bivariate KDE based on the oversmoothing factor introduced by Terrell (1990).
}
\usage{
OS(data, nstar = NULL, scaler = NA)
}
\arguments{
  \item{data}{
An object of type \code{\link{data.frame}}, \code{\link{list}}, \code{\link{matrix}}, or \code{\link[spatstat]{ppp}} giving the observed data from which we wish to calculate the OS bandwidth. See `Details' for further information.
}
  \item{nstar}{
A single numeric, positive value to use in place of the number of observations \emph{n} in the OS formula. If \code{NULL} (default), \emph{n} will simply be the number of observations in \code{data}.
}
  \item{scaler}{
A single numeric, positive value to use for transforming the result with respect to the scale of the recorded data. If \code{NA} (default), the scaling value is set as the mean of the interquartile ranges (IQR) of the \emph{x} and \emph{y} data values divided by 1.34 (Gaussian IQR). This approach was used in Davies and Hazelton (2010).
}
}

  
\details{
This function calculates a smoothing bandwidth for kernel density estimates of bivariate data, following the maximal smoothing priciple of Terrell (1990). If the \code{data} argument is a \code{\link{data.frame}} or a \code{\link{matrix}}, this must have exactly two columns containing the \emph{x} (\code{[,1]}) and \emph{y} (\code{[,2]}) data values. Should \code{data} be a \code{\link{list}}, this must have two vector components of equal length named \code{x} and \code{y}. Alternatively, \code{data} may be an object of class \code{\link[spatstat]{ppp}} (see \code{\link[spatstat]{ppp.object}}).
}

\value{
A single numeric value of the estimated bandwidth.
}

\references{
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.\cr\cr
Terrell, G.R. (1990), The maximal smoothing principle in density estimation, \emph{Journal of the American Statistical Association}, \bold{85}, 470-477.
}
\author{
T.M. Davies
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{

data(PBC)
PBC.casedata <- split(PBC)[[1]]
PBC.controldata <- split(PBC)[[2]]

pbc.h <- OS(PBC, nstar = sqrt(PBC.casedata$n * PBC.controldata$n))


##Scaling bandwidth for case data using standard deviations
sd.x <- sd(PBC.casedata$x)
sd.y <- sd(PBC.casedata$y)

OS(data = PBC.casedata, scaler = mean(sd.x, sd.y))

}

