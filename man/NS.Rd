\name{NS}
\Rdversion{1.1}
\alias{NS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normal scale rule for bivariate KDE bandwidths
}
\description{
Provides the (isotropic) optimal bandwidth for a bivariate normal density based on a simple expression.
}
\usage{
NS(data, nstar = NULL, scaler = NA)
}
\arguments{
  \item{data}{
An object of type \code{\link{data.frame}}, \code{\link{list}}, \code{\link{matrix}}, or \code{\link[spatstat]{ppp}} giving the observed data from which we wish to calculate the NS bandwidth. See `Details' for further information.
}
  \item{nstar}{
A single numeric, positive value to use in place of the number of observations \emph{n} in the NS formula. If \code{NULL} (default), \emph{n} will simply be the number of observations in \code{data}.
}
  \item{scaler}{
A single numeric, positive value to use for transforming the result with respect to the scale of the recorded data (i.e. a scalar representation of the standard deviation of the data). If \code{NA} (default), the scaling value is set as the mean of the interquartile ranges (IQR) of the \emph{x} and \emph{y} data values divided by 1.34 (Gaussian IQR).
}
}

  
\details{
This function calculates a smoothing bandwidth for kernel density estimates of 2-dimensional data: the optimal value which would minimise the asymptotic mean integrated squared error of the bivariate normal density function, assuming the standard Gaussian kernel function. See Wand and Jones (1995) for example. If the \code{data} argument is a \code{\link{data.frame}} or a \code{\link{matrix}}, this must have exactly two columns containing the \emph{x} (\code{[,1]}) and \emph{y} (\code{[,2]}) data values. Should \code{data} be a \code{\link{list}}, this must have two vector components of equal length named \code{x} and \code{y}. Alternatively, \code{data} may be an object of class \code{\link[spatstat]{ppp}} (see \code{\link[spatstat]{ppp.object}}).
}

\value{
A single numeric value of the estimated bandwidth.
}

\references{
Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
}
\author{
T.M. Davies
}

\section{Warning}{
The NS bandwidth is an approximation, and assumes \emph{that the target density is bivariate normal}. This is considered rare in e.g. epidemiological applications. Nevertheless, it remains a quick and easy `rule-of-thumb' method with which one may obtain a smoothing parameter in general applications.
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{

data(PBC)
PBC.casedata <- split(PBC)[[1]]
PBC.controldata <- split(PBC)[[2]]

pbc.h <- NS(PBC, nstar = sqrt(PBC.casedata$n * PBC.controldata$n))


##Scaling bandwidth for case data using standard deviations
sd.x <- sd(PBC.casedata$x)
sd.y <- sd(PBC.casedata$y)

NS(data = PBC.casedata, scaler = mean(sd.x, sd.y))

}