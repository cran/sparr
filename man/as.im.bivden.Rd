\name{as.im.bivden}
\Rdversion{1.1}
\alias{as.im.bivden}
\alias{as.im.rrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converting a \code{sparr} bivariate kernel density estimate or relative risk surface object into a \code{\link[spatstat]{spatstat}} pixel image.
}
\description{
\code{as.im} methods for classes \code{"bivden"} and \code{"rrs"}
}
\usage{
\method{as.im}{bivden}(X, \dots)
\method{as.im}{rrs}(X, \dots)
}
\arguments{
  \item{X}{
An object of class \code{"bivden"} resulting from a call to \code{\link{bivariate.density}}, or an object of class \code{"rss"} resulting from a call to \code{\link{risk}}.
}
  \item{...}{
Ignored.
}
}

\value{
An object of class \code{\link[spatstat]{im}} corresponding to the supplied argument. Additional return information originally part of \code{X} is lost.
}

%\references{
%}

\author{
T.M. Davies
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[spatstat]{im}}, \code{\link[spatstat]{as.im}}
}
\examples{
data(chorley)

ch.bivden <- bivariate.density(chorley, ID = "lung", adaptive = FALSE, pilotH = 1.5)
ch.im <- as.im(ch.bivden)
summary(ch.im)

}

