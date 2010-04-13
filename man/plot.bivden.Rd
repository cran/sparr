\name{plot.bivden}
\Rdversion{1.1}
\alias{plot.bivden}
\alias{plot.rrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting a bivariate kernel density estimate object
}
\description{
\code{plot} methods for classes \code{"bivden"} and \code{"rrs"}
}
\usage{
\method{plot}{bivden}(x, \dots, display = c("heat", "contour", "persp", "3d"),
 show.WIN = TRUE)
\method{plot}{rrs}(x, \dots, display = c("heat", "contour", "persp", "3d"),
 show.WIN = TRUE)
}
\arguments{
  \item{x}{
An object of class \code{"bivden"} resulting from a call to \code{\link{bivariate.density}}, or an object of class \code{"rss"} resulting from a call to \code{\link{risk}}.
}
  \item{...}{
Additional graphical parameters to be passed to the relevant plot command depending on the value of \code{display}.
}
  \item{display}{
One of four possible character strings indicating the kind of plot desired (see 'Details'). Defaults to \code{"heat"}.
}
  \item{show.WIN}{
Boolean. Whether or not to draw the study region as an aesthetic enhancement to the plot of the density/risk surface. Defaults to \code{TRUE}.
}
}

  
\details{
There are currently four implemented plot types to visualise the estimated density or risk function. \code{"heat"} selects a heatplot, \code{"contour"} is simply a contour plot and \code{"persp"} creates a perspective plot. Selection of \code{"3d"} uses functions from the \code{\link[rgl]{rgl}} package to open an RGL graphics device and creates a 3-dimensional surface which the user can interact with using the mouse. To use \code{\dots} to improve the appearance of the four possible plot types \code{"heat"}, \code{"contour"}, \code{"persp"} and \code{"3d"}, the reader is highly recommended to consult the relevant documentation in the help pages \code{\link[spatstat]{plot.im}}, \code{\link{contour}}, \code{\link{persp}} and \code{\link[rgl]{persp3d}} respectively.
}

\value{
Plots to the relevant graphics device.
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
\code{\link{bivariate.density}, \link{risk}, \link{plot.default}, \link[spatstat]{plot.im}, \link{contour}, \link{persp}, \link[rgl]{persp3d}, \link{par}, \link[rgl]{par3d}}
}
\examples{
## see Examples in documentation for functions 'bivariate.density',
##  'risk' and 'tolerance'.
}

