\name{KBivQ}
\Rdversion{1.1}
\alias{KBivQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard bivariate quartic (biweight) kernel
}
\description{
Evaluates the standard bivariate quartic (biweight) kernel function at specified values, for either the spherical or product derivation of the function.
}
\usage{
KBivQ(X,type="spher")
}
\arguments{
  \item{X}{
A numeric vector of length 2 or a data frame with 2 columns.
}
  \item{type}{
A character string.
\tabular{lr}{
  \code{"spher"} \tab (default) selects spherical method of calculating the bivariate quartic kernel function\cr
  \code{"prod"} \tab uses the product approach to calculating the function
}
}
}
  
\details{
If \code{X} is a vector of length 2, then the two components \code{X[1]} and \code{X[2]} are taken to be the x and y coordinates respectively. For multiple evaluations at differing coordinates, \code{X} must be a data frame with \code{X[,1]} and \code{X[,2]} as the corresponding pairs of x and y coordinates respectively.\cr\cr

Unlike the bivariate Gaussian kernel, it is necessary to specify the method of extending the univariate quartic kernel to the bivariate case; this can be done in two different ways, one way resulting in a slightly different kernel to the other. An explanation of these `spherical' and `product' approaches is given in Wand and Jones (1995).
}
\value{
A single numeric value if \code{X} is a vector, or \code{nrow(X)} values if \code{X} is a data frame, giving the result of the standard bivariate quartic kernel at the specified coordinate(s) for the elected function derivation type.
}
\references{
Wand, M.P. and Jones, C.M. (1995), \emph{Kernel Smoothing}, Chapman & Hall, London.
}
\author{
T.M. Davies
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
KBivQ(c(0.1,0.4))

x <- y <- seq(-0.9,0.9,length=50)
z.spher <- KBivQ(data.frame(cbind(sort(rep(x,50)), rep(y,50))))
z.prod <- KBivQ(data.frame(cbind(sort(rep(x,50)), rep(y,50))), "prod")

par(mfrow=c(1,2))
contour(x, y, matrix(z.spher, 50, 50, byrow = TRUE),
 main = "bivariate quartic kernel (spherical)")
contour(x, y, matrix(z.prod, 50, 50, byrow=TRUE),
 main = "bivariate quartic kernel (product)")

}


