\name{KBivN}
\Rdversion{1.1}
\alias{KBivN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard bivariate normal kernel
}
\description{
Evaluates the standard bivariate normal (Gaussian) kernel function at specified values.
}
\usage{
KBivN(X)
}
\arguments{
  \item{X}{
A numeric vector of length 2 or a data frame with 2 columns.
}
}
\details{
If \code{X} is a vector of length 2, then the two components \code{X[1]} and \code{X[2]} are taken to be the x and y coordinates respectively. For multiple evaluations at differing coordinates, \code{X} must be a data frame with \code{X[,1]} and \code{X[,2]} as the corresponding pairs of x and y coordinates respectively.
}
\value{
A single numeric value if \code{X} is a vector, or \code{nrow(X)} values if \code{X} is a data frame, giving the result of the standard bivariate normal kernel at the specified coordinate(s).
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
KBivN(c(0.1,0.4))

x <- y <- seq(-4,4,length=50)
z <- KBivN(data.frame(cbind(sort(rep(x,50)),rep(y,50))))
persp(x,y,matrix(z,50,50,byrow=TRUE),main="bivariate Gaussian kernel",
 phi=30,theta=-30)
}

