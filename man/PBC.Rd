\name{PBC}
\Rdversion{1.1}
\alias{PBC}
\docType{data}
\title{
Primary biliary cirrhosis data
}
\description{
Data of the locations of 761 cases of primary biliary cirrhosis in several adjacent health regions of north-eastern England, along with 3020 controls representing the at-risk population, collected between 1987 and 1994. These data were first presented and analysed by Prince et al. (2001); subsequent analysis of these data in the spirit of \code{\link{sparr}} was performed in Davies and Hazelton (2010). Also included is the polygonal study region.
}
\usage{data(PBC)}
\format{\code{PBC} is a named list with two elements.
\describe{
\item{\code{data}}{an object of type \code{\link{data.frame}} with 3781 rows and 3 columns.
\tabular{rlll}{
\code{PBC$data[,1]} \tab \code{x} \tab Easting coordinate (m) \cr
\code{PBC$data[,2]} \tab \code{y} \tab Northing coordinate (m) \cr
\code{PBC$data[,3]} \tab \code{ID} \tab Dichotomous indicator. 1: case, 0: control
}
}
%\item{data}{ a data frame with 3781 rows and 3 columns \code{x}, \code{y}, and \code{ID}, giving the Easting and Northing coordinates in meters and case(1)/control(0) designation respectively.}
\item{\code{owin}}{ an object of class \code{\link[spatstat]{owin}} descibing the geographical study region in the form of a polygon with 115 vertices.}
}
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Prince et al. (2001), The geographical distribution of primary biliary cirrhosis in a well-defined cohort, \emph{Hepatology}, \bold{34}, 1083-1088.\cr
}
\section{Acknowledgements}{
The authors thank Prof. Peter Diggle at Lancaster University (\url{http://www.lancs.ac.uk/staff/diggle/}) for providing access to these data.
}
\references{
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, in press [accepted 20 May 2010].
}
\examples{
data(PBC)
plot(PBC$owin,main="PBC liver disease data")

points(PBC$data$x[PBC$data$ID==0],
 PBC$data$y[PBC$data$ID==0]) #controls
points(PBC$data$x[PBC$data$ID==1],
 PBC$data$y[PBC$data$ID==1], col = 2) #cases

legend("bottomleft", legend = c("case", "control"),
 pch = c(1, 1), col = c(2, 1))
}

\keyword{datasets}
