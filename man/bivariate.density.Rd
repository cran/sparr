\name{bivariate.density}
\Rdversion{1.1}
\alias{bivariate.density}
\alias{bivden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bivariate kernel density estimates
}
\description{
Provides an adaptive or fixed bandwidth kernel density estimate of bivariate data.
}
\usage{
bivariate.density(data, ID = NULL, pilotH, globalH = pilotH,
		  adaptive = TRUE, edgeCorrect = TRUE, res = 50, WIN = NULL,
		  counts = NULL, intensity = FALSE, xrange = NULL,
		  yrange = NULL, trim = 5, gamma = NULL, atExtraCoords = NULL,
		  use.ppp.methods = TRUE, comment = TRUE)
}
\arguments{
  \item{data}{
An object of type \code{\link{data.frame}}, \code{\link{list}}, \code{\link{matrix}}, or \code{\link[spatstat]{ppp}} giving the observed data from which we wish to calculate the density estimate. Optional ID information (e.g. a dichotomous indicator for cases and controls) may also be provided in these four data structures. See `Details' for further information on how to properly specify each one.
}
  \item{ID}{
If \code{data} is a data structure with a third component/column indicating case (1) or control (0) status, \code{ID} must specify which of these groups we wish to estimate a density for. If \code{ID} is \code{NULL} (default), a density is estimated for all present observations, regardless of any status information.
}
  \item{pilotH}{
A single numeric, positive `smoothing parameter' or `bandwidth'. When \code{adaptive} is \code{TRUE} (default), this value is taken to be the pilot bandwidth, used to construct the bivariate pilot density required for adaptive smoothing (see `Details'). For a fixed bandwidth kernel density estimate, \code{pilotH} simply represents the fixed amount of smoothing. Currently, all smoothing is isotropic in nature.
}
  \item{globalH}{
A single numeric, positive smoothing multiplier referred to as the global bandwidth, used to calculate the adaptive bandwidths (see `Details'). When \code{adaptive} is \code{TRUE}, this defaults to be the same as the pilot bandwidth. Ignored for a fixed density estimate.
}
  \item{adaptive}{
Boolean. Whether or not to produce an adaptive (variable bandwidth) density estimate, with the alternative being a fixed bandwith density estimate. Defaults to \code{TRUE}.
}
  \item{edgeCorrect}{
Boolean. Whether or not to perform edge-correction on the density estimate according to the methods demonstrated by Diggle (1985) (fixed bandwidth) and Marshall and Hazelton (2010) (adaptive). This can have a noticable effect on computation time in some cases. Defaults to \code{TRUE}. When \code{adaptive = TRUE}, the fixed-bandwidth pilot density is also edge-corrected according to \code{edgeCorrect}.
}
  \item{res}{
A single, numeric, positive integer indicating the square root of the desired resolution of the evaluation grid. That is, each of the evaluation grid axes will have length \code{res}. Currently, only \code{res*res} grids are supported. Defaults to 50 for computational reasons.
}
  \item{WIN}{
A polygonal object of class \code{\link[spatstat]{owin}} from the package \code{\link[spatstat]{spatstat}} giving the study region or `window'. All functions in the package \code{\link{sparr}} that require knowledge of the specific study region make use of this class; no other method of defining the study region is currently supported. If no window is supplied (default), the function defines (and returns) it's own rectangular \code{owin} based on \code{xrange} and \code{yrange}. Ignored if \code{data} is an object of type \code{\link[spatstat]{ppp}}.
}
  \item{counts}{
To perform binned kernel estimation, a numeric, positive, integer vector of giving counts associated with each observed coordinate in \code{data}, if \code{data} contains unique observations. If \code{NULL} (default), the function assumes each coordinate in \code{data} corresponds to one observation at that point. Should the \code{data} being supplied to \code{bivariate.density} contain duplicated coordinates, the function computes the counts vector internally (overriding any supplied value for \code{counts}), issues a warning, and continues with binned estimation. Non-integer values are rounded to the nearest integer.
}
  \item{intensity}{
A boolean value indicating whether or not to return an intensity (interpreted as the the expected number of observations per unit area and integrating to the number of observations in the study region) function, rather than a density (integrating to one). Defaults to \code{FALSE}.
}
  \item{xrange}{
Required only when no study region is supplied (\code{WIN = NULL}) and \code{data} is not an object of class \code{\link[spatstat]{ppp}}, and ignored otherwise. A vector of length 2 giving the upper and lower limits of the estimation interval for the x axis, in which case an evenly spaced set of values of length \code{res} is generated.
}
  \item{yrange}{
As above, but for the y axis.
}
  \item{trim}{
A numeric value (defaulting to \code{5}) that prevents excessively large bandwidths in adaptive smoothing by trimming the originally computed bandwidths \bold{h} by \code{trim} times \code{median(}\bold{h}\code{)}. A value of \code{NA} or a negative numeric value requests no trimming. Ignored when \code{adaptive} is \code{FALSE}.
}
  \item{gamma}{
An optional positive numeric value to use in place of \code{gamma} for adaptive bandwidth calculation (see `Details'). For adaptive relative risk estimation, this value can sensibly be chosen as common for both case and control densities (such as the \code{gamma} value from the adaptive density estimate of the `pooled' (full) dataset) - see Davies and Hazelton (2010). If nothing is supplied (default), this value is computed from the data being used to estimate the density in the defined fashion (again, see `Details'). Ignored for fixed bandwidth estimation.
}
  \item{atExtraCoords}{
It can occasionally be useful to retrieve the values of the estimated density at specific coordinates that are not the specific observations or the exact grid coordinates, for further analysis or plotting. \code{atExtraCoords} allows the user to specify an additional object of type \code{\link{data.frame}} with 2 colums giving the x \code{atExtraCoords[,1]} and y \code{atExtraCoords[,2]} coordinates at which to calculate and return the estimated density and other statistics (see `Value').
}
  \item{use.ppp.methods}{
Boolean. Whether or not to switch to using methods defined for objects of class \code{\link[spatstat]{ppp.object}} from the package \code{\link[spatstat]{spatstat}} to estimate the density. This approach is much, much faster than forcing \code{bivariate.density} to do the explicit calculations (due to implementation of a Fast Fourier Transform; see \code{\link[spatstat]{density.ppp}}) and is highly recommended for large datasets. To further reduce computation time in the adaptive case when \code{use.ppp.methods = TRUE}, the variable edge-correction factors are calculated using the integer percentiles of the varying bandwidths. Defaults to \code{TRUE}.
}
  \item{comment}{
Boolean. Whether or not to print function progress (including starting and ending times) during execution. Defaults to \code{TRUE}.
}
}

  
\details{
This function calculates an adaptive or fixed bandwidth bivariate kernel density estimate, using the bivariate Gaussian kernel. Abramson's method is used for adaptive smoothing (Abramson, 1982). Suppose our \code{data} argumnent is a \code{\link{data.frame}} or \code{\link{matrix}}. Then for each observation \code{data[i,1:2]} (\emph{i} = 1, 2, ... \emph{n}), the bandwidth \code{h[i]} is given by 
\cr\cr
\code{h[i]}=\code{globalH} / ( \emph{w}(\code{data[i,1:2]}; \code{pilotH})^(1/2)*\code{gamma} )
\cr\cr
where \emph{w} is the fixed bandwidth pilot density constructed with bandwidth \code{pilotH} and the scaling parameter \code{gamma} is the geometric mean of the \emph{w}^(-1/2) values. A detailed discussion on this construction is given in Silverman (1986).

If the \code{data} argument is a \code{\link{data.frame}} or a \code{\link{matrix}}, this must have exactly two columns containing the x (\code{[,1]}) and y (\code{[,2]}) data values, or exactly three columns with the third (rightmost) column giving ID information by way of a numeric, dichotomous indicator. Should \code{data} be a \code{\link{list}}, this must have two vector components of equal length named \code{x} and \code{y}. The user may specify a third component with the name \code{ID} giving the vector of corresponding ID information (must be of equal length to \code{x} and \code{y}). Alternatively, \code{data} may be an object of class \code{\link[spatstat]{ppp}} (see \code{\link[spatstat]{ppp.object}}). ID information can be stored in such an object through the argument \code{marks}. If \code{data} is a \code{\link[spatstat]{ppp}} object, the value of \code{window} of this object overrides the value of the argument \code{WIN} above.

}

\value{
An object of class \code{"bivden"}. This is effectively a list with the following components:
\item{Zm}{ a numeric matrix giving the value of the estimated (edge-corrected if elected) density at each of the coordinates of the grid. Values corresponding to points on the grid that fall outside the study region \code{WIN} are set to \code{NA}}
\item{X}{ a the sequence of values that were used as x grid coordinates. Will have length \code{res}}
\item{Y}{ a the sequence of values that were used as y grid coordinates. Will have length \code{res}}
\item{kType}{ the kernel function used in estimation. Currently fixed at \code{"gaus"}}
\item{h}{ a numeric vector with length equal to the number of observations, giving the bandwidths assigned to each observation in the order they appeared in \code{data}. For a fixed bandwidth estimate, this will simply be the identical value passed to and returned as \code{pilotH}}
\item{pilotH}{ the pilot or fixed bandwidth depending on whether adaptive smoothing is employed or not, respectively}
\item{globalH}{ the global bandwidth \code{globalH} if adaptive smoothing is employed, \code{NA} for fixed smoothing}
\item{hypoH}{ the matrix of `hypothetical' bandwidths (with element placement corresponding to \code{Zm}) for each coordinate of the evaluation grid. That is, these values are the bandwidths at that grid coordinate if, hypothetically, there was an observation there (along with the original data). These are used for edge-correction in adaptive densities (Marshall and Hazelton, 2010). Will be \code{NA} for fixed bandwidth estimates}
\item{zSpec}{ a numeric vector with length equal to the number of observations used, giving the values of the density at the specific coordinates of the observations. Order corresponds to the order of the observations in \code{data}}
\item{zExtra}{ as \code{zSpec} for the observations in \code{atExtraCoords}, \code{NA} if \code{atExtraCoords} is not supplied}
\item{WIN}{ the object of class \code{\link{owin}} used as the study region}
\item{qhz}{ a numeric matrix of the edge-correction factors for the entire evaluation grid (with placement corresponding to \code{Zm}. If \code{edgeCorrect = FALSE}, all edge correction factors are set to and returned as \code{1}}
\item{qhzSpec}{ edge-correction factors for the individual observations; order corresponding to \code{data}}
\item{qhzExtra}{ as \code{qhzSpec} for the observations in \code{atExtraCoords}; \code{NA} if \code{atExtraCoords} is not supplied}
\item{pilotvals}{ the values of the pilot density used to compute the adaptive bandwidths. Order corresponds to the order of the observations in \code{data}. \code{NULL} when \code{adaptive = FALSE}}
\item{gamma}{ the value of \code{gamma} that was passed to the function, or the geometric mean term of the reciprocal of the square root of the pilot density values used to scale the adaptive bandwidths if \code{gamma} is not supplied. \code{NULL} when \code{adaptive = FALSE}}
\item{counts}{ the counts vector used in estimation of the density/intensity. If all values in \code{data} were unique and \code{counts = NULL}, the returned \code{counts} will be a vector of ones equal to the number of coordinates in \code{data}}
\item{data}{ a two-column numeric data frame giving the observations in the originally supplied \code{data} that were used for the density estimation. If \code{data} originally contained duplicated coordinates, the returned \code{data} will contain only the unique coordinates, and should be viewed with respect to the returned value of \code{counts}}


}
\references{
Abramson, I. (1982). On bandwidth variation in kernel estimates --- a square root law, \emph{Annals of Statistics}, \bold{10}(4), 1217-1223.\cr\cr
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.\cr\cr
Diggle, P.J. (1985), A kernel method for smoothing point process data, \emph{Journal of the Royal Statistical Society, Series C}, \bold{34}(2), 138-147.\cr\cr
Marshall, J.C. and Hazelton, M.L. (2010) Boundary kernels for adaptive density estimators on regions with irregular boundaries, \emph{Journal of Multivariate Analysis}, \bold{101}, 949-963.\cr\cr
Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis}, Chapman & Hall, New York. \cr\cr
%Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
}
\author{
T.M. Davies
}
%\note{
%Future versions of this function will aim to further minimise the computational requirements.
%}

\section{Warning}{
Explicit calculation of bivariate kernel density estimates is computationally expensive. The decision to produce adaptive over fixed bandwidth estimates, the size of the dataset, the evaluation grid resolution specified by \code{res}, the complexity of the study region and electing to edge-correct all have a direct impact upon the time it will take to estimate the density. Keeping \code{use.ppp.methods = TRUE} can drastically reduce this computational cost at the expense a degree of accuracy that is generally considered negligible for most practical purposes.
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
##Chorley-Ribble laryngeal cancer data ('spatstat' library)
data(chorley)

ch.lar.density <- bivariate.density(data = chorley, ID = "larynx",
 pilotH = 1.5, adaptive = FALSE)

plot(ch.lar.density, col = "lightblue", phi = 30, theta = -30,
 ticktype = "detailed", main = "chorley.larynx", display = "persp")

\dontrun{

##PBC liver disease data
data(PBC)
pbc.adaptive.density <- bivariate.density(data = PBC$data, ID = 1,
 pilotH = 350, WIN = PBC$owin)

#3D plot - may need to adjust size of RGL device. Hold left click
# to pan, hold right to zoom.
plot(pbc.adaptive.density, display = "3d", col = heat.colors(20),
 main = "Density of PBC in north-east England", aspect = 1:2)
}
}
