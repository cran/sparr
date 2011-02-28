\name{risk}
\Rdversion{1.1}

\alias{risk}
\alias{rrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bivariate relative risk function
}
\description{
Estimates a relative risk function based on the ratio of two bivariate kernel density estimates over identical grids and regions. In geographical epidemiology, the two densities would represent a set of disease cases (numerator) and a sample of controls illustrating the at-risk population (denominator).
}
\usage{
risk(f, g, delta = 0, log = TRUE, h = NULL, adaptive = FALSE, res = 50,
 WIN = NULL, tolerate = FALSE, plotit = TRUE, comment = TRUE)
}
\arguments{
  \item{f}{
Either a pre-calculated object of class \code{"bivden"} representing the `case' density estimate, or an object of type \code{\link{data.frame}}, \code{\link{list}}, \code{\link{matrix}}, or \code{\link[spatstat]{ppp}} giving the observed case data. If this raw data is provided, a kernel density estimate is computed internally, with certain options available to the user in \code{\link{bivariate.density}} chosen/calculated automatically. See `Details' for further information.
}
  \item{g}{
As for argument \code{f}, but for the controls. Whatever the type, the class of \code{g} must match that of \code{f}.
}
  \item{delta}{
A single numeric scaling parameter used for an optional additive constant to the densities; occasionally used for risk surface construction (see `Details'). A negative or zero value for \code{delta} requests no additive constant (default).
}
  \item{log}{
Boolean. Whether or not to return the (natural) log-transformed relative risk function as recommended by Kelsall and Diggle (1995a). Defaults to \code{TRUE} with the alternative being the raw density ratio.
}
  \item{h}{
Ignored if \code{f} and \code{g} are already \code{"bivden"} objects. An optional numeric vector of length 1 OR 2, giving the global bandwidth(s) for internal estimation of the case and control densities if \code{adaptive = TRUE}, or the fixed bandwidth(s) if \code{adapative = FALSE}. When \code{h} is a single numeric value, this is elected as the common global/fixed bandwidth for case and control densities. When \code{h} has length 2, the values \code{h[1]} and \code{h[2]} are assigned as the case and control global/fixed bandwidths respectively. By default, a value of \code{h = NULL} tells the function to use the global/fixed smoothing parameters as outlined in `Details' below. Note that for adaptive estimation, this argument does not affect calculation of the pilot bandwiths.
}
  \item{adaptive}{
Ignored if \code{f} and \code{g} are already \code{"bivden"} objects. A boolean value specifying whether or not to employ adaptive smoothing for internally estimating the densities. A value of \code{FALSE} (default) elects use of fixed-bandwidth estimates.
}
  \item{res}{
Ignored if \code{f} and \code{g} are already \code{"bivden"} objects. A numeric value giving the desired resolution (of one side) of the evaluation grid. Higher values increase resolution at the expense of computational efficiency. Defaults to a 50 by 50 grid.
}
  \item{WIN}{
Ignored if \code{f} and \code{g} are already \code{"bivden"} objects OR objects of class \code{\link[spatstat]{ppp}} (in which case the study region is set to the value of the resident \code{window} component). A polygonal object of class \code{\link[spatstat]{owin}} giving the relevant study region in which the \code{f} and \code{g} data was collected.
}
  \item{tolerate}{
Ignored if \code{f} and \code{g} are already \code{"bivden"} objects. A boolean value specifying whether or not to calculate a corresponding asymptotic p-value surface (for tolerance contours) for the estimated relative risk function. If \code{TRUE}, the p-value surface tests for elevated risk only (equivalent to setting \code{test = "greater"} in \code{\link{tolerance}}) and is evaluated over a maximum grid resolution of 50 by 50. Defaults to \code{FALSE} for computational reasons.
}
  \item{plotit}{
Boolean. If \code{TRUE} (default), a heatplot of the estimated relative risk function is produced. If \code{tolerate = TRUE}, asymptotic tolerance contours are automatically added to the plot at a significance level of 5\%.

}
  \item{comment}{
Ignored if \code{f} and \code{g} are already \code{"bivden"} objects. Boolean. Whether or not to print function progress (including starting and ending date-times) during execution. Defaults to \code{TRUE}.
}
}
\details{
This function estimates a relative risk function via the density ratio method using fixed or adaptive bandwidth bivariate kernel density estimates. Both densities must be estimated using the same evaluation grid (and the same study window) in \code{\link{bivariate.density}}. In geographical epidemiology, the argument \code{f} represents the spatial distribution of the disease cases, and \code{g} the at-risk (control) population.

The option to supply the raw case and control data is available. If this is done, the function runs \code{\link{bivariate.density}} internally, abstracting certain decisions about the density estimation away from the user. If the user sets \code{adaptive = TRUE} (and \code{h} remains at \code{NULL}), the smoothing parameters are calculated as per the approach taken in Davies and Hazelton (2010): a common global bandwidth using the pooled data from \code{\link{OS}}, and differing pilot bandwidths using \code{\link{CV.sm}} on the case and control data separately. The scaling parameter \code{gamma} is common for the case and control density estimates, set as the \code{gamma} component of the pooled estimate. If a fixed relative risk is desired (\code{adaptive = FALSE}) and no specific bandwidths are given via the argument \code{h}, the case and control densities share a common bandwidth computed from the pooled data using \code{\link{OS}}. In supplying raw data to \code{\link{risk}}, the user must also specify an evaluation grid resolution (defaulting to 50 by 50) and the study region \code{WIN} (unless \code{f} and \code{g} are objects of class \code{\link[spatstat]{ppp}}, in which case the resident \code{window} component overrides \code{WIN}). All other arguments are set to their defaults as in \code{\link{bivariate.density}}.

If more flexibility is required for estimation of the case and control densities, the user must supply `pre-calculated' objects of class \code{"bivden"} (from \code{\link{bivariate.density}}) as the \code{f} and \code{g} arguments. This drastically reduces the running time of a call to \code{\link{risk}} (as the density estimation step is already complete). However, the option of internally computing the asymptotic p-value surfaces (via the argument \code{tolerate}) is unavailable in this case; the user must run the \code{\link{tolerance}} function separately if tolerance contours are desired.

The relative risk function is defined here as the ratio of the `case' density to the `control' (Bithell, 1990; 1991). Using kernel density estimation to model these densities (Diggle, 1985), we obtain a workable estimate thereof. This function defines the risk function \emph{r} in the following fashion:
\cr\cr
\emph{r}\code{ = (f + delta*max(g))/(g + delta*max(g))}
\cr\cr
Note the (optional) additive constants defined by \code{delta} times the maximum of each of the densities in the numerator and denominator respectively (see Bowman and Azzalini, 1997).

The log-risk function \emph{rho}, given by \emph{rho} = log[\emph{r}], is argued to be preferable in practice as it imparts a sense of symmetry in the way the case and control densities are treated (Kelsall and Diggle, 1995a;b). The option of log-transforming the returned risk function is therefore selected by default.
}

\value{
An object of class \code{"rrs"}. This is a marked list with the following components:
\item{rsM}{a numeric \code{res*res} matrix (where \code{res} is the grid resolution as specified in the calls to \code{\link{bivariate.density}} for calculation of \code{f} and \code{g}) giving the values of the risk surface over the evaluation grid. Values corresponding to grid coordinates outside the study region are assigned \code{NA}}
\item{f}{the object of class \code{"bivden"} used as the numerator or `case' density estimate}
\item{g}{the object of class \code{"bivden"} used as the denominator or `control' density estimate}
\item{log}{whether or not the returned risk function is on the log-scale}
\item{pooled}{the object of class \code{"bivden"} (based on the pooled data) calculated internally if \code{f} and \code{g} were raw data arguments, \code{NA} otherwise}
\item{P}{a numeric 50 by 50 matrix of the asymptotic p-value surface if \code{tolerate = TRUE} and \code{f} and \code{g} were raw data arguments, \code{NA} otherwise}
}

\references{
Bithell, J.F. (1990), An application of density estimation to geographical epidemiology, \emph{Statistics in Medicine}, \bold{9}, 691-701.\cr\cr
Bithell, J.F. (1991), Estimation of relative risk functions, \emph{Statistics in Medicine}, \bold{10}, 1745-1751.\cr\cr
Bowman, A.W. and Azzalini A. (1997), \emph{Applied Smoothing Techniques for Data Analysis: The Kernel Approach with S-Plus Illustrations}, Oxford University Press Inc., New York.\cr\cr
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.\cr\cr
Diggle, P.J. (1985), A kernel method for smoothing point process data, \emph{Journal of the Royal Statistical Society Series C}, \bold{34}(2), 138-147.\cr\cr
Kelsall, J.E. and Diggle, P.J. (1995a), Kernel estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.\cr\cr
Kelsall, J.E. and Diggle, P.J. (1995b), Non-parametric estimation of spatial variation in relative risk, \emph{Statistics in Medicine}, \bold{14}, 2335-2342.
}
\author{
T.M. Davies, M.L. Hazelton and J.C. Marshall
}
%\note{
%%  ~~further notes~~
%}

\section{Warning}{
If raw data is supplied to \code{\link{risk}}, as opposed to previously computed objects of class \code{"bivden"}, the running time of this function will be greater. This is particularly the case if the user has also selected \code{tolerate = TRUE}. In the same fashion as \code{\link{bivariate.density}} and \code{\link{tolerance}}, setting \code{comment = TRUE} can keep the user appraised of the function progress during run-time.
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(PBC)
PBC.casedata <- PBC$data[PBC$data[,3]==1,1:2]
PBC.controldata <- PBC$data[PBC$data[,3]==0,1:2]

pbc.h <- OS(PBC$data[,1:2],
 nstar = sqrt(nrow(PBC.casedata) * nrow(PBC.controldata)))

pbc.pool <- bivariate.density(data = PBC$data[,1:2], pilotH = pbc.h,
 adaptive = FALSE, WIN = PBC$owin)
pbc.case <- bivariate.density(data = PBC.casedata, pilotH = pbc.h,
 adaptive = FALSE, WIN = PBC$owin)
pbc.con <- bivariate.density(data = PBC.controldata, pilotH = pbc.h,
 adaptive = FALSE, WIN = PBC$owin)

pbc.rrs <- risk(f = pbc.case, g = pbc.con, plotit = FALSE)
summary(pbc.rrs)

}

