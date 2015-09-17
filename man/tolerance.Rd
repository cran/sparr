\name{tolerance}
\Rdversion{1.1}

\alias{tolerance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic p-value surfaces
}
\description{
Calculates pointwise p-values based on asymptotic theory or Monte-Carlo (MC) permutations describing the extremity of risk over a given fixed or adaptive kernel-smoothed relative risk function.
}
\usage{
tolerance(rs, test = "upper", method = "ASY", 
    pooled = NULL, symchoice = NULL, hpsim = NULL, 
    h0sim = NULL, reduce = 1, ITER = 1000, 
    exactL2 = TRUE, comment = TRUE)
}
\arguments{
  \item{rs}{
An object of class \code{"rrs"} resulting from a call to \code{\link{risk}}, giving the fixed or adaptive kernel-smoothed risk function.
}
  \item{test}{
A character string indicating the kind of test desired to yield the p-values. Must be one of \code{"upper"} (default - performs upper tailed tests examining heighted risk `hotspots'), \code{"lower"} (lower tailed tests examining `troughs') or \code{"double"} (double-sided tests). See `Details' for further information.
}
  \item{method}{
A character string, either "ASY" (default) or "MC" indicating which method to use for calculating the p-value surface (asymptotic and Monte-Carlo approaches respectively). The MC approach is far more computationally expensive than the asymptotic method (see `Warnings'). See `Details' for more information on the Monte-Carlo implementation.
}
  \item{pooled}{
Required if \code{method = "ASY"}, ignored otherwise. An object of class \code{"bivden"} resulting from a call to \code{\link{bivariate.density}} (or the component \code{pooled} from \code{rs} if it was created using raw data arguments) representing a density estimate based on the `pooled' dataset of both `case' and `control' points. If separate from \code{rs}, this pooled density estimate must follow the same smoothing approach, evaluation grid and study region window as the densities used to create \code{rs}.
}
  \item{symchoice}{
Currently only implemented for \code{method  = "MC"}; ignored otherwise (asymptotic version still in development). A character string denoting the density to use as the common pilot estimate for the symmetric adaptive estimator (Davies et al. 2015). If \code{symchoice = "f"} the case density estimate is used, if \code{symchoice = "g"} it's the control, and if \code{symchoice = "pool"} then a pooled estimate is used. If nothing is specified, i.e. \code{symchoice} is set at the default \code{NULL}, then the function implements the adaptive asymmetric estimator for the iterations. See `Details'.
}
  \item{hpsim}{
Used to specify the pilot bandwidth(s) for the Monte-Carlo simulated tolerance contours. Can be \code{NULL} (default), in which case the pilot bandwidths used at each iteration are the same values used for the original density estimates used in creation of \code{rs}; a numeric vector of length 1 (pilot bandwidths for case and control densities are the same) or 2 (in the order of \code{c(}pilot f\code{,}pilot g\code{)}); or a function (see `Details' for the latter). If the risk function in \code{rs} is a fixed-bandwidth estimate, then \code{hpsim} is used to directly specify the case and control density estimate bandwidths at each iteration. Ignored if \code{method = "ASY"}. 
}
  \item{h0sim}{
As for \code{hpsim}, but for the global bandwidths used at each iteration when computing Monte-Carlo tolerance contours for adaptive estimates in \code{rs}. Not used if \code{rs} is a fixed-bandwidth estimate.
}
  \item{reduce}{
A numeric value greater than zero and less than or equal to one giving the user the option to reduce the resolution of the evaluation grid for the pointwise p-values by specifying a proportion of the size of the evaluation grid for the original density estimates. For example, if the case and control \code{"bivden"} objects were calculated using \code{res = 100} and \code{tolerance} was called with \code{reduce = 0.5}, the p-value surface will be evaluated over a 50 by 50 grid. A non-integer value resulting from use of \code{reduce} will be ceilinged.
}
  \item{ITER}{
An integer value specifying the number of iterations to be used if \code{method = "MC"} (defaulting to 1000). Non-integer numeric values are rounded. Ignored when \code{method = "ASY"}.
}
  \item{exactL2}{
Ignored if \code{rs} (and \code{pooled}) are fixed-bandwidth density estimates, or if \code{method = "MC"}. A boolean value indicating whether or not to separately calculate the `L2' integral components for adaptive tolerance contours. A value of \code{FALSE} will approximate these components based on the `K2' integrals for faster execution (depending on the size of the evaluation grid, this improvement may be small) at the expense of a small degree of accuracy. Defaults to \code{TRUE}. See the reference for adaptive p-value surfaces in `Details' for definitions of these integral components.
}
  \item{comment}{
Boolean. Whether or not to print function progress (including starting and ending times) during execution. Defaults to \code{TRUE}.
}
}
  
\details{
This function implements developments in Hazelton and Davies (2009) (fixed) and Davies and Hazelton (2010) (adaptive) to compute pointwise p-value surfaces based on asymptotic theory of kernel-smoothed relative risk surfaces. Alternatively, the user may elect to calculate the p-value surfaces using Monte-Carlo methods (see Kelsall and Diggle, 1995). Superimposing upon a plot of the risk surface contours of these p-values at given significance levels (i.e. `tolerance contours') can be an informative way of exploring the statistical significance of the extremity of risk across the defined study region. The asymptotic approach to the p-value calculation is advantageous over a Monte-Carlo method, which can lead to excessive computation time for adaptive risk surfaces and large datasets. See the aforementioned references for further comments.

Choosing different options for the argument \code{test} simply manipulates the `direction' of the p-values. That is, plotting tolerance contours at a significance level of 0.05 for a p-value surface calculated with \code{test = "double"} is equivalent to plotting tolerance contours at significance levels of 0.025 and 0.975 for \code{test = "upper"}.

Implementation of the Monte-Carlo contours for the fixed-bandwidth estimator simply involves random allocation of case/control marks and re-estimation of the risk surface \code{ITER} times, against which the original estimate is compared. The bandwidth(s) for case and control densities in the permuted estimates are controlled by \code{hpsim}. If your risk surface is adaptive, \code{hpsim} is used to control the pilot bandwidths, \code{h0sim} the global bandwidths. In particular, for the adaptive symmetric estimator (Davies et al. 2015), it is assumed that the original estimate was itself calculated as a symmetric estimate via use of the \code{pdef} argument. The \code{symchoice} argument governs the specific permuted data set to use for the pilot estimate, and if \code{hpsim} is \code{NULL}, the pilot bandwidth thereof is taken from the stored \code{pdef} object in the original estimate. An error will occur if you attempt to set \code{symchoice} with an \code{rs} argument in this function that does not contain density estimates with present \code{pdef} objects of class \code{"bivden"}. See the help file for \code{\link{bivariate.density}} for details on using the \code{pdef} argument.

In addition to the usage noted above, you may define either \code{hpsim} and/or \code{h0sim} as functions which re-calculate the case and control pilot (or fixed) bandwidth(s) and the global bandwidth(s) at each iteration, based on the data set of the permuted case/control marks. If so, these must strictly be functions that take the case data as the first argument and the control data as the second argument, each as a two-column matrix of the x-y coordinates. The function must strictly return a numeric vector of length 1 or 2; these entries to be assigned to the relevant density estimates as per the usage notes on supply of a numeric vector for \code{hpsim}. Take care -- warnings will be issued if, for example, you specify a \code{hpsim} function that returns two numeric values, but your \code{rs} object is a symmetric-adaptive estimate (in which case it only makes sense to yield one pilot bandwidth)!
}

\value{
A list with four components:
\item{X}{the equally spaced sequence of length \code{ceiling(reduce*res)} giving the evaluation locations on the x axis (where \code{res} is the grid resolution as specified in the calls to \code{\link{bivariate.density}} for calculation of the densities for \code{rs} and \code{pooled})}
\item{Y}{as above, for the y axis}
\item{Z}{a numeric \code{ceiling(reduce*res)*ceiling(reduce*res)} matrix giving the values of the risk surface over the evaluation grid. Values corresponding to grid coordinates outside the study region are assigned \code{NA}. If \code{method = "MC"}, this will be a single value of \code{NA}}
\item{P}{a \code{ceiling(reduce*res)*ceiling(reduce*res)} matrix giving the p-values corresponding to the evaluation grid in light of the elected \code{test}}

}
\references{
Kelsall, J.E. and Diggle, P.J. (1995), Kernel estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.\cr\cr
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.\cr\cr
Davies, T.M., Jones, K. and Hazelton, M.L. (2015), Symmetric adaptive smoothing regimens for estimation of the spatial relative risk function, \emph{Submitted for publication}.\cr\cr
Hazelton, M.L. and Davies, T.M. (2009), Inference based on kernel estimates of the relative risk function in geographical epidemiology, \emph{Biometrical Journal}, \bold{51}(1), 98-109.
}
\author{
T.M. Davies and M.L. Hazelton
}
%\note{
%%  ~~further notes~~
%}

\section{Warning}{
Though far computationally intensive than calculation of Monte-Carlo p-value surfaces, the asymptotic p-value surfaces (particularly for adaptive relative risk surfaces) can still take some time to complete. The argument of \code{reduce} provides an option to reduce this computation time by decreasing the resolution of the evaluation grid. However, the accuracy and appearance of the resulting tolerance contours can be severely degraded if \code{reduce} is assigned too small a value. Care must therefore be taken and consideration given to the resolution of the original evaluation grid when altering \code{reduce} from its default value. For most practical purposes, we have found a value of \code{reduce} resulting in evaluation of a p-value surface of size 50 by 50 is adequate.

The MC approach is currently the only way to obtain tolerance contours for the adaptive-symmetric estimator. Given the computational expense, especially for adaptive estimates, it's recommended you do some preliminary runs with a small \code{ITER} value and/or make use of \code{reduce} to get acceptable running times.
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
data(chorley)
ch.h <- LSCV.density(chorley, hlim = c(0.1, 2))

ch.pool <- bivariate.density(data = chorley, pilotH = ch.h,
 adaptive = FALSE)
ch.case <- bivariate.density(data = chorley, ID = "larynx", pilotH = ch.h,
 adaptive = FALSE)
ch.con <- bivariate.density(data = chorley, ID = "lung", pilotH = ch.h,
 adaptive = FALSE)

##Compute log-risk surface and asymptotic p-value surfaces
ch.rrs <- risk(f = ch.case, g = ch.con) 
ch.tol <- tolerance(rs = ch.rrs, pooled = ch.pool)
contour(ch.tol$X, ch.tol$Y, ch.tol$P, levels = 0.05, add = TRUE)



data(PBC)
PBC.casedata <- split(PBC)[[1]]
PBC.controldata <- split(PBC)[[2]]

pbc.rrs.rawdata <- risk(f = PBC.casedata, g = PBC.controldata,
 adaptive = TRUE, tolerate = TRUE)
 
plot(pbc.rrs.rawdata, display = "3d", aspect = 1:2, col = heat.colors(12)[12:1],
 tolerance.matrix = pbc.rrs.rawdata$P, tol.opt = list(col = "white", raise = 0.03))
 
}

}

