\name{tolerance}
\Rdversion{1.1}

\alias{tolerance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic p-value surfaces
}
\description{
Calculates pointwise p-values based on asymptotic theory describing the extremity of risk over a given fixed or adaptive kernel-smoothed relative risk function.
}
\usage{
tolerance(rs, pooled, test = "greater", reduce = 1, exactL2 = TRUE, comment = TRUE)
}
\arguments{
  \item{rs}{
An object of class \code{"rrs"} resulting from a call to \code{\link{risk}}, giving the fixed or adaptive kernel-smoothed risk function.
}
  \item{pooled}{
An object of class \code{"bivden"} resulting from a call to \code{\link{bivariate.density}} (or the component \code{pooled} from \code{rs} if it was created using raw data arguments) representing a density estimate based on the `pooled' dataset of both `case' and `control' points. If separate from \code{rs}, this pooled density estimate must follow the same smoothing approach, evaluation grid and study region window as the densities used to create \code{rs}.
}
  \item{test}{
A character string indicating the kind of test desired to yield the p-values. Must be one of \code{"greater"} (default - performs upper tailed tests examining heighted risk `hotspots'), \code{"less"} (lower tailed tests examining `troughs') or \code{"both"} (double-sided tests).
}
  \item{reduce}{
A numeric value greater than zero and less than or equal to one giving the user the option to reduce the resolution of the evaluation grid for the pointwise p-values by specifying a proportion of the size of the evaluation grid for the original density estimates. For example, if the case and control \code{"bivden"} objects were calculated using \code{res = 100} and \code{tolerance} was called with \code{reduce = 0.5}, the p-value surface will be evaluated over a 50 by 50 grid. A non-integer value resulting from use of \code{reduce} will be ceilinged.
}
  \item{exactL2}{
Ignored if \code{rs} (and \code{pooled}) are fixed-bandwidth density estimates. A boolean value indicating whether or not to separately calculate the `L2' integral components for adaptive tolerance contours. A value of \code{FALSE} will approximate these components based on the `K2' integrals for faster execution (depending on the size of the evaluation grid, this improvement may be small) at the expense of a small degree of accuracy. Defaults to \code{TRUE}. See the reference for adaptive p-value surfaces in `Details' for definitions of these integral components.
}
  \item{comment}{
Boolean. Whether or not to print function progress (including starting and ending times) during execution. Defaults to \code{TRUE}.
}
}
  
\details{
This function implements developments in Hazelton and Davies (2009) (fixed) and Davies and Hazelton (2010) (adaptive) to compute pointwise p-value surfaces based on asymptotic theory of kernel-smoothed relative risk surfaces. Superimposing upon a plot of the risk surface contours of these p-values at given empirical significance levels (i.e. `tolerance contours') can be an informative way of exploring the statistical significance of the extremity of risk across the defined study region. The asymptotic approach to the p-value calculation is advantageous over a Monte-Carlo method, which would lead to excessive computation time for adaptive risk surfaces. See the aforementioned references for further comments.

}

\value{
A list with four components:
\item{X}{the equally spaced sequence of length \code{ceiling(reduce*res)} giving the evaluation locations on the x axis (where \code{res} is the grid resolution as specified in the calls to \code{\link{bivariate.density}} for calculation of the densities for \code{rs} and \code{pooled})}
\item{Y}{as above, for the y axis}
\item{Z}{a numeric \code{ceiling(reduce*res)*ceiling(reduce*res)} matrix giving the values of the risk surface over the evaluation grid. Values corresponding to grid coordinates outside the study region are assigned \code{NA}}
\item{P}{a \code{ceiling(reduce*res)*ceiling(reduce*res)} matrix giving the p-values corresponding to \code{Z} in light of the elected \code{test}}

}
\references{
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Submitted for publication}.\cr\cr
Hazelton, M.L. and Davies, T.M. (2009), Inference based on kernel estimates of the relative risk function in geographical epidemiology, \emph{Biometrical Journal}, \bold{51}(1), 98-109.
}
\author{
T.M. Davies and M.L. Hazelton
}
%\note{
%%  ~~further notes~~
%}

\section{Warning}{
Though far less expensive computationally than calculation of Monte-Carlo p-value surfaces, these asymptotic p-value surfaces (particularly for adaptive relative risk surfaces) can still take some time to complete. The argument of \code{reduce} provides an option to reduce this computation time by decreasing the resolution of the evaluation grid. However, the accuracy and appearance of the resulting tolerance contours can be severely degraded if \code{reduce} is assigned too small a value. Care must therfore be taken and consideration given to the resolution of the original evaluation grid when altering \code{reduce} from its default value. For most practical purposes, we have found a value of \code{reduce} resulting in evaluation of a p-value surface of size 50 by 50 is adequate.

}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

data(chorley)
ch.h <- CV.sm(chorley)

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


\dontrun{
data(PBC)
PBC.casedata <- PBC$data[PBC$data[,3]==1,1:2]
PBC.controldata <- PBC$data[PBC$data[,3]==0,1:2]

pbc.rrs.rawdata <- risk(f = PBC.casedata, g = PBC.controldata,
 WIN = PBC$owin, adaptive = TRUE, tolerate = TRUE)
}

}

