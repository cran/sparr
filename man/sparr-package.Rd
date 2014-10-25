\name{sparr-package}
\Rdversion{1.1}
\alias{sparr-package}
\alias{sparr}
\docType{package}
\title{
The sparr Package: SPAtial Relative Risk
}
\description{
Provides functions to estimate fixed and adaptive kernel-smoothed relative risk surfaces via the density-ratio method and perform subsequent inference.
}
\details{
\tabular{ll}{
Package: \tab sparr\cr
Version: \tab 0.3-6\cr
Date: \tab 2014-10-25\cr
License: \tab GPL (>= 2) \cr
}

Kernel smoothing, and the flexibility afforded by this methodology, provides an attractive approach to estimating complex probability density functions. This is particularly of interest when exploring problems in geographical epidemiology, the study of disease dispersion throughout some spatial region, given a population. The so-called `relative risk surface', constructed as a ratio of estimated case to control densities (Bithell, 1990; 1991), describes the variation in the `risk' of the disease, given the underlying at-risk population. This is a technique that has been applied successfully for mainly exploratory purposes in a number of different examples (see for example Sabel et al., 2000; Prince et al., 2001; Wheeler, 2007).

This package provides functions for bivariate kernel density estimation (KDE), implementing both fixed and `variable' or `adaptive' (Abramson, 1982) smoothing parameter options (see the function documentation for more information). A selection of bandwidth calculators for bivariate KDE and the relative risk function are provided, including one based on the maximal smoothing principle (Terrell, 1990), and others involving a leave-one-out least-squares cross-validation (see below). In addition, the ability to construct asymptotically derived p-value surfaces (`tolerance' contours of which signal statistically significant sub-regions of `extremity' in a risk surface - Hazelton and Davies, 2009; Davies and Hazelton, 2010), as well as some flexible visualisation tools, are provided.

The content of \code{sparr} can be broken up as follows:\cr

\emph{Datasets}\cr
\code{\link{PBC}} a case/control planar point pattern (\code{\link[spatstat]{ppp}}) concerning liver disease in northern England. Also available is the case/control dataset \code{\link[spatstat]{chorley}} of the \code{\link[spatstat]{spatstat}} package, which concerns the distribution of laryngeal cancer in an area of Lancashire, England.\cr\cr
\emph{Bandwidth calculators}\cr
\code{\link{OS}} estimation of an isotropic smoothing parameter for bivariate KDE, based on the oversmoothing principle introduced by Terrell (1990).\cr
\code{\link{NS}} estimation of an isotropic smoothing parameter for bivariate KDE, based on the optimal value for a normal density (bivariate normal scale rule - see e.g. Wand and Jones, 1995).\cr
\code{\link{LSCV.density}} a least-squares cross-validated (LSCV) estimate of an isotropic bandwidth for bivariate KDE (see e.g. Bowman and Azzalini, 1997).\cr
\code{\link{LSCV.risk}} a least-squares cross-validated (LSCV) estimate of a jointly optimal, common isotropic case-control bandwidth for the kernel-smoothed risk function (see Kelsall and Diggle, 1995a;b and Hazelton, 2008).\cr\cr
\emph{Bivariate functions}\cr
\code{\link{KBivN}} bivariate normal (Gaussian) kernel\cr
\code{\link{KBivQ}} bivariate quartic (biweight) kernel\cr
\code{\link{bivariate.density}} kernel density estimate of bivariate data; fixed or adaptive smoothing\cr\cr
\emph{Relative risk and p-value surfaces}\cr
\code{\link{risk}} estimation of a (log) relative risk function\cr
\code{\link{tolerance}} calculation of asymptotic p-value surface\cr\cr
\emph{Printing and summarising objects}\cr
\code{S3} methods (\code{\link{print.bivden}}, \code{\link{print.rrs}}, \code{\link{summary.bivden}} and \code{\link{summary.rrs}}) are available for the bivariate density and risk function objects.\cr\cr
\emph{Visualisation}\cr
Most applications of the relative risk function in practice require plotting the relative risk within the study region (especially for an inspection of tolerance contours). To this end, \code{\link{sparr}} provides a number of different ways to achieve attractive and flexible visualisation. The user may produce a heat plot, a perspective plot, a contour plot, or an interactive 3D perspective plot (that the user can pan around and zoom - courtesy of the powerful \code{\link[rgl]{rgl}} package; see below) for either an estimated relative risk function or a bivariate density estimate. These capabilities are available through \code{S3} support of the \code{plot} function; see\cr
\code{\link{plot.bivden}} for visualising a single bivariate density estimate from \code{\link{bivariate.density}}, and\cr
\code{\link{plot.rrs}} for visualisation of an estimated relative risk function from \code{\link{risk}}.
}

\section{Dependencies}{
The \code{sparr} package depends upon/imports some other important contributions to CRAN in order to operate; their uses here are indicated:\cr\cr
\code{\link[spatstat]{spatstat}} - Fast-fourier transform assistance with fixed and adaptive density estimation, as well as region handling; see Baddeley and Turner (2005).\cr
\code{\link[rgl]{rgl}} - Interactive 3D plotting of densities and surfaces; see Adler and Murdoch (2009).\cr
\code{MASS} - Utility support for internal functions; see Venables and Ripley (2002).
}
\author{
T.M. Davies\cr
Dept. of Mathematics & Statistics, University of Otago, Dunedin, New Zealand;\cr
M.L. Hazelton and J.C. Marshall\cr
Institute of Fundamental Sciences - Statistics, Massey University, Palmerston North, New Zealand.\cr

Maintainer: T.M.D. \email{tdavies@maths.otago.ac.nz}\cr
Feedback welcomed.
}
\section{Citation}{
To cite use of \code{sparr} in publications, the user may refer to the following work:\cr
Davies, T.M., Hazelton, M.L. and Marshall, J.C. (2011), \code{sparr}: Analyzing spatial relative risk using fixed and adaptive kernel density estimation in \code{R}, \emph{Journal of Statistical Software} \bold{39}(1), 1-14.

}
\references{
Abramson, I. (1982), On bandwidth variation in kernel estimates --- a square root law, \emph{Annals of Statistics}, \bold{10}(4), 1217-1223.\cr
Adler, D. and Murdoch, D. (2009), rgl: 3D visualization device system (OpenGL). R package version 0.87; URL: http://CRAN.R-project.org/package=rgl\cr
Baddeley, A. and Turner, R. (2005), Spatstat: an R package for analyzing spatial point patterns, \emph{Journal of Statistical Software}, \bold{12}(6), 1-42.\cr
Bithell, J.F. (1990), An application of density estimation to geographical epidemiology, \emph{Statistics in Medicine}, \bold{9}, 691-701.\cr
Bithell, J.F. (1991), Estimation of relative risk function,. \emph{Statistics in Medicine}, \bold{10}, 1745-1751.\cr
Bowman, A.W. and Azzalini, A. (1997), \emph{Applied Smoothing Techniques for Data Analysis: The Kernel Approach with S-Plus Illustrations.} Oxford University Press Inc., New York. ISBN 0-19-852396-3.\cr
%Bowman, A.W. and Azzalini, A. (2010), R package 'sm': nonparametric smoothing methods (version 2.2-4), URL: http://www.stats.gla.ac.uk/~adrian/sm; http://azzalini.stat.unipd.it/Book_sm\cr
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.\cr
Hazelton, M. L. (2008), Letter to the editor: Kernel estimation of risk surfaces without the need for edge correction, \emph{Statistics in Medicine}, \bold{27}, 2269-2272.\cr
Hazelton, M.L. and Davies, T.M. (2009), Inference based on kernel estimates of the relative risk function in geographical epidemiology, \emph{Biometrical Journal}, \bold{51}(1), 98-109.\cr
Kelsall, J.E. and Diggle, P.J. (1995a), Kernel estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.\cr
Kelsall, J.E. and Diggle, P.J. (1995b), Non-parametric estimation of spatial variation in relative risk, \emph{Statistics in Medicine}, \bold{14}, 2335-2342.\cr
Prince, M. I., Chetwynd, A., Diggle, P. J., Jarner, M., Metcalf, J. V. and James, O. F. W. (2001), The geographical distribution of primary biliary cirrhosis in a well-defined cohort, \emph{Hepatology} \bold{34}, 1083-1088.\cr
Sabel, C. E., Gatrell, A. C., Loytonenc, M., Maasiltad, P. and Jokelainene, M. (2000), Modelling exposure opportunitites: estimating relative risk for motor disease in Finland, \emph{Social Science & Medicine} \bold{50}, 1121-1137.\cr
Terrell, G.R. (1990), The maximal smoothing principle in density estimation, \emph{Journal of the American Statistical Association}, \bold{85}, 470-477.\cr
Venables, W. N. and Ripley, B. D. (2002). \emph{Modern Applied Statistics with S}, Fourth Edition, Springer, New York.\cr
Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.\cr
Wheeler, D. C. (2007), A comparison of spatial clustering and cluster detection techniques for childhood leukemia incidence in Ohio, 1996-2003, \emph{International Journal of Health Geographics}, \bold{6}(13).
}

\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{

%}
