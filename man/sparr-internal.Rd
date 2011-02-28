\name{sparr-internal}
\Rdversion{1.1}
\alias{sparr-internal}
\alias{AREA_BASIC}
\alias{areaLsq}
\alias{assignColors}
\alias{compute.coord}
\alias{dupli.data.frame}
\alias{getNearest}
\alias{getQhz_Adaptive}
\alias{getQhz_Fixed}
\alias{identical_windows}
\alias{idQuan}
\alias{KBivN_OPT}
\alias{KERNEL}
\alias{KSPEC}
\alias{Lsq}
\alias{nearestBC_Vals}
\alias{run_ppp}
\alias{translateCoord}
\title{
Internal sparr functions
}
\description{
Internal sparr functions, usually not to be called by the user.
}
\section{Usage}{

\code{AREA_BASIC(ndim = 2, lower, upper, minpts = 100, functn, ...)}\cr
\code{areaLsq(xh,WIN,iter=1000)}\cr
\code{assignColors(cols,MAT)}\cr
\code{compute.coord(coord,data,h,n,WIN,type,q.opt=1,forAdapt=F,counts)}\cr
\code{dupli.data.frame(data,w,comment)}\cr
\code{getNearest(obs,gridx,gridy,WIN,anypoint = FALSE)}\cr
\code{getQhz_Adaptive(ARG,kType,WIN,both=T,onlyk2=F)}\cr
\code{getQhz_Fixed(Xseq,Yseq,kType,WIN,h,both=T,onlyk2=F)}\cr
\code{identical_windows(w1,w2)}\cr
\code{idQuan(h,q)}\cr
\code{KBivN_OPT(x,WIN=NULL,mu=c(0,0),Sigma=matrix(c(1,0,0,1),2,2),sq=F,forSimul=F,y=NULL)}\cr
\code{KERNEL(X,type,dim=2,Qcomp="spher")}\cr
\code{KSPEC(location,data,h,type,counts)}\cr
\code{Lsq(x,uh,WIN)}\cr
\code{nearestBC_Vals(framevals,WIN)}\cr
\code{rsmc(object,pooled,ITER,corrgrid,comment)}\cr
\code{run_ppp(data,xy,h,WIN,counts)}\cr
\code{tol3d(X,Y,TOL,RHO,levels=c(0.05),raise=0.1,col="white",lwd=2,lty=1)}\cr
\code{tol3d.component(locations,gridx,gridy,gridz,raise)}\cr
\code{translateCoord(gridCoord,mvmat,WIN)}\cr

}
\keyword{internal}