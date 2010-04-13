run_ppp <- function(data,xy,h,WIN){
	ppp.full <- ppp(x=data[,1],y=data[,2],window=WIN,check=F)
	ppp.den.QEX <- density(x=ppp.full,sigma=h,xy=xy,weights=rep(1/ppp.full$n,ppp.full$n),spill=1)
	return(ppp.den.QEX)
}