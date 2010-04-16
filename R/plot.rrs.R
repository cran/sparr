plot.rrs <- function(x, ..., display = c("heat","contour","persp","3d"), show.WIN = TRUE){
	extras <- list(...)
	if(is.null(extras$xlab)) extras$xlab <- ""
	if(is.null(extras$ylab)) extras$ylab <- ""
	if(is.null(extras$zlab)&&(display=="persp"||display=="3d")) extras$zlab <- ""
	if(is.null(extras$main)) extras$main <- ""
	if(length(display)==4||display=="heat"){
		image <- im(t(x$rsM),xcol=x$f$X,yrow=x$f$Y)[x$f$WIN,drop=F]
		if(is.null(list(...)$axes)) plot(image,...,axes=T)
		else plot(image,...)
		if(show.WIN) plot(x$f$WIN,add=T,lwd=2)
	} else if(display=="contour"){
		contour(x$f$X,x$f$Y,x$rsM,...)
		if(show.WIN) plot(x$f$WIN,add=T,lwd=2)
	} else if(display=="persp"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				zfacet <- x$rsM[-1,-1]+x$rsM[-1,-ncol(x$rsM)]+x$rsM[-nrow(x$rsM),-1]+x$rsM[-nrow(x$rsM),-ncol(x$rsM)]
				facetcol <- cut(zfacet,length(extras$col))
				extras$col <- extras$col[facetcol]
			}
		}
		do.call("persp",c(list(x=x$f$X,y=x$f$Y,z=x$rsM), extras))
		if(show.WIN) lines(x$f$WIN,lwd=2)
	} else if(display=="3d"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				extras$col <- assignColors(extras$col,as.matrix(im(t(x$rsM),xcol=x$f$X,yrow=x$f$Y)))
			}
		}
		do.call("persp3d", c(list(x=x$f$X,y=x$f$Y,z=x$rsM), extras)) 
		if(show.WIN){
			gridLocs <- apply(as.data.frame(vertices(x$f$WIN)),1,getNearest,gridx=sort(rep(x$f$X,length(x$f$X))),gridy=rep(x$f$Y,length(x$f$Y)),WIN=x$f$WIN)
			lines3d(c(vertices(x$f$WIN)$x,vertices(x$f$WIN)$x[1]),c(vertices(x$f$WIN)$y,vertices(x$f$WIN)$y[1]),c(as.vector(t(x$rsM))[gridLocs],as.vector(t(x$rsM))[gridLocs][1]),lwd=4)
		}
	}
}
