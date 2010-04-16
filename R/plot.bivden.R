plot.bivden <- function(x, ..., display = c("heat","contour","persp","3d"), show.WIN = TRUE){
	extras <- list(...)
	if(is.null(extras$xlab)) extras$xlab <- ""
	if(is.null(extras$ylab)) extras$ylab <- ""
	if(is.null(extras$zlab)&&(display=="persp"||display=="3d")) extras$zlab <- ""
	if(is.null(extras$main)) extras$main <- ""
	if(length(display)==4||display=="heat"){
		image <- im(t(x$Zm),xcol=x$X,yrow=x$Y)[x$WIN,drop=F]
		if(is.null(list(...)$axes)) plot(image,...,axes=T)
		else plot(image,...)
		if(show.WIN) plot(x$WIN,add=T,lwd=2)
	} else if(display=="contour"){
		contour(x$X,x$Y,x$Zm,...)
		if(show.WIN) plot(x$WIN,add=T,lwd=2)
	} else if(display=="persp"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				zfacet <- x$Zm[-1,-1]+x$Zm[-1,-ncol(x$Zm)]+x$Zm[-nrow(x$Zm),-1]+x$Zm[-nrow(x$Zm),-ncol(x$Zm)]
				facetcol <- cut(zfacet,length(extras$col))
				extras$col <- extras$col[facetcol]
			}
		}
		do.call("persp",c(list(x=x$X,y=x$Y,z=x$Zm), extras))

	} else if(display=="3d"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				extras$col <- assignColors(extras$col,as.matrix(im(t(x$Zm),xcol=x$X,yrow=x$Y)))
			}
		}
		do.call("persp3d", c(list(x=x$X,y=x$Y,z=x$Zm), extras)) 
		if(show.WIN){
			gridLocs <- apply(as.data.frame(vertices(x$WIN)),1,getNearest,gridx=sort(rep(x$X,length(x$X))),gridy=rep(x$Y,length(x$Y)),WIN=x$WIN)
			lines3d(c(vertices(x$WIN)$x,vertices(x$WIN)$x[1]),c(vertices(x$WIN)$y,vertices(x$WIN)$y[1]),c(as.vector(t(x$Zm))[gridLocs],as.vector(t(x$Zm))[gridLocs][1]),lwd=4)
		}
	}
}
