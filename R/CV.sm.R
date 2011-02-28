CV.sm <- function(data){
	if(class(data)=="data.frame"){
		if(ncol(data)!=2) stop("data.frame must have exactly two columns")
		if(nrow(data)<10) warning("less than 10 observations!")
	} else if(class(data)=="list"){
		if(is.null(data$x)||is.null(data$y)) stop("data list must have two components named 'x' and 'y'")
		if(length(data$x)!=length(data$y)) stop("data components 'x' and 'y' are of unequal lengths")
		if(length(data$x)<10) warning("less than 10 observations!")
		data <- as.data.frame(data)
	} else if(class(data)=="matrix"){
		if(ncol(data)!=2) stop("data matrix must have exactly two columns")
		if(nrow(data)<10) warning("less than 10 observations!")
		data <- as.data.frame(data)
	} else if(class(data)=="ppp"){
		if(is.null(data$x)||is.null(data$y)) stop("data ppp.object must have two non-empty components named 'x' and 'y'")
		if(length(data$x)<10) warning("less than 10 observations!")
		data <- data.frame(cbind(data$x,data$y))
	} else {
		stop("data must be an object of type 'data.frame', 'list', 'matrix', or 'ppp'")
	}
	
	result <- h.select(x=as.matrix(data),method="cv",structure.2d="common")[1]
	return(result)
}