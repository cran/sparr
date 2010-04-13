summary.bivden <- function(object, ...){

	print.bivden(x=object)
	
	if(length(object)==7) cat("Computed using ppp methods of package 'spatstat'\n")

	cat("Evaluated over",nrow(object$Zm),"by",ncol(object$Zm),"rectangular grid.\nDefined study region is a polygon with",length(vertices(object$WIN)$x),"vertices.\n\n")
	
	cat("Estimated density description\n")
	print(summary(as.vector(object$Zm)))

}
	