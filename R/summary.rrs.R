summary.rrs <- function(object, ...){
	if(!object$log) cat("Relative risk function.\n\n")
	else cat("Log-Relative risk function.\n\n")
	
	cat("Surface (Z) summary:\n")
	print(summary(as.vector(object$rsM)))
	
	cat("\n--Numerator (case) density--\n")
	summary.bivden(object$f)
	cat("\n--Denominator (control) density--\n")
	summary.bivden(object$g)
	
}
	