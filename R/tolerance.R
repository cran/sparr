tolerance <- function(rs, pooled, test = "greater", reduce = 1, exactL2 = TRUE, comment = TRUE){ 
	if(comment) print(date())
	
	if(class(rs)!="rrs") stop("'rs' must be of class 'rrs'")
	if(class(pooled)!="bivden") stop("'pooled' must be an object of class 'bivden'")
	
	if(!all(c(length(rs$f$zVec)==length(pooled$zVec),length(rs$g$zVec)==length(pooled$zVec)))){
		stop("'pooled' appears to have been estimated using a different evaluation grid to that used in 'rs' - evaluation grids must be identical")
	}
	if(!all(c(identical_windows(rs$f$WIN,pooled$WIN),identical_windows(rs$g$WIN,pooled$WIN)))){
		stop("'pooled$WIN' does not appear to match study region window used in 'rs' - study regions must be identical")
	}
	if(length(pooled$hypoH)!=length(rs$f$hypoH)) stop("smoothing approach (fixed or adaptive) of 'pooled' must match approach used in 'rs'")
	
	if(all(c("greater","less","both")!=test)) stop("'test' must be one of 'greater', 'less' or 'both'")
	
	if(reduce<=0) stop("'reduce' must be greater than zero or less than or equal to one")
	if(reduce>1) stop("'reduce' must be greater than zero or less than or equal to one")

	edgef <- range(as.vector(rs$f$qhz),na.rm=T)[1]!=range(as.vector(rs$f$qhz),na.rm=T)[2]
	edgep <- range(as.vector(pooled$qhz),na.rm=T)[1]!=range(as.vector(pooled$qhz),na.rm=T)[2]
	
	if((edgef+edgep==1)) stop("edge-correction is inconsistent. all densities must be either edge-corrected or not.")
	
	
	adaptive <- (length(pooled$hypoH)>1)
	fvec <- as.vector(t(rs$f$Zm))
	gvec <- as.vector(t(rs$g$Zm))

	xr <- range(rs$f$X)
	yr <- range(rs$f$Y)
	gsize <- ceiling(reduce*length(rs$f$X))

	if(gsize < 10) warning("given reduction results in smaller p-value grid than 10x10")

	grx <- sort(rep(seq(xr[1],xr[2],length=gsize),gsize))
	gry <- rep(seq(yr[1],yr[2],length=gsize),gsize)
		
	if(reduce==1){
		corrGridSpec <- 1:(length(rs$f$X)^2)
	} else {
		corrGridSpec <- apply(as.matrix(data.frame(cbind(grx,gry))),1,getNearest,gridx=sort(rep(rs$f$X,length(rs$f$X))),gridy=rep(rs$f$Y,length(rs$f$Y)),WIN=rs$f$WIN,anypoint=T)
	}
	
	datarange.list <- list(x=seq(xr[1],xr[2],length=gsize),y=seq(yr[1],yr[2],length=gsize))
	
	if(edgep){
		if(adaptive){
			if(comment) cat("\n--Adaptive-bandwidth asymptotics--\n")
			if(comment) cat("calculating integrals K2...\n")
		
			#if(comment) cat("--f--\n")
			#fk2 <- apply(as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$f$hypoH))[corrGridSpec]))),1,getQhz_Adaptive,kType="gaus",WIN=rs$f$WIN,both=F,onlyk2=T)[2,]
			#if(comment) cat("--g--\n")
			#gk2 <- apply(as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$g$hypoH))[corrGridSpec]))),1,getQhz_Adaptive,kType="gaus",WIN=rs$g$WIN,both=F,onlyk2=T)[2,]
			
			if(comment) cat("--f--\n")
			hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
			corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
			qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$f$data,xy=datarange.list,WIN=rs$f$WIN)
			fk2 <- rep(-1,gsize*gsize)
			fk2[is.na(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec])] <- NA
			for(i in 1:length(hypoQuan)) fk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
			fk2 <- (1/(4*pi))*fk2 
			
			if(comment) cat("--g--\n")
			hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
			corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
			qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$g$data,xy=datarange.list,WIN=rs$f$WIN)
			gk2 <- rep(-1,gsize*gsize)
			gk2[is.na(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec])] <- NA
			for(i in 1:length(hypoQuan)) gk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
			gk2 <- (1/(4*pi))*gk2 
			
			if(exactL2){
				xrL <- sort(rep(seq(-4,4,length=10),10))
				yrL <- rep(seq(-4,4,length=10),10)
				grL <- matrix(c(xrL,yrL),100,2)
				Lsq_gr <- apply(grL,1,Lsq,uh=c(0,0,1),WIN=NULL)
				
				
				if(comment) cat("calculating integrals L2...\n--f--\n")
				#S1rzK <- (1/(as.vector(t(rs$f$qhz))[corrGridSpec]^2))*(2*fk2 + 0.25*apply(as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$f$hypoH))[corrGridSpec]))),1,areaLsq,WIN=rs$f$WIN,iter=4))
				coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$f$hypoH))[corrGridSpec])))
				fL2 <- c()
				for(i in 1:nrow(coords)){
					if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
						tempLsq <- NA
					} else {
						#h <<- coords[i,3]
						temp.xrL <- xrL*coords[i,3]+coords[i,1]
						temp.yrL <- yrL*coords[i,3]+coords[i,2]
						tempLsq1 <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
						temp.ygap <- temp.yrL[2]-temp.yrL[1]
						tempLsq <- sum((tempLsq1/(coords[i,3]^4))*temp.ygap^2)
					}
					fL2 <- append(fL2,tempLsq)
				}		
				S1rzK <- (1/(as.vector(t(rs$f$qhz))[corrGridSpec]^2))*(2*fk2 + 0.25*fL2)
				
				
				if(comment) cat("--g--\n\n")
				#S2rzK <- (1/(as.vector(t(rs$g$qhz))[corrGridSpec]^2))*(2*gk2 + 0.25*apply(as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$g$hypoH))[corrGridSpec]))),1,areaLsq,WIN=rs$g$WIN,iter=4))
				coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$g$hypoH))[corrGridSpec])))
				gL2 <- c()
				for(i in 1:nrow(coords)){
					if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
						tempLsq <- NA
					} else {
						temp.xrL <- xrL*coords[i,3]+coords[i,1]
						temp.yrL <- yrL*coords[i,3]+coords[i,2]
						tempLsq <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
						temp.ygap <- temp.yrL[2]-temp.yrL[1]
						tempLsq <- sum((tempLsq/(coords[i,3]^4))*temp.ygap^2)
					}
					gL2 <- append(gL2,tempLsq)
				}
				S2rzK <- (1/(as.vector(t(rs$g$qhz))[corrGridSpec]^2))*(2*gk2 + 0.25*gL2)
			} else {
				if(comment) cat("calculating integrals L2...\n--f--\n")
				S1rzK <- (1/(as.vector(t(rs$f$qhz))[corrGridSpec]^2))*(2*fk2 + .5*fk2)
				if(comment) cat("--g--\n\n")
				S2rzK <- (1/(as.vector(t(rs$g$qhz))[corrGridSpec]^2))*(2*gk2 + .5*gk2)
			}
			
			denominator <- sqrt(((S1rzK*rs$f$gamma^2)/(nrow(rs$f$data)*rs$f$globalH^2))+((S2rzK*rs$g$gamma^2)/(nrow(rs$g$data)*rs$g$globalH^2)))
		} else {
			if(comment) cat("\n--Fixed-bandwidth asymptotics--\n")
			if(comment) cat("calculating integrals K2...")
		
			#k2fix <- getQhz_Fixed(grx,gry,"gaus",pooled$WIN,pooled$pilotH,F,T)$qhz_sq
			k2fix <- (1/(4*pi))*as.vector(run_ppp(data=pooled$data,xy=datarange.list,h=(sqrt(0.5*pooled$pilotH^2)),WIN=pooled$WIN)$edg$v)
			
			if(comment) cat("done.\n\n")
			h <- unique(pooled$h)
			RrzK <- k2fix/(as.vector(t(pooled$qhz))[corrGridSpec]^2)
			denominator <- sqrt(RrzK*(nrow(rs$f$data)^(-1)+nrow(rs$g$data)^(-1)))/(h*sqrt(as.vector(t(pooled$Zm)))[corrGridSpec])
		}
	} else {
		#warning("densities have not been edge-corrected. additional computation time required for components.")
		
		if(adaptive){
			if(comment) cat("\n--Adaptive-bandwidth asymptotics--\n")
			if(comment) cat("calculating integrals K and K2...\n")
		
			if(comment) cat("--f--\n")
			#fk <- apply(as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$f$hypoH))[corrGridSpec]))),1,getQhz_Adaptive,kType="gaus",WIN=rs$f$WIN,both=T)
			hypoQuan <- unique(quantile(as.vector(t(rs$f$hypoH))[corrGridSpec],(1:100)/100,na.rm=T))
			corrQuan <- apply(as.matrix(as.vector(t(rs$f$hypoH))[corrGridSpec]),1,idQuan,q=hypoQuan)
			qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$f$data,xy=datarange.list,WIN=rs$f$WIN)
			fk <- rep(-1,gsize*gsize)
			fk[is.na(as.vector(t(rs$f$hypoH))[corrGridSpec])] <- NA
			for(i in 1:length(hypoQuan)) fk[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
			
			hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
			corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
			qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$f$data,xy=datarange.list,WIN=rs$f$WIN)
			fk2 <- rep(-1,gsize*gsize)
			fk2[is.na(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec])] <- NA
			for(i in 1:length(hypoQuan)) fk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
			fk2 <- (1/(4*pi))*fk2 
			
			
			if(comment) cat("--g--\n")
			#gk <- apply(as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$g$hypoH))[corrGridSpec]))),1,getQhz_Adaptive,kType="gaus",WIN=rs$g$WIN,both=T)
			hypoQuan <- unique(quantile(as.vector(t(rs$g$hypoH))[corrGridSpec],(1:100)/100,na.rm=T))
			corrQuan <- apply(as.matrix(as.vector(t(rs$g$hypoH))[corrGridSpec]),1,idQuan,q=hypoQuan)
			qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$g$data,xy=datarange.list,WIN=rs$g$WIN)
			gk <- rep(-1,gsize*gsize)
			gk[is.na(as.vector(t(rs$g$hypoH))[corrGridSpec])] <- NA
			for(i in 1:length(hypoQuan)) gk[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
			
			hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
			corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
			qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$g$data,xy=datarange.list,WIN=rs$g$WIN)
			gk2 <- rep(-1,gsize*gsize)
			gk2[is.na(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec])] <- NA
			for(i in 1:length(hypoQuan)) gk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
			gk2 <- (1/(4*pi))*gk2 
			
			if(exactL2){
				xrL <- sort(rep(seq(-4,4,length=10),10))
				yrL <- rep(seq(-4,4,length=10),10)
				grL <- matrix(c(xrL,yrL),100,2)
				Lsq_gr <- apply(grL,1,Lsq,uh=c(0,0,1),WIN=NULL)
				
				if(comment) cat("calculating integrals L2...\n--f--\n")
				coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$f$hypoH))[corrGridSpec])))
				fL2 <- c()
				for(i in 1:nrow(coords)){
					if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
						tempLsq <- NA
					} else {
						#h <<- coords[i,3]
						temp.xrL <- xrL*coords[i,3]+coords[i,1]
						temp.yrL <- yrL*coords[i,3]+coords[i,2]
						tempLsq1 <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
						temp.ygap <- temp.yrL[2]-temp.yrL[1]
						tempLsq <- sum((tempLsq1/(coords[i,3]^4))*temp.ygap^2)
					}
					fL2 <- append(fL2,tempLsq)
				}		
				S1rzK <- (1/(fk^2))*(2*fk2 + 0.25*fL2)
				
				
				if(comment) cat("--g--\n\n")
				coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$g$hypoH))[corrGridSpec])))
				gL2 <- c()
				for(i in 1:nrow(coords)){
					if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
						tempLsq <- NA
					} else {
						temp.xrL <- xrL*coords[i,3]+coords[i,1]
						temp.yrL <- yrL*coords[i,3]+coords[i,2]
						tempLsq <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
						temp.ygap <- temp.yrL[2]-temp.yrL[1]
						tempLsq <- sum((tempLsq/(coords[i,3]^4))*temp.ygap^2)
					}
					gL2 <- append(gL2,tempLsq)
				}
				S2rzK <- (1/(gk^2))*(2*gk2 + 0.25*gL2)
			} else {
				if(comment) cat("calculating integrals L2...\n--f--\n")
				S1rzK <- (1/(fk^2))*(2*fk2 + .5*fk2)
				if(comment) cat("--g--\n\n")
				S2rzK <- (1/(gk^2))*(2*gk2 + .5*gk2)
			}
			
			denominator <- sqrt(((S1rzK*rs$f$gamma^2)/(nrow(rs$f$data)*rs$f$globalH^2))+((S2rzK*rs$g$gamma^2)/(nrow(rs$g$data)*rs$g$globalH^2)))
		} else {
			if(comment) cat("\n--Fixed-bandwidth asymptotics--\n")
			if(comment) cat("calculating integrals K and K2...")
		
			#k2fix <- getQhz_Fixed(grx,gry,"gaus",pooled$WIN,pooled$pilotH,T,F)
			kfix <- as.vector(run_ppp(data=pooled$data,xy=datarange.list,h=pooled$pilotH,WIN=pooled$WIN)$edg$v)
			k2fix <- (1/(4*pi))*as.vector(run_ppp(data=pooled$data,xy=datarange.list,h=(sqrt(0.5*pooled$pilotH^2)),WIN=pooled$WIN)$edg$v)
			
			if(comment) cat("done.\n\n")
			
			RrzK <- k2fix/(kfix^2)
			denominator <- sqrt(RrzK*(unique(rs$f$h)^(-2)*nrow(rs$f$data)^(-1)+unique(rs$g$h)^(-2)*nrow(rs$g$data)^(-1)))/(sqrt(as.vector(t(pooled$Zm)))[corrGridSpec])
		}
	}
	
	if(rs$log) {
		numerator <- as.vector(t(rs$rsM))[corrGridSpec]
	} else {
		numerator <- as.vector(t(rs$rsM))[corrGridSpec]-1
	}
	Zstandard <- numerator/denominator
	
	if(test=="greater"){
		P <- pnorm(Zstandard,lower.tail=F)
	} else if (test=="less"){
		P <- pnorm(Zstandard,lower.tail=T)
	} else {
		P <- 2*pnorm(abs(Zstandard),lower.tail=F)
	}
	if(comment) print(date())
    return(list(X=seq(xr[1],xr[2],length=gsize),Y=seq(yr[1],yr[2],length=gsize),Z=matrix(Zstandard,gsize,gsize,byrow=T),P=matrix(P,gsize,gsize,byrow=T)))
}
