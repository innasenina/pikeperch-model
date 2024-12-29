########  FUNCTION SECTION #########

# a function to get the n-th element of each item in the list
# can be achieved also with as.numeric(as.data.frame(lst)[n,]) 
# if each element has the same dimension
get.nth <- function(lst, n) as.numeric(unlist(sapply(lst, `[`,n)))

fnorm<-function(x,mu,sig){
	
	return(exp(-((x-mu)/sig)^2))	  
}

ricker <- function(N1,N2,a1,a2,b) (a1*N1+a2*N2)*exp(-b*(N1+N2))

r.model<-function(a1, a2, b, sigS, sigT, N1, N2, S, T,S0=11,T0=10){

	N0<-(a1*N1+a2*N2)*exp(-b*(N1+N2))*fnorm(S,S0,sigS)*fnorm(T,T0,sigT)

	return(N0)
}
  
pp.model<-function(x){

	for (i in 1:length(x)) 
	assign(par.names[i],x[i])
  
	N0<-array(NA,nY); N1<-N0; N2<-N0
	 
	N0[1]<-r.model(a1,a2,b,sigS,sigT,N1.0,N2.0,S[1],T[1],S0,T0) + INTRO
	N1[1]<-p0*N0.0*(1-h1[1])
	N2[1]<-(p1*N1.0+p2*N2.0)*(1-h2[1])

	for (n in 2:nY){	  

	N0[n]<-r.model(a1,a2,b,sigS,sigT,N1[n-1],N2[n-1],S[n],T[n],S0,T0)
	
		if (n<5) N0[n] <- N0[n]+INTRO

	N1[n]<-p0*N0[n-1]*(1-h1[n])

	N2[n]<-(p1*N1[n-1]+p2*N2[n-1])*(1-h2[n])
	}

	return(list(recruits=N0, yearlings=N1, adults=N2))		
}

L <- function(x){

	mod<-pp.model(x)
	N0<-mod$recruits; N1<-mod$yearlings; N2<-mod$adults
	sum((N0-N0.obs)^2)+sum((N1-N1.obs)^2)+sum((N2-N2.obs)^2)
}

get.N1.from.N2 <- function(N2t,N2t1,h2t1,p1,p2){
	return((N2t1/(1-h2t1) - p2*N2t)/p1)
}

randn<-function(xmin, xmax) xmin+abs(rnorm(length(xmin),0,0.3*(xmax-xmin)))

lgrey <- rgb(.8,.8,.8)

rcoef <- function(obs, mod) round(cor(obs,mod,use="complete.obs"),2)
nrmse <- function(obs, mod) round(sqrt(sum((obs-mod)^2,na.rm=TRUE))/(mean(obs,na.rm=TRUE)*sqrt(length(which(!is.na(obs) & !is.na(mod))))),2)

nice.plot<-function(x,y,x.at=pretty(x),y.at=pretty(y),
			my.title,add.legend=FALSE,
			legend.text="Temperature",
			title.cex=2.5,labs.cex=1.5,
			my.ylab="",...){

	#y.max<-10*ceiling(.1*max(y,na.rm=T))

	h.lines <- y.at
	v.lines <- x.at
   
	par(las=1,mar=c(2,5,2,1),mgp=c(3,.5,0),...)

	plot(x,y,type="n",bty='n',xaxt='n',yaxt='n',
	 xlab='',ylab='',ylim=range(y.at),lwd=2); 
	axis(1,col="white"); axis(2,col="white");
	
	abline(h=h.lines,lwd=0.75,col=lgrey)
	abline(v=v.lines,lwd=0.75,col=lgrey)
	
	lines(x,y,lty=1,...); 
	points(x,y,lwd=1,cex=1.25,bg='white')

	if (add.legend)
	legend("topright",legend.text,
		   lty=1,col=col,lwd=lwd,
		   bty='n')
	par(las=0)
	mtext(my.title,3,line=0,cex=title.cex)
	mtext(my.ylab,2,line=2.5,cex=labs.cex);
}

add.metrics.text <- function(r,err){
	legend("topright",paste0("(r = ",r,"; NRMSE = ",err,")"),bty='n',cex=2.25)
}

nice.plot2<-function(Y, N.obs,N.est,my.title,plot.legend=FALSE){

	y.max<-10*ceiling(.1*max(N.obs,na.rm=T))
	if (max(N.obs,na.rm=T)<5) y.max <- 5

	h.lines<-seq(0,y.max,ifelse(y.max>=80,20,
				ifelse(y.max>5 & y.max<=10,2,
					   ifelse(y.max<=5,1,10))))	
   
	plot(Y,N.obs,type="l",bty='n',xaxt='n',yaxt='n',ylim=c(0,y.max),lwd=2); 
	
	lines(Y,N.est,lty=1,lwd=2.0,col=mycol); 
	axis(1,col="white"); axis(2,col="white");
	abline(h=h.lines,lty=3,lwd=0.75)
	if (plot.legend)
	legend("topright",c("data","model"),
		   lty=c(1,1),col=c(1,mycol),
		   lwd=c(2,2),bty='n',ncol=2)
	mtext(my.title,3,line=0,cex=2)
}

get.quantiles<-function(mat,pmin,pmax){

	nt<-ncol(mat)	  
	qnt.min<-array(NA,nt)	  
	qnt.max<-array(NA,nt)	  
	for (n in 1:nt){
	qnt.min[n]<-quantile(mat[,n],pmin,na.rm=TRUE)
	qnt.max[n]<-quantile(mat[,n],pmax,na.rm=TRUE)
	}
	return(list(qmin=qnt.min,qmax=qnt.max))
}

#Get average variable over an array of lists
#with recruits, yearlings and adults
lst.ave <- function(arr.lst,var){

	nbr <- length(arr.lst)
	nbt <- length(arr.lst[[1]][[var]])
	res <- array(0,nbt)
	for (n in 1:nbr)
	res <- res + arr.lst[[n]][[var]]/nbr

	return(res)
}

#Get average variable over an array of lists
#with recruits, yearlings and adults
alst2mat <- function(arr.lst,var){

	nbr <- length(arr.lst)
	nbt <- length(arr.lst[[1]][[var]])
	mat <- array(0,c(nbr,nbt))
	for (n in 1:nbr)
	mat[n,] <- arr.lst[[n]][[var]]

	return(mat)
}

colmove <- function(color1,color2,percentage=10) colorRampPalette(c(color1, color2))(100)[percentage]
darken  <- function(cols,p=10) sapply(cols,colmove,color2=1,percentage=p)

#ty is time as years; 
#ts1 is the one-dimensional time series, e.g. based on estimated pars 
#ts2 is two-dimensional time series, i.e., including scenarios/replications
nice.plot.ts.polygon<-function(ty,ts1,ts2,var.name="",var.units="Nb",
							   ylims=c(0,120),val0=NA,mar=c(4,7,4,3),
							   col='grey80',lwd=3,add=FALSE,...){
	if (!is.na(val0)){
		ts1 <- c(val0,ts1)
		ts2 <- t(sapply(seq_along(1:dim(ts2)[1]),function(i) c(val0,ts2[i,])))
	}

	if (!add){
		par(las=1,cex.lab=2,cex.axis=2,mar=mar)#,xpd='NA')
		plot(ty,ts1,type="n",xlab="",bty="l",ylab="",ylim=ylims,xaxt='n')
		if (is.na(val0)) myticks <- seq(ty[1],max(ty),10)
		else myticks <- seq(ty[2],max(ty),10)
		axis(1,at=myticks)

		par(cex.main=3)
		title(var.name,font=2,family="sans",cex=3)
	
		par(las=0)
		mtext(var.units,2,line=4.0,cex=cex.ylab)	
	}
	res <- get.quantiles(ts2,0.05,0.95)
	nt<-length(ty)

	polygon(c(ty,ty[nt:1]),c(res$qmin,res$qmax[nt:1]),
		col=col,border='NA',...)

	if (!add) lines(ty,ts1,lwd=2)
	if (add)  lines(ty,ts1,lwd=lwd,col=darken(col,40))
}	

plot.var <- function(yy,var.obs,var.est,var.est.mc,
					 var.cor=NA,var.err=NA,
					 var.name="Catch",var.units="Thousand tons",
					 var.max=7,show.scores=TRUE){
	
	nice.plot.ts.polygon(yy,var.est,var.est.mc,var.name,
						 var.units,ylims=c(0,var.max))
	lines(yy,var.obs); 
	points(yy,var.obs,pch=21,cex=1.5,bg="white")
	if (show.scores)
		add.metrics.text(var.cor,var.err)
}

get.probabilities <- function(N.sc,all=TRUE){

	# To compute probability of survival over all nbm x nmc runs
	last.adult.vals <- sapply(seq_along(1:(nbm*nmc)),
							  function(x)N.sc[[x]]$adults[10])
		
	p.stock <- round(mean(last.adult.vals>adu.cr),4)

	p.trend <- NA
	p.both  <- NA
	if (all){
		xn <- (nY.cc-4):nY.cc; 
		# To evaluate the trend over last five years
		a.trend <- sapply(seq_along(1:(nbm*nmc)),
						  function(x) {
								  reg <- lm(N.sc[[x]]$adults[xn]~xn); 
								  return(reg$coefficients[2])
						  })
	
		p.trend <- round(mean(a.trend>0),4)
		p.both <- round(mean(last.adult.vals>adu.cr & a.trend>0),4)
	}
	
	return(c(p.stock,p.trend,p.both))
}


##### END OF FUNCTION SECTION ######

