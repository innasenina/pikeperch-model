# The pike-perch (M1) model published in Tyutyunov et al., EcoMod (2002).
 
#********** PARAMETER SECTION **********#

#1. Model parameters

#Initial conditions for model variables in 1949
Nn.1949 <- list(recruits=25, yearlings=19.309, adults=36.5951)

#Fixed or variables b?
b.fixed <- TRUE
#Optimal value of intraspecific competition parameter
#obtained by 7-parameter optimization with jitter run strategy:
fixed.b <- .03949 

#Survival of adults (not estimated)
p2<-0.622

#Optimal values of salinity and temperature
S0<-11; T0<-10

#Fishing effort: same for yearlings and adults? 
h.unique <- TRUE
#Fishing effort during 1981-2018 when neither C.obs nor N1.obs and N2.obs are available
h.missing <- 0.05
#conversion factor from milion fish to thous.tons
n2w <- 1.285714
#Assumed catch composition
c1 <- 0.475; 
#mean weights of yearlings and adults
w1 <- 0.5; w2 <- (n2w-c1*w1)/(1-c1) #w2=2.5 in EcoMod paper
if (h.unique) w2 <- 1.5

#If account for unreported and illegal fishing, set C.iuu to the parts of legal catch.
C.iuu <- 0
if (C.iuu>0) h.approx <- TRUE
fishing.ban.year <- 2017

#These are parameter names used in the code
par.names <- c("a1","a2","b","sigS","sigT","p0","p1")

#2. Optimization and bootstrap parameters
#number of jitter runs:
nbj <- 1	
#parameter bounds
a1.b <- c(0,20); a2.b <- c(0,20)
b.b  <- c(0,1);  p0.b <- c(0,1); p1.b <- c(0,1)
sigS.b <- c(0.1,3); sigT.b <- c(0.1,3); 

#number of bootstrap realisations
nbs <- 10000
if (nbs>100 & bootstrap) verbose <- FALSE

#where bootstrap output is written
bootstrap.file  <- paste0(wdir,"res/bootstrap-out.txt")

#3. CC scenario parameters

#Local control parameter: means running additional to base scenarios
run.scenarios <- scenario

#Number of MCMC x CC runs
#If nmc <= 2, it's treated as no MCMC, but only Effort runs
#Meaningless to put 1 as in this case only h=0 will be used
nmc <- 10000
# number of CC models
nbm <- 200;
#Set lower abundance value for stock restoration
adu.cr   <- 1.38 # mln according to Zhivoglyadov&Loukyanov-2018

#Data files
file.in <- paste0(wdir,"data/pike-perch-data.csv")
file.cc <- paste0(wdir,"data/base-scenario-2021-2030.csv")

#***** END OF PARAMETER SECTION ********#


#********** EXECUTION SECTION **********#
if (run.scenarios)
	write.prob <- TRUE		

if (!run.scenarios){	
	N0.intro <- 0
	T.offset <- 0
	S.offset <- 0
}
#Local model parameter used by function pp.model
INTRO <- N0.intro 

if (b.fixed) {
	par.names <- par.names[-3]
	b <- fixed.b
}
nbp <- length(par.names)
p.bound <- sapply(mget(paste0(par.names,".b")),list)

dat<-read.csv(file.in,sep=";",nrows=32)
  
Y <- dat$Year; 
nY<- length(Y)
S <- dat$Salinity
T <- dat$Temperature
C <- dat$Catch*(1+C.iuu)
N0.obs <- dat$N0.obs 
N1.obs <- dat$N1.obs
N2.obs <- dat$N2.obs
N0.0 <- Nn.1949$recruits; 
N1.0 <- Nn.1949$yearlings; 
N2.0 <- Nn.1949$adults

if (h.unique){
	h1 <- C/(N1.obs+N2.obs);
	h2 <- h1
} else {
	h1 <- c1*C/N1.obs
	h2 <- (1-c1)*C/N2.obs
}

if (do.optim){
  
	message("ENTERING PARAMETER ESTIMATION...")	
	Lmin <- 1e5	  
	for (n in 1:nbj){

	x.lower <- get.nth(p.bound,1)
	x.upper <- get.nth(p.bound,2)

	x.init<-randn(x.lower,x.upper)


	res<-optim(x.init,L,method = "L-BFGS-B",
		   hessian=TRUE,control=list(maxit=300),
		   lower=x.lower,upper=x.upper)
	L.n<-L(res$par)

	print(c(n,L(x.init),L.n,Lmin))
	  
	if (L.n<Lmin){
		Lmin<-L.n
		xmin<-res$par
		Hmin<-res$hessian
		print(res)
	}
	} #Optimization done
	
	print(Lmin); 
	print(xmin)
	#Errors and correlations:
	Cov<-chol2inv(chol(Hmin))
	sigmas<-sqrt(diag(Cov))
	Cor <- cov2cor(Cov)
	cor <- Cor; diag(cor) <- NA; cor[upper.tri(cor)] <- NA; 
	if (any(abs(cor)>sqrt(0.5),na.rm=TRUE)){
		message("!!!Found pairs of highly correlated parameters!!! 
		They are: ")
	print(which(abs(cor)>sqrt(0.5),arr.ind=TRUE))
	}

	#Compute an error for the parameter estimates
	out<-cbind(xmin,sqrt(diag(solve(Hmin))))
	colnames(out)<-c("Estimate","Std.Error")
	print(out)

	x.est <- xmin
}
#Done fitting

if (!do.optim){

	#Optimal solution, L = 17261.7
	x.est <- c(6.7745815,12.5939987,1.2636094,1.9210239,0.5451080,0.4712031)
	if (!fixed.b){ # L = 17261.87 (Biophysics)
		x.est <- c(6.82171026,12.66743120,0.03967993,1.26311400,1.92381950,0.54442711,0.47134960)
	}
}

mod <- pp.model(x.est)
N0 <- mod$recruits; 
N1 <- mod$yearlings; 
N2 <- mod$adults
C.est <- h1*w1*N1+h2*w2*N2

if (!run.scenarios){
	r.N0<-round(cor(N0.obs,N0),2)
	r.N1<-round(cor(N1.obs,N1),2)
	r.N2<-round(cor(N2.obs,N2),2)
	r.N1.N2<-round(cor(N1+N2,N1.obs+N2.obs),2)
	r.C<-round(cor(n2w*C,C.est),2)

	message("Fit metrics with estimated parameters:")
	message("Pearson correlations: \nr(N0) = ",r.N0,"; r(N1) = ", 
			r.N1,"; r(N2) = ", r.N2, "; r(C) = ", r.C)
	message("NRMSE: \nnrmse(N0) = ",nrmse(N0.obs,N0),"; nrmse(N1) = ", 
			nrmse(N1.obs,N1),"; nrmse(N2) = ", nrmse(N2.obs,N2), 
			"; nrmse(C) = ", nrmse(n2w*C,C.est))
}

if (do.optim & plot.figures){
	# Plotting the results with estimated parameters
	png(paste0(wdir,"figs/fit/pp-model-fit.png"),800,1000)
	par(mfcol=c(3,1),mar=c(2,3,2,2),cex=1.25,las=1,mgp=c(3,0.5,0))
	nice.plot2(Y,N0.obs,N0,paste0("Juveniles (r =",r.N0,")"))
	nice.plot2(Y,N1.obs,N1,paste0("Yearlings (r =",r.N1,")"))
	nice.plot2(Y,N2.obs,N2,paste0("Adults (r = ",r.N2,", r(N1,N2) = ", r.N1.N2,")"))
	dev.off()
}
  
if (bootstrap){

	message("ENTERING BOOTSTRAP, will be running ",nbs," realisations...")	

	write(c("N","Lopt",par.names,paste0("err.",par.names)),
	  bootstrap.file,sep=";",ncolumns=14)
	
	x.bts<-array(NA,c(nbs,length(x.est)))
	x.err<-x.est; x.err[]<-0		 

	mod<-pp.model(x.est)
	N0<-mod$recruits; N1<-mod$yearlings; N2<-mod$adults;

	#residuals
	resid.N0<-N0.obs-N0
	resid.N1<-N1.obs-N1
	resid.N2<-N2.obs-N2

	if (!verbose) pb<-txtProgressBar(style=3)

	for (i in 1:nbs){
		if (!verbose) setTxtProgressBar(pb,i/nbs)			
		#Warning: modifying global variables, 
		#restore true obs once bootstrap is done!
		my.sample<-sample(1:nY,replace=TRUE)  
		N0.obs<-N0+sign(rnorm(1))*resid.N0[my.sample]
		my.sample<-sample(1:nY,replace=TRUE)  
		N1.obs<-N1+sign(rnorm(1))*resid.N1[my.sample]
		my.sample<-sample(1:nY,replace=TRUE)  
		N2.obs<-N2+sign(rnorm(1))*resid.N2[my.sample]

		#Now, run optimizations with perturbed data
		Lmin<-1e5	  
		for (n in 1:nbj){

			x.lower <- get.nth(p.bound,1)
			x.upper <- get.nth(p.bound,2)
			x.init<-x.est*(1+rnorm(1,0,0.1))


			res<-optim(x.init,L,method = "L-BFGS-B",
				   hessian=TRUE,control=list(maxit=300),
				   lower=x.lower,upper=x.upper)
			L.n<-L(res$par)
			if (verbose)
				print(c(n,L(x.init),L.n))

			if (L.n<Lmin){
				Lmin<-L.n
				xmin<-res$par
				Hmin<-res$hessian
			}
		}
		if (verbose){
			print(Lmin)
			print(xmin)
		}

		x.bts[i,]<-xmin
		if (i==1){
			x.err  <- sqrt(diag(solve(Hmin)))
			x.mean <- xmin
		}
		if (i>=2){
			x.err  <- apply(x.bts,2,sd,na.rm=TRUE)	
			x.mean <- apply(x.bts,2,mean,na.rm=TRUE)
		}

		if (verbose){
			message("Sample #",i,"...")
			print(x.mean)
			print(x.err)
		}
		write.table(cbind(i,round(Lmin,2),t(x.mean),t(x.err)),
				file=bootstrap.file,append=TRUE,sep=";",quote=FALSE,
				col.names=FALSE,row.names=FALSE)
	}
	if (!verbose) close(pb)


	message("Final bootstrap statistics:")
	message("Mean parameters:  ", paste(round(x.mean,2),collapse=" "))
	message("Parameter errors: ", paste(round(x.err,2),collapse=" "))
	message("Relative errors00:  ", paste(round(x.err/x.mean,2),collapse=" "))
  
	#Restoring observations after bootstrap
	N0.obs <- dat$N0.obs
	N1.obs <- dat$N1.obs
	N2.obs <- dat$N2.obs
}

#Read bootstrap file. It must exist, otherwise, run bootstrap first
#Here using the file 'bootstrap.txt' pointing to the file with largest 
#number of realisations, otherwise take the most recent bootstrap.file
bfile <- gsub("-out","",bootstrap.file)
if (!file.exists(bfile))
	bfile <- bootstrap.file
if (!file.exists(bfile))
	stop("Error: no bootstrap file found. Run the bootstrap first. Exit now...")

message("Reading parameters with uncertainties from file ",bfile)
out  <- read.csv(bfile,sep=";")
nbr  <- nrow(out)

if (!bootstrap){
	#Use previously computed parameters
	x.mean <- out[nbr,3:8]
	x.err <- out[nbr,9:14]
}

if (plot.figures){
	#Plot bootstrap results
	png(paste0(wdir,"figs/fit/bootstrapping.png"),400,1200)
	par(mfcol=c(nbp,1),mar=c(2,4,2,1))
	for (par in par.names){
		ylim <- range(out[,par],x.est[which(par.names==par)])
		plot(1:nbr,out[,par],type="l",ylim=ylim,ylab=par);
		abline(h=x.est[which(par.names==par)],lty=2,col=2)
	}
	dev.off()
}

if (forecast | validate | run.scenarios){

	if (validate)
		message("RUNNING MODEL VALIDATION WITH EFFORT AND PARAMETER MCMC...")	
	# read all data, including verification period, until 2020
	dat.all <- read.csv(file.in,sep=";")
	ind.20  <- which(dat.all$Year<=2020)
	dat.all <- dat.all[ind.20,]

	Y <- dat.all$Year; 
	S <- dat.all$Salinity; 
	T <- dat.all$Temperature; 
	#C.obs  <- dat.all$Catch; # <- catch data collected before November 2024
	C.obs  <- (1+C.iuu)*dat.all$Catch.upd2024; # <- added new catch data in November 2024
	N0.obs <- dat.all$N0.obs; 
	N1.obs <- dat.all$N1.obs; 
	N2.obs <- dat.all$N2.obs

	N2.Belousov<-dat.all$N2.B2004
	N2.Podoinitsin<-dat.all$N2.Pod
	N2.Protokoly<-dat.all$N2.Prot

	nY <- length(Y)
	ind <- which(Y>1981)
	if (h.unique){
		h  <- C.obs/(N1.obs+N2.obs);
		h[ind] <- h.missing
		h1 <- h; h2 <- h
	} else {
		h1 <- c1*C.obs/N1.obs
		h2 <- (1-c1)*C.obs/N2.obs
		h1[ind] <- h.missing
		h2[ind] <- h.missing
	}	

	h.min <- min(h2)
	h.max <- max(h2)
	
	h1[which(Y>fishing.ban.year)] <- 0.0; 
	h2[which(Y>fishing.ban.year)] <- 0.0; 
	h1.ref <- h1; h2.ref <- h2
	# re-run the model with all data, without any scenarios
	INTRO <- 0
	mod <- pp.model(x.est)
#	mod <- pp.model(as.numeric(x.mean))
print(x.mean)
	N0 <- mod$recruits; 
	N1 <- mod$yearlings; 
	N2 <- mod$adults

	#Run MCMC with parameter uncertainty estimated by bootstrapping
	arr   <- array(0,c(nmc,nY));
	if (nmc == 2) {
		h.2vals <- c(h.min,h.max)
	}
	x.mc  <- array(0,c(nmc,nbp))
	N.mc  <- list(recruits=arr, yearlings=arr, adults=arr)
	C.est.mc <- arr
	for (n in 1:(nmc)){
		if (nmc<=2) {
			x <- x.est
			h1[ind] <- h.2vals[n]
			h2[ind] <- h.2vals[n]
		}
		if (nmc>2){
			x <- as.numeric(x.mean + sapply(x.err,function(x)
											sign(rnorm(1))*randn(0,x)))
			h1[ind] <- randn(h.min,h.max)
			h2[ind] <- randn(h.min,h.max)
		}
		res <- pp.model(x)
	
		N.mc$recruits[n,] <- res$recruits
		N.mc$yearlings[n,] <- res$yearlings
		N.mc$adults[n,] <- res$adults
		
		C.est.mc[n,] <- h1*w1*N.mc$yearlings[n,]+h2*w2*N.mc$adults[n,] # in thous tons

		#remember x.mc for CC simulations
		x.mc[n,] <- x
	}
	N0.mean <- apply(N.mc$recruits,2,mean)
	N1.mean <- apply(N.mc$yearlings,2,mean)
	N2.mean <- apply(N.mc$adults,2,mean)
	C.est <- apply(C.est.mc,2,mean)
	if (validate){

		#Model validation without CC but with effort scenarios
		#and parameters MCMC using bootstrap-estimated means and errors 			
		# Evaluate the model fit to data
		r.N0 <- rcoef(N0.obs,N0)
		r.N1 <- rcoef(N1.obs,N1)
		r.N2 <- rcoef(N2.obs,N2)
		r.C  <- round(cor(n2w*C.obs,C.est,use="complete.obs"),2)
		#Same for average stats
		r.N0m <- rcoef(N0.obs,N0.mean)
		r.N1m <- rcoef(N1.obs,N1.mean)
		r.N2m <- rcoef(N2.obs,N2.mean)	
		
		e.N0 <- nrmse(N0.obs,N0); e.N1 <- nrmse(N1.obs,N1)
		e.N2 <- nrmse(N2.obs,N2); e.C  <- nrmse(n2w*C.obs,C.est)

		e.N0m <- nrmse(N0.obs,N0.mean); 
		e.N1m <- nrmse(N1.obs,N1.mean)
		e.N2m <- nrmse(N2.obs,N2.mean); 

		message("Pearson correlations: \nr(N0) = ",r.N0m,"(",r.N0,"); r(N1) = ", 
				r.N1m,"(",r.N1,"); r(N2) = ", r.N2m, "(",r.N2,"); r(C) = ", r.C)
		message("NRMSE: \nnrmse(N0) = ",e.N0m,"(",e.N0,"); nrmse(N1) = ", 
				e.N1m,"(",e.N1,"); nrmse(N2) = ", e.N2m, "(",e.N2,
				"); nrmse(C) = ", e.C)

		#plots will be done whatever the plot.figures is
		val.file <- paste0(wdir,"figs/fit/pp-model-validation.png")
		message("Plotting validation to file ",val.file)

		png(val.file,800,1200)
		layout(1:4)
		
		plot.var(Y,N0.obs,N0,N.mc$recruits,r.N0,e.N0,"Recruits",
				 var.units=expression(italic(N)[0]~(10^6)),120)
		lines(Y,N0.mean,lwd=2,col=2)
		plot.var(Y,N1.obs,N1,N.mc$yearlings,r.N1,e.N1,"Yearlings",
				 var.units=expression(italic(N)[1]~(10^6)),60)
		lines(Y,N1.mean,lwd=2,col=2)
		plot.var(Y,N2.obs,N2,N.mc$adults,r.N2,e.N2,"Adults",
				 var.units=expression(italic(N)[2]~(10^6)),50)
		lines(Y,N2.mean,lwd=2,col=2)
		abline(h=adu.cr,lty=2)	
		plot.var(Y,n2w*C.obs,C.est,C.est.mc,r.C,e.C)

		dev.off()

		fig3.file <- paste0(wdir,"figs/paper/Fig-3.pdf")
		message("Plotting validation to file ",fig3.file)
		pdf(fig3.file,15,6); 
		save.cex <- cex.ylab; cex.ylab <- 2.2
		plot.var(Y,n2w*C.obs,C.est,C.est.mc,r.C,e.C,show.scores=FALSE)
		dev.off()
		cex.ylab <- save.cex

		if (opt!="all")
			return("Done estimation and validation.")
	}
	
	message("RUNNING MODEL WITH FORECAST, INCLUDING CC SCENARIOS...")	

	# Initial conditions in 2020:
	i.2020 <- which(Y==2020)
	N0.0<-N0[i.2020]+INTRO; 
	N1.0<-N1[i.2020]; 
	N2.0<-N2[i.2020]
  
	#Remember historical time vector and environment
	Y.hist <- Y
	T.hist <- T
	S.hist <- S

	#CLIMATE CHANGE SIMULATIONS
	# read CC data
	dat.cc <- read.csv(file.cc,sep=";")
	# read years from the first model (N=0)
	yy.cc  <- dat.cc$Year[which(dat.cc$N==0)]

	T.cc <- matrix(dat.cc$Temperature,ncol=nbm)
	S.cc <- matrix(dat.cc$Salinity,ncol=nbm)

	#CC scenarios with MCMC with model parameters

	#to show progress bar:
	if (nmc>=100) show.progress <- TRUE
	message("1. CCMC with base CC scenarios...")

	nY.cc <- length(yy.cc)
	N.cc  <- list()
	# Compute model with CC data:
	INTRO <- N0.intro
	message("INTRO = ",N0.intro,
			"; Toffset = ",T.offset,
			"; Soffset = ",S.offset)
	if (show.progress) pb<-txtProgressBar(style=3)
	Y <- yy.cc; 
	nY <- nY.cc
	for (m in 1:nbm){
		if (show.progress) setTxtProgressBar(pb,m/nbm)			
		T <- T.cc[,m]-T.offset
		S <- S.cc[,m]-S.offset
		#Run CC with MCMC parameters, so we now have CCMC
		for (n in 1:nmc)
			N.cc[[m + (n-1)*nbm]] <- pp.model(x.mc[n,])
	}
	
	if (show.progress) close(pb)

	if (!run.scenarios){	
		#Add two more scenarios: 
		#1) with optimistic scenario for environment
		#2) with (1)+larval introduction
		T.offset <- 2; S.offset <- 3; N.cc2 <- N.cc; N.cc3 <- N.cc
		message("\n2. CCMC with additional scenarios...")
		message("INTRO = (0,10)",
				"; Toffset = ",T.offset,
				"; Soffset = ",S.offset)
		if (show.progress) pb<-txtProgressBar(style=3)
		for (m in 1:nbm){
			if (show.progress) setTxtProgressBar(pb,m/nbm)			

			INTRO <- 0; 
			T <- T.cc[,m]-T.offset
			S <- S.cc[,m]-S.offset
			for (n in 1:nmc)
			N.cc2[[m + (n-1)*nbm]] <- pp.model(x.mc[n,])

			INTRO <- 10
			for (n in 1:nmc)
				N.cc3[[m + (n-1)*nbm]] <- pp.model(x.mc[n,])
		}
		if (show.progress) close(pb)
	}
	# Finally, concatenate historical and the mean CC variable:
	Y  <- c(Y.hist,yy.cc)
	nY <- length(Y.hist) + nY.cc
	N0.cc.ave <- lst.ave(N.cc,"recruits");  N0 <- c(N0,N0.cc.ave) 
	N1.cc.ave <- lst.ave(N.cc,"yearlings"); N1 <- c(N1,N1.cc.ave) 
	N2.cc.ave <- lst.ave(N.cc,"adults");	N2 <- c(N2,N2.cc.ave) 

	# Just for visualization, update Nobs vectors
	N0.obs <- c(N0.obs,rep(NA,nY.cc))
	N1.obs <- c(N1.obs,rep(NA,nY.cc))
	N2.obs <- c(N2.obs,rep(NA,nY.cc))
	# Also, add 0 to N.mc for 2021-2030
	# so it will have the full length
	# but will be empty on polygon plot
	arr   <- array(0,c(nmc,nY.cc));
	N.mc$recruits  <- cbind(N.mc$recruits,arr)
	N.mc$yearlings <- cbind(N.mc$yearlings,arr)
	N.mc$adults	<- cbind(N.mc$adults,arr)


	prob <- get.probabilities(N.cc,!run.scenarios)
	if (!run.scenarios){
		message("Survival probabilities with base scenario\n", 
				"(above 1.38M, positive trend over last 5 years, both):")
		print(prob)
		message("Survival probabilities with T and S offsets:")
		prob.cc2 <- get.probabilities(N.cc2)
		print(prob.cc2)
		prob.cc3 <- get.probabilities(N.cc3)
		message("Survival probabilities with T and S offsets and introduction of 10M:")
		print(prob.cc3)
		
	}
	if (run.scenarios){
		message("p = ",prob[1])			
		write(cbind(N0.intro, T.offset, S.offset, round(prob[1],2)),
		  "res/probabilities.txt",append=TRUE,sep=";")
	}
  
	#***************#
	#  CC FIGURES   #
	#***************#
	if (plot.figures){	
		# Figure 1: a) T and b) S for 1950-2030	
		op <- par(no.readonly=TRUE)
	
		T.all <- c(T.hist,array(NA,nY.cc))
		S.all <- c(S.hist,array(NA,nY.cc))
		Y.cc2 <- c(2020,yy.cc)
		T.cc2 <- rbind(rep(T.all[which(Y==2020)],nbm),T.cc)
		S.cc2 <- rbind(rep(S.all[which(Y==2020)],nbm),S.cc)

		pdf("figs/paper/Fig-1a.pdf",10,5)
		nice.plot(Y,T.all,x.at=seq(1950,2030,10),y.at=seq(10,16,1),
			  cex=1,pch=22,bg="white",lwd=lwd.af,labs.cex=2,cex.axis=1.5,
			  my.title="",my.ylab=expression(italic(T)~("\u00B0"*C)))
		for (m in 1:nbm) lines(Y.cc2,T.cc2[,m],col="grey50",lwd=.5*lwd.af)
		lines(Y.cc2,rowMeans(T.cc2),lwd=lwd.af); 
		dev.off()

		cairo_pdf("figs/paper/Fig-1b.pdf",10,5) #pdf device does not print permiles
		nice.plot(Y,S.all,x.at=seq(1950,2030,10),y.at=seq(9,16,1),
			  cex=1,pch=22,bg="white",lwd=lwd.af,labs.cex=2,cex.axis=1.5,
			  my.title="",my.ylab=expression(italic(S)~("\u{2030}")))
		for (m in 1:nbm) lines(Y.cc2,S.cc2[,m],col="grey50",lwd=.5*lwd.af)
		lines(Y.cc2,rowMeans(S.cc2),lwd=lwd.af); 
		dev.off()
	
		par(op)

		# Figure 2: N0, N1 and N2
		# Base scenario with lines, optimistic and release scenarios as side boxplots
		#Note, here N0 includes the average of CC model runs
		#Also, let's add 1949 to start from the same point
		Y <- c(1949,Y); nY <- nY+1; 
		N0.49 <- as.numeric(Nn.1949$recruits); 
		N1.49 <- as.numeric(Nn.1949$yearlings); 
		N2.49 <- as.numeric(Nn.1949$adults); 
		N0.obs <- c(N0.49,N0.obs); N0 <- c(N0.mean,rep(NA,nY.cc))
		N1.obs <- c(N1.49,N1.obs); N1 <- c(N1.mean,rep(NA,nY.cc))
		N2.obs <- c(N2.49,N2.obs); N2 <- c(N2.mean,rep(NA,nY.cc))
		#and CC projections to start from 2020 
		yy.cc <- c(2020,yy.cc)

		pdf(paste0(wdir,"figs/paper/Fig-2.pdf"),9.5,12.5)
		layout(1:3,heights=c(1,1,1.2))
		nice.plot.ts.polygon(Y,N0,N.mc$recruits,"Recruits",val0=N0.49,
							 var.units=expression(italic(N)[0]~(10^6)))
		lines(Y,N0.obs,lwd=lwd.pp); points(Y,N0.obs,pch=21,cex=1.5,bg="white")
		nice.plot.ts.polygon(yy.cc,lst.ave(N.cc3,"recruits"),alst2mat(N.cc3,"recruits"),
							 val0=N0[i.2020],col=cc.cols[3],add=TRUE)
		nice.plot.ts.polygon(yy.cc,lst.ave(N.cc2,"recruits"),alst2mat(N.cc2,"recruits"),
							 val0=N0[i.2020],col=cc.cols[2],add=TRUE)
		nice.plot.ts.polygon(yy.cc,N0.cc.ave,alst2mat(N.cc,"recruits"),
							 val0=N0[i.2020],add=TRUE,col=cc.cols[1])
		text(x=2030,y=-1,expression(italic(T)[0]~italic(S)[0]),pos=4,cex=2,xpd=NA)
		text(x=2030,y=10,expression(italic(T)[2]~italic(S)[3]),pos=4,cex=2,xpd=NA)
		text(x=2030,y=20,expression(italic(T)[2]~italic(S)[3]^'+'),pos=4,cex=2,xpd=NA)
		

		nice.plot.ts.polygon(Y,N1,N.mc$yearlings,"Yearlings",val0=N1.49,
							 var.units=expression(italic(N)[1]~(10^6)),ylims=c(0,60))
		lines(Y,N1.obs,lwd=lwd.pp); points(Y,N1.obs,pch=21,cex=1.5,bg="white")
		nice.plot.ts.polygon(yy.cc,lst.ave(N.cc3,"yearlings"),alst2mat(N.cc3,"yearlings"),
							 val0=N1[i.2020],col=cc.cols[3],add=TRUE)
		nice.plot.ts.polygon(yy.cc,lst.ave(N.cc2,"yearlings"),alst2mat(N.cc2,"yearlings"),
							 val0=N1[i.2020],col=cc.cols[2],add=TRUE)
		nice.plot.ts.polygon(yy.cc,N1.cc.ave,alst2mat(N.cc,"yearlings"),
							 val0=N1[i.2020],add=TRUE,col=cc.cols[1])
		text(x=2030,y=-1,expression(italic(T)[0]~italic(S)[0]),pos=4,cex=2,xpd=NA)
		text(x=2030,y=5,expression(italic(T)[2]~italic(S)[3]),pos=4,cex=2,xpd=NA)
		text(x=2030,y=10,expression(italic(T)[2]~italic(S)[3]^'+'),pos=4,cex=2,xpd=NA)
		
	
		nice.plot.ts.polygon(Y,N2,N.mc$adults,"Adults",val0=N2.49,
							 var.units=expression(italic(N)[2]~(10^6)),
							 ylims=c(0,50),mar=c(7,7,4,3))
		lines(Y,N2.obs,lwd=lwd.pp); points(Y,N2.obs,pch=21,cex=1.5,bg="white")
		lines(Y.hist,N2.Podoinitsin); points(Y.hist,N2.Podoinitsin,cex=1.5,pch=20)
		lines(Y.hist,N2.Belousov);	points(Y.hist,N2.Belousov,pch=25,cex=1)
		lines(Y.hist,N2.Protokoly);   points(Y.hist,N2.Protokoly,pch=24,cex=1)

		nice.plot.ts.polygon(yy.cc,lst.ave(N.cc3,"adults"),alst2mat(N.cc3,"adults"),
							 val0=N2[i.2020],col=cc.cols[3],add=TRUE)
		nice.plot.ts.polygon(yy.cc,lst.ave(N.cc2,"adults"),alst2mat(N.cc2,"adults"),
							 val0=N2[i.2020],col=cc.cols[2],add=TRUE)
		nice.plot.ts.polygon(yy.cc,N2.cc.ave,alst2mat(N.cc,"adults"),
							 val0=N2[i.2020],add=TRUE,col=cc.cols[1])
		text(x=2030,y=-1,expression(italic(T)[0]~italic(S)[0]),pos=4,cex=2,xpd=NA)
		text(x=2030,y=4,expression(italic(T)[2]~italic(S)[3]),pos=4,cex=2,xpd=NA)
		text(x=2030,y=14,expression(italic(T)[2]~italic(S)[3]^'+'),pos=4,cex=2,xpd=NA)


		abline(h=adu.cr,lty=2)
		axis(1,at=c(1950,1981),labels=rep("",2),tcl=0.3,line=2.5); 
		mtext("Identification",1,line=4,at=1965,cex=1.5)
		axis(1,at=c(1982,2020),labels=rep("",2),tcl=0.3,line=2.5); 
		mtext("Validation",1,line=4,at=2000,cex=1.5)
		axis(1,at=c(2021,2030),labels=rep("",2),tcl=0.3,line=2.5); 
		mtext("Projection",1,line=4,at=2025,cex=1.5)
		dev.off()
   
	}
}
if (!run.scenarios) message("ALL done.")


