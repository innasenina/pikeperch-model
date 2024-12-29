args <- commandArgs(TRUE)
wdir <- args[1] #model code directory
opt  <- args[2] #execution option


#CONTROL FLAGS
do.optim <- FALSE
estimate <- FALSE
bootstrap<- FALSE
validate <- FALSE
forecast <- FALSE
scenario <- FALSE

if (opt == "")
	do.optim <- TRUE		

if (opt == "bootstrap")
	bootstrap <- TRUE

if (opt == "validate"){
	do.optim <- TRUE
	validate <- TRUE
}

if (opt == "forecast")
	forecast <- TRUE

if (opt == "scenario")
	scenario <- TRUE

if (opt == "estimate"){
	do.optim <- TRUE
	bootstrap<- TRUE
}

if (opt == "historical"){
	do.optim  <- TRUE		
	bootstrap <- TRUE
	validate  <- TRUE
}

if (opt == "all"){
	do.optim  <- TRUE		
	bootstrap <- TRUE
	validate  <- TRUE
	forecast  <- TRUE
}

setwd(wdir)

#Initialize graphical and misc parameters
verbose  <- TRUE
show.progress <- FALSE
plot.figures  <- TRUE
mycol  <- 4
lwd.pp <- 1.5
lwd.af <- 3
cex.ylab <- 1.5
#cc.cols <- c(rgb(1,.5,.5),rgb(.5,.5,1),rgb(.5,1,0.5))
cc.cols <- c(rgb(1,.25,.25),rgb(.45,.65,.9,.7),rgb(.5,.75,.5))

if (scenario) plot.figures <- FALSE

source("R/funcs.R")  

figs.dir <- paste0("figs/",c("fit","cc","paper"))
for (n in 1:length(figs.dir)){
	if (!dir.exists(figs.dir[n])) 
		dir.create(figs.dir[n],recursive=TRUE)
}

if (scenario) source("R/run-and-get-prob.R")
if (!scenario) 
	#main file with execution section
	source("R/pp-model.R")


