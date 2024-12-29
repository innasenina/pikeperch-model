#To run this script with different projection parameters:

write(c("N0.intro","T.offset","S.offset","p"),"R/res/probabilities.txt",ncol=4,sep=";");

Tvals <- seq(0,3,1) 
Svals <- seq(0,3,1) 
for (N0.intro in c(0,10))
    for (T.offset in Tvals)
		for (S.offset in Svals)
			source("R/pp-model.R")
		


#Then read the file and re-arrange it into a Table
#Run	 |  Tref  |  Tref-2  |  Tref-3
#-------------------------------------
#Sref	 |   p1	...
#Sref-2  | ...
#Sref-3	 |

out <- read.csv("R/res/probabilities.txt",sep=";")
N0.int <- unique(out$N0.intro)

for (n in 1:2){
    ind <- which(out$N0.intro==N0.int[n])
    res <- out[ind,]
    T.off <- res$T.offset; 
    S.off <- res$S.offset; 
    Tu <- sort(unique(T.off))
    Su <- sort(unique(S.off))
    nT <- length(Tu)
    nS <- length(Su)
    p <- res$p

    tab <- matrix(0,nrow=nS,ncol=nT,
		  dimnames=list(rownames=paste0("S-",Su),
				colnames=paste0("T-",Tu)))

    for (i in 1:nS)
	for (j in 1:nT)
	    tab[i,j] <- p[which(S.off==Su[i] & T.off==Tu[j])]
    print(tab)
}
message("ALL done.")

		 
