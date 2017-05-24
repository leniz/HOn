rm(list=ls())


library(mutoss)
library(knitr)
source('~/R/codeANDdata_pruebas/hierBH.R')
source('~/R/codeANDdata_pruebas/functions.R')
source('~/R/codeANDdata_pruebas/hierasymptkappa.R')

# Setting ----------------------------------------------------------------------------
# bs = "100HxFam" # number of setting 
CONFIG = 1  # 1 or 2, structure of the families
n.supfam <- 3
m <- n.supfam*4*100
k <- 12 #48  # k <- 12
alpha = 0.05
kappa =0.001

# x all null families, the first one containing only 100 true null hyp.
pi = array(1/k,k)
pinull <- array(.8,k)
pinull[c(1:2)] <- .1
# pinull[c(3:4)] <- .5
B <- 100
mu<-c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5) 


# FDR computation  ------------------------------------------------
# set.seed()
source('~/R/codeANDdata_pruebas/fdr_computation.R')

# results --
fam <- data.frame(unlist(famFDR, recursive = FALSE))
names(fam) <- c(array("mu_1",dim=3),array("mu_2",dim=3),array("mu_3",dim=3),
                array("mu_4",dim=3),array("mu_5",dim=3),array("mu_6",dim=3),
                array("mu_7",dim=3),array("mu_8",dim=3),array("mu_9",dim=3),
                array("mu_10",dim=3),array("mu_11",dim=3))
nfam <- dim.data.frame(fam) [1]
z <- cbind(array("fam",dim=nfam),1:nfam)
row.names(fam) <- paste(z[,1],z[,2])

# Standard Deviation
SEs = data.frame(mu,SEofdrHO, SEofdrHO3, SEofdrHBH)
names(SEs) <- c("mu","HO","HO3","HBH")

settings = list(n.supfam = n.supfam, k = k, alpha = alpha, kappa = kappa, pi = pi, 
                pinull = pinull, B = B, mu = mu)
save(settings,fam,famFDR,ofdrHBH,ofdrHO,ofdrHO3,
     powerHBH,powerHO,powerHO3,
     SEofdrHBH,SEofdrHO,SEofdrHO3,
     file= paste("~/R/codeANDdata_pruebas/Results/FDR config_",bs, ".RData",sep=""))
# names(fam) <- array(rbind(mu,mu,mu), dim = c(1,33))
# fam <- rbind(fam, apply(fam,2,mean)) # mean FDR
# write(print(xtable(fam[, 10:21], digits = 5),scalebox = 0.7), 
#       file = "~/R/codeANDdata_pruebas/Results/table config_pba.tex")

# PLOT FDR --
pdf(file = paste("~/R/codeANDdata_pruebas/Results/FDR conf_",bs,".pdf",sep=""))  # open PDF
plot(mu,ofdrHO,main="",xlab=expression(mu),ylab="gFDR",type="l",ylim=c(0,.1))
points(mu,ofdrHO,pch=1)
lines(mu,ofdrHO3,lty=1, col="red")
points(mu,ofdrHO3,pch=2, col="red")
# lines(mu,ofdrHO3a,lty=1, col="orange")
# points(mu,ofdrHO3a,pch=5, col="orange")
lines(mu,ofdrHBH,lty=1, col="blue")
points(mu,ofdrHBH,pch=4, col="blue")
lines(mu, array(0.05,length(mu)),lty = 1, col ="green")

legend("topright",
       c(expression(varphi^{HO}), expression(varphi^{HO3}), expression(varphi^{Bog})),
       lty=c(1,1,1),
       pch=c(1,2,4),
       merge=TRUE, horiz=TRUE)
dev.off()

# PLOT POWER --
pdf(file = paste("~/R/codeANDdata_pruebas/Results/power conf_",bs,".pdf",sep=""))  # open PDF
plot(mu,powerHO,main="",xlab=expression(mu),ylab="power",type="l",ylim=c(0,1))
points(mu,powerHO,pch=1)
lines(mu,powerHO3,lty=1, col="red")
points(mu,powerHO3,pch=2, col="red")
# lines(mu,ofdrHO3a,lty=1, col="orange")
# points(mu,ofdrHO3a,pch=5, col="orange")
lines(mu,powerHBH,lty=1, col="blue")
points(mu,powerHBH,pch=4, col="blue")


legend("topleft",
       c(expression(varphi^{HO}), expression(varphi^{HO3}), expression(varphi^{Bog})),
       lty=c(1,1,1),
       pch=c(1,2,4),
       merge=TRUE)
dev.off()

#=================================================================================#
