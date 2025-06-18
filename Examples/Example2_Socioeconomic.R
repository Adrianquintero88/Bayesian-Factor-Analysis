
library(rjags)
library(dclone)  # To run MCMC in
library(snow)    # several cores
library(factor.switching) #Post-processing for loadings

setwd("F:/Documentos/PhD/Factorial_analysis/code/JEBS/Supplementary material")

### Call functions
source("func_rotate.R")

### Creating multi-cluster for JAGS
cl <- makeCluster(3, type="SOCK")
tmp <- clusterEvalQ(cl, library(dclone))

### Loading data
Y <- read.table("Data_SES.txt", header = TRUE)
# colnames(Y) <- c('Father_Edu', 'Mother_Edu', 'Internet', 'Books', 'Computer',
#                  'Milk', 'Meat', 'Fruits')
# Y <- Y[,c('Mother_Edu', 'Father_Edu', 'Internet', 'Books', 'Computer','Milk', 'Meat', 'Fruits')]
# write.table(Y, file = 'Data_SES.txt', row.names = FALSE)
p <- ncol(Y)
n <- nrow(Y)

### Standardize data
for(i in 1:p)
  Y[,i] <- (Y[,i]-mean(Y[,i]))/sd(Y[,i])

###
### IDIFA m=4
###

m <- 4
F <- matrix(rnorm(n*3), ncol=3)

inits_f <- function(){
  list(tau=rep(10,p), gamma=c(1,1,1,0), nu=c(1,1,1,0), eta=cbind(F,1))
}

data_f <- list(N=n, y=Y, nresp=p, zero=rep(0,m), ident=diag(m), m=m)

parJagsModel(cl, name="mcmcM", file="JAGS_IDIFA.txt", data=data_f, n.chains=3, inits=inits_f)
parUpdate(cl, "mcmcM", n.iter=5000)
model_s <- parCodaSamples(cl, model="mcmcM", n.iter=10000, thin=10, 
                          variable.names=c("dgamma","gamma","sig2","lam2","remove"))

###
### Results
###

### Selecting the number of factors
slabs <- unlist(model_s[,"dgamma"])
remiter <- unlist(model_s[,"remove"]) 
freqm <- table(slabs[remiter==0]) #Remove iterations with non-identified models
freqm/sum(freqm) #Posterior probabilities

#Factor loadings estimates with post-processing
lam_est <- fun_rotate(model_s, 3)

summ.mat <- summary(lam_est)
facload <- summ.mat$statistics[,1]
facload <- matrix(facload, ncol = 3, byrow = TRUE)
facload

