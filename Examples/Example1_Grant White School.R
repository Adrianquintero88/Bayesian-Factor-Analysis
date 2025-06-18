
library(rjags)
library(dclone)  # To run MCMC in
library(snow)    # several cores
library(lavaan) #Library with Grant-White school data set 
library(factor.switching) #Post-processing for loadings

### Call functions
source("func_rotate.R") #Function for post-processing

#Data set
data(HolzingerSwineford1939)

### Creating clusters to run a chain in each core for the MCMC
cl <- makeCluster(3, type="SOCK")
tmp <- clusterEvalQ(cl, library(dclone))

Y <- HolzingerSwineford1939[,7:15]

m <- 5
p <- ncol(Y)
n <- nrow(Y)

for(i in 1:p)
  Y[,i] <- (Y[,i]-mean(Y[,i]))/sd(Y[,i])

data_f <- list(N=n, y=Y, nresp=p, zero=rep(0,m), ident=diag(m), m=m)

###
### IDIFA method m=5
###

parJagsModel(cl, name="mod", 
             file=paste0("JAGS_IDIFA.txt"), data=data_f, n.chains=3)
parUpdate(cl, "mod", n.iter=5000, thin=1)
model_s <- parCodaSamples(cl, model="mod", n.iter=10000, thin=10,
                          variable.names=c("dgamma","gamma","sig2","lam2","remove"))

### Selecting the number of factors
slabs <- unlist(model_s[,"dgamma"])
remiter <- unlist(model_s[,"remove"]) 
freqm <- table(slabs[remiter==0]) #Remove iterations with non-identified models
freqm/sum(freqm) #Posterior probabilities

#Factor loadings estimates with post-processing
lam_est <- fun_rotate(model_s, 3) 

### Quality of posterior sample for factor loadings
postSum <- summary(lam_est)
gelman.diag(lam_est, multivariate = FALSE)
MCSE <- batchSE(lam_est, batchSize=100)
sd <- postSum$statistics[,2]
100*MCSE/sd #All of them are below 5%
sort(effectiveSize(lam_est)) #Effective sample size
traceplot(lam_est[,"LambdaV6_2"])
traceplot(lam_est[,"LambdaV7_2"])

summ.mat <- summary(lam_est)
facload <- summ.mat$statistics[,1]
facload <- matrix(facload, ncol = 3, byrow = TRUE)

# Varimax rotation of factor loadings matrix
facload <- facload %*% varimax(facload)$rotmat
facload
