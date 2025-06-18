
library(rjags)
library(dclone) # To run MCMC in
library(snow)   # several cores

### Creating clusters to run a chain in each core
detectCores()
cl <- makeCluster(3, type="SOCK")
tmp <- clusterEvalQ(cl, library(dclone))

###
### Creating data
###

delta <- cbind(c(0.87, 0.00, 0.58, 0.53, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
               c(0.00, 0.67, 0.50, 0.00, 0.45, 0.64, 0.00, 0.40, 0.00, 0.00),
               c(0.00, 0.00, 0.60, 0.00, 0.00, 0.00, 0.72, 0.00, 0.50, 0.83))

Sigma <- diag(c(0.2431, 0.5511, 0.0536, 0.7191, 0.7975, 0.5904, 0.4816, 
                0.8400,0.7500, 0.3111))

n <- 100
p <- 10
k <- ncol(delta)

m <- 4 #maximum number of factors

nchain <- 3 #Number of MCMC chains

set.seed(15+5)
F <- matrix(rnorm(n*k), ncol=k)

Y <- matrix(0, nrow=n, ncol=p)
for(i in 1:n)
{
  set.seed(15+10*i)
  Y[i,] <- delta%*%F[i,] + rnorm(p, sd=sqrt(diag(Sigma)))
}

###
### BF via path sampling
###

inits_f <- function(){  list(tau=rep(10,p))}

u <- matrix(0,nrow=10, ncol=m) #Components U_bar to compute BF

for(d in 1:m)
  for(t in 1:10) #Grid points
  {
    data_f <- list(N=n, y=Y, nresp=p, zero=rep(0,d), ident=diag(d), m=d, 
                   t=(t-1)/9)
    
    parJagsModel(cl, name="mod", 
                 file=paste0("JAGS_PATH_p10_m",d-1,"tom",d,".txt"), 
                 data=data_f, n.chains=nchain, inits=inits_f)
    parUpdate(cl, "mod", n.iter=1000, thin=1)
    model_s <- parCodaSamples(cl, model="mod", variable.names=c("ubar"), 
                              n.iter=5000, thin=1)
    
    u[t,d] <- mean(unlist(model_s[,"ubar"]))
  }

logBF <- rep(0,m)
for(d in 1:m)
  for(t in 1:9)
    logBF[d] <- logBF[d] + 0.1/2*(u[t+1,d]+u[t,d])

logBF_0 <- rep(0,m)
for(d in 1:m)
  logBF_0[d] <- sum(logBF[1:d])

BF_0 <- exp(c(0,logBF_0))

BF_0 #Evidence for dimensions 0, 1, 2, ..., m
which.max(BF_0)-1 #Selected number of factors
  
