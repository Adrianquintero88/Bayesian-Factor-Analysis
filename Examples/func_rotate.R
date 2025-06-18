
fun_rotate <- function(model_s, nfac)
{
  nchains <- length(model_s)
  niter <- length(model_s[[1]][,"dgamma"])

  lam_sel <- vector('list',length=nchains)
  indCol <- vector('list',length=nchains)
  slabs <- vector('list',length=nchains)
  
  #Selecting iterations with selected number of factors
  for(chain in 1:nchains)
  {
    slabs[[chain]] <- as.matrix(model_s[[chain]][,"dgamma"])
    
    lam_selii <- array(0, dim = c(p,m,niter))
    for(ii in 1:p)
      for(jj in 1:m)
        lam_selii[ii, jj,] <- 
      as.matrix(model_s[[chain]][,paste0("lam2[",ii,",",jj,"]")])
    
    lam_sel[[chain]] <- lam_selii[,,slabs[[chain]]==nfac]
    
    indColii <- matrix(0, nrow = niter, ncol = m)
    for(jj in 1:m)
      indColii[,jj] <- as.matrix(model_s[[chain]][,paste0("gamma[", jj, "]")])
    
    indCol[[chain]] <- indColii[slabs[[chain]]==nfac,]
    
  }
  
  #Factor loadings matrix across iterations
  lam_mat <- vector('list',length=nchains)
  iters <- rep(0, nchains)
  for(chain in 1:nchains)
  {
    iters[chain] <- dim(lam_sel[[chain]])[3]
    lam_red <- matrix(0, nrow = iters[chain], ncol = nfac*p)
    for(kk in 1:iters[chain])
    {
      count_col <- 0
      for(jj in 1:m)
      {
        if(indCol[[chain]][kk,jj] == 1)
        {
          count_col <- count_col + 1
          lam_red[kk,(p*(count_col-1)+1):(p*count_col)] <- 
            lam_sel[[chain]][, jj, kk]
        }
      }
    }
    lam_mat[[chain]] <- lam_red
    colnames(lam_mat[[chain]]) <- paste(paste0('LambdaV', rep(1:p, nfac)), 
                                        sort(rep(1:nfac, p)), sep = '_')
  }
  
  #Factor loadings matrix across iterations
  lam_ord <- vector('list',length=nchains)
  for(chain in 1:nchains)
  {
    lam_ord[[chain]] <- matrix(0, nrow = nrow(lam_mat[[chain]]), 
                               ncol = ncol(lam_mat[[chain]]))
    colnames(lam_ord[[chain]]) <- paste(paste0('LambdaV', sort(rep(1:p, nfac))), 
                                        1:nfac, sep = '_')
    for(jj in colnames(lam_ord[[chain]]))
      lam_ord[[chain]][, jj] <- lam_mat[[chain]][, jj]
    
    #Taking same number of MCMC iterations for all chains
    lam_ord[[chain]] <- lam_ord[[chain]][1:min(iters), ]
  }
  
  orderPost <- vector('list',length=nchains)
  for(chain in 1:nchains){
    orderPost[[chain]] <- rsp_exact(lambda_mcmc = lam_ord[[chain]])
  }

  #Same solution across chains: signs change
  orderPost <- compareMultipleChains(rspObjectList=orderPost)
  
  return(orderPost)
  
}
