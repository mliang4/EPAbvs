# Wrapper function for the core EPA variable selection function in cpp
bvsEPA = function(
  iter = NULL,
  n_i = NULL,
  n_betn = NULL,
  n_shuff = NULL,
  flex_covas = "default", # flexible covariates (trivial cluster allowed)
  clus_covas = "default", # flexible covariates for clustering update
  between_subj = NULL, # flexible subjects (to perform between steps)
  Y = NULL,
  X = NULL,
  T = NULL,
  Ustar = NULL,
  Ustar_dems = NULL,
  dist = NULL,
  betainit = NULL, 
  xiinit = NULL, 
  sinit = NULL,
  gammainit = NULL,
  alphainit = NULL,
  deltainit = NULL,
  gammadeltainit = NULL,
  etainit = NULL,
  sigmainit = NULL,
  muinit = NULL,
  lambda2init = NULL,
  tau2init = NULL,
  a_alp = 5,
  b_alp = 1,
  a_del = 2,
  b_del = 2,
  a_gam = 1,
  b_gam = 1,
  a_et = 1,
  b_et = 1,
  simfunc = 2,
  model = 'EPA',
  seed = 123,
  zval = 0.995
){
  
  set.seed(seed)
  
  N = length(n_i)
  P = ncol(X)
  U <- matrix( rep( T, ncol( X) ), ncol = ncol( X) )
  
  if( is.null(n_betn) ){ n_betn = 1 }
  if( is.null(n_shuff)){ n_shuff = 5 }
  if( length(flex_covas) > 0 && flex_covas == "default" ){ 
    flex_covas = c(1:(P-1)) 
  }
  
  if( length(clus_covas) > 0 && clus_covas == "default" ){ 
    clus_covas = c(1:(P-1)) 
  }
  
  if( is.null(between_subj) ){ 
    between_subj = c(0:(N-1)) 
  }
  
  subject = NULL
  for(i in 1:N){
    subject = c(subject,rep(i-1,n_i[i]))
  }
  
  subject_dems = 0
  for(i in (1:N)-1){
    subject_dems = c(subject_dems, subject_dems[i]+length(which(subject==i)))
  }
  subject_dems = c(subject_dems,length(subject))
  
  onesP = rep(1,P);ones3P = rep(1,3*P)
  
  # Matrix of barX. Rows indexed by ij. Columns indexed by 2P ( x1u, x1, x2u, x2,...xPu, xP )
  Xbar <- numeric()
  for( p in 1:ncol( U ) ){
    tmp <- X[ , p]*U[ , p ]
    Xbar <- cbind( Xbar, tmp, X[ , p ] )
  }
  
  #### make Ustar and Ustar_dem if not provided
  
  # Ustar - Matrix of spline functions. Rows indexed by ij. Columns indexed by sum r_p over p. ( Should be a list if r_p != r_p' )
  # Ustar_dems - Input vector of dimensions for each spline function. Element is equal to starting indicies for corresponding p. Last term is number of columns. Ustar_dems[ 0 ] = 0 and length = P + 1
  # Allows for only one input for U
  
  if( is.null(Ustar) && is.null(Ustar_dems) ){
    if( ncol( as.matrix( U ) ) == 1 ){
      Ustar <- numeric()
      Ustar_dems <- c( 0 )
      tmp <- sm( U, rankZ = zval )
      for( p in 1:P ){
        Ustar <- cbind( Ustar, tmp )
        Ustar_dems <- c( Ustar_dems, ( Ustar_dems[ p ] + ncol( tmp ) ) )
      }
      U <- matrix( rep( U, P ), ncol = P , nrow = length( Y ) )
    }else{
      Ustar <- numeric()
      Ustar_dems <- c( 0 )
      for( p in 1:ncol( U ) ){
        tmp <- sm( U[ , p ], rankZ = zval )
        Ustar <- cbind( Ustar, tmp )
        Ustar_dems <- c( Ustar_dems, ( Ustar_dems[ p ] + ncol( tmp ) ) )
      }
    }
  }
  
  if( is.null(betainit) ){ betainit = matrix(0, nrow = N, ncol = 3*P) }
  if( is.null(xiinit) ){ xiinit = matrix(data = 1, nrow = N, ncol = ncol(Ustar)) }
  if( is.null(sinit) ){ sinit = matrix(data = 0, nrow = N, ncol = P) }
  if( is.null(gammainit) ){ 
    gammainit = matrix(data = 0, nrow = N, ncol = P)
    gammainit[,1] = 1
  }
  
  winit <- rep(0.1, nrow(X))
  if( is.null(alphainit) ){ alphainit <-  rep(1, P) }
  if( is.null(deltainit) ){ deltainit <-  rep(0, P) }
  if( is.null(gammadeltainit) ){ gammadeltainit <- rep(0, P) }
  if( is.null(etainit) ){ etainit <- rep(0, P) }
  if( is.null(sigmainit) ){ sigmainit <- seq(0, N-1)  }
  if( is.null(muinit) ){ muinit <- rep(1,ncol(Ustar)) }
  
  if( is.null(lambda2init) ){ lambda2init = matrix(rep(ones3P, each = N),nrow = N) } # N by 3P all-one matrix 
  if( is.null(tau2init) ){ tau2init = onesP } # P-dim all-one vector
  nulambda2init = matrix(rep(ones3P, each = N),nrow = N)
  nutau2init = onesP
  
  rtime = Sys.time()
  
  if(model == 'EPA'){
    out = bvsEPAcpp(
      iterations = iter,
      n_clust_betn = n_betn,
      n_shuffle = n_shuff, 
      between_covas = flex_covas,
      cluster_covas = clus_covas,
      subject = subject,
      subject_dems = subject_dems,
      Y = Y,
      Xbar = as.matrix(Xbar),
      Ustar = as.matrix(Ustar),
      Ustar_dems = Ustar_dems,
      dist = as.matrix(dist),
      subj_clust_free = between_subj,
      beta_init = betainit, 
      xi_init = xiinit, 
      s_init = sinit,
      gamma_init = gammainit,
      w = winit, 
      alpha = alphainit,
      delta = deltainit,
      gamma_delta = gammadeltainit,
      eta = etainit,
      sigma = sigmainit,
      mu = muinit,
      lambda2_init = lambda2init,
      tau2_init = tau2init,
      nu_lambda_init = nulambda2init,
      nu_tau_init = nutau2init,
      a_alpha = a_alp,
      b_alpha = b_alp,
      a_delta = a_del,
      b_delta = b_del,
      a_gamma = a_gam,
      b_gamma = b_gam,
      a_eta = a_et,
      b_eta = b_et,
      ftype = simfunc
    )
    
    names( out ) <- c( "beta", "xi", "s", "gamma", "lambda2", "tau2", "nu_lambda", "nu_tau", "alpha", "delta", "eta", "sigma" )
  } else if(model == 'DP'){
    out = bvsPYcpp(
      iterations = iter,
      n_clust_betn = n_betn,
      n_shuffle = n_shuff, 
      between_covas = flex_covas,
      cluster_covas = clus_covas,
      subject = subject,
      subject_dems = subject_dems,
      Y = Y,
      Xbar = as.matrix(Xbar),
      Ustar = as.matrix(Ustar),
      Ustar_dems = Ustar_dems,
      subj_clust_free = between_subj,
      beta_init = betainit, 
      xi_init = xiinit, 
      s_init = sinit,
      gamma_init = gammainit,
      w = winit, 
      alpha = alphainit,
      delta = deltainit,
      mu = muinit,
      lambda2_init = lambda2init,
      tau2_init = tau2init,
      nu_lambda_init = nulambda2init,
      nu_tau_init = nutau2init,
      a_gamma = a_gam,
      b_gamma = b_gam,
      a_alpha = a_alp,
      b_alpha = b_alp,
      a_delta = a_del,
      b_delta = b_del,
      DP = TRUE
    )
    
    names( out ) = c( "beta", "xi", "s", "gamma", "lambda2", "tau2", "nu_lambda", "nu_tau", "alpha", "delta" )

  } else if(model == 'PY'){
    
    out = bvsPYcpp(
      iterations = iter,
      n_clust_betn = n_betn,
      n_shuffle = n_shuff, 
      between_covas = flex_covas,
      subject = subject,
      subject_dems = subject_dems,
      Y = Y,
      Xbar = as.matrix(Xbar),
      Ustar = as.matrix(Ustar),
      Ustar_dems = Ustar_dems,
      subj_clust_free = between_subj,
      beta_init = betainit, 
      xi_init = xiinit, 
      s_init = sinit,
      gamma_init = gammainit,
      w = winit, 
      alpha = alphainit,
      delta = deltainit,
      mu = muinit,
      lambda2_init = lambda2init,
      tau2_init = tau2init,
      nu_lambda_init = nulambda2init,
      nu_tau_init = nutau2init,
      a_gamma = a_gam,
      b_gamma = b_gam,
      a_alpha = a_alp,
      b_alpha = b_alp,
      a_delta = a_del,
      b_delta = b_del,
      DP = FALSE
    )
    
    names( out ) = c( "beta", "xi", "s", "gamma", "lambda2", "tau2", "nu_lambda", "nu_tau", "alpha", "delta" )
    
  }
  
  
  mcmc_time = Sys.time()-rtime
  
  
  return( list( mcmc = out, mcmc_time = mcmc_time, T = U, Tstar = Ustar, Tstar_dems = Ustar_dems, iter = iter, n_i = n_i ) )
  
}

# Evaluate clustering and clustering performance (if truth provided)
clustering = function(
  clust_samples = NULL,
  sel_samples = NULL, # clustering based on active samples (gamma=1) only if provided
  true_clust = NULL,
  loss = "VI.lb" # one of "binder.psm", "VI.lb" (default), "omARI.approx"
){
  
  library(salso)
  library(mcclust)
  
  N = dim(clust_samples)[1]; P = dim(clust_samples)[2]; iter = dim(clust_samples)[3]

  if( !is.null(sel_samples) ){ ## if clustering is inferred on active samples only
    
    gaSub = apply(sel_samples,c(1,2),mean) ## subject-level MPPI after burn-in
    clus = list() ## clustering with active sample pairs only
    set.seed(123)
    for(p in 1:P){
      clus[[p]] = list()
      subSel = which(gaSub[,p] > 0.5)
      if(length(subSel) == 0){
        next
      }else if(length(subSel) == 1){
        clus[[p]][[1]] = 1
        clus[[p]][[2]] = subSel
        next
      }

      psm = diag(length(subSel))
      for(i in c(1:(nrow(psm)-1))){
        for(j in c((i+1):nrow(psm))){
          subi = subSel[i]; subi_ga = sel_samples[subi, p, ]
          subj = subSel[j]; subj_ga = sel_samples[subj, p, ]
          clusi = clust_samples[subi, p, ]
          clusj = clust_samples[subj, p, ]
          both_active = which(subi_ga * subj_ga == 1)
          denom = length(both_active)
          same_clus = which(clusi[both_active] == clusj[both_active])
          numer = length(same_clus)
          psm[i,j] = psm[j,i] = numer/denom
        }
      }
      clus[[p]][[1]] = salso(psm, loss=loss, maxZealousAttempts = 0, probSequentialAllocation = 1)
      clus[[p]][[2]] = subSel
    }

    s = matrix(data = 0, nrow = N, ncol = P)
    for(p in 1:P){
      if(length(clus[[p]]) == 0){next}
      s[clus[[p]][[2]],p] = clus[[p]][[1]]
    }
    
  }else{ ## clustering is inferred on all samples
    
    estClus = NULL
    for(p in 1:P){
      bnmatrix = clust_samples[, p, ]
      estClus = rbind(estClus, salso(t(bnmatrix)))
    }
    s = t(estClus)
    
  } ## end clustering
  
  ## evaluate performance (if provided)
  clus_perf = NULL
  if( !is.null(true_clust)){
    for(p in 1:P){
      eval_VI = vi.dist(true_clust[,p], s[,p])
      eval_arandi = arandi(true_clust[,p], s[,p], adjust = TRUE)
      clus_perf = rbind(clus_perf,c( p, eval_VI, eval_arandi ))
    }
    colnames(clus_perf) = c("p", "VI", "ARI")
  }
  
  return( list(clus = s, clus_perf = clus_perf) )
}

# Evaluate selection and selection performance (if truth provided)

selPerf = function(estimate, truth){
  Gindex1 = which(truth==1)
  Gindex0 = which(truth==0)
  pos = which(estimate==1)
  neg = which(estimate==0)
  
  FP = sum(pos %in% Gindex0)
  FN = sum(neg %in% Gindex1)
  TP = sum(pos %in% Gindex1)
  TN = sum(neg %in% Gindex0)
  
  # fpr
  fpr = FP/(FP+TN)
  # fnr
  fnr = FN/(FN+TP)
  
  tpr = 1 - fnr
  tnr = 1 - fpr
  
  # Matthews correlation coefficient
  denom = sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  
  mcc = (TP*TN - FP*FN)/denom
  
  return(list(TPR = tpr, TNR = tnr, FPR = fpr, FNR = fnr, MCC = mcc))
}

gSelect2d = function(
  sel_samples, # N by P by iter matrix
  truth, # N by P matrix
  threshold = 0.5
){
  
  P = dim(truth)[2]
  selePerf = NULL
  mppi = apply(sel_samples,c(1,2),mean)
  sout = mppi > threshold
  
  if( !is.null(truth) ){
    for(p in 1:P){
      temp = as.data.frame(selPerf(sout[,p], truth[,p]))
      selePerf = rbind(selePerf, temp)
      rownames(selePerf)[p] = p
    }
  }

  return(list( g_sel = sout, g_mppi = mppi, g_perf = selePerf))

}


lSelect = function(bvsEPA_out = NULL, s = NULL, burnin = NULL, plot = FALSE){
  
  bebe = bvsEPA_out$mcmc$beta
  gaga = bvsEPA_out$mcmc$gamma
  iter = dim(bebe)[3]
  bdf = lSel = NULL
  
  for(sf in 1:5){
    
    sfbeta = (sf-1)*3+(1:3)
    uniques = unique(s[,sf])
    
    for(us in uniques){
      
      if(us == 0){next}
      
      slist = which(s[,sf] == us) # all subjects in cluster us
      blist = NULL
      
      for(subj in slist){
        
        gtemp = gaga[subj,sf,(burnin+1):iter]
        
        for(i in 1:length(gtemp)){
          if(gtemp[i] == 0){
            gtemp[i] = NA
          }
        }
        
        gtemp = rbind(gtemp,gtemp,gtemp)
        
        btemp = bebe[subj,sfbeta,(burnin+1):iter]*gtemp
        blist = rbind(blist, as.numeric(apply(btemp, 1, quantile, c(0.025,0.5,0.975), na.rm = T)))
        
      } ## end for slist
      
      blist = colMeans(blist)
      lSel = rbind(lSel, c(sf, us, blist[1]*blist[3] >= 0, blist[4]*blist[6] >= 0, blist[7]*blist[9] >= 0))
      bout = rbind(blist[1:3], blist[4:6], blist[7:9])
      title = as.data.frame(cbind(rep(paste0("Smooth function",sf, "\nCluster", us),3), c("Nlin", "Linr", "Main")))
      bout = cbind(title, bout)
      bdf = rbind(bdf, bout)
    }
  }
  colnames(bdf) = c("Cluster", "Type", "LL", "LogOddsRatio", "UL")
  colnames(lSel) = c("Smooth function", "Cluster", "Non-linear", "Linear", "Main")
  
  if(plot){
    print(
      ggplot(bdf, aes(x = Cluster, y = LogOddsRatio, colour = Type)) +
        geom_errorbar(width = .2, aes(ymax = UL, ymin = LL), position = position_dodge(width = .3), size = 1) +
        geom_point(position = position_dodge(0.3), aes(shape=Type), size = 2) + 
        geom_hline(yintercept= 0, linetype="dashed", size = 1)
    )
  }
  
  return(lSel)
}


# Plot each estimated smooth function (with truth included if provided)

plotSf <- function( 
  bvsEPA_out = NULL, 
  s = NULL, 
  sf = NULL,
  sfClus = NULL,
  plot_subj = NULL, 
  burnin = NULL, 
  xlab = NULL, 
  ylab = NULL, 
  main = NULL, # Title of the figure
  # exp = FALSE, # False: to plot the curves in log scale, True: original scale
  cl = NULL, # Customized color scheme for the curves
  rm_legend = F,  # True: remove legend in the plot
  truth = NULL # A function: if provided the truth curve will be added to the figure
){
  
  if( !is.null(sfClus) && !is.null(plot_subj) ){
    stop("Can only plot the smooth function for a single subject or for a cluster, please remove either 'sfClus' or 'plot_subj'")
  }

  iter = bvsEPA_out$iter
  U = bvsEPA_out$T
  Ustar = bvsEPA_out$Tstar
  Ustar_dems = bvsEPA_out$Tstar_dems
  if( is.null(burnin) ){ burnin = 0 }
  n_i = bvsEPA_out$n_i
  
  bebe = bvsEPA_out$mcmc$beta
  gaga = bvsEPA_out$mcmc$gamma
  xixi = bvsEPA_out$mcmc$xi
  
  if( is.null(sfClus) ){
    sfClus = unique(s[,sf])
  }
  
  rstart = c(1,cumsum(n_i)+1)
  rstart = rstart[-length(rstart)] ## starting index for each subject
  rend = cumsum(n_i) ## ending index for each subject
  
  if( is.null( xlab ) ){
    xlab <- "Time"
  }
  # if( is.null( ylab ) & exp ){
  #   ylab <- "Odds Ratio"
  # }
  # if( is.null( ylab ) & !exp ){
  #   ylab <- "Log odds Ratio"
  # }
  if( is.null( ylab ) ){
    ylab <- "Log odds Ratio"
  }
  if( is.null( main ) ){
    main <- paste0("Smooth Function ", sf)
  }
  yint = 0
  # if( exp ){ yint = 1 }
  
  nd1 = NULL
  
  if( !is.null(plot_subj) ){ ## plot single subject

    subj = plot_subj

    xtmp = Ustar_dems[sf+1] - Ustar_dems[sf]
    bind = ((sf-1)*3+1):(sf*3)
    xind = ((sf-1)*xtmp+1):(sf*xtmp)
    
    post.UL = NULL
    
    post.beta = bebe[subj,bind,]
    post.gamma = gaga[subj,sf,]
    post.xi = xixi[subj,xind,]
    
    nobs = length(rstart[subj]:rend[subj])
    
    post.smline = matrix(ncol = nobs, nrow = iter - burnin)
    for(i in (burnin+1):iter){
      if(post.gamma[i] != 0){
        post.smline[i-burnin,] = post.beta[1,i] * Ustar[(rstart[subj]:rend[subj]),xind]%*%(post.xi[,i]) +
          post.beta[2,i] * U[(rstart[subj]:rend[subj]),sf] + 
          post.beta[3,i]
      }
    }

    post.UL = post.smline
    
    x = U[rstart[subj]:rend[subj]]
    post.ul = apply(post.UL,2,quantile,c(0.025,0.5,0.975),na.rm = T)
    y = post.ul[2,]
    loci = post.ul[1,]
    hici = post.ul[3,]
    Curve = paste0("Subject", plot_subj)
    nd1 = rbind(nd1, cbind(x,y,loci,hici,Curve))
    
    nd1 = as.data.frame(nd1)
    nd1[,1:4] = as.numeric(as.matrix(nd1[,1:4]))
    colnames(nd1)[2:4] = c('y','loci','hici')
    
    p <- ggplot(nd1, aes(x = x, y = y, colour = Curve)) +
      geom_smooth( aes(ymin = loci, ymax = hici, fill = Curve ), stat = "identity", alpha = 0.25) +
      xlab(xlab) +
      ylab(ylab) + 
      geom_hline(yintercept=c(yint), linetype="dotted", size=1) + 
      ggtitle(paste0("Time-varying effect for ",main)) +
      theme(plot.title = element_text(hjust = 0.5)) + guides(fill = "none") + labs( colour = "" ) 
    
  }else{ ## plot the average for entire cluster
    
    for(sfc in sfClus){
      sfClusSubj = which(s[,sf] == sfc)
      
      if(sfc == 0){
        next
      }
      
      xtmp = Ustar_dems[sf+1] - Ustar_dems[sf]
      bind = ((sf-1)*3+1):(sf*3)
      xind = ((sf-1)*xtmp+1):(sf*xtmp)
      
      for(subj in sfClusSubj){
        post.beta = bebe[subj,bind,]
        post.gamma = gaga[subj,sf,]
        post.xi = xixi[subj,xind,]
        
        nobs = length(rstart[subj]:rend[subj])
        
        post.smline = matrix(ncol = nobs, nrow = iter - burnin)
        for(i in (burnin+1):iter){
          if(post.gamma[i] != 0){
            post.smline[i-burnin,] = post.beta[1,i] * Ustar[(rstart[subj]:rend[subj]),xind]%*%(post.xi[,i]) +
              post.beta[2,i] * U[(rstart[subj]:rend[subj]),sf] + 
              post.beta[3,i]
          }
        }
        # post.UL = rbind(post.UL, colMeans(post.smline, na.rm = T))
        x = U[rstart[subj]:rend[subj]]
        post.ul = apply(post.smline,2,quantile,c(0.025,0.5,0.975), na.rm = T)
        y = post.ul[2,]
        loci = post.ul[1,]
        hici = post.ul[3,]
        Cluster = paste0("Cluster", sfc)
        nd1 = rbind(nd1, cbind(x,y,loci,hici,Cluster))
        
      }
      
    }
    
    nd1 = as.data.frame(nd1)
    nd1[,1:4] = as.numeric(as.matrix(nd1[,1:4]))
    colnames(nd1)[2:4] = c('y','loci','hici')
    
    p <- ggplot(nd1, aes(x = x, y = y)) +
      geom_smooth( aes(ymin = loci, ymax = hici, fill = Cluster, colour = Cluster ), stat = "identity", alpha = 0.25) +
      xlab(xlab) +
      ylab(ylab) + 
      geom_hline(yintercept=c(yint), linetype="dotted", size=1) + 
      ggtitle(paste0("Time-varying effect for ",main)) +
      theme(plot.title = element_text(hjust = 0.5)) + guides(fill = "none") + labs( colour = "" ) 
    
  }
  
  if( !is.null(cl) ){ p <- p + scale_fill_manual(values=c(cl)) + scale_color_manual(values = c(cl)) }
  if( rm_legend ){ p <- p + theme(legend.position = "none")  }
  if( !is.null(truth) ){ p <- p + stat_function(fun = truth, aes(colour = "Truth"))   }
  
  print(p)

}