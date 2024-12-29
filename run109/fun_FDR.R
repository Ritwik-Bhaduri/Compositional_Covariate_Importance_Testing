myBH <- function(x, q){
  x_sorted = sort(x,decreasing = FALSE) 
  cutoffs = (1:length(x))*q/length(x)
  n_rejections = ifelse(any(x_sorted <= cutoffs), max(which(x_sorted <= cutoffs)), 0)
  rejections = which(x <= sort(x)[n_rejections])
  return(list("rejections" = rejections, "n_rejections" = length(rejections)))
}

Storey <- function(x, q, lambda){
  pi_0_hat <- min((1+n-sum(x<lambda))/((1-lambda)*length(x)),1)
  x_sorted = sort(x,decreasing = FALSE) 
  cutoffs = (1:length(x))*q/length(x)
  n_rejections = ifelse(any((pi_0_hat * x_sorted) <= cutoffs), 
                        max(which((pi_0_hat * x_sorted) <= cutoffs)), 0)
  rejections = which(x <= sort(x)[n_rejections])
  return(list("rejections" = rejections, "n_rejections" = length(rejections)))
}

FDR_control <- function(p, alpha = 0.1, m_1, method = "max-BH", p_filter_tolerance = 0.9, sparsity_cutoff_m_value){  
  # m_1 is required for oracle methods
  # p_filter_tolerance - 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable
  if(length(p) == 2){
    p_univariate = p$p_univariate
    p_univariate = p_univariate[1:sparsity_cutoff_m_value]
  }
  p = p$p_bivariate
  p = p[1:sparsity_cutoff_m_value, 1:sparsity_cutoff_m_value]
  m = nrow(p)
  
  if(method == "max-BH"){
    library(sgof)
    diag(p) = 0
    p_i = apply(p, 2, max)
    rejections = myBH(p_i, q = alpha)$rejections
    # rejections = which(p_i <= sort(p_i)[BH(p_i, alpha = alpha)$Rejections])
    remove(p_i)
  }
  
  if(method == "LOO-BH-BH"){
    diag(p) = 1
    n_rejections = apply(p,2, function(p_i) myBH(p_i, q = alpha)$n_rejections)
    null_var = sample(which(n_rejections==min(n_rejections)),1)
    p_i = p[,null_var]
    rejections = myBH(x = p_i, q = alpha)$rejections
    remove(p_i)
  }
  
  if(method == "LKO-Rand-BH"){
    diag(p) = 1
    k=ceiling(log(alpha, base = 0.5))+1
    selected_cols = sort(sample(1:m, size = k))
    p_i = apply(p[,selected_cols], 1, max)
    rejections = myBH(p_i, q = alpha)$rejections
    remove(p_i)
  }
  
  if(grepl("LKO-Rand-BH-changed-alpha",method, fixed = TRUE)){
    small_alpha = as.numeric(tail(strsplit(method, split='-', fixed=TRUE)[[1]],1))
    
    diag(p) = 1
    k=ceiling(log(2*alpha-small_alpha, base = 0.5))+1
    selected_cols = sort(sample(1:m, size = k))
    p_i = apply(p[,selected_cols], 1, max)
    rejections = myBH(p_i, q = small_alpha)$rejections
    remove(p_i)
  }
  if(method == "max-BH-old"){
    diag(p) = 0
    p_i = apply(p, 2, max)
    bh = sgof::BH(p_i, alpha = alpha)
    # rejections = which(p_i <= sort(p_i)[sum(bh$Adjusted.pvalues <= alpha)])
    rejections = which(bh$Adjusted.pvalues < alpha)
    remove(p_i)
  }
  
  if(method == "p-filter"){ # note that this is not an oracle method
    diag(p) = 1
    p_vec = as.vector(p)
    m = nrow(p)
    group_p = as.vector(unlist(sapply(1:m, function(j) rep(j, m))))
    groups = cbind(1:m^2, group_p)
    alphas = c(0.1, 0.1) # fist is singleton group fdr, second is our variable wise group
    Discoveries = pfilter(p_vec,alphas,groups)
    rejections = which(apply(matrix(Discoveries, nrow = m), 1, 
                             function(x) sum(x) >= (m-1)*p_filter_tolerance))
    remove(p_vec)
  }
  
  if(method == "univariate-BH"){
    if(all(is.na(p_univariate))) return(integer(0))
    
    # library(sgof)
    # bh = sgof::BH(p_univariate, alpha = alpha)
    # rejections = which(p_univariate <= sort(p_univariate)[sum(bh$Adjusted.pvalues <= alpha)])
    
    rejections = myBH(x = p_univariate, q = alpha)$rejections
  }
  if(grepl("DHH_storey_",method, fixed = TRUE)){
    lambda = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    diag(p) = 0
    p_i = apply(p, 2, max)
    
    cutoffs = rep(0.6, length(p_i))
    accepted = which(p_i > cutoffs)
    p_new = p_i[-accepted]/cutoffs[-accepted]
    # rejections_raw = myBH(p_new, q = alpha)$rejections
    rejections_raw = Storey(p_new, q = alpha, lambda = lambda)$rejections
    rejections = (1:length(p_i))[-accepted][rejections_raw]
    remove(p_i)
  }
  
  if(method == "DHH_BH"){
    diag(p) = 0
    p_i = apply(p, 2, max)
    
    cutoffs = rep(0.6, length(p_i))
    accepted = which(p_i > cutoffs)
    p_new = p_i[-accepted]/cutoffs[-accepted]
    rejections_raw = myBH(p_new, q = alpha)$rejections
    # rejections_raw = Storey(p_new, q = alpha, lambda = 0.6)$rejections
    rejections = (1:length(p_i))[-accepted][rejections_raw]
    remove(p_i)
  }
  if(method == "DHH_adaptive"){
    diag(p) = 0
    p_i = apply(p, 2, max)
    
    tau_array = seq(0.5,0.99, length.out = 50)
    x_of_tau <- function(tau=0.8,p_i){
      x_array = sort(p_i)/tau
      # x_array = x_array[which(x_array<=alpha)]
      S_tau = length(which(p_i<=tau))
      S_taux = sapply(x_array, function(x) length(which(p_i<=tau*x)))
      x_best = x_array[max(which((S_tau*x_array)/sapply(S_taux, function(y) max(y,1)) <= alpha))]
      return(x_best)
    }
    x_tau = sapply(tau_array, function(tau) x_of_tau(tau = tau, p_i))
    tau_best = tau_array[which.max(sapply(1:length(tau_array), function(i) length(which(p_i<=tau_array[i]*x_tau[i]))))]
    
    cutoffs = rep(tau_best, length(p_i))
    accepted = which(p_i > cutoffs)
    p_new = p_i[-accepted]/cutoffs[-accepted]
    rejections_raw = myBH(p_new, q = alpha)$rejections
    rejections = (1:length(p_i))[-accepted][rejections_raw]
    remove(p_i)
  }
  
  Simes <- function(x,r){
    x_ord = sort(x)
    # return(min(sapply(1:(n-r+1), function(i) ((n-r+1)/i)*x_ord[r-1+i])))
    return(min(sapply(1:((m-1)-r+1), function(i) (((m-1)-r+1)/i)*x_ord[r-1+i])))
  }
    
  
  if(grepl("Simes_oracle_",method, fixed = TRUE)){
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    p = sapply(1:nrow(p), function(i) p[i, -i])
    p_simes = apply(p, 2, function(x) Simes(x,r))
    rejections = myBH(p_simes, q = alpha)$rejections
  }
  
  if(method == "adafilter"){
    # library(devtools)
    # install_github("jingshuw/adaFilter")
    # library("adaFilter") 
    p = sapply(1:nrow(p), function(i) p[i, -i])
    result <- adafilter_BH(p, r = m-1, alpha = alpha)
    rejections = which(result$decision == 1)
  }
  
  if(grepl("adafilter_oracle_",method, fixed = TRUE)){
    # library(devtools)
    # install_github("jingshuw/adaFilter")
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    
    # library("adaFilter") 
    p = sapply(1:nrow(p), function(i) p[i, -i])
    result <- adafilter_BH(p = p, r = r, alpha = alpha)
    rejections = which(result$decision == 1)
  }
  
  if(grepl("bonferroni_oracle_",method, fixed = TRUE)){
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    p = sapply(1:nrow(p), function(i) p[i, -i])
    
    p_r_B = apply(p, 2, function(p_i) min((m-r) * sort(p_i)[r],1))
    rejections = myBH(x = p_r_B, q = alpha)$rejections
  }
  
  
  if(method == "oracle fdp control - independent"){
    m_0 = m-m_1
    diag(p) = 0
    p_i = apply(p, 2, max)
    temp = c()
    for(t in sort(c(p_i, 1-p_i))){
      # assuming all independent (best case)
      low_tail = pbeta(t, m_0-1, 1, log.p = TRUE)
      c = qbeta(low_tail, m_0-1, 1, log.p = TRUE, lower.tail = FALSE)
      fdp_hat = (1+sum(p_i > c))/max(1, sum(p_i < t))
      
      temp = c(temp, fdp_hat<alpha)
    }
    cut_off = sort(c(p_i, 1-p_i))[which(temp==0)[1]]
    rejections = which(p_i < cut_off)
  }
  
  if(method == "oracle fdp control - worst case"){ # note that this is not an oracle method, in fact it same as adapt without covariate
    diag(p) = 0
    p_i = apply(p, 2, max)
    temp = c()
    for(t in sort(c(p_i, 1-p_i))){
      # assuming all dependent (worst case)
      fdp_hat = (1+sum(p_i > (1-t)))/max(1, sum(p_i < t))
      temp = c(temp, fdp_hat<alpha)
    }
    cut_off = sort(c(p_i, 1-p_i))[which(temp==0)[1]]
    rejections = which(p_i < cut_off)
  }
  
  
  return(rejections)
}

pfilter = function(P,alphas,groups){
  # P in [0,1]^n = vector of p-values
  # alphas in [0,1]^M = vector of target FDR levels
  # groups is a n-by-M matrix; 
  #	groups[i,m] = which group does P[i] belong to,
  #		for the m-th grouping
  
  n = length(P)
  M = length(alphas)
  G = apply(groups,2,max) # G[m] = # groups, for grouping m
  
  
  Simes = list()
  for(m in 1:M){
    Simes[[m]]=rep(0,G[m])
    for(g in 1:G[m]){
      group = which(groups[,m]==g)
      Simes[[m]][g]=min(sort(P[group])*length(group)/(1:length(group)))
    }
  }
  
  
  # initialize
  thresh = alphas
  Sh = 1:n
  for(m in 1:M){
    pass_Simes_m = which(is.element(groups[,m],which(Simes[[m]]<=thresh[m])))
    Sh = intersect(Sh,pass_Simes_m)
  }
  done = FALSE
  
  
  while(!done){
    thresh_old = thresh
    for(m in 1:M){
      # which groups, for the m-th grouping, 
      #	have any potential discoveries?
      Shm = sort(unique(groups[Sh,m]))
      
      # run BH on Simes[[m]], constraining to Shm
      Pvals_m = rep(1.01,G[m]); # >1 for groups not in Dm
      Pvals_m[Shm] = Simes[[m]][Shm]
      khatm = max(0,which(sort(Pvals_m)<=(1:G[m])/G[m]*alphas[m]))
      thresh[m] = khatm/G[m]*alphas[m]
      Sh = intersect(Sh,
                     which(is.element(groups[,m],which(Pvals_m<=thresh[m]))))
    }
    if(all(thresh_old==thresh)){done = TRUE}
  }
  
  Sh_temp = Sh;
  Sh = rep(0,n); Sh[Sh_temp] = 1
  Sh
  
}

adafilter_BH <- function(p, r, alpha = 0.1){
  n = nrow(p)
  m = ncol(p)
  f = apply(p, 2, function(x) (n-r+1)*sort(x)[r-1])
  s = apply(p, 2, function(x) (n-r+1)*sort(x)[r])
  order_s = order(s)
  sorted_s = sort(s)
  m_j_af = sapply(sorted_s, function(x) sum(f < x))
  s_h.m_h.h = sorted_s * m_j_af / (1:m)
  p_adjusted = sapply(1:m, function(j) min(min(s_h.m_h.h[j:m]),1))
  n_rejections = sum(p_adjusted < alpha)
  if(n_rejections > 0){
    decision = (s <= sort(s)[n_rejections])*1
  }else{
    decision = rep(0,m)
  }
  return(data.frame("decision" = decision, "adjusted.p" = p_adjusted, 
                    selection.p = s, filter.p = f))
}
