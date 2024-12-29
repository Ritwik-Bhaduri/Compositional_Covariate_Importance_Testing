Bonferroni <- function(x, q){
  rejections = which(x <= q/length(x))
  return(list("rejections" = rejections, "n_rejections" = length(rejections)))
}

Holm <- function(x, q){
  x_sorted = sort(x,decreasing = FALSE) 
  cutoffs = q/(length(x):1)
  t = (x_sorted < cutoffs)*1
  
  if(all(t == 1)){
    n_rejections = length(x)
  }else if(t[1]==0){
    n_rejections = 0
  }else{
    n_rejections = min(which(diff(t) != 0))
  }
  
  rejections = which(x <= sort(x)[n_rejections])
  return(list("rejections" = rejections, "n_rejections" = length(rejections)))
}

Simes <- function(x,r){
  x_ord = sort(x)
  n = length(x)
  return(min(sapply(1:(n-r+1), function(i) ((n-r+1)/i)*x_ord[r-1+i])))
}

Bonferroni_pval <- function(x,r){
  n = length(x)
  return(min((n-r+1) * sort(x)[r],1))
}

myBH <- function(x, q){
  x_sorted = sort(x,decreasing = FALSE) 
  cutoffs = (1:length(x))*q/length(x)
  n_rejections = ifelse(any(x_sorted <= cutoffs), max(which(x_sorted <= cutoffs)), 0)
  rejections = which(x <= sort(x)[n_rejections])
  return(list("rejections" = rejections, "n_rejections" = length(rejections)))
}

FWER_control <- function(p, alpha = 0.1, m_1, method = "max-BH", p_filter_tolerance = 0.9, sparsity_cutoff_m_value = NULL){  
  # m_1 is required for oracle methods
  # p_filter_tolerance - 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable
  if(length(p) == 2){
    p_univariate = p$p_univariate
    if(!is.null(sparsity_cutoff_m_value)) p_univariate = p_univariate[1:sparsity_cutoff_m_value]
  }
  p = p$p_bivariate
  if(!is.null(sparsity_cutoff_m_value)) p = p[1:sparsity_cutoff_m_value, 1:sparsity_cutoff_m_value]
  m = nrow(p)
  
  if(method == "max-Bonferroni"){
    diag(p) = 0
    p_i = apply(p, 2, max)
    rejections = Bonferroni(p_i, q = alpha)$rejections
  }
  
  if(method == "max-Holm"){
    diag(p) = 0
    p_i = apply(p, 2, max)
    rejections = Bonferroni(p_i, q = alpha)$rejections
  }
  
  # if(grepl("Simes_holm_oracle_",method, fixed = TRUE)){
  #   r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
  #   p = sapply(1:nrow(p), function(i) p[i, -i])
  #   p_simes = apply(p, 2, function(x) Simes(x,r))
  #   rejections = Holm(p_simes, q = alpha)$rejections
  # }
  
  
  if(grepl("adaptive_Simes_holm_oracle",method, fixed = TRUE)){
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    p = sapply(1:nrow(p), function(i) p[i, -i])
    p_original = p
    m = ncol(p); n = nrow(p)
    
    rejections = c()
    for(i in n:r){
      if(length(rejections) != 0){
        p = p_original[-rejections, ]
      }else{
        p = p_original
      }
      
      p_simes = apply(p, 2, function(x) Simes(x,r=r))
      if(i == n){
        if(min(p_simes) < alpha/m){
          rejections = c(rejections, which.min(p_simes))
        } else{
          break
        }
      }
      
      if(i<n){
        if(min(p_simes[-rejections]) < alpha/(m-(n-i))){
          t = which.min(p_simes[-rejections])
          new_rejection = (1:m)[-rejections][min(t)]
          rejections = c(rejections, new_rejection)
        } else{
          break
        }
      }
      
    }
    return(rejections)
  }
  
  if(grepl("adaptive_bonferroni_holm_oracle",method, fixed = TRUE)){
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    p = sapply(1:nrow(p), function(i) p[i, -i])
    p_original = p
    m = ncol(p); n = nrow(p)
    
    rejections = c()
    for(i in n:r){
      if(length(rejections) != 0){
        p = p_original[-rejections, ]
      }else{
        p = p_original
      }
      
      # p_simes = apply(p, 2, function(x) Simes_pval(x,r=r))
      p_bonferroni = apply(p, 2, function(x) Bonferroni_pval(x, r=r))
      if(i == n){
        if(min(p_bonferroni) < alpha/m){
          rejections = c(rejections, which.min(p_bonferroni))
        } else{
          break
        }
      }
      
      if(i<n){
        if(min(p_bonferroni[-rejections]) < alpha/(m-(n-i))){
          t = which.min(p_bonferroni[-rejections])
          new_rejection = (1:m)[-rejections][min(t)]
          rejections = c(rejections, new_rejection)
        } else{
          break
        }
      }
      
    }
    return(rejections)
  }
  
  if(method== "adaptive-Holm"){
    p = sapply(1:nrow(p), function(i) p[i, -i])
    p = sapply(1:ncol(p), function(i) sort(p[, i], decreasing = FALSE))
    m = ncol(p); n = nrow(p)
    rejections = c()
    for(i in n:1){
      p_simes = apply(p, 2, function(x) Simes(x,r=i))
      
      if(i == n){
        if(min(p_simes) < alpha/m){
          rejections = c(rejections, which.min(p_simes))
        } else{
          break
        }
      }
      
      if(i<n){
        if(min(p_simes[-rejections]) < alpha/(m-(n-i))){
          t = which.min(p_simes[-rejections])
          new_rejection = (1:m)[-rejections][min(t)]
          rejections = c(rejections, new_rejection)
        } else{
          break
        }
      }
      
    }
    return(rejections)
  }
  if(method== "adaptive-Holm-aggressive"){
    # p = sapply(1:nrow(p), function(i) p[i, -i])
    # p = sapply(1:ncol(p), function(i) sort(p[, i], decreasing = FALSE))
    m = ncol(p); n = nrow(p)-1
    rejections = c()
    for(i in n:1){
      if(i == n){
        p_temp = sapply(1:nrow(p), function(i) p[i, -i])
        p_simes = apply(p_temp, 2, function(x) Simes(x,r=i))
        if(min(p_simes) < alpha/m){
          rejections = c(rejections, which.min(p_simes))
        } else{
          return(rejections)
        }
      }
      
      if(i<n){
        p_temp = p[-rejections, -rejections]
        p_temp = sapply(1:nrow(p_temp), function(j) p_temp[j, -j])
        p_simes = apply(p_temp, 2, function(x) Simes(x,r=i))
        if(min(p_simes) < alpha/(m-(n-i))){
          t = which.min(p_simes)
          new_rejection = (1:m)[-rejections][min(t)]
          rejections = c(rejections, new_rejection)
        } else{
          return(rejections)
        }
      }
    }
    return(rejections)
  }
  
  
  if(method == "adafilter_bonferroni"){
    p = sapply(1:nrow(p), function(i) p[i, -i])
    result <- adafilter_bonferroni(p, r = m-1, alpha = alpha)
    rejections = which(result$decision == 1)
  }
  
  if(grepl("adafilter_bonferroni_oracle_",method, fixed = TRUE)){
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    p = sapply(1:nrow(p), function(i) p[i, -i])
    result <- adafilter_bonferroni(p = p, r = r, alpha = alpha)
    rejections = which(result$decision == 1)
  }
  
  if(method == "univariate-Bonferroni"){
    if(all(is.na(p_univariate))) return(integer(0))
    p_univariate = p_univariate[!is.na(p_univariate)]
    rejections = Bonferroni(x = p_univariate, q = alpha)$rejections
  }
  
  if(method == "univariate-Holm"){
    if(all(is.na(p_univariate))) return(integer(0))
    p_univariate = p_univariate[!is.na(p_univariate)]
    rejections = Holm(x = p_univariate, q = alpha)$rejections
  }
  
  if(method == "max-BH"){
    diag(p) = 0
    p_i = apply(p, 2, max)
    rejections = myBH(p_i, q = alpha)$rejections
  }
  
  if(method == "adafilter"){
    p = sapply(1:nrow(p), function(i) p[i, -i])
    result <- my_adafilter(p, r = m-1, alpha = alpha)
    rejections = which(result$decision == 1)
  }
  
  return(rejections)
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

adafilter_bonferroni <- function(p, r, alpha = 0.1){
  n = nrow(p)
  m = ncol(p)
  f = apply(p, 2, function(x) (n-r+1)*sort(x)[r-1])
  s = apply(p, 2, function(x) (n-r+1)*sort(x)[r])
  order_s = order(s)
  sorted_s = sort(s)
  m_j_af = sapply(sorted_s, function(x) sum(f < x))
  s_h.m_h.h = sorted_s * m_j_af / (1:m)
  p_adjusted = s_h.m_h.h
  n_rejections = sum(p_adjusted < alpha)
  if(n_rejections > 0){
    decision = (s <= sort(s)[n_rejections])*1
  }else{
    decision = rep(0,m)
  }
  return(data.frame("decision" = decision, "adjusted.p" = p_adjusted, 
                    selection.p = s, filter.p = f))
}
