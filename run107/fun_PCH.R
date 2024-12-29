# for conditioning method

PCH_pval <- function(p, alpha = 0.1, m_1, method = NULL, sparsity_cutoff_m_value, j = NULL){  
  # m_1 is required for oracle methods
  
  if(length(p) == 2){
    p_univariate = p$p_univariate[1:sparsity_cutoff_m_value]
  }
  p = p$p_bivariate
  m_original = nrow(p)
  p = p[1:sparsity_cutoff_m_value, 1:sparsity_cutoff_m_value]
  m = nrow(p)
  
  Simes <- function(x,r){
    x_ord = sort(x)
    m=length(x)
    return(min(sapply(1:(m-r+1), function(i) ((m-r+1)/i)*x_ord[r-1+i])))
  }
    
  
  if(grepl("Simes_oracle_",method, fixed = TRUE)){
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    if(r == m_original) r = m-1
    p_j = p[-j,j] 
    pval = Simes(p_j,r)
  }
  
  if(grepl("bonferroni_oracle_",method, fixed = TRUE)){
    r = as.numeric(tail(strsplit(method, split='_', fixed=TRUE)[[1]],1))
    if(r == m_original) r = m-1
    p_j = p[-j,j]
    m = length(p_j)
    pval = min((m-r+1) * sort(p_j)[r],1)
  }
  
  if(method == "LOO"){
    LOO_var = sample(1:m, 1)
    pval = p[LOO_var,j] 
  }
  
  return(pval)
}
