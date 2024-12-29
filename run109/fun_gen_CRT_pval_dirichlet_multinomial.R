# remove(list = ls())
library(MASS)
library(lars)
library(glmnet)
library("CompQuadForm")
library(gtools)
library(CCA)
library(MCMCprecision)

gen_data = function(n = 100, m = 10, model = "dirichlet_multinomial", var.model = "autocorrelation", 
                    rho = 0.6, m_1, alpha_vec = NULL, snr = 10, non_nulls, ...){
  if(model %in% c("dirichlet_multinomial")){
    N_multinomial = 200
    theta = rdirichlet(n, alpha_vec)
    X = t(sapply(1:n, function(i) rmultinom(1,N_multinomial,prob = theta[i,])))
    model.X = list("model" = model, "alpha_vec" = alpha_vec, "N_multinomial" = N_multinomial)
  }
  if(model == "dirichlet"){
    X = rdirichlet(n, alpha_vec)+1e-300
    model.X = list("model" = model, "alpha_vec" = alpha_vec)
  }
  
  m_0 = m-m_1
  beta = rep(0,m)
  beta[non_nulls] = rnorm(m_1, mean = 0, sd = snr)
  noise = rnorm(n, 0, 1)
  Y = log(1+abs(X)) %*% as.matrix(beta) + noise
  return(list(X = X, Y = Y, model.X = model.X))
}

CRT_p_val = function(data, method = "lasso", rho = 0.6, CRT.replications = 100, 
                     CRT.method_array = c("d0_CRT"), suppress.messages = TRUE, seed = NULL, L = 1, ...){
  start_time = Sys.time()
  set.seed(seed)
  
  check_CRT.method = function(methods){
    return(length(intersect(methods, CRT.method_array)) > 0)
  }
  
  r.squared = function(lin_model, Y){ # get r squared from .lm.fit 
    SSE <- sum((lin_model$residuals)^2)
    SST <- sum((Y - mean(Y))^2)
    return(1- (SSE / SST))
  }
  my_rdirichlet <- function(n, x){
    ans = rdirichlet(n,x)
    ans[is.na(ans)]=0
    return(ans)
  }
  
  X = data$X
  Y = data$Y
  model.X = data$model.X
  n = length(Y)
  m = ncol(X)
  alpha_vec = model.X$alpha_vec
  if(model.X$model == "dirichlet_multinomial") N_multinomial = model.X$N_multinomial
  
  # computation for dirichlet multinomial
  pmf_DM = function(x, i=1,j=2,alpha_vec){
    N = sum(x)
    log_prob = lfactorial(N) - sum(lfactorial(x)) + lgamma(sum(alpha_vec)) - sum(lgamma(alpha_vec)) + sum(lgamma(alpha_vec+x)) - lgamma(sum(alpha_vec + x))
    return(exp(log_prob))
  }
  prob_vector_DM = function(x, i, j, alpha_vec){
    candidates = 0:(x[i]+x[j])
    prob_vector = candidates*0
    for(k in 1:length(candidates)){
      candidate = candidates[k]
      y = x
      y[i] = candidate; y[j] = x[i]+x[j]-candidate
      prob_vector[k] = pmf_DM(y, i=i, j=j, alpha_vec = alpha_vec)
    }
    prob_vector = exp(prob_vector)/sum(exp(prob_vector))
    return(prob_vector)
  }
  X_i_sampler = function(sample_size, X,i,j,alpha_vec){
    X_i = matrix(0, nrow = nrow(X), ncol = sample_size)
    for(row in 1:nrow(X)){
      prob_vector = prob_vector_DM(X[row,], i=i, j=j, alpha_vec = alpha_vec)
      if(is.nan(sum(prob_vector))) {print(row); print(X[row,]); saveRDS(X[row, ],file="temp.rds")}
      X_i[row,] = sample(0:(X[row,i]+X[row,j]), size = sample_size, replace = TRUE, prob = prob_vector)
    }
    return(X_i)
  }
  
  p = rep(list(list("p_bivariate" = diag(m)*0, "p_univariate" = rep(0,m))), each = length(CRT.method_array))
  names(p) <- CRT.method_array
  
  for(j in 2:m){
    for(i in 1:(j-1)){
      set.seed(seed + i+j*(j-1))
      
      Z = X[, -c(i,j)]
      
      lasso_cv = cv.glmnet(log(1+abs(Z)), Y, alpha = 1, nfolds = 5)
      lambda_best <- lasso_cv$lambda.min
      d_y = predict(lasso_cv, s = lasso_cv$lambda.min , newx = log(1+abs(Z)))
      
      Y_res = Y-d_y
      
      if(model.X$model == "dirichlet_multinomial"){
        X_i_samples = X_i_sampler(sample_size = CRT.replications, X=X, i=i, j=j, alpha_vec = alpha_vec)
        X_j_samples = rowSums(X[,c(i,j)]) - X_i_samples
      }else if(model.X$model == "dirichlet"){
        alpha_vec = model.X$alpha_vec
        draw.X = function() sweep(my_rdirichlet(n, alpha_vec[c(i,j)]), MARGIN=1, rowSums(X[, c(i,j)]), `*`)
        mu_X_ij = sweep(matrix(rep((alpha_vec[c(i,j)]/sum(alpha_vec[c(i,j)])), each = n), nrow = n),
                        MARGIN=1, 1:n, `*`)
      }
      
      # mm = as.vector(X_i_samples)
      # par(mfrow=c(2,1))
      # hist(mm/rowSums(X[,c(i,j)]), main = "Sampling from density", breaks = seq(0,1,length.out=30), freq = FALSE)
      # hist(X[,i]/rowSums(X[,c(i,j)]), main = "true data",breaks = seq(0,1,length.out=30), freq = FALSE)
      
      
      T_crt = rep(list(rep(NA, CRT.replications+1)), length(CRT.method_array))
      # notation: last index in each array will be data
      names(T_crt) <- CRT.method_array
      
      for(k in 1:(1+CRT.replications)){
        if(k == (1+CRT.replications)){
          X_ij = X[, c(i,j)]
        }else{
          if(model.X$model == "dirichlet_multinomial") X_ij = cbind(X_i_samples[,k], X_j_samples[,k])
          if(model.X$model == "dirichlet") X_ij = draw.X()+1e-300
        }
        
        if(check_CRT.method(c("d_CRT_with_r_squared", "d_CRT_with_max_coeff"))){
          lin_model = .lm.fit(x = cbind(rep(1, nrow(X_ij)), log(1+abs(X_ij))), y = Y_res)
          if(check_CRT.method("d_CRT_with_max_coeff")) T_crt[["d_CRT_with_max_coeff"]][k] = max(abs(lin_model$coefficients[-1]))
          if(check_CRT.method("d_CRT_with_r_squared")) T_crt[["d_CRT_with_r_squared"]][k] = r.squared(lin_model, Y_res)
        }
        
        if(check_CRT.method("d_CRT_with_canonical_correlation")){
          cc_results = cancor(log(1+abs(X_ij)),(Y-d_y))
          CC1_X <- log(1+abs(as.matrix(X_ij))) %*% cc_results$xcoef[, 1]
          CC1_Y <- as.matrix((Y-d_y)) %*% cc_results$ycoef[, 1]
          CC1_X_centered = CC1_X- mean(CC1_X)
          CC1_Y_centered = CC1_Y- mean(CC1_Y)
          T_crt[["d_CRT_with_canonical_correlation"]][k] = abs(sum(CC1_Y_centered * CC1_X_centered) / sum(CC1_X_centered^2))
        }
        
        if(check_CRT.method("trial")){
          lin_model = .lm.fit(x = cbind(rep(1, nrow(X_ij)), log(1+abs(X_ij)), rowSums(log(1+abs(Z)))), y = Y_res)
          T_crt[["trial"]][k] = r.squared(lin_model, Y_res)
        }
      }
      
      for(CRT.method in CRT.method_array){
        # if(CRT.method == "d_CRT_with_r_squared") print((1 + sum(T_crt[[CRT.method]][-(CRT.replications+1)] >= T_crt[[CRT.method]][CRT.replications+1]))/(1+CRT.replications))
        p[[CRT.method]]$p_bivariate[i,j] = p[[CRT.method]]$p_bivariate[j,i] =
          (1 + sum(T_crt[[CRT.method]][-(CRT.replications+1)] >= T_crt[[CRT.method]][CRT.replications+1]))/(1+CRT.replications)
      }
      cat("(",i,",",j,"), seed =", seed + i+j*(j-1), ", p[", j, ",",i, "] = ", p[[CRT.method]]$p_bivariate[i,j],"\n")
    }      
  }
  
  if(model.X$model %in% c("dirichlet_multinomial")){
    j = sample(1:m, size = 1)
    for(CRT.method in CRT.method_array){
      p[[CRT.method]]$p_univariate = as.vector(p[[CRT.method]]$p_bivariate[,j])
      p[[CRT.method]]$p_univariate[j] = 1
    }
  }
  
  return(p)
  
}