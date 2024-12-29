# remove(list = ls())
library(MASS)
library(lars)
library(glmnet)
library("CompQuadForm")
library(gtools)
library(CCA)
library(MCMCprecision)

gen_data = function(n = 100, m = 10, model = "normal", var.model = "autocorrelation", 
                    rho = 0.6, m_1, alpha_vec = NULL, snr = 10, ...){
  if(model %in% c("normal", "normal_spiked")){
    if(model == "normal") 
      mu = matrix(0, nrow = n, ncol = m)
    if(model == "normal_spiked")
      mu = rbind(diag(m)*exp(5), matrix(0, nrow = n-m, ncol = m))
    
    if(var.model == "independent") Sigma = diag(m) # only works for n>= m
    if(var.model == "equicorrelation") Sigma = matrix(rho, nrow = m, ncol = m) + (1-rho) * diag(m)
    if(var.model == "autocorrelation") Sigma =  rho^abs(row(diag(m))-col(diag(m)))
    
    X_centered = mvrnorm(n, mu = rep(0, m), Sigma = Sigma)+1e-300
    X = X_centered + mu
    model.X = list("model" = model, "var.model" = var.model, "mu" = mu, "Sigma" = Sigma, "rho" = rho)
  }
  
  if(model == "dirichlet"){
    X = rdirichlet(n, alpha_vec)+1e-300
    model.X = list("model" = model, "alpha_vec" = alpha_vec)
  }
  
  if(model == "dirichlet_estimated"){
    X = rdirichlet(n, alpha_vec)+1e-300
    estimated_alpha = fit_dirichlet(X, maxit = 1e+06, abstol = min(alpha_vec) * 0.005)$alpha
    model.X = list("model" = model, "alpha_vec" = estimated_alpha)
  }
  
  if(model == "beta"){
    X = matrix(rbeta(m*n, alpha, alpha*(m-1)),nrow = n)
    model.X = list("model" = model, "alpha" = alpha)
  } 
  
  m_0 = m-m_1
  beta = c(rnorm(m_1, mean = 0, sd = snr), rep(0, m_0))
  noise = rnorm(n, 0, 1)
  Y = log(abs(X)) %*% as.matrix(beta) + noise
  plot(Y, log(abs(X)) %*% as.matrix(beta))
  
  return(list(X = X, Y = Y, model.X = model.X))
}

CRT_p_val = function(data, method = "lasso", rho = 0.6, CRT.replications = 100, 
                     CRT.method_array = c("d0_CRT"), suppress.messages = TRUE, seed = NULL,...){
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
  
  if(model.X$model %in% c("normal", "normal_spiked")){
    mu = model.X$mu
    Sigma = model.X$Sigma
  }
  
  if(model.X$model == "beta") alpha = model.X$alpha
  
  p = rep(list(list("p_bivariate" = diag(m)*0, "p_univariate" = rep(0,m))), each = length(CRT.method_array))
  names(p) <- CRT.method_array
  
  for(j in 2:m){
    for(i in 1:(j-1)){
      set.seed(seed + i+j*(j-1))
      cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)))
      
      Z = X[, -c(i,j)]
      
      lasso_cv = cv.glmnet(log(abs(Z)), Y, alpha = 1, nfolds = 5)
      lambda_best <- lasso_cv$lambda.min
      d_y = predict(lasso_cv, s = lasso_cv$lambda.min , newx = log(abs(Z)))
      
      Y_res = Y-d_y
      
      if(model.X$model %in% c("normal", "normal_spiked")){
        Sigma_11 = Sigma[c(i,j), c(i,j)]
        Sigma_12 = Sigma[c(i,j), -c(i,j)]
        Sigma_21 = Sigma[-c(i,j), c(i,j)]
        Sigma_22 = Sigma[-c(i,j), -c(i,j)]
        Sigma_22_inv = solve(Sigma_22)
        
        mu_X_ij = mu[, c(i,j)] + t(Sigma_12 %*% Sigma_22_inv %*% t(Z - mu[, -c(i,j)]))
        Sigma_X_ij = Sigma_11 - Sigma_12 %*% Sigma_22_inv %*% Sigma_21
        
        # draw.X = function() t(apply(mu_X_ij, 1, function(mu) mvrnorm(n = 1, mu = mu, Sigma = Sigma_X_ij)))
        draw.X = function() mu_X_ij + mvrnorm(n = n, mu = c(0,0), Sigma = Sigma_X_ij) ### correct this
        
      }
      
      if(model.X$model %in% c("dirichlet", "dirichlet_estimated")){
        alpha_vec = model.X$alpha_vec
        draw.X = function() sweep(my_rdirichlet(n, alpha_vec[c(i,j)]), MARGIN=1, rowSums(X[, c(i,j)]), `*`)
        mu_X_ij = sweep(matrix(rep((alpha_vec[c(i,j)]/sum(alpha_vec[c(i,j)])), each = n), nrow = n),
                        MARGIN=1, 1:n, `*`)
      }
      
      if(model.X$model == "beta"){
        draw.X = function() matrix(rbeta(2*n, alpha, alpha*(m-1)),nrow = n)
      }
      
      T_crt = rep(list(rep(NA, CRT.replications+1)), length(CRT.method_array))
      # notation: last index in each array will be data
      names(T_crt) <- CRT.method_array
      
      for(k in 1:(1+CRT.replications)){
        if(k == (1+CRT.replications)){
          X_ij = X[, c(i,j)]+1e-300
        }else{X_ij = draw.X()+1e-300}
        # if(k==5) print(c(i,j,k))
        # if(all(c(i,j,k)==c(2,5,5))) print(X_ij[32,])
        # if(i==2){if(k==5) print(rowSums(X_ij))}
        
        if(check_CRT.method(c("d_CRT_with_r_squared", "d_CRT_with_max_coeff"))){
          # print(log(abs(X_ij+1)))
          lin_model = .lm.fit(x = cbind(rep(1, nrow(X_ij)), log(abs(X_ij))), y = Y_res)
          if(check_CRT.method("d_CRT_with_max_coeff")) T_crt[["d_CRT_with_max_coeff"]][k] = max(abs(lin_model$coefficients[-1]))
          if(check_CRT.method("d_CRT_with_r_squared")) T_crt[["d_CRT_with_r_squared"]][k] = r.squared(lin_model, Y_res)
          # if(check_CRT.method("d_CRT_with_r_squared")) T_crt[["d_CRT_with_r_squared"]][k] =
          # summary(lm(Y_res ~ log(abs(X_ij)) %*% lin_model$coefficients))$r.squared
        }
        
        if(check_CRT.method("d_CRT_with_canonical_correlation")){
          cc_results = cancor(log(abs(X_ij)),(Y-d_y))
          CC1_X <- log(abs(as.matrix(X_ij))) %*% cc_results$xcoef[, 1]
          CC1_Y <- as.matrix((Y-d_y)) %*% cc_results$ycoef[, 1]
          CC1_X_centered = CC1_X- mean(CC1_X)
          CC1_Y_centered = CC1_Y- mean(CC1_Y)
          T_crt[["d_CRT_with_canonical_correlation"]][k] = abs(sum(CC1_Y_centered * CC1_X_centered) / sum(CC1_X_centered^2))
        }
        
        if(check_CRT.method("trial")){
          lin_model = .lm.fit(x = cbind(rep(1, nrow(X_ij)), log(abs(X_ij)), rowSums(log(abs(Z)))), y = Y_res)
          T_crt[["trial"]][k] = r.squared(lin_model, Y_res)
        }
        
      }
      
      for(CRT.method in CRT.method_array){
        # if(CRT.method == "d_CRT_with_r_squared") print((1 + sum(T_crt[[CRT.method]][-(CRT.replications+1)] >= T_crt[[CRT.method]][CRT.replications+1]))/(1+CRT.replications))
        p[[CRT.method]]$p_bivariate[i,j] = p[[CRT.method]]$p_bivariate[j,i] =
          (1 + sum(T_crt[[CRT.method]][-(CRT.replications+1)] >= T_crt[[CRT.method]][CRT.replications+1]))/(1+CRT.replications)
      }
      cat(", p[", j, ",",i, "] = ", p[[CRT.method]]$p_bivariate[i,j],"\n")

    }
  }
  
  if(!(model.X$model %in% c("dirichlet","dirichlet_estimated"))){
    for(i in 1:m){
      if(suppress.messages == FALSE) print(paste0("Univariate test: ", i, " out of ", m))
      Z = X[, -i]
      lasso_cv = cv.glmnet(log(abs(Z)), Y, alpha = 1, nfolds = 5)
      lambda_best <- lasso_cv$lambda.min
      d_y = predict(lasso_cv, s = lasso_cv$lambda.min , newx = log(abs(Z)))
      Y_res = Y-d_y
      
      if(model.X$model %in% c("normal", "normal_spiked")){
        Sigma_11 = Sigma[i, i]
        Sigma_12 = Sigma[i, -i]
        Sigma_21 = Sigma[-i, i]
        Sigma_22 = Sigma[-i, -i]
        Sigma_22_inv = solve(Sigma_22)
        
        mu_X_i = mu[,i] + t(matrix(Sigma_12)) %*% Sigma_22_inv %*% t(Z - mu[,-i])
        Sigma_X_i = Sigma_11 - t(matrix(Sigma_12)) %*% Sigma_22_inv %*% matrix(Sigma_21)
        
        draw.X = function() t(mu_X_i+ rnorm(n = n, mean = 0, sd = sqrt(Sigma_X_i)))
      }
      
      if(model.X$model == "beta"){
        draw.X = function() matrix(rbeta(n,alpha, alpha*(m-1)), nrow=n)
      }
      
      
      T_crt = rep(list(rep(NA, CRT.replications+1)), length(CRT.method_array))
      # notation: last index in each array will be data
      names(T_crt) <- CRT.method_array
      
      for(k in 1:(1+CRT.replications)){
        if(k == (1+CRT.replications)){
          X_i = matrix(X[, i])
        }else{X_i = draw.X()+1e-300}
        
        if(check_CRT.method(c("d_CRT_with_r_squared", "d_CRT", "d_CRT_with_max_coeff"))){
          lin_model = .lm.fit(x = cbind(rep(1, nrow(X_i)), log(abs(X_i))), y = Y_res)
          if(check_CRT.method("d_CRT")) T_crt[["d_CRT"]][k] = abs(lin_model$coefficients)[2]
          if(check_CRT.method("d_CRT_with_max_coeff")) T_crt[["d_CRT_with_max_coeff"]][k] = max(abs(lin_model$coefficients[-1]))
          if(check_CRT.method("d_CRT_with_r_squared")) T_crt[["d_CRT_with_r_squared"]][k] = r.squared(lin_model, Y_res)
          
        }
        
        if(check_CRT.method("d_CRT_with_canonical_correlation")){
          cc_results = cancor(log(abs(X_i)),(Y-d_y))
          CC1_X <- log(abs(as.matrix(X_i))) %*% cc_results$xcoef[, 1]
          CC1_Y <- as.matrix((Y-d_y)) %*% cc_results$ycoef[, 1]
          CC1_X_centered = CC1_X- mean(CC1_X)
          CC1_Y_centered = CC1_Y- mean(CC1_Y)
          T_crt[["d_CRT_with_canonical_correlation"]][k] = abs(sum(CC1_Y_centered * CC1_X_centered) / sum(CC1_X_centered^2))
        }
        
        if(check_CRT.method("trial")){
          T_crt[["trial"]][k] = 0
        }
      }
      
      for(CRT.method in CRT.method_array){
        p[[CRT.method]]$p_univariate[i] = 
          (1 + sum(T_crt[[CRT.method]][-(CRT.replications+1)] >= T_crt[[CRT.method]][CRT.replications+1]))/(1+CRT.replications)
      }
      
    }
  }
  
  if(model.X$model %in% c("dirichlet", "dirichlet_estimated")){
    j = sample(1:m, size = 1)
    for(CRT.method in CRT.method_array){
      p[[CRT.method]]$p_univariate = as.vector(p[[CRT.method]]$p_bivariate[,j])
      p[[CRT.method]]$p_univariate[j] = 1
    }
  }
  
  return(p)
  
}