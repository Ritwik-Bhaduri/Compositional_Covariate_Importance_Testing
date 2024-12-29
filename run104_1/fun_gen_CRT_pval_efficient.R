# remove(list = ls())
library(MASS)
library(lars)
library(glmnet)
library("CompQuadForm")
library(gtools)
library(CCA)

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
  
  
  if(model == "beta"){
    X = matrix(rbeta(m*n, alpha, alpha*(m-1)),nrow = n)
    model.X = list("model" = model, "alpha" = alpha)
  } 
  
  m_0 = m-m_1
  beta = c(rnorm(m_1, mean = 0, sd = snr), rep(0, m_0))
  noise = rnorm(n, 0, 1)
  Y = log(abs(X)) %*% as.matrix(beta) + noise
  
  return(list(X = X, Y = Y, model.X = model.X))
}

CRT_p_val_efficient = function(data, method = "lasso", rho = 0.6, CRT.replications = 100, 
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
  
  # Preliminary screening using LASSO
  lasso_cv = cv.glmnet(log(abs(X)), Y, alpha = 1, nfolds = 5)
  # lambda_best <- lasso_cv$lambda.min
  # lasso_model <- glmnet(log(abs(X)), Y, alpha = 1, lambda = lambda_best)
  # coefficients <- coef(lasso_model, s = lambda_best)
  lambda_1se <- lasso_cv$lambda.1se
  lasso_model <- glmnet(log(abs(X)), Y, alpha = 1, lambda = lambda_1se)
  coefficients <- coef(lasso_model, s = lambda_1se)
  non_zero_coef <- coefficients[-1, , drop = FALSE] != 0
  screened_indices <- which(non_zero_coef) # These indices are selected by a LASSO of Y on X
  
  cat(paste0("screened indices by LASSO Y|log(X): "))
  # print(coefficients)
  
  cat(screened_indices)
  cat("\n")
  
  p = rep(list(list("p_bivariate" = matrix(1,nrow=m,ncol=m), "p_univariate" = rep(0,m))), each = length(CRT.method_array))
  names(p) <- CRT.method_array
  
  ##################### Parallel computing ###############################################################
  lapply(c("pbapply", "parallelly", "foreach", "doParallel"), require, character.only = TRUE)
  writeLines(as.character(Sys.time()),paste0("result files/parallel_output_", current_seed+base_seed,".txt"))
  
  test_col_length = 10
  randomly_selected_rows = sample(1:m, size = test_col_length, replace = FALSE)
  cat(paste0("\nrandomly selected indices by LASSO Y|log(X): "))
  cat(randomly_selected_rows)
  cat("\n\n")

  p_value_mat_test = matrix(1, nrow  = test_col_length, ncol = m)
  for(j in screened_indices){
    for(i in randomly_selected_rows){
      if(i==j){
        p_value_mat_test[which(randomly_selected_rows == i), j]=1
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)), ", p[", i, ",",j, "] = ", 1,"\n")
      }else{
        set.seed(seed + i+j*(j-1))
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)))
        
        Z = X[, -c(i,j)]
        X_i_plus_j = rowSums(X[, c(i,j)])
        
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
        
        if(model.X$model == "dirichlet"){
          alpha_vec = model.X$alpha_vec
          
          draw.X = function(){
            X_i = rbeta(n,alpha_vec[i], alpha_vec[j])
            X_i[is.na(X_i)]=0
            X_i = X_i * X_i_plus_j
            X_j = X_i_plus_j - X_i
            return(cbind(X_i, X_j))
          }
          
          # draw.X = function() sweep(my_rdirichlet(n, alpha_vec[c(i,j)]), MARGIN=1, rowSums(X[, c(i,j)]), `*`)
          mu_X_ij = sweep(matrix(rep((alpha_vec[c(i,j)]/sum(alpha_vec[c(i,j)])), each = n), nrow = n),
                          MARGIN=1, 1:n, `*`)
        }
        
        if(model.X$model == "beta"){
          draw.X = function() matrix(rbeta(2*n, alpha, alpha*(m-1)),nrow = n)
        }
        
        T_crt = rep(list(rep(NA, CRT.replications+1)), length(CRT.method_array))
        # notation: last index in each array will be data
        names(T_crt) <- CRT.method_array
        
        flag = 0 # flag is 1 if we think P_ij is null after 100 iterations
        for(k in 1:(1+CRT.replications)){
          if(k == 1){
            X_ij = X[, c(i,j)]+1e-300
          }else{X_ij = draw.X()+1e-300}
          
          if(check_CRT.method(c("d_CRT_with_r_squared", "d_CRT_with_max_coeff"))){
            lin_model = .lm.fit(x = cbind(rep(1, nrow(X_ij)), log(abs(X_ij))), y = Y_res)
            if(check_CRT.method("d_CRT_with_r_squared")) T_crt[["d_CRT_with_r_squared"]][k] = r.squared(lin_model, Y_res)
          }
          
          if(k==101){
            value = (1 + sum(T_crt[["d_CRT_with_r_squared"]][2:101] >= T_crt[["d_CRT_with_r_squared"]][1]))/(101)
            if(value > 1/10){
              flag = 1
              break
            }
          }
        }
        
        if(flag == 1){
          value = 1
        }else{
          value = (1 + sum(T_crt[["d_CRT_with_r_squared"]][-1] >= T_crt[["d_CRT_with_r_squared"]][1]))/(1+CRT.replications)
        }
        p_value_mat_test[which(randomly_selected_rows == i), j] = value
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)), ", p[", i, ",",j, "] = ", value,"\n")
      }
    }
  }
  cat(paste0("filtering step done", "\n"))
  
  filtering_function = function(x,r=round(test_col_length/2)){ # we are using simes as the filtering function with half the number of generated p-vals per column as r
    x_ord = sort(x); m=length(x)
    return(min(sapply(1:((m-1)-r+1), function(i) (((m-1)-r+1)/i)*x_ord[r-1+i])))
  }
  filtered_indices = (1:m)[(apply(p_value_mat_test,MARGIN = 2, filtering_function) < 0.1)] 
  
  cat(paste0("test indices: "))
  cat(randomly_selected_rows)
  cat("\n")
  
  
  cat(paste0("Filtered indices: "))
  cat(filtered_indices)
  cat("\n")
  
  # stop("filtering done!!!!")
  
  n.cores <- parallelly::availableCores()
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "FORK",
    outfile=paste0("result files/parallel_output_", current_seed+base_seed,".txt")
  )
  doParallel::registerDoParallel(cl = my.cluster)
  writeLines(as.character(Sys.time()),paste0("result files/parallel_output_", current_seed+base_seed,".txt"))
  
  for(j in filtered_indices){
    for(i in 1:m){
      if(i==j){
        p[["d_CRT_with_r_squared"]]$p_bivariate[i,j] = 1
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)), ", p[", i, ",",j, "] = ", 1,"\n")
      }else{
        set.seed(seed + i+j*(j-1))
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)))
        
        Z = X[, -c(i,j)]
        X_i_plus_j = rowSums(X[, c(i,j)])
        
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
        
        if(model.X$model == "dirichlet"){
          alpha_vec = model.X$alpha_vec
          
          draw.X = function(){
            X_i = rbeta(n,alpha_vec[i], alpha_vec[j])
            X_i[is.na(X_i)]=0
            X_i = X_i * X_i_plus_j
            X_j = X_i_plus_j - X_i
            return(cbind(X_i, X_j))
          }
          
          # draw.X = function() sweep(my_rdirichlet(n, alpha_vec[c(i,j)]), MARGIN=1, rowSums(X[, c(i,j)]), `*`)
          mu_X_ij = sweep(matrix(rep((alpha_vec[c(i,j)]/sum(alpha_vec[c(i,j)])), each = n), nrow = n),
                          MARGIN=1, 1:n, `*`)
        }
        
        if(model.X$model == "beta"){
          draw.X = function() matrix(rbeta(2*n, alpha, alpha*(m-1)),nrow = n)
        }
        
        T_crt = rep(list(rep(NA, CRT.replications+1)), length(CRT.method_array))
        # notation: last index in each array will be data
        names(T_crt) <- CRT.method_array
        
        flag = 0 # flag is 1 if we think P_ij is null after 100 iterations
        for(k in 1:(1+CRT.replications)){
          if(k == 1){
            X_ij = X[, c(i,j)]+1e-300
          }else{X_ij = draw.X()+1e-300}
          
          if(check_CRT.method(c("d_CRT_with_r_squared", "d_CRT_with_max_coeff"))){
            lin_model = .lm.fit(x = cbind(rep(1, nrow(X_ij)), log(abs(X_ij))), y = Y_res)
            if(check_CRT.method("d_CRT_with_max_coeff")) T_crt[["d_CRT_with_max_coeff"]][k] = max(abs(lin_model$coefficients[-1]))
            if(check_CRT.method("d_CRT_with_r_squared")) T_crt[["d_CRT_with_r_squared"]][k] = r.squared(lin_model, Y_res)
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
          
          if(k==101){
            value = (1 + sum(T_crt[["d_CRT_with_r_squared"]][2:101] >= T_crt[["d_CRT_with_r_squared"]][1]))/(101)
            if(value > 1/10){
              flag = 1
              break
            }
          }
        }
        
        if(flag == 1){
          value = 1
        }else{
          value = (1 + sum(T_crt[["d_CRT_with_r_squared"]][-1] >= T_crt[["d_CRT_with_r_squared"]][1]))/(1+CRT.replications)
        }
        p[["d_CRT_with_r_squared"]]$p_bivariate[i,j] = value
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)), ", p[", i, ",",j, "] = ", value,"\n")
      }
    }
  }
  
  if(model.X$model != "dirichlet"){
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
  
  if(model.X$model == "dirichlet"){
    j = sample(1:m, size = 1)
    for(CRT.method in CRT.method_array){
      p[[CRT.method]]$p_univariate = as.vector(p[[CRT.method]]$p_bivariate[,j])
      p[[CRT.method]]$p_univariate[j] = 1
    }
  }
  
  return(p)
  
}