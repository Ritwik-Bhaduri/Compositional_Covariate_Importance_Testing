# remove(list = ls())
library(MASS)
library(lars)
library(glmnet)
library("CompQuadForm")
library(gtools)
library(CCA)
library(MCMCprecision)

gen_data = function(n = 100, m = 10, model = "softmax_log_normal", var.model = "autocorrelation", 
                    rho = 0.6, m_1, alpha_vec = NULL, snr = 10, ...){
  if(model %in% c("softmax_log_normal")){
    mu = rep(0,m)
    mu_matrix <- matrix(rep(mu, n), nrow = n, byrow = TRUE)
    
    if(var.model == "independent") Sigma = diag(m) # only works for n>= m
    if(var.model == "equicorrelation") Sigma = matrix(rho, nrow = m, ncol = m) + (1-rho) * diag(m)
    if(var.model == "autocorrelation") Sigma =  rho^abs(row(diag(m))-col(diag(m)))
    
    X_normal = MASS::mvrnorm(n, mu = mu, Sigma = Sigma)
    X_sums = rowSums(exp(X_normal))
    X = t(sapply(1:n, function(i) exp(X_normal[i,])/X_sums[i]))
    Sigma_inv = solve(Sigma)
    det_Sigma_inv = det(Sigma_inv)
    
    model.X = list("model" = model, "var.model" = var.model, "mu" = mu, "Sigma" = Sigma, "Sigma_inv" = Sigma_inv, "det_Sigma_inv" = det_Sigma_inv)
  }
  
  m_0 = m-m_1
  beta = c(snr, rnorm(m_1-1, mean = 0, sd = 1), rep(0, m_0))
  noise = rnorm(n, 0, 1)
  Y = log(abs(X)) %*% as.matrix(beta) + noise
  plot(Y, log(abs(X)) %*% as.matrix(beta))
  
  return(list(X = X, Y = Y, model.X = model.X))
}

CRT_p_val = function(data, method = "lasso", rho = 0.6, CRT.replications = 100, 
                     CRT.method_array = c("d0_CRT"), suppress.messages = TRUE, seed = NULL, L = 5, null_j=NULL, ...){
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
  
  X = data$X
  Y = data$Y
  model.X = data$model.X
  n = length(Y)
  m = ncol(X)
  mu = model.X$mu
  Sigma = model.X$Sigma
  Sigma_inv = model.X$Sigma_inv
  det_Sigma_inv = model.X$det_Sigma_inv
  
  # precomputations for log density
  one_m = matrix(rep(1, m))
  E_m = diag(m) # I am not sure what E_m is. Have to ask Siyuan
  F_0 = E_m - 1/m * one_m %*% t(one_m)
  omega = Sigma_inv
  v = (omega - sum(rowSums(omega) * colSums(omega)) / sum(omega))
  log_det_Sigma_inv = log(det_Sigma_inv)
  log_sum_omega = log(sum(omega))
  
  log_density = function(x, mu, Sigma_inv = Sigma_inv, det_Sigma_inv = det_Sigma_inv, 
                         F_01, u_2, u2_t_v_u2, v_u2, F_01_t_v_F_01, F_01_t_v_u2, index_1, index_2) {
    log_x = log(x)
    u_1 = F_01 %*% matrix((log_x-mu)[c(index_1, index_2)])
    u = u_1 + u_2
    Q_1 = t(log_x[c(index_1, index_2)]) %*% F_01_t_v_F_01 %*% log_x[c(index_1, index_2)] + 2 * t(log_x[c(index_1, index_2)]) %*% F_01_t_v_u2 + u2_t_v_u2
    log_density = -(m-1)/2 * log(2*pi) + 1/2 * log_det_Sigma_inv - 1/2 * log_sum_omega - sum(log_x) - 1/2 * Q_1
    return(log_density)
  }
  
  MH_independent = function(x, index_1, index_2, n_samples, proposal_alpha = c(3,3), reversed = FALSE, 
                            F_01, u_2, u2_t_v_u2, v_u2, F_01_t_v_F_01, F_01_t_v_u2){
    y = rep(NULL, n_samples)
    a = 0; b = x[index_1] + x[index_2]
    for(i in 1:n_samples){
      old_value = ifelse(i == 1, x[index_1], y[i-1])
      temp = a + rbeta(1, proposal_alpha[1], proposal_alpha[2])*(b-a)
      vec_old = x; vec_old[index_1] = old_value; vec_old[index_2] = b - old_value
      vec_new = x; vec_new[index_1] = temp;   vec_new[index_2] = b - temp
      if(reversed == FALSE){
        acceptance_prob = min(1,exp(log_density(vec_new, mu, Sigma_inv, det_Sigma_inv, F_01 = F_01, u_2 = u_2, u2_t_v_u2 = u2_t_v_u2, v_u2 = v_u2, F_01_t_v_F_01 = F_01_t_v_F_01, F_01_t_v_u2 = F_01_t_v_u2, index_1 = index_1, index_2 = index_2) - 
                                      log_density(vec_old, mu, Sigma_inv, det_Sigma_inv, F_01 = F_01, u_2 = u_2, u2_t_v_u2 = u2_t_v_u2, v_u2 = v_u2, F_01_t_v_F_01 = F_01_t_v_F_01, F_01_t_v_u2 = F_01_t_v_u2, index_1 = index_1, index_2 = index_2)) * 
                                dbeta((old_value-a)/(b-a), proposal_alpha[1], proposal_alpha[2]) / 
                                dbeta((temp-a)/(b-a), proposal_alpha[1], proposal_alpha[2])) # independent random_walk
      }else if(reversed == TRUE){
        acceptance_prob = min(1,exp(log_density(vec_old, mu, Sigma_inv, det_Sigma_inv, F_01 = F_01, u_2 = u_2, u2_t_v_u2 = u2_t_v_u2, v_u2 = v_u2, F_01_t_v_F_01 = F_01_t_v_F_01, F_01_t_v_u2 = F_01_t_v_u2, index_1 = index_1, index_2 = index_2) - 
                                      log_density(vec_new, mu, Sigma_inv, det_Sigma_inv, F_01 = F_01, u_2 = u_2, u2_t_v_u2 = u2_t_v_u2, v_u2 = v_u2, F_01_t_v_F_01 = F_01_t_v_F_01, F_01_t_v_u2 = F_01_t_v_u2, index_1 = index_1, index_2 = index_2)) * 
                                dbeta((temp-a)/(b-a), proposal_alpha[1], proposal_alpha[2]) / 
                                dbeta((old_value-a)/(b-a), proposal_alpha[1], proposal_alpha[2])) # inverted independent random_walk
      }
      y[i] = ifelse(runif(1) < acceptance_prob, temp, old_value)
    }
    return(y)
  }
  
  permuted_serial_sampler_one_index = function(x, index_1, index_2, L, M, proposal_alpha = c(3,3)){ 
    #L is the thinning rate of the Markov chain. So we accept only the Lth entries 
    #M is the number of CRT resamples that we want to draw
    pi_permutation = sample(0:M, replace = FALSE)
    output = rep(NULL, M+1)
    m_star = which(pi_permutation == 0)
    useful_samples_1 = c(); useful_samples_2 = c()
    
    # precomputations for density
    F_01 = F_0[,c(index_1, index_2)]; F_02 = F_0[, -c(index_1, index_2)]
    u_2 = F_02 %*% matrix((log(x)-mu)[-c(index_1, index_2)])
    u2_t_v_u2 = t(u_2) %*% v %*% u_2  
    v_u2 = v %*% u_2
    F_01_t_v_F_01 = t(F_01) %*% v %*% F_01
    F_01_t_v_u2 = t(F_01) %*% v_u2
    
    if(m_star > 1){
      chain_1 = MH_independent(x, index_1 = index_1, index_2 = index_2, n_samples = L*(m_star-1), proposal_alpha = proposal_alpha, reversed = FALSE,
                               F_01 = F_01, u_2 = u_2, u2_t_v_u2 = u2_t_v_u2, v_u2 = v_u2, F_01_t_v_F_01 = F_01_t_v_F_01, F_01_t_v_u2 = F_01_t_v_u2)
      output[pi_permutation[(m_star-1):1]] = chain_1[(1 : (m_star-1)) * L]
    }
    if(m_star < M + 1){
      chain_2 = MH_independent(x, index_1 = index_1, index_2 = index_2, n_samples = L*(M-m_star+1), proposal_alpha = proposal_alpha, reversed = TRUE,
                               F_01 = F_01, u_2 = u_2, u2_t_v_u2 = u2_t_v_u2, v_u2 = v_u2, F_01_t_v_F_01 = F_01_t_v_F_01, F_01_t_v_u2 = F_01_t_v_u2)
      output[pi_permutation[(m_star+1):(M+1)]] = chain_2[(1 : (M-m_star+1)) * L]
    }
    output[0] = x[index_1]
    return(output)
  }
  
  permuted_serial_sampler = function(X, index_1, index_2, L, M, proposal_alpha = c(3,3)){ 
    #L is the thinning rate of the Markov chain. So we accept only the Lth entries 
    #M is the number of CRT resamples that we want to draw
    output_matrix = matrix(0, nrow = nrow(X), ncol = M)
    for(i in 1:nrow(X)){
      output_matrix[i, ] = permuted_serial_sampler_one_index(x=X[i,], index_1 = index_1, index_2 = index_2, L = L, M = M, proposal_alpha = proposal_alpha)
    }
    return(output_matrix)  
  }
  
  p = rep(list(list("p_bivariate" = diag(m)*0, "p_univariate" = rep(0,m))), each = length(CRT.method_array))
  names(p) <- CRT.method_array
  
  for(j in c(1,null_j)){
    for(i in 1:m){
      set.seed(seed + i+j*(j-1))
      cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)))
      
      Z = X[, -c(i,j)]
      
      lasso_cv = cv.glmnet(log(abs(Z)), Y, alpha = 1, nfolds = 5)
      lambda_best <- lasso_cv$lambda.min
      d_y = predict(lasso_cv, s = lasso_cv$lambda.min , newx = log(abs(Z)))
      
      Y_res = Y-d_y
      
      X_i_samples = permuted_serial_sampler(X, index_1 = i, index_2 = j, L = L, M = CRT.replications, proposal_alpha = c(3,3))
      X_j_samples = rowSums(X[,c(i,j)]) - X_i_samples
      T_crt = rep(list(rep(NA, CRT.replications+1)), length(CRT.method_array))
      # notation: last index in each array will be data
      names(T_crt) <- CRT.method_array
      
      for(k in 1:(1+CRT.replications)){
        if(k == (1+CRT.replications)){
          X_ij = X[, c(i,j)]
        }else{
          X_ij = cbind(X_i_samples[,k], X_j_samples[,k])
        }
        
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
        
      }
      
      for(CRT.method in CRT.method_array){
        p[[CRT.method]]$p_bivariate[i,j] = (1 + sum(T_crt[[CRT.method]][-(CRT.replications+1)] >= T_crt[[CRT.method]][CRT.replications+1]))/(1+CRT.replications)
      }
      cat(paste0(", p[", i, ",",j, "] = ", p[[CRT.method]]$p_bivariate[i,j],"\n"))
    }
  }
  
  if(model.X$model %in% c("softmax_log_normal")){
    j = sample(1:m, size = 1)
    for(CRT.method in CRT.method_array){
      p[[CRT.method]]$p_univariate = as.vector(p[[CRT.method]]$p_bivariate[,j])
      p[[CRT.method]]$p_univariate[j] = 1
    }
  }
  
  return(p)
  
}