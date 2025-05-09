# remove(list = ls())
library(MASS)
library(lars)
library(glmnet)
library("CompQuadForm")
library(gtools)
library(CCA)
library(randomForest)
library(nnet)
library(MCMCprecision)

library(nnet)

CRT_p_val = function(data, rho = 0.6, CRT.replications = 100, 
                     CRT.method_array = c("d0_CRT"), suppress.messages = TRUE, seed = NULL, L = 1, ...){
  start_time = Sys.time()
  set.seed(seed)
  
  check_CRT.method = function(methods){
    return(length(intersect(methods, CRT.method_array)) > 0)
  }
  
  r.squared = function(lin_model, y){ # get r squared from .lm.fit 
    SSE <- sum((lin_model$residuals)^2)
    SST <- sum(y^2)  # sum((y - mean(y))^2)
    return(1- (SSE / SST))
  }
  
  X = data$X
  totals = data$totals
  Y = data$Y
  model.X = data$model.X
  n = length(Y)
  m = ncol(X)
  
  ## two-step
  
  log_with_offset <- function(x) {
    if (is.matrix(x) || is.data.frame(x))
      return(apply(x, 2, log_with_offset))
    
    # Otherwise assume vector:
    if(min(x) == 0){
      offset = 1e-5
    }else{
      offset <- min(x[x > 0], na.rm = TRUE) / 2
    }
    
    return(log(x + offset))
  }
  
  X_j_sampler_two_step <- function(sample_size, X, i, j, model.X) { 
    # models log(X_j + offset) accurately
    X_j <- matrix(NA, nrow = nrow(X), ncol = sample_size)
    
    total_ij = X[, i] + X[, j]
    
    X_prop = X / totals
    
    classes = rep("between", length(total_ij))
    classes[X[,j] == total_ij] = "1"
    classes[X[,j] == 0] = "0"
    classes = factor(classes, levels = c("0", "between", 1))
    
    predictors = log_with_offset(X[, -c(i,j)] / totals)
    
    category_model = multinom(classes ~ predictors, trace = FALSE)
    category_probs = category_model$fitted.values
    
    levs <- unique(classes)
    
    if(length(unique(classes)) == 2) {
      # the category_probs vector is P(Y == levs[2]); P(Y == levs[1]) = 1 – that
      category_probs_matrix = matrix(0, nrow = length(total_ij), ncol = 3)
      colnames(category_probs_matrix) = levels(classes)
      category_probs_matrix[,levs[2]] = as.numeric(category_probs)
      category_probs_matrix[,levs[1]] = 1 - as.numeric(category_probs)
    }else{
      category_probs_matrix = category_probs
      colnames(category_probs_matrix) = levels(classes)
    }
    
    # predicting the values between 0 and 1
    
    if(sum(classes == "between") > 1){
      target_lm = log_with_offset(X[classes == "between",j] / totals[classes == "between"])
      predictors_lm = log_with_offset(X[classes == "between", -c(i,j)] / totals[classes == "between"])
      combined_data = data.frame(cbind(target_lm, predictors_lm))
      colnames(combined_data) = c("target", colnames(predictors_lm))
      
      reg_model = lm(target ~ ., data = combined_data)
      
      reg_prediction_data <- data.frame(log_with_offset(X[, -c(i,j)] / totals))
      colnames(reg_prediction_data) = colnames(predictors_lm)
      
      reg_predictions = predict(reg_model, reg_prediction_data)
      reg_residuals = reg_model$residuals
    }else if(sum(classes == "between") == 1){
      reg_predictions = target_lm = 
        rep(log_with_offset(X[classes == "between",j] / totals[classes == "between"]), nrow(X))
      reg_residuals = 0
    }
    
    classes_sampled = apply(category_probs_matrix, 1, function(p) {
      sample(levels(classes), size = sample_size, prob = p, replace = TRUE)
    })
    
    for (s in 1:sample_size) {
      X_j[classes_sampled[s,] == "1", s] = total_ij[classes_sampled[s,] == "1"]
      X_j[classes_sampled[s,] == "0", s] = 0
      
      if(any(classes_sampled[s,] == "between")){
        is_continuous_sampled = which(classes_sampled[s,] == "between")
        continuous_predictions_sampled = reg_predictions[is_continuous_sampled]
        continuous_res_sampled = sample(reg_residuals, length(is_continuous_sampled), replace = TRUE)
        
        log_X_prop_j_with_offset_sampled = pmin(0, continuous_predictions_sampled + continuous_res_sampled)
        
        offset = exp(min(target_lm)) / 2
        X_prop_j_sampled = pmax(0, exp(log_X_prop_j_with_offset_sampled) - offset)
        
        X_j[is_continuous_sampled, s] = X_prop_j_sampled * total_ij[classes_sampled[s,] == "between"]
      }
    }
    
    return(X_j)
  }
  
  p = rep(list(list("p_bivariate" = diag(m), "p_univariate" = rep(1,m))), each = length(CRT.method_array))
  names(p) <- CRT.method_array
  
  #################### Parallel computing ###############################################################
  lapply(c("pbapply", "parallelly", "foreach", "doParallel"), require, character.only = TRUE)
  n.cores <- parallelly::availableCores()
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "FORK",
    outfile=paste0("parallel_output_", seed,".txt")
  )
  doParallel::registerDoParallel(cl = my.cluster)
  writeLines(as.character(Sys.time()),paste0("parallel_output_", seed,".txt"))
  
  p_bivariate_flat <- foreach(j = 1:(m-1), .combine = rbind, .inorder = FALSE) %:%
    foreach(i = j:m, .combine = rbind, .inorder = FALSE) %dopar%{
      if(i == j){
        value = 1
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)), ", p[", i, ",",j, "] = ", value,"\n")
        c('i' = i, 'j' = j, 'val' = value)
      }else{
        set.seed(seed + i+j*(j-1))
        
        Z = X[, -c(i,j)] / totals
        log_with_offset_Z = log_with_offset(Z)
        
        # multinomial logit
        multi_logit = multinom(as.factor(Y) ~ log_with_offset_Z, trace = FALSE)
        d_y = multi_logit$fitted.values
        
        Y_mat = matrix(0, nrow = length(Y), ncol = length(unique(Y)))
        for(row in 1:length(Y)){
          Y_mat[row, Y[row]] = 1
        }
        Y_res = Y_mat - d_y
        
        X_j_samples = X_j_sampler_two_step(sample_size = CRT.replications, X=X, i=i, j=j, model.X=model.X)
        X_i_samples = rowSums(X[,c(i,j)]) - X_j_samples
        
        T_crt = rep(list(rep(NA, CRT.replications+1)), length(CRT.method_array))
        # notation: last index in each array will be data
        names(T_crt) <- CRT.method_array
        
        for(k in 1:(1+CRT.replications)){
          if(k == (1+CRT.replications)){
            X_ij = X[, c(i,j)]
          }else{
            X_ij = cbind(X_i_samples[,k], X_j_samples[,k])
          }
          
          canonical_r2 = function(X, Y) {
            cca_result = cancor(X, Y)
            return(cca_result$cor[1]^2)  # Squared first canonical correlation
          }
          sum_r2 = function(X, Y) {
            return(sum(sapply(1:ncol(Y), function(j) r.squared(.lm.fit(x = X, y = Y[,j]), Y[,j]))))
          }
          
          if(check_CRT.method("trial")){
            rowsums_X_ij = rowSums(X)
            normalized_X_ij = X_ij / rowsums_X_ij
            normalized_X_ij[rowsums_X_ij == 0, ] = c(0,0)
            
            # T_crt[["trial"]][k] = sum_r2(matrix(log_with_offset(normalized_X_ij[,2])), Y_res)
            T_crt[["trial"]][k] = sum_r2(log_with_offset(normalized_X_ij), Y_res)
          }
        }
        value = (1 + sum(T_crt[["trial"]][-(CRT.replications+1)] >= T_crt[["trial"]][CRT.replications+1]))/(1+CRT.replications)
        cat(paste0("(",i,",",j,", seed =", seed + i+j*(j-1)), ", p[", j, ",",i, "] = ", value,"\n")
        c('i' = i, 'j' = j, 'val' = value)
      }
    }
  # }
  
  
  parallel::stopCluster(my.cluster)
  
  for(l in 1:nrow(p_bivariate_flat)){
    p[["trial"]]$p_bivariate[p_bivariate_flat[l,1], p_bivariate_flat[l,2]] = 
      p[["trial"]]$p_bivariate[p_bivariate_flat[l,2], p_bivariate_flat[l,1]] =
      p_bivariate_flat[l,3]
  }
  
  return(p)
}