remove(list = ls())
versionid = 107
setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid))

source("fun_gen_CRT_pval_dirichlet_multinomial.R")
source("fun_PCH.R")

taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# taskid = 50
set.seed(taskid)
cat("taskid = ", taskid)

##### Fixed Parameters ############################################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data

# sparsity_cutoff_m_values = c(13,14)
# oracle_r = c(m_1+1,m_1+2)

sparsity_cutoff_m_values = c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
oracle_r = c(m_1+1,m_1+2,10, 20, 25, m)

## special parameters ##############################################################################
rev_logit = function(x,a=1,b=0,mu=m/2,slope=5){ 
  # starts at a, ends at b, the center is at mu, less the slope variable steeper the slope
  return((a-b)/(1+exp((x-mu)/slope))+b)
}
x_array = seq(1,m,length.out = m)
y_array = sapply(x_array, rev_logit)
# plot(x_array, y_array)

dirichlet_parameters = sapply(x_array, function(x) rev_logit(x=x, a=1,b=0))
non_null_probs = sapply(x_array, function(x) rev_logit(x=x, a=1,b=1))

####################################################################################################

alpha = 0.05 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
simulation.iter = 10 # number of replications of the simulations
target_simulation.iter = 500
CRT.method_array = c("d_CRT_with_r_squared")
model_values = c("dirichlet_multinomial")

####### Methods ###############################################################################
pch_methods = sapply(oracle_r, function(x) paste0("Simes_oracle_",x))
pch_methods = c(pch_methods, sapply(oracle_r, function(x) paste0("bonferroni_oracle_",x)))
pch_methods = c(pch_methods, "LOO")

###############################################################################################

base_seed = taskid*target_simulation.iter + versionid*simulation.iter*1000
print(paste0("base_seed = ", base_seed))

##### Variable Parameters arrays ##################################################### 
model_array = rep("dirichlet_multinomial", ceiling(target_simulation.iter/simulation.iter))
seed_array = (1:ceiling(target_simulation.iter/simulation.iter)) * simulation.iter
var_params = data.frame("model" = model_array, "seed" = seed_array)

##### Variable Parameters for this run ############################################### 
snr = 0.85
model = model_array[taskid]
current_seed = var_params[taskid, "seed"]
print(paste0("snr = ", snr, ", model = ",model, ", seed = ", current_seed))

######################################################################################  
simulate_fdr_fwer_crt_unit <- function(non_nulls, alpha, m, n, m_1, snr = 1, suppress.messages = TRUE, model = 'normal', 
                                       pch_methods=NULL, CRT.method_array=CRT.method_array, seed = NULL, ...){
  start.time = Sys.time()
  nulls = (1:m)[-non_nulls]
  
  data = gen_data(n = n, m = m, m_1 = m_1, model = model, snr = snr, alpha_vec = dirichlet_parameters, non_nulls=non_nulls)
  
  null_j = sample((m_1+1):min(sparsity_cutoff_m_values), size = 1)
  
  p = CRT_p_val(data = data, method = "lasso", CRT.replications = CRT.replications,CRT.method_array = CRT.method_array,seed = seed, null_j=null_j)
  
  remove(data)
  
  result = rep(list(NULL),length(sparsity_cutoff_m_values))
  for(i in 1:length(sparsity_cutoff_m_values)){
    sparsity_cutoff_m_value = sparsity_cutoff_m_values[i]
    
    pval_non_null = lapply(CRT.method_array, function(CRT.method)
      lapply(pch_methods, function(pch_method) PCH_pval(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = pch_method, sparsity_cutoff_m_value = sparsity_cutoff_m_value, j = 1)))
    pval_null = lapply(CRT.method_array, function(CRT.method)
      lapply(pch_methods, function(pch_method) PCH_pval(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = pch_method, sparsity_cutoff_m_value = sparsity_cutoff_m_value, j = null_j)))
    
    m_0 = m-m_1
    end.time = Sys.time()  
    
    if(suppress.messages == FALSE) print(end.time - start.time)
    
    result[[i]] = lapply(1:length(CRT.method_array), function(i) 
      list("pval_non_null" = unlist(pval_non_null[[i]]),
           "pval_null" = unlist(pval_null[[i]]),
           "rejection_non_null" = (unlist(pval_non_null[[i]]) <= alpha),
           "rejection_null" = (unlist(pval_null[[i]]) <= alpha), 
           "pch_methods" = pch_methods, 
           "sparsity_cutoff_m_value" = sparsity_cutoff_m_value,
           "null_index" = null_j,
           "p" = p[[i]], 
           "runtime" = end.time - start.time))
    names(result[[i]]) <- CRT.method_array
  }
  return(result)
}

simulate_fdr_fwer_crt <- function(simulation.iter = 2, alpha, m, n, m_1, snr = 1,model = 'normal', 
                                  pch_methods=pch_methods, CRT.method_array=CRT.method_array,
                                  non_null_probs = non_null_probs,...){
  result = rep(list(NULL),simulation.iter)
  for(i in 1:simulation.iter){
    seed_for_this_iteration = i + current_seed * simulation.iter + base_seed
    cat(paste0("\niteration ", i, ", base_seed = ", base_seed, ", current_seed = ", current_seed, ", seed for this iteration = ", seed_for_this_iteration, "\n"))
    set.seed(seed_for_this_iteration)
    non_nulls = sample(1:m, size = m_1, prob = non_null_probs, replace = FALSE)
    result[[i]] = simulate_fdr_fwer_crt_unit(non_nulls=non_nulls, 
                                             alpha = alpha, m = m, n = n, m_1 = m_1, snr = snr, model = model, 
                                             pch_methods=pch_methods, 
                                             CRT.method_array=CRT.method_array, seed = seed_for_this_iteration, ...)
  }
  return(result)
}

start_time = Sys.time()
result = simulate_fdr_fwer_crt(simulation.iter = simulation.iter, alpha = alpha,
                               m = m, n = n, m_1 = m_1, snr = snr,
                               p_filter_tolerance = p_filter_tolerance,
                               model = model,pch_methods=pch_methods,
                               CRT.method_array=CRT.method_array, non_null_probs = non_null_probs)
saveRDS(object = result, file = paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
                                       ",alpha=", alpha, ",CRT.replications=", CRT.replications,
                                       ",simulation.iter=", simulation.iter,
                                       ",snr=", snr, ",seed=", current_seed))
total_time = Sys.time()-start_time
print(paste0("Total time: ", round(total_time,2), " ",paste0(units(total_time))))

length(result[[1]])
result[[1]][[1]][[1]]$pval_non_null
result[[1]][[1]][[1]]$pch_methods

