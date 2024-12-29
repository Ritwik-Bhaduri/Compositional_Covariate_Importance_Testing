remove(list = ls())
versionid = 106
setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid))

source("fun_gen_CRT_pval.R")
source("fun_PCH.R")

taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# taskid = 30
set.seed(taskid)
cat("taskid = ", taskid)

##### Fixed Parameters ############################################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
non_nulls = 1:m_1
alpha = 0.05 # level of significance
CRT.replications = 1500 # number of knockoff variables generated
simulation.iter = 100 # number of replications of the simulations
target_simulation.iter = 5000

#### Variable parameters #############################################################
CRT.method_array = c("d_CRT_with_r_squared")
model_values = c("dirichlet")
# snr_values = c(0.01,0.02,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5)
snr_values = c(0.01,0.02,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5)

####### Methods ###############################################################################
oracle_r = c(m_1+1,m_1+2,round(m/4), round(m/2),m-1)

pch_methods = sapply(oracle_r, function(x) paste0("Simes_oracle_",x))
pch_methods = c(pch_methods, sapply(oracle_r, function(x) paste0("bonferroni_oracle_",x)))
pch_methods = c(pch_methods, "LOO")
###############################################################################################

base_seed = taskid*simulation.iter + versionid*1e4
print(paste0("base_seed = ", base_seed))

##### Variable Parameters arrays ##################################################### 
snr_array = rep(rep(snr_values, times = length(model_values)), each = ceiling(target_simulation.iter/simulation.iter))
model_array = as.vector(sapply(model_values, function(x) rep(x, each = length(snr_array)/length(model_values)))) # model for generating CRT resamples
seed_array = rep(1:ceiling(target_simulation.iter/simulation.iter), times = length(snr_values))
var_params = data.frame("snr" = snr_array, "model" = model_array, "seed" = seed_array)

##### Variable Parameters for this run ############################################### 
snr = snr_array[taskid]
model = model_array[taskid]
current_seed = var_params[taskid, "seed"]
print(paste0("model = ",model, ", snr = ", snr, ", seed = ", current_seed))

######################################################################################  
simulate_pch_crt_unit <- function(non_nulls=non_nulls, nulls = nulls, alpha, m, n, m_1, snr = 1, suppress.messages = TRUE, model = 'normal', 
                                       pch_methods=NULL, CRT.method_array=CRT.method_array, seed = NULL, ...){
  start.time = Sys.time()
  
  data = gen_data(n = n, m = m, m_1 = m_1, model = model, snr = snr, alpha_vec = rep(2, m))
  null_j = sample((m_1+1):m, size = 1)
  
  p = CRT_p_val(data = data, method = "lasso", CRT.replications = CRT.replications,CRT.method_array = CRT.method_array,seed = seed, null_j = null_j)
  remove(data)
  
  pval_non_null = lapply(CRT.method_array, function(CRT.method)
    lapply(pch_methods, function(pch_method) PCH_pval(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = pch_method, j = 1)))
  
  pval_null = lapply(CRT.method_array, function(CRT.method)
    lapply(pch_methods, function(pch_method) PCH_pval(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = pch_method, j = null_j)))
  
  m_0 = m-m_1
  end.time = Sys.time()  
  
  if(suppress.messages == FALSE) print(end.time - start.time)
  
  result = lapply(1:length(CRT.method_array), function(i) 
    list("pval_non_null" = unlist(pval_non_null[[i]]),
         "pval_null" = unlist(pval_null[[i]]),
         "rejection_non_null" = (unlist(pval_non_null[[i]]) <= alpha),
         "rejection_null" = (unlist(pval_null[[i]]) <= alpha), 
         "pch_methods" = pch_methods, 
         "null_index" = null_j,
         "p" = p[[i]], 
         "runtime" = end.time - start.time))
  names(result) <- CRT.method_array
  
  return(result)
  
}

simulate_pch_crt <- function(non_nulls=non_nulls,nulls=nulls, simulation.iter = 2,
                                  alpha, m, n, m_1, snr = 1,model = 'normal', pch_methods=pch_methods, 
                                  CRT.method_array=CRT.method_array,...){
  result = rep(list(NULL),simulation.iter)
  for(i in 1:simulation.iter){
    seed_for_this_iteration = i + current_seed * simulation.iter + base_seed
    cat(paste0("\niteration ", i, ", base_seed = ", base_seed, ", current_seed = ", current_seed, ", seed for this iteration = ", seed_for_this_iteration, "\n"))
    set.seed(seed_for_this_iteration)
    result[[i]] = simulate_pch_crt_unit(non_nulls=non_nulls, nulls = nulls, 
                                             alpha = alpha, m = m, n = n, m_1 = m_1, snr = snr, model = model, 
                                             pch_methods=pch_methods, 
                                             CRT.method_array=CRT.method_array, seed = seed_for_this_iteration, ...)
  }
  return(result)
}

start_time = Sys.time()
result = simulate_pch_crt(simulation.iter = simulation.iter, alpha = alpha,
                               m = m, n = n, m_1 = m_1, snr = snr,
                               model = model,pch_methods=pch_methods,
                               CRT.method_array=CRT.method_array)

saveRDS(object = result, file = paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
                                       ",alpha=", alpha, ",CRT.replications=", CRT.replications,
                                       ",simulation.iter=", simulation.iter,
                                       ",snr=", snr, ",seed=", current_seed))

total_time = Sys.time()-start_time
print(paste0("Total time: ", round(total_time,2), " ",paste0(units(total_time))))
