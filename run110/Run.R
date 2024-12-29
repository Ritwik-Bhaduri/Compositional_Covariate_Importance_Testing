remove(list = ls())
versionid = 110
setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid))

# source("fun_gen_CRT_pval_parallel.R")
source("fun_gen_CRT_pval.R")
source("fun_FDR.R")
source("fun_FWER.R")

taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# taskid = 3000
set.seed(taskid)
cat("taskid = ", taskid)

##### Fixed Parameters ############################################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
non_nulls = 1:m_1
alpha = 0.1 # target fdr level
CRT.replications = 25000 # number of knockoff variables generated
target_simulation.iter = 500
simulation.iter = 1 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

#### Variable parameters #############################################################
CRT.method_array = c("d_CRT_with_r_squared")
model_values = c("dirichlet")
# snr_values = c(0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
snr_values = c(0.01,0.02,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5)

####### Methods ###############################################################################
fdr_methods = c("univariate-BH", "max-BH", "adafilter","LOO-BH-BH", "DHH_BH")

DHH_storey_lambda = c(0.3,0.5,0.7)
fdr_methods = c(fdr_methods, sapply(DHH_storey_lambda,
                                    function(x) paste0("DHH_storey_",x)))
oracle_r = sort(union(c(m_1+1,m_1+2,round(m/4), round(m/2),m-1), c(11,20,30,40,50,60,70,80,90)))

# fdr_methods = c(fdr_methods, sapply(oracle_r,
#                                     function(x) paste0("adafilter_oracle_",x)))
fdr_methods = c(fdr_methods, sapply(oracle_r, 
                                    function(x) paste0("Simes_oracle_",x)))
fdr_methods = c(fdr_methods, sapply(oracle_r,
                                    function(x) paste0("bonferroni_oracle_",x)))
fwer_methods = c("max-Bonferroni","max-Holm","univariate-Bonferroni","univariate-Holm",
                 "adaptive-Holm")
fwer_methods = c(fwer_methods, sapply(oracle_r, 
                                      function(x) paste0("adaptive_Simes_holm_oracle_",x)))
fwer_methods = c(fwer_methods, sapply(oracle_r,
                                      function(x) paste0("adaptive_bonferroni_holm_oracle_",x)))
# fwer_methods = c(fwer_methods, sapply(oracle_r, 
#                                       function(x) paste0("adafilter_bonferroni_oracle_",x)))

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
simulate_fdr_fwer_crt_unit <- function(non_nulls=non_nulls, nulls = nulls, alpha, m, n, m_1, snr = 1, suppress.messages = TRUE, model = 'normal', 
                                       fdr_methods=NULL, fwer_methods=NULL, CRT.method_array=CRT.method_array, seed = NULL, ...){
  start.time = Sys.time()
  
  data = gen_data(n = n, m = m, m_1 = m_1, model = model, snr = snr, alpha_vec = rep(2, m))
  p = CRT_p_val(data = data, method = "lasso", CRT.replications = CRT.replications,CRT.method_array = CRT.method_array,seed = seed)
  # saveRDS(p, paste0("result files/p_", model, "_", snr, "_.rds"))
  # p = readRDS(paste0("result files/p_", model, "_", snr, "_.rds"))
  remove(data)
  
  # if(model == 'dirichlet') fdr_methods = fdr_methods[fdr_methods != "univariate-BH"]
  rejections_fdr = lapply(CRT.method_array, function(CRT.method)
    lapply(fdr_methods, function(fdr_method) FDR_control(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = fdr_method)))
  rejections_fwer = lapply(CRT.method_array, function(CRT.method)
    lapply(fwer_methods, function(fwer_method) FWER_control(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = fwer_method)))
  
  m_0 = m-m_1
  end.time = Sys.time()  
  
  if(suppress.messages == FALSE) print(end.time - start.time)
  
  result_fdr = lapply(1:length(CRT.method_array), function(i) 
    list("power" = sapply(rejections_fdr[[i]], function(x) length(intersect(x, 1:m_1))/m_1),
         "fdr" = sapply(rejections_fdr[[i]], function(x) ifelse(length(x)>0,length(intersect(x,m_1 + 1:m_0))/length(x),0)),
         "rejections" = rejections_fdr[[i]], 
         "fdr_methods" = fdr_methods, 
         "p" = p[[i]], 
         "runtime" = end.time - start.time))
  names(result_fdr) <- CRT.method_array
  
  result_fwer = lapply(1:length(CRT.method_array), function(i) 
    list("power" = sapply(rejections_fwer[[i]], function(x) length(intersect(x, 1:m_1))/m_1),
         "fwer" = sapply(rejections_fwer[[i]], function(x) ifelse(length(intersect(x,m_1 + 1:m_0))>0, 1,0)),
         "rejections" = rejections_fwer[[i]], 
         "fwer_methods" = fwer_methods, 
         "p" = p[[i]], 
         "runtime" = end.time - start.time))
  names(result_fwer) <- CRT.method_array
  return(list(result_fdr, result_fwer))
  
}

simulate_fdr_fwer_crt <- function(non_nulls=non_nulls,nulls=nulls, simulation.iter = 2,
                                  alpha, m, n, m_1, snr = 1,model = 'normal', fdr_methods=fdr_methods, 
                                  fwer_methods=fwer_methods, CRT.method_array=CRT.method_array,...){
  result = rep(list(NULL),simulation.iter)
  for(i in 1:simulation.iter){
    seed_for_this_iteration = i + current_seed * simulation.iter + base_seed
    cat(paste0("\niteration ", i, ", base_seed = ", base_seed, ", current_seed = ", current_seed, ", seed for this iteration = ", seed_for_this_iteration, "\n"))
    set.seed(seed_for_this_iteration)
    result[[i]] = simulate_fdr_fwer_crt_unit(non_nulls=non_nulls, nulls = nulls, 
                                             alpha = alpha, m = m, n = n, m_1 = m_1, snr = snr, model = model, 
                                             fdr_methods=fdr_methods, fwer_methods=fwer_methods, 
                                             CRT.method_array=CRT.method_array, seed = seed_for_this_iteration, ...)
  }
  return(result)
}

start_time = Sys.time()
result = simulate_fdr_fwer_crt(simulation.iter = simulation.iter, alpha = alpha,
                               m = m, n = n, m_1 = m_1, snr = snr,
                               p_filter_tolerance = p_filter_tolerance,
                               model = model,fdr_methods=fdr_methods,fwer_methods=fwer_methods,
                               CRT.method_array=CRT.method_array)

fdr <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[1]][[1]]$fdr_methods), function(j)
  round(mean(sapply(1:simulation.iter, function(i)  result[[i]][[1]][[CRT.method]]$fdr[j])),2)))
power_fdr <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[1]][[1]]$fdr_methods), function(j)
  round(mean(sapply(1:simulation.iter, function(i)  result[[i]][[1]][[CRT.method]]$power[j])),2)))
rownames(fdr) = rownames(power_fdr) = result[[1]][[1]][[1]]$fdr_methods

fwer <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[2]][[1]]$fwer_methods), function(j)
  round(mean(sapply(1:simulation.iter, function(i)  result[[i]][[2]][[CRT.method]]$fwer[j])),2)))
power_fwer <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[2]][[1]]$fwer_methods), function(j)
  round(mean(sapply(1:simulation.iter, function(i)  result[[i]][[2]][[CRT.method]]$power[j])),2)))
rownames(fwer) = rownames(power_fwer) = result[[1]][[2]][[1]]$fwer_methods

saveRDS(object = result, file = paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
                                       ",alpha=", alpha, ",CRT.replications=", CRT.replications,
                                       ",simulation.iter=", simulation.iter,
                                       ",snr=", snr,
                                       ",p_filter_tolerance=", p_filter_tolerance, ",seed=", current_seed))

total_time = Sys.time()-start_time
print(paste0("Total time: ", round(total_time,2), " ",paste0(units(total_time))))

# print(fdr)
print(power_fdr)

print(paste0("Total time: ", round(total_time,2), " ",paste0(units(total_time))))

# print(round(result[[1]][[1]]$d_CRT_with_r_squared$p$p_bivariate,1))
# pheatmap::pheatmap(result[[1]][[1]]$d_CRT_with_r_squared$p$p_bivariate, cluster_cols = FALSE, cluster_rows = FALSE)
