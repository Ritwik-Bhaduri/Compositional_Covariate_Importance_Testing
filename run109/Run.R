remove(list = ls())
versionid = 109
setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid))

source("fun_gen_CRT_pval_dirichlet_multinomial.R")
source("fun_FDR.R")
source("fun_FWER.R")

taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# taskid = 50
set.seed(taskid)
cat("taskid = ", taskid)

##### Fixed Parameters ############################################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data

# sparsity_cutoff_m_values = c(12,14)
# oracle_r = c(m_1+1,m_1+2,10)

sparsity_cutoff_m_values = c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
oracle_r = c(m_1+1,m_1+2,10, 20, 25)

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

alpha = 0.1 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
target_simulation.iter = 500
simulation.iter = 1 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable
CRT.method_array = c("d_CRT_with_r_squared")
model_values = c("dirichlet_multinomial")

####### Methods ###############################################################################
fdr_methods = c("univariate-BH", "max-BH", "adafilter","LOO-BH-BH", "DHH_BH")

DHH_storey_lambda = c(0.3,0.5,0.7)
fdr_methods = c(fdr_methods, sapply(DHH_storey_lambda,
                                    function(x) paste0("DHH_storey_",x)))
fdr_methods = c(fdr_methods, sapply(oracle_r,
                                    function(x) paste0("adafilter_oracle_",x)))
fdr_methods = c(fdr_methods, sapply(oracle_r, 
                                    function(x) paste0("Simes_oracle_",x)))
fdr_methods = c(fdr_methods, sapply(oracle_r,
                                    function(x) paste0("bonferroni_oracle_",x)))
fwer_methods = c("max-Bonferroni","max-Holm","univariate-Bonferroni","univariate-Holm",
                 "adaptive-Holm", "adaptive-Holm-aggressive")
fwer_methods = c(fwer_methods, sapply(oracle_r, 
                                      function(x) paste0("Simes_holm_oracle_",x)))
fwer_methods = c(fwer_methods, sapply(oracle_r, 
                                      function(x) paste0("adafilter_bonferroni_oracle_",x)))

###############################################################################################

base_seed = taskid*target_simulation.iter + versionid*simulation.iter*1000
print(paste0("base_seed = ", base_seed))

##### Variable Parameters arrays ##################################################### 
model_array = rep("dirichlet_multinomial", ceiling(target_simulation.iter/simulation.iter))
seed_array = rep(1:ceiling(target_simulation.iter/simulation.iter))
var_params = data.frame("model" = model_array, "seed" = seed_array)

##### Variable Parameters for this run ############################################### 
snr = 1
model = model_array[taskid]
current_seed = var_params[taskid, "seed"]
print(paste0("snr = ", snr, ", model = ",model, ", seed = ", current_seed))

######################################################################################  
simulate_fdr_fwer_crt_unit <- function(non_nulls, alpha, m, n, m_1, snr = 1, suppress.messages = TRUE, model = 'normal', 
                                       fdr_methods=NULL, fwer_methods=NULL, CRT.method_array=CRT.method_array, seed = NULL, ...){
  start.time = Sys.time()
  nulls = (1:m)[-non_nulls]
  
  data = gen_data(n = n, m = m, m_1 = m_1, model = model, snr = snr, alpha_vec = dirichlet_parameters, non_nulls=non_nulls)
  p = CRT_p_val(data = data, method = "lasso", CRT.replications = CRT.replications,CRT.method_array = CRT.method_array,seed = seed)
  remove(data)
  
  for(i in 1:length(sparsity_cutoff_m_values)){
    sparsity_cutoff_m_value = sparsity_cutoff_m_values[i]
    rejections_fdr = lapply(CRT.method_array, function(CRT.method)
      lapply(fdr_methods, function(fdr_method) FDR_control(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = fdr_method, 
                                                           sparsity_cutoff_m_value=sparsity_cutoff_m_value)))
    rejections_fwer = lapply(CRT.method_array, function(CRT.method)
      lapply(fwer_methods, function(fwer_method) FWER_control(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = fwer_method,
                                                              sparsity_cutoff_m_value=sparsity_cutoff_m_value)))
    m_0 = m-m_1
    end.time = Sys.time()  
    
    if(suppress.messages == FALSE) print(end.time - start.time)
    result_fdr_temp = lapply(1:length(CRT.method_array), function(i) 
      list("power" = sapply(rejections_fdr[[i]], function(x) length(intersect(x, non_nulls))/m_1),
           "fdr" = sapply(rejections_fdr[[i]], function(x) ifelse(length(x)>0,length(intersect(x,nulls))/length(x),0)),
           "rejections" = rejections_fdr[[i]], 
           "non_nulls" = non_nulls, 
           "fdr_methods" = fdr_methods, 
           "sparsity_cutoff_m_value" = sparsity_cutoff_m_value,
           "p" = p[[i]], 
           "runtime" = end.time - start.time))
    names(result_fdr_temp) <- CRT.method_array
    
    result_fwer_temp = lapply(1:length(CRT.method_array), function(i) 
      list("power" = sapply(rejections_fwer[[i]], function(x) length(intersect(x, non_nulls))/m_1),
           "fwer" = sapply(rejections_fwer[[i]], function(x) ifelse(length(intersect(x,nulls))>0, 1,0)),
           "rejections" = rejections_fwer[[i]], 
           "non_nulls" = non_nulls, 
           "fwer_methods" = fwer_methods, 
           "sparsity_cutoff_m_value" = sparsity_cutoff_m_value,
           "p" = p[[i]], 
           "runtime" = end.time - start.time))
    names(result_fwer_temp) <- CRT.method_array
    
    if(i==1){
      result_fdr  = result_fdr_temp
      result_fwer = result_fwer_temp
    }else{
      result_fdr = c(result_fdr, result_fdr_temp)
      result_fwer = c(result_fwer, result_fwer_temp)
    }
  }
  
  return(list(result_fdr, result_fwer))
  
}

simulate_fdr_fwer_crt <- function(simulation.iter = 2, alpha, m, n, m_1, snr = 1,model = 'normal', 
                                  fdr_methods=fdr_methods, fwer_methods=fwer_methods, CRT.method_array=CRT.method_array,
                                  non_null_probs = non_null_probs,...){
  result = rep(list(NULL),simulation.iter)
  for(i in 1:simulation.iter){
    seed_for_this_iteration = i + current_seed * simulation.iter + base_seed
    cat(paste0("\niteration ", i, ", base_seed = ", base_seed, ", current_seed = ", current_seed, ", seed for this iteration = ", seed_for_this_iteration, "\n"))
    set.seed(seed_for_this_iteration)
    non_nulls = sample(1:m, size = m_1, prob = non_null_probs, replace = FALSE)
    result[[i]] = simulate_fdr_fwer_crt_unit(non_nulls=non_nulls, 
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
                               CRT.method_array=CRT.method_array, non_null_probs = non_null_probs)
length(result[[1]][[1]])

fdr <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[1]][[1]]$fdr_methods), function(j)
  round(mean(sapply(1:length(sparsity_cutoff_m_values), function(i)  result[[1]][[1]][[i]]$fdr[j])),2)))
power_fdr <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[1]][[1]]$fdr_methods), function(j)
  round(mean(sapply(1:length(sparsity_cutoff_m_values), function(i)  result[[1]][[1]][[i]]$power[j])),2)))
rownames(fdr) = rownames(power_fdr) = result[[1]][[1]][[1]]$fdr_methods

fwer <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[2]][[1]]$fwer_methods), function(j)
  round(mean(sapply(1:length(sparsity_cutoff_m_values), function(i)  result[[1]][[2]][[i]]$fwer[j])),2)))
power_fwer <- sapply(CRT.method_array, function(CRT.method) sapply(1:length(result[[1]][[2]][[1]]$fwer_methods), function(j)
  round(mean(sapply(1:length(sparsity_cutoff_m_values), function(i)  result[[1]][[2]][[i]]$power[j])),2)))
rownames(fwer) = rownames(power_fwer) = result[[1]][[2]][[1]]$fwer_methods

saveRDS(object = result, file = paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
                                       ",alpha=", alpha, ",CRT.replications=", CRT.replications,
                                       ",simulation.iter=", simulation.iter,
                                       ",snr=", snr,
                                       ",p_filter_tolerance=", p_filter_tolerance, ",seed=", current_seed))

total_time = Sys.time()-start_time
print(paste0("Total time: ", round(total_time,2), " ",paste0(units(total_time))))

print(fdr)
print(power_fdr)

print(paste0("Total time: ", round(total_time,2), " ",paste0(units(total_time))))

# print(round(result[[1]][[1]]$d_CRT_with_r_squared$p$p_bivariate,1))
# pheatmap::pheatmap(result[[1]][[1]]$d_CRT_with_r_squared$p$p_bivariate, cluster_cols = FALSE, cluster_rows = FALSE)

