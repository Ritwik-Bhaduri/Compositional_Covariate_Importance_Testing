remove(list = ls())
versionid = 114

setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid))
source("fun_FDR.R")
source("fun_FWER.R")
##### Fixed Parameters ############################################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
sparsity_cutoff_m_values = c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
oracle_r = c(m_1+1,m_1+2,10, 20, 25)

alpha = 0.1 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
target_simulation.iter = 500
simulation.iter = 1 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable
CRT.method_array = c("d_CRT_with_r_squared")
model_values = c("dirichlet_multinomial")
snr = 1
model = model_values[1]
#################################################### 
fdr_methods = c("univariate-BH", "max-BH", "LOO-BH-BH", "DHH_BH")

DHH_storey_lambda = c(0.3,0.5,0.7)
fdr_methods = c(fdr_methods, sapply(DHH_storey_lambda,
                                    function(x) paste0("DHH_storey_",x)))
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

#################################################### 
file_list = list.files(path="../run109/result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]
file_list <- file_list[grep(paste0("simulation.iter=",simulation.iter,","), file_list)]
sink("/dev/null")
lapply(c("parallelly", "foreach", "doParallel"), require, character.only = TRUE)
sink()

n.cores <- parallelly::availableCores()
my.cluster <- parallel::makeCluster(
  n.cores-1, 
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)

foreach(file_index = 1:length(file_list), .packages = c('methods'), .combine = 'list') %dopar% {
# for(file_index in 1:length(file_list)){
  cat(paste("file index=",file_index, "out of ", length(file_list), "\n"), file = "progress.txt")
  file_name = file_list[file_index]
  result_original = readRDS(paste0("../run109/result files/", file_name)) 
  result_new = list()
  for(simulation_index in 1:length(result_original)){
    p = lapply(1:length(sparsity_cutoff_m_values), function(l) result_original[[simulation_index]][[1]][[l]]$p)
    names(p) = sparsity_cutoff_m_values
    
    rejections_fdr = lapply(sparsity_cutoff_m_values, function(sparsity_cutoff_m)
      lapply(fdr_methods, function(fdr_method) FDR_control(p[[as.character(sparsity_cutoff_m)]], alpha = alpha, m_1 = m_1, method = fdr_method, 
                                                           sparsity_cutoff_m_value=sparsity_cutoff_m)))
    rejections_fwer = lapply(sparsity_cutoff_m_values, function(sparsity_cutoff_m)
      lapply(fwer_methods, function(fwer_method) FWER_control(p[[as.character(sparsity_cutoff_m)]], alpha = alpha, m_1 = m_1, method = fwer_method, 
                                                              sparsity_cutoff_m_value=sparsity_cutoff_m)))
    
    m_0 = m-m_1
    non_nulls = result_original[[simulation_index]][[1]][[1]]$non_nulls
    result_fdr_temp = lapply(1:length(sparsity_cutoff_m_values), function(l) 
      list("non_nulls" = non_nulls, 
           "power" = sapply(rejections_fdr[[l]], function(x) length(intersect(x, non_nulls))/m_1),
           "fdr" = sapply(rejections_fdr[[l]], function(x) ifelse(length(x)>0,length(intersect(x,(1:m)[-non_nulls]))/length(x),0)),
           "rejections" = rejections_fdr[[l]], 
           "fdr_methods" = fdr_methods, 
           "p" = p[[l]], 
           "runtime" =  result_original[[simulation_index]][[1]][[l]]$runtime))
    names(result_fdr_temp) <- sparsity_cutoff_m_values
    
    result_fwer_temp = lapply(1:length(sparsity_cutoff_m_values), function(l) 
      list("non_nulls" = non_nulls,
           "power" = sapply(rejections_fwer[[l]], function(x) length(intersect(x, non_nulls))/m_1),
           "fwer" = sapply(rejections_fwer[[l]], function(x) ifelse(length(intersect(x,(1:m)[-non_nulls]))>0, 1,0)),
           "rejections" = rejections_fwer[[l]], 
           "fwer_methods" = fwer_methods, 
           "p" = p[[l]], 
           "runtime" =  result_original[[simulation_index]][[2]][[l]]$runtime))
    names(result_fwer_temp) <- sparsity_cutoff_m_values
    
    result = list(list(result_fdr_temp, result_fwer_temp))
    result_new[[simulation_index]] = result
  }
  saveRDS(object = result_new, file = paste0("result files from pval/Result-", "model=", model,
                  ",n=", n, ",m=", m, ",m_1=", m_1,
                  ",alpha=", alpha, ",CRT.replications=", CRT.replications,
                  ",simulation.iter=", simulation.iter,
                  ",snr=", snr,
                  ",p_filter_tolerance=", p_filter_tolerance,",seed=",strsplit(file_name, "seed=")[[1]][2]))
  cat(paste("file index=",file_index, "out of ", length(file_list), "\n"), file = "progress.txt")
  return(paste0("Completed: ", file_index))
}
stopCluster(my.cluster)

