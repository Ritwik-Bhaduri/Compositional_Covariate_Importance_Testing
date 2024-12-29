remove(list = ls())
versionid = 92

setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid, "_1"))
source("fun_FDR.R")
source("fun_FWER.R")
##### Fixed Parameters ############################################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
alpha = 0.1 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
target_simulation.iter = 250
simulation.iter = 1 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

CRT.method_array = c("d_CRT_with_r_squared")
snr_values = c(0.01,0.02,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5)

model_array = c("dirichlet", "normal")
snr_list = list(snr_values)

#################################################### 
fdr_methods = c("univariate-BH", "max-BH", "adafilter","LOO-BH-BH", "DHH_BH")

DHH_storey_lambda = c(0.3,0.5,0.7)
fdr_methods = c(fdr_methods, sapply(DHH_storey_lambda,
                                    function(x) paste0("DHH_storey_",x)))
oracle_r = c(m_1+1,m_1+2,round(m/4), round(m/2),m-1)

fdr_methods = c(fdr_methods, sapply(oracle_r,
                                    function(x) paste0("adafilter_oracle_",x)))
fdr_methods = c(fdr_methods, sapply(oracle_r, 
                                    function(x) paste0("Simes_oracle_",x)))
fdr_methods = c(fdr_methods, sapply(oracle_r,
                                    function(x) paste0("bonferroni_oracle_",x)))
fwer_methods = c("max-Bonferroni","max-Holm","univariate-Bonferroni","univariate-Holm",
                 "adaptive-Holm")
fwer_methods = c(fwer_methods, sapply(oracle_r, 
                                      function(x) paste0("adaptive_Simes_holm_oracle_",x)))
fwer_methods = c(fwer_methods, sapply(oracle_r, 
                                      function(x) paste0("adafilter_bonferroni_oracle_",x)))

#################################################### 

snr_array_union = sort(unique(unlist(snr_list)))
result_file_name = function(snr_scale_factor, model, seed){
  paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
         ",alpha=", alpha, ",CRT.replications=", CRT.replications,
         ",simulation.iter=", simulation.iter,
         ",snr=", snr_scale_factor,
         ",p_filter_tolerance=", p_filter_tolerance,",seed=", seed)
}

sink("/dev/null")
lapply(c("parallelly", "foreach", "doParallel"), require, character.only = TRUE)
sink()

n.cores <- parallelly::availableCores()
my.cluster <- parallel::makeCluster(
  n.cores-1, 
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
# Create a function to update progress in a file
update_progress <- function(progress) {
  progress_file <- "progress.txt"
  writeLines(paste("Percent completed =", progress, "%"), progress_file)
}

system.time({
  for(i in 1:length(model_array)){
    for(j in 1:length(snr_array_union)){
      # The '%dopar%' operator indicates the loop's body is executed in parallel
      results <- foreach(k = 1:target_simulation.iter, .packages = c('methods'), .export = c('FDR_control', 'FWER_control', 'result_file_name', 'update_progress'), .combine = 'list') %dopar% {
        result = readRDS(result_file_name(snr_array_union[j], model = model_array[i],seed=k))
        p = lapply(1:length(CRT.method_array), function(l) result[[1]][[1]][[l]]$p)
        names(p) = CRT.method_array
        
        rejections_fdr = lapply(CRT.method_array, function(CRT.method)
          lapply(fdr_methods, function(fdr_method) FDR_control(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = fdr_method)))
        rejections_fwer = lapply(CRT.method_array, function(CRT.method)
          lapply(fwer_methods, function(fwer_method) FWER_control(p[[CRT.method]], alpha = alpha, m_1 = m_1, method = fwer_method)))
        
        m_0 = m-m_1
        non_nulls = 1:m_1 #result[[1]][[1]][[1]]$non_nulls
        result_fdr_temp = lapply(1:length(CRT.method_array), function(l) 
          list("non_nulls" = non_nulls, 
               "power" = sapply(rejections_fdr[[l]], function(x) length(intersect(x, non_nulls))/m_1),
               "fdr" = sapply(rejections_fdr[[l]], function(x) ifelse(length(x)>0,length(intersect(x,(1:m)[-non_nulls]))/length(x),0)),
               "rejections" = rejections_fdr[[l]], 
               "fdr_methods" = fdr_methods, 
               "p" = p[[l]], 
               "runtime" =  result[[1]][[1]][[l]]$runtime))
        names(result_fdr_temp) <- CRT.method_array
        
        result_fwer_temp = lapply(1:length(CRT.method_array), function(l) 
          list("non_nulls" = non_nulls,
               "power" = sapply(rejections_fwer[[l]], function(x) length(intersect(x, non_nulls))/m_1),
               "fwer" = sapply(rejections_fwer[[l]], function(x) ifelse(length(intersect(x,(1:m)[-non_nulls]))>0, 1,0)),
               "rejections" = rejections_fwer[[l]], 
               "fwer_methods" = fwer_methods, 
               "p" = p[[l]], 
               "runtime" =  result[[1]][[l]]$runtime))
        names(result_fwer_temp) <- CRT.method_array
        
        # Construct the result to be returned to the master process
        result = list(list(result_fdr_temp, result_fwer_temp))
        
        saveRDS(object = result, file = paste0("result files from pval/Result-", "model=", model_array[i],
                                               ",n=", n, ",m=", m, ",m_1=", m_1,
                                               ",alpha=", alpha, ",CRT.replications=", CRT.replications,
                                               ",simulation.iter=", simulation.iter,
                                               ",snr=", snr_array_union[j],
                                               ",p_filter_tolerance=", p_filter_tolerance,",seed=",k))
        
        # Return a status or filename to indicate completion
        cat(paste("i=",i,", j=",j, ", k=",k, "\n"), file = "progress.txt")
        
        return(paste0("Completed: ", k))
      }

      # # Update the progress
      # progress <- (j+(i-1)*length(snr_array_union)) / (length(snr_array_union)*length(model_array)) * 100
      # update_progress(progress)
    }
  }
})

stopCluster(my.cluster)
