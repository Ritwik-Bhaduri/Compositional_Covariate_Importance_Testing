remove(list = ls())
library(ggplot2); library(stringr)
library(ggpubr)
library(dplyr)
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run104_1/")

##### Parameters ################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
alpha = 0.1 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
simulation.iter = 1 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

CRT.method_array = c("d_CRT_with_r_squared")
model_array = c('dirichlet') # model for generating X
# snr_list = list(c( 0.01,0.02,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5))
snr_list = list(c(0.01,0.5,1,2,5))
CI_coverage = 90
#################################################### 

snr_array_union = sort(unique(unlist(snr_list)))
result_file_name = function(snr_scale_factor, model,seed){
  paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
         ",alpha=", alpha, ",CRT.replications=", CRT.replications,
         ",simulation.iter=", simulation.iter,
         ",snr=", snr_scale_factor,
         ",p_filter_tolerance=", p_filter_tolerance,",seed=",seed)
}
df <- data.frame(matrix(NA, nrow = 2*length(model_array)*length(snr_array_union)*simulation.iter*length(CRT.method_array), ncol = 7))
sink("/dev/null")
lapply(c("parallelly", "foreach", "doParallel"), require, character.only = TRUE)
sink()

n.cores <- parallelly::availableCores()
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
i=j=k=1;CRT.method = CRT.method_array[1]
doParallel::registerDoParallel(cl = my.cluster)

file_list = list.files(path="./result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]

df <- foreach(k = 1:length(file_list)) %:%
  foreach(i = 1:simulation.iter) %dopar%{
    file = file_list[k]
    result = readRDS(paste0("result files/", file))
    rejections_fdr_cbpch  = list(result[[i]][[1]][[1]][[CRT.method]]$rejections[[17]], result[[i]][[2]][[1]][[CRT.method]]$rejections[[17]])
    rejections_fwer_cbpch = list(result[[i]][[1]][[2]][[CRT.method]]$rejections[[9]], result[[i]][[2]][[2]][[CRT.method]]$rejections[[9]])
    cat("k=", k,file = "progress.txt")
    list(rejections_fdr_cbpch, rejections_fwer_cbpch)
  }
parallel::stopCluster(my.cluster)
df_backup = df

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(union(a,b))
  if(union==0) return(1)
  return(intersection/union)
}

jaccard_fdr = jaccard_fwer = c()
for(i in 1:length(df)){
  rejections_fdr_fast = df[[i]][[1]][[1]][[1]]
  rejections_fdr_slow = df[[i]][[1]][[1]][[2]]
  rejections_fwer_fast = df[[i]][[1]][[2]][[1]]
  rejections_fwer_slow = df[[i]][[1]][[2]][[2]]
  jaccard_fdr[i] = jaccard(rejections_fdr_fast, rejections_fdr_slow)
  jaccard_fwer[i] = jaccard(rejections_fwer_fast, rejections_fwer_slow)
}
mean(jaccard_fdr)
mean(jaccard_fwer)

plot(jaccard_fdr)
plot(jaccard_fwer)

