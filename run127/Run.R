remove(list = ls())
versionid = 127
setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid))

source("fun_gen_CRT_pval_parallel_temp.R")
source("fun_FDR.R")
source("fun_FWER.R")

taskid = 1
set.seed(taskid)

######## data #######################################################################
checkData <- function(compositions,totals,labels){
  X = as.matrix(compositions)
  reads = as.matrix(totals)
  q = ncol(X)
  counts = matrix(NA,nrow(X),ncol(X))
  for (j in 1:q){
    counts[,j] = round(X[,j]*reads)
  }
  
  idx0 = which(apply(counts,2,sum)<=5)
  X.counts = counts[,apply(counts,2,sum)>5]
  idx = which(apply(X.counts,1,sum)<10)
  X.counts = X.counts[-idx,]
  Y = as.matrix(labels)
  Y = Y[-idx]
  X = X[-idx,]
  outputs = list(compositions=X,counts = X.counts, labels = Y)
  return(outputs)
}

compositions = read.csv('./data/HMPdataL2.txt',header=FALSE)
totals = read.table("./data/HMPdataL2reads.txt")[[1]]
labels = read.table("./data/HMPlabels.txt")[[1]]
HMPdata = checkData(compositions,totals,labels)

X = HMPdata$counts
totals = rowSums(X)
Y = HMPdata$labels

n = nrow(X)
m = ncol(X)

## reorder X according to sparsity
col_order = order(colSums(X > 0), decreasing = TRUE)
X = X[, col_order]
colnames(X) = sapply(1:ncol(X), function(i) paste0("X", i)) 
X_prop = X / totals

##### Fixed Parameters ############################################################### 
sparsity_cutoff_m_values = sapply(c(0.01, 0.02, 0.03), function(sparsity) sum(colMeans(X>0) >= sparsity))
oracle_r_list = list(c(5, 7, 13, 19), c(7, 13, 16), c(7, 13, 14))
# oracle_r_list = list(c(4))

alpha = 0.1 # target fdr level
CRT.replications = 2500 # number of knockoff variables generated
CRT.method_array = c("trial")
model = "two-step"

data = list(X = X, X_prop = X_prop, totals = totals, Y = Y,
            model.X = list("model" = model, "step1" = "logistic", "step2" = "logistic-lm"))

start.time = Sys.time()
p = CRT_p_val(data = data, CRT.replications = CRT.replications, 
              CRT.method_array = CRT.method_array, seed = taskid)
end.time = Sys.time()
end.time-start.time

####################################################################################################
library(tibble); library(dplyr); library(tidyr); library(knitr)

res <- tibble(r=integer(), cutoff=integer(), result=character())
for(r in 5:23) {
  for(cutoff in c(14,20,23)) {
    p_temp <- p$trial$p_bivariate[1:cutoff, 1:cutoff]
    p_no_diag <- matrix(NA, nrow(p_temp)-1, ncol(p_temp))
    for(j in seq_len(ncol(p_temp))) p_no_diag[,j] <- p_temp[-j, j]
    p_simes <- apply(p_no_diag, 2, function(x) Simes(x, r))
    
    rejections_raw = myBH(p_simes, q = alpha)$rejections # FDR results
    # rejections_raw = Holm(p_simes, q = alpha)$rejections # FWER results
    rejections_formatted <- sort(col_order[rejections_raw])
    res <- add_row(res, r=r, cutoff=cutoff, result=if(length(rejections_formatted)) paste(rejections_formatted, collapse=", ") else "")
  }
}

df_wide <- res %>%
  pivot_wider(names_from=cutoff,
              values_from=result,
              names_prefix="cutoff_") %>%
  arrange(r)

View(df_wide)

kable(df_wide,format="latex",booktabs=TRUE,align="c|c|c|c",
      caption="Results of varying $r$ and sparsity cutoff",
      label="tab:results",
      col.names=c("$r$","cutoff = 14","cutoff = 20","cutoff = 23"))
