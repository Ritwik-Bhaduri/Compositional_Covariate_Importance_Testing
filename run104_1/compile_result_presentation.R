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
simulation.iter = 500 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

CRT.method_array = c("d_CRT_with_r_squared")
model_array = c('dirichlet') # model for generating X
snr_list = list(c( 0.01,0.02,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5))
# snr_list = list(c(0.01,0.5,1,2,5))
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

df <- foreach(k = 1:length(file_list), .combine = 'rbind') %dopar%{
  file = file_list[k]
  result = readRDS(paste0("result files/", file))
  d_fdr = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                     'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                     'control_methods' = result[[1]][[1]][[1]][[CRT.method]]$fdr_methods,
                     'CRT_methods' = CRT.method,
                     'error'  = result[[1]][[1]][[1]][[CRT.method]]$fdr,
                     'power' = result[[1]][[1]][[1]][[CRT.method]]$power, 
                     'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                     'error_type' = 'fdr')
  d_fwer = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                      'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                      'control_methods' = result[[1]][[1]][[2]][[CRT.method]]$fwer_methods,
                      'CRT_methods' = CRT.method,
                      'error'  = result[[1]][[1]][[2]][[CRT.method]]$fwer,
                      'power' = result[[1]][[1]][[2]][[CRT.method]]$power, 
                      'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                      'error_type' = 'fwer')
  cat("i=",i,", j=", j,", k=", k,file = "progress.txt")
  rbind(d_fdr,d_fwer)
}
parallel::stopCluster(my.cluster)
df_backup = df
df = df_backup

df[df$control_methods=="adaptive-Holm", "control_methods"] = "Max-Holm"
df[df$control_methods=="univariate-BH", "control_methods"] = "Univariate-BH"
df[df$control_methods=="univariate-Holm", "control_methods"] = "Univariate-Holm"
df[df$control_methods=="Simes_oracle_11", "control_methods"] = "CB-PCH-oracle"
df[df$control_methods=="Simes_oracle_50", "control_methods"] = "CB-PCH"
df[df$control_methods=="adaptive_Simes_holm_oracle_11", "control_methods"] = "CB-PCH-oracle"
df[df$control_methods=="adaptive_Simes_holm_oracle_50", "control_methods"] = "CB-PCH"
df[df$control_methods=="adafilter_oracle_11", "control_methods"] = "adafilter-BH-10"
df[df$control_methods=="adafilter_oracle_50", "control_methods"] = "adafilter-BH-50"
df[df$control_methods=="adafilter_bonferroni_oracle_11", "control_methods"] = "adafilter-bonferroni-10  "
df[df$control_methods=="adafilter_bonferroni_oracle_50", "control_methods"] = "adafilter-bonferroni-50  "

df_speedup = df_speedup_backup =df

####################################################################################################
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run104_1/") # for the slow version

n.cores <- parallelly::availableCores()
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
i=j=k=1;CRT.method = CRT.method_array[1]
doParallel::registerDoParallel(cl = my.cluster)

file_list = list.files(path="./result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]

df <- foreach(k = 1:length(file_list), .combine = 'rbind') %dopar%{
  file = file_list[k]
  result = readRDS(paste0("result files/", file))
  d_fdr = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                     'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                     'control_methods' = result[[1]][[2]][[1]][[CRT.method]]$fdr_methods,
                     'CRT_methods' = CRT.method,
                     'error'  = result[[1]][[2]][[1]][[CRT.method]]$fdr,
                     'power' = result[[1]][[2]][[1]][[CRT.method]]$power, 
                     'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                     'error_type' = 'fdr')
  d_fwer = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                      'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                      'control_methods' = result[[1]][[2]][[2]][[CRT.method]]$fwer_methods,
                      'CRT_methods' = CRT.method,
                      'error'  = result[[1]][[2]][[2]][[CRT.method]]$fwer,
                      'power' = result[[1]][[2]][[2]][[CRT.method]]$power, 
                      'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                      'error_type' = 'fwer')
  cat("i=",i,", j=", j,", k=", k,file = "progress.txt")
  rbind(d_fdr,d_fwer)
}
parallel::stopCluster(my.cluster)
df_backup2 = df

df[df$control_methods=="adaptive-Holm", "control_methods"] = "Max-Holm"
df[df$control_methods=="univariate-BH", "control_methods"] = "Univariate-BH"
df[df$control_methods=="univariate-Holm", "control_methods"] = "Univariate-Holm"
df[df$control_methods=="Simes_oracle_11", "control_methods"] = "CB-PCH-oracle"
df[df$control_methods=="Simes_oracle_50", "control_methods"] = "CB-PCH"
df[df$control_methods=="adaptive_Simes_holm_oracle_11", "control_methods"] = "CB-PCH-oracle"
df[df$control_methods=="adaptive_Simes_holm_oracle_50", "control_methods"] = "CB-PCH"
df[df$control_methods=="adafilter_oracle_11", "control_methods"] = "adafilter-BH-10"
df[df$control_methods=="adafilter_oracle_50", "control_methods"] = "adafilter-BH-50"
df[df$control_methods=="adafilter_bonferroni_oracle_11", "control_methods"] = "adafilter-bonferroni-10  "
df[df$control_methods=="adafilter_bonferroni_oracle_50", "control_methods"] = "adafilter-bonferroni-50  "

df_slow = df_slow_backup =df
############################################################################################################
df_slow <- df_slow %>%
  mutate(
    Method = case_when(
      str_detect(control_methods, "CB-PCH-oracle") & error_type=="fdr" ~ "BCP(s+1)-BH",
      str_detect(control_methods, "CB-PCH-oracle") & error_type=="fwer" ~ "BCP(s+1)-Holm",
      str_detect(control_methods, "CB-PCH") & error_type=="fdr" ~ "BCP(p/2)-BH",
      str_detect(control_methods, "CB-PCH") & error_type=="fwer" ~ "BCP(p/2)-Holm",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Method %in% c("BCP(p/2)-BH","BCP(s+1)-BH","BCP(p/2)-Holm","BCP(s+1)-Holm"))
df_slow$control_methods = factor(df_slow$Method, levels = c("BCP(p/2)-BH","BCP(p/2)-Holm", "BCP(s+1)-BH", "BCP(s+1)-Holm"))

df_speedup <- df_speedup %>%
  mutate(
    Method = case_when(
      str_detect(control_methods, "CB-PCH-oracle") & error_type=="fdr" ~ "BCP(s+1)-BH",
      str_detect(control_methods, "CB-PCH-oracle") & error_type=="fwer" ~ "BCP(s+1)-Holm",
      str_detect(control_methods, "CB-PCH") & error_type=="fdr" ~ "BCP(p/2)-BH",
      str_detect(control_methods, "CB-PCH") & error_type=="fwer" ~ "BCP(p/2)-Holm",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Method %in% c("BCP(p/2)-BH","BCP(s+1)-BH","BCP(p/2)-Holm","BCP(s+1)-Holm"))
df_speedup$control_methods = factor(df_speedup$Method, levels = c("BCP(p/2)-BH","BCP(p/2)-Holm", "BCP(s+1)-BH", "BCP(s+1)-Holm"))

##########################################################################################################

fdr_method_shortlist = c("BCP(p/2)-BH","BCP(s+1)-BH")
fwer_method_shortlist = c("BCP(p/2)-Holm","BCP(s+1)-Holm") 
df_speedup$control_methods = factor(sapply(df_speedup$control_methods, function(x) paste0(x, "_fast")), 
                                    levels = c("BCP(p/2)-BH_fast","BCP(p/2)-Holm_fast", "BCP(s+1)-BH_fast", "BCP(s+1)-Holm_fast"))
df = rbind(df_speedup, df_slow)
df$control_methods = factor(df$control_methods, levels = c("BCP(p/2)-BH_fast","BCP(p/2)-BH", "BCP(s+1)-BH_fast", "BCP(s+1)-BH",
                                                           "BCP(p/2)-Holm_fast","BCP(p/2)-Holm", "BCP(s+1)-Holm_fast", "BCP(s+1)-Holm"))
df_backup = df
cbPalette <- c("#4767CD", "#FE0000")[c(1,1,2,2)]
fdr_method_shortlist = c("BCP(p/2)-BH","BCP(s+1)-BH", "BCP(p/2)-BH_fast","BCP(s+1)-BH_fast")
fwer_method_shortlist = c("BCP(p/2)-Holm","BCP(s+1)-Holm", "BCP(p/2)-Holm_fast","BCP(s+1)-Holm_fast") 

model_name = "dirichlet"
snr_array = snr_list[[1]]

convert_to_expression <- function(label) {
  # Replace 's+1' with 's+1' etc.
  label <- gsub("p_([0-9]+)", "p[\\1]", label)
  parse(text=label)
}

my_plot <- function(name = "ERROR", model_name = "dirichlet", error_type_name = 'fdr'){
  df = subset(df, snr %in% snr_array &  model %in% model_name & error_type %in% error_type_name)
  df = subset(df, CRT_methods %in% CRT.method_array)
  
  df.summary <- subset(df, model == model_name & CRT_methods == "d_CRT_with_r_squared" & error_type %in% error_type_name) %>%
    group_by(snr, control_methods) %>%
    summarise(
      upper_error = min(mean(error) + 1.96 * sd(error)/sqrt(simulation.iter), rep(1,length(snr))),
      lower_error = max(mean(error) - 1.96 * sd(error)/sqrt(simulation.iter), rep(0,length(snr))),
      error = mean(error),
      upper_power = min(mean(power) + 1.96 * sd(power)/sqrt(simulation.iter), rep(1,length(snr))),
      lower_power = max(mean(power) - 1.96 * sd(power)/sqrt(simulation.iter), rep(0,length(snr))),
      power = mean(power)
    )
  df.summary = subset(df.summary, snr %in% snr_array)
  if(error_type_name == 'fdr') df.summary = subset(df.summary, control_methods %in% fdr_method_shortlist)
  if(error_type_name == 'fwer') df.summary = subset(df.summary, control_methods %in% fwer_method_shortlist)
  
  n_methods = ifelse(error_type_name == 'fdr', length(fdr_method_shortlist), length(fwer_method_shortlist))
  
  p = ggplot(df.summary, aes_string("snr", tolower(name))) +
    geom_line(aes(group = control_methods, color = control_methods, linetype = control_methods), data = df.summary) +
    geom_errorbar(aes(ymin = get(paste0("lower_", tolower(name))), ymax = get(paste0("upper_", tolower(name))), 
                      color = control_methods), width = 0, show.legend = FALSE) +
    theme_classic() + 
    labs(x = "SNR", y = ifelse(name == 'Power', 'Average Power', toupper(error_type_name)), 
         color = "Methods ", linetype = "Methods ")+
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text.align = 0)
  if(name == 'ERROR') p=p+coord_cartesian(ylim = c(0, 6*alpha)) + geom_hline(aes(yintercept = 0.1), linetype = "dashed") 
  if(name == 'Power') p=p+coord_cartesian(ylim = c(0, 1))

  unique_methods <- as.character(unique(df.summary$control_methods))
  labels_expr <- sapply(unique_methods, convert_to_expression)
  p = p + scale_color_manual(values = cbPalette, labels = labels_expr) +
    scale_linetype_manual(values = c("solid","dashed")[c(1,2,1,2)],
                          labels = labels_expr)
  return(p)
  
}

setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run104_1/")

fdr_plot_1 = my_plot(name = "ERROR", model_name = "dirichlet", error_type_name = 'fdr')
fdr_plot_2 = my_plot(name = "Power", model_name = "dirichlet", error_type_name = 'fdr')
fdr_plot = ggarrange(fdr_plot_1, fdr_plot_2, common.legend = TRUE, legend = "right")
fdr_plot
ggsave(paste0("plots/presentation/","fdr_dirichlet.png"), 
       plot=fdr_plot, device = "png", width = 7, height = 3, units = "in", bg="white",dpi = 600)
ggsave("plots/presentation/fdr_dirichlet.pdf", plot = fdr_plot, device = "pdf", width = 7, height = 3, units = "in")

df = df_backup
fwer_plot_1 = my_plot(name = "ERROR", model_name = "dirichlet", error_type_name = 'fwer')
fwer_plot_2 = my_plot(name = "Power", model_name = "dirichlet", error_type_name = 'fwer')
fwer_plot = ggarrange(fwer_plot_1, fwer_plot_2, common.legend = TRUE, legend = "right")
fwer_plot 
ggsave(paste0("plots/presentation/","fwer_dirichlet.png"), 
       plot=fwer_plot, device = "png", width = 7, height = 3, units = "in", bg="white",dpi = 600)
ggsave("plots/presentation/fwer_dirichlet.pdf", plot = fwer_plot, device = "pdf", width = 7, height = 3, units = "in")

