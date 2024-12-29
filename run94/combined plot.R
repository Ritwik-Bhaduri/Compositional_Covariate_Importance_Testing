remove(list = ls())
library(ggplot2); library(stringr)
library(ggpubr)
library(dplyr)
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run94/")

##### Parameters ################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
alpha = 0.1 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
target_simulation.iter = 500
simulation.iter = 25 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

CRT.method_array = c("d_CRT_with_r_squared"); CRT.method=CRT.method_array[1]
model_array = c('dirichlet') # model for generating X
snr_list = list(c( 0.01,0.02,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5))
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

file_list = list.files(path="./result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]

result_list=list()
for(k in 1:length(file_list)){
  file = file_list[k]
  result = readRDS(paste0("result files/", file))
  d_fdr = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                     'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                     'control_methods' = result[[1]][[1]][[CRT.method]]$fdr_methods,
                     'CRT_methods' = CRT.method,
                     'error'  = result[[1]][[1]][[CRT.method]]$fdr,
                     'power' = result[[1]][[1]][[CRT.method]]$power, 
                     'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                     'error_type' = 'fdr')
  d_fwer = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                      'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                      'control_methods' = result[[1]][[2]][[CRT.method]]$fwer_methods,
                      'CRT_methods' = CRT.method,
                      'error'  = result[[1]][[2]][[CRT.method]]$fwer,
                      'power' = result[[1]][[2]][[CRT.method]]$power, 
                      'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                      'error_type' = 'fwer')
  cat(paste0(k, " files read out of ", length(file_list)) ,"\r")
  result_list[[k]]= rbind(d_fdr,d_fwer)
}
df = do.call(rbind, result_list)
df_backup2 = df

df[df$control_methods=="univariate-BH", "control_methods"] = "Univariate-BH"
df[df$control_methods=="univariate-Holm", "control_methods"] = "Univariate-Holm"
df[df$control_methods=="Simes_oracle_11", "control_methods"] = "BCP(s+1)-BH-estimated"
df[df$control_methods=="Simes_oracle_50", "control_methods"] = "BCP(p/2)-BH-estimated"
df[df$control_methods=="Simes_oracle_99", "control_methods"] = "BCP(p-1)-BH-estimated"
df[df$control_methods=="Simes_holm_oracle_11", "control_methods"] = "BCP(s+1)-Holm-estimated"
df[df$control_methods=="Simes_holm_oracle_50", "control_methods"] = "BCP(p/2)-Holm-estimated"
df[df$control_methods=="adaptive-Holm", "control_methods"] = "BCP(p-1)-Holm-estimated"
df[df$control_methods=="adafilter_oracle_11", "control_methods"] = "adafilter-BH-10"
df[df$control_methods=="adafilter_oracle_50", "control_methods"] = "adafilter-BH-50"
df[df$control_methods=="adafilter_oracle_99", "control_methods"] = "adafilter-BH-99"
df[df$control_methods=="adafilter_bonferroni_oracle_11", "control_methods"] = "adafilter-bonferroni-10  "
df[df$control_methods=="adafilter_bonferroni_oracle_50", "control_methods"] = "adafilter-bonferroni-50  "

df_estimated = df_estimated_backup =df

####################################################################################################
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run103/") # for the non estimated version

file_list = list.files(path="result files/")
file_list <- file_list[grep("^Result-model=dirichlet", file_list)]
extract_seed <- function(s) {
  matches <- regmatches(s, regexec("seed=(\\d+)", s))
  as.numeric(matches[[1]][2]) # Extract the numeric part of the match
}
seed_values <- sapply(file_list, extract_seed)
file_list <- file_list[seed_values <= target_simulation.iter/simulation.iter]

result_list=list()
for(k in 1:length(file_list)){
  file = file_list[k]
  result = readRDS(paste0("result files/", file))
  result_fields = c('model', 'snr','control_methods','error','power','m_1','error_type')
  d_fdr = d_fwer = setNames(data.frame(matrix(ncol = length(result_fields), nrow = 0)), result_fields)
  for(i in 1:simulation.iter){
    d_fdr_temp = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                            'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                            'control_methods' = result[[i]][[1]][[CRT.method]]$fdr_methods,
                            'CRT_methods' = CRT.method,
                            'error'  = result[[i]][[1]][[CRT.method]]$fdr,
                            'power' = result[[i]][[1]][[CRT.method]]$power, 
                            'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                            'error_type' = 'fdr')
    d_fwer_temp = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                             'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                             'control_methods' = result[[i]][[2]][[CRT.method]]$fwer_methods,
                             'CRT_methods' = CRT.method,
                             'error'  = result[[i]][[2]][[CRT.method]]$fwer,
                             'power' = result[[i]][[2]][[CRT.method]]$power, 
                             'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                             'error_type' = 'fwer')
    d_fdr <- rbind(d_fdr, d_fdr_temp)
    d_fwer <- rbind(d_fwer, d_fwer_temp)
  }
  cat(paste0(k, " files read out of ", length(file_list)) ,"\r")
  result_list[[k]]= rbind(d_fdr,d_fwer)
}

df = do.call(rbind, result_list)
df_backup2 = df

df[df$control_methods=="univariate-BH", "control_methods"] = "Univariate-BH"
df[df$control_methods=="Simes_oracle_11", "control_methods"] = "BCP(s+1)-BH"
df[df$control_methods=="Simes_oracle_50", "control_methods"] = "BCP(p/2)-BH"
df[df$control_methods=="max-BH", "control_methods"] = "BCP(p-1)-BH"
df[df$control_methods=="adaptive_Simes_holm_oracle_11", "control_methods"] = "BCP(s+1)-Holm"
df[df$control_methods=="adaptive_Simes_holm_oracle_50", "control_methods"] = "BCP(p/2)-Holm"
df[df$control_methods=="adaptive-Holm", "control_methods"] = "BCP(p-1)-Holm"
df[df$control_methods=="adafilter_oracle_11", "control_methods"] = "adafilter-BH-10"
df[df$control_methods=="adafilter_oracle_50", "control_methods"] = "adafilter-BH-50"
df[df$control_methods=="adafilter_bonferroni_oracle_11", "control_methods"] = "adafilter-bonferroni-10  "
df[df$control_methods=="adafilter_bonferroni_oracle_50", "control_methods"] = "adafilter-bonferroni-50  "

df_non_estimated = df_non_estimated_backup = df

##########################################################################################################
df_estimated = df_estimated_backup
df_non_estimated = df_non_estimated_backup

fdr_method_shortlist = c("BCP(p-1)-BH", "BCP(p/2)-BH", "BCP(s+1)-BH")
fwer_method_shortlist = c("BCP(p-1)-Holm", "BCP(p/2)-Holm", "BCP(s+1)-Holm")
method_shortlist = c(fdr_method_shortlist, fwer_method_shortlist)
method_estimated_shortlist = sapply(method_shortlist, function(x) paste0(x,"-estimated"))


df_non_estimated <- df_non_estimated %>%
  rename(Method = control_methods) %>%
  filter(Method %in% method_shortlist) %>%
  mutate(Method = factor(Method, levels = method_shortlist))

df_estimated <- df_estimated %>%
  rename(Method = control_methods) %>%
  filter(Method %in% method_estimated_shortlist) %>%
  mutate(Method = factor(Method, levels = method_estimated_shortlist))

##########################################################################################################
df = rbind(df_estimated, df_non_estimated)
df_backup = df
cbPalette <- rep(c("#47B9CD", "#FA9400","#FE0000"),times = 2)
fdr_method_shortlist = c(fdr_method_shortlist, sapply(fdr_method_shortlist, function(x) paste0(x, "-estimated")))
fwer_method_shortlist = c(fwer_method_shortlist, sapply(fwer_method_shortlist, function(x) paste0(x, "-estimated")))

model_name = "dirichlet_estimated"
snr_array = snr_list[[1]]

convert_to_expression <- function(label) {
  # Replace 's+1' with 's+1' etc.
  label <- gsub("p_([0-9]+)", "p[\\1]", label)
  parse(text=label)
}

my_plot <- function(name = "ERROR", model_name = "dirichlet", error_type_name = 'fdr'){
  df = subset(df, snr %in% snr_array & error_type %in% error_type_name)
  df = subset(df, CRT_methods %in% CRT.method_array)
  
  df.summary <- subset(df, CRT_methods == "d_CRT_with_r_squared" & error_type %in% error_type_name) %>%
    group_by(snr, Method) %>%
    summarise(
      upper_error = min(mean(error) + 1.96 * sd(error)/sqrt(simulation.iter), rep(1,length(snr))),
      lower_error = max(mean(error) - 1.96 * sd(error)/sqrt(simulation.iter), rep(0,length(snr))),
      error = mean(error),
      upper_power = min(mean(power) + 1.96 * sd(power)/sqrt(simulation.iter), rep(1,length(snr))),
      lower_power = max(mean(power) - 1.96 * sd(power)/sqrt(simulation.iter), rep(0,length(snr))),
      power = mean(power)
    )
  df.summary = subset(df.summary, snr %in% snr_array)
  if(error_type_name == 'fdr') df.summary = subset(df.summary, Method %in% fdr_method_shortlist)
  if(error_type_name == 'fwer') df.summary = subset(df.summary, Method %in% fwer_method_shortlist)
  
  n_methods = ifelse(error_type_name == 'fdr', length(fdr_method_shortlist), length(fwer_method_shortlist))
  
  p = ggplot(df.summary, aes_string("snr", tolower(name))) +
    geom_line(aes(group = Method, color = Method, linetype = Method), data = df.summary) +
    geom_errorbar(aes(ymin = get(paste0("lower_", tolower(name))), ymax = get(paste0("upper_", tolower(name))), 
                      color = Method), width = 0, show.legend = FALSE) +
    theme_classic() + 
    labs(x = "SNR", y = ifelse(name == 'Power', 'Power', toupper(error_type_name)), 
         color = "Methods ", linetype = "Methods ")+
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text.align = 0)
  if(name == 'ERROR') p=p+coord_cartesian(ylim = c(0, 6*alpha)) + geom_hline(aes(yintercept = 0.1), linetype = "dashed") 
  if(name == 'Power') p=p+coord_cartesian(ylim = c(0, 1))
  
  unique_methods <- as.character(unique(df.summary$Method))
  labels_expr <- sapply(unique_methods, convert_to_expression)
  p = p + scale_color_manual(values = cbPalette, labels = labels_expr) +
    scale_linetype_manual(values = rep(c("solid","dashed"),each=3),
                          labels = labels_expr)
  return(p)
  
}

setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run94/")

fdr_plot_1 = my_plot(name = "ERROR", model_name = "dirichlet_estimated", error_type_name = 'fdr')
fdr_plot_2 = my_plot(name = "Power", model_name = "dirichlet_estimated", error_type_name = 'fdr')
fdr_plot = ggarrange(fdr_plot_1, fdr_plot_2, common.legend = TRUE, legend = "right")
fdr_plot
ggsave(paste0("plots/presentation/","fdr_dirichlet_estimated.png"), 
       plot=fdr_plot, device = "png", width = 7, height = 3, units = "in", bg="white",dpi = 600)
ggsave("plots/presentation/fdr_dirichlet_estimated.pdf", plot = fdr_plot, device = "pdf", width = 7, height = 3, units = "in")

df = df_backup
fwer_plot_1 = my_plot(name = "ERROR", model_name = "dirichlet_estimated", error_type_name = 'fwer')
fwer_plot_2 = my_plot(name = "Power", model_name = "dirichlet_estimated", error_type_name = 'fwer')
fwer_plot = ggarrange(fwer_plot_1, fwer_plot_2, common.legend = TRUE, legend = "right")
fwer_plot 
ggsave(paste0("plots/presentation/fwer_dirichlet_estimated.png"), 
       plot=fwer_plot, device = "png", width = 7, height = 3, units = "in", bg="white",dpi = 600)
ggsave("plots/presentation/fwer_dirichlet_estimated.pdf", plot = fwer_plot, device = "pdf", width = 7, height = 3, units = "in")

