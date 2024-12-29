remove(list = ls())
library(ggplot2); library(stringr)
library(ggpubr)
library(dplyr)
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run114/")

##### Parameters ################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
alpha = 0.1 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
simulation.iter = 500 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

CRT.method_array = c("d_CRT_with_r_squared")
model_array = c('dirichlet_multinomial') # model for generating X

snr_val = 1
sparsity_cutoff_m_values = c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)

CI_coverage = 90
#################################################### 

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

file_list = list.files(path="./result files from pval/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]

df <- foreach(k = 1:length(file_list), .combine = 'rbind') %dopar%{
  file = file_list[k]
  result = readRDS(paste0("./result files from pval/", file))
  for(i in 1:length(sparsity_cutoff_m_values)){
    sparsity_cutoff_m = sparsity_cutoff_m_values[i]
    cat(paste(i/length(file_list), "\n", "\r"), file = "progress.txt")
    d_fdr_temp = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2],
                            'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                            'sparsity_cutoff_m' = sparsity_cutoff_m,
                            'control_methods' = result[[1]][[1]][[1]][[i]]$fdr_methods,
                            'CRT_methods' = CRT.method,
                            'error'  = result[[1]][[1]][[1]][[i]]$fdr,
                            'power' = result[[1]][[1]][[1]][[i]]$power,
                            'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                            'error_type' = 'fdr')
    d_fwer_temp = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2],
                             'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                             'sparsity_cutoff_m' = sparsity_cutoff_m,
                             'control_methods' = result[[1]][[1]][[2]][[i]]$fwer_methods,
                             'CRT_methods' = CRT.method,
                             'error'  = result[[1]][[1]][[2]][[i]]$fwer,
                             'power' = result[[1]][[1]][[2]][[i]]$power,
                             'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2],
                             'error_type' = 'fwer')
    if(i==1){
      d_fdr  = d_fdr_temp
      d_fwer = d_fwer_temp
    }else{
      d_fdr = rbind(d_fdr, d_fdr_temp)
      d_fwer = rbind(d_fwer, d_fwer_temp)
    }
  }
  rbind(d_fdr,d_fwer)
}
parallel::stopCluster(my.cluster)
df_backup = df
dim(df)
unique(df$control_methods)
df=df_backup


# View(df_backup)
df <- df %>%
  mutate(
    Method = case_when(
      str_detect(control_methods, "max-BH") & error_type=="fdr" ~ "BCP(p-1)-BH",
      str_detect(control_methods, "Simes_oracle_25") & error_type=="fdr" ~ "BCP(p/4)-Simes-BH",
      str_detect(control_methods, "Simes_oracle_11") & error_type=="fdr" ~ "BCP(s+1)-Simes-BH",
      
      str_detect(control_methods, "bonferroni_oracle_25") & error_type=="fdr" ~ "BCP(p/4)-Bonferroni-BH",
      str_detect(control_methods, "bonferroni_oracle_11") & error_type=="fdr" ~ "BCP(s+1)-Bonferroni-BH",
      
      str_detect(control_methods, "adaptive-Holm") & error_type=="fwer" ~ "BCP(p-1)-Holm",
      str_detect(control_methods, "adaptive_Simes_holm_oracle_25") & error_type=="fwer" ~ "BCP(p/4)-Simes-Holm",
      str_detect(control_methods, "adaptive_Simes_holm_oracle_11") & error_type=="fwer" ~ "BCP(s+1)-Simes-Holm",
      
      str_detect(control_methods, "adaptive_bonferroni_holm_oracle_25") & error_type=="fwer" ~ "BCP(p/4)-Bonferroni-Holm",
      str_detect(control_methods, "adaptive_bonferroni_holm_oracle_11") & error_type=="fwer" ~ "BCP(s+1)-Bonferroni-Holm",
      
      TRUE ~ "Other"
    )
  ) %>%
  filter(Method %in% c("BCP(p-1)-BH","BCP(p/4)-Simes-BH","BCP(s+1)-Simes-BH",
                       "BCP(p/4)-Bonferroni-BH","BCP(s+1)-Bonferroni-BH",
                       "BCP(p-1)-Holm","BCP(p/4)-Simes-Holm","BCP(s+1)-Simes-Holm",
                       "BCP(p/4)-Bonferroni-Holm","BCP(s+1)-Bonferroni-Holm"))
df$control_methods = factor(df$Method, levels = c("BCP(p-1)-BH","BCP(p/4)-Simes-BH","BCP(s+1)-Simes-BH",
                            "BCP(p/4)-Bonferroni-BH","BCP(s+1)-Bonferroni-BH",
                            "BCP(p-1)-Holm","BCP(p/4)-Simes-Holm","BCP(s+1)-Simes-Holm",
                            "BCP(p/4)-Bonferroni-Holm","BCP(s+1)-Bonferroni-Holm"))
df$Method = df$control_methods

convert_to_expression <- function(label) {
  if (label == "BCP(p-1)-BH") label = expression("BCP(|" * italic(D) * "|-1)-BH")
  if (label == "BCP(p/4)-Simes-BH") label = expression("BCP(p/4)-Simes-BH")
  if (label == "BCP(s+1)-Simes-BH") label = expression("BCP(s+1)-Simes-BH")
  if (label == "BCP(p/4)-Bonferroni-BH") label = expression("BCP(p/4)-Bonferroni-BH")
  if (label == "BCP(s+1)-Bonferroni-BH") label = expression("BCP(s+1)-Bonferroni-BH")
  if (label == "BCP(p-1)-Holm") label = expression("BCP(|" * italic(D) * "|-1)-Holm")
  if (label == "BCP(p/4)-Simes-Holm") label = expression("BCP(p/4)-Simes-Holm")
  if (label == "BCP(s+1)-Simes-Holm") label = expression("BCP(s+1)-Simes-Holm")
  if (label == "BCP(p/4)-Bonferroni-Holm") label = expression("BCP(p/4)-Bonferroni-Holm")
  if (label == "BCP(s+1)-Bonferroni-Holm") label = expression("BCP(s+1)-Bonferroni-Holm")
  return(label)
}


cbPalette <- c("#8A2BE2", "#FA9400", "#FE0000", "#4767CD", "#47B9CD")
df$cov_conditioned_on = m - df$sparsity_cutoff_m

fdr_method_shortlist <- levels(df$control_methods)[grepl("BH", levels(df$control_methods))]
fwer_method_shortlist <- levels(df$control_methods)[grepl("Holm", levels(df$control_methods))]

my_plot <- function(name = "ERROR", model_name = "dirichlet", error_type_name = 'fdr'){
  df = subset(df, sparsity_cutoff_m %in% sparsity_cutoff_m_values & model %in% model_name & error_type %in% error_type_name & snr == snr_val)
  df = subset(df, CRT_methods %in% CRT.method_array)
  df.summary <- subset(df, model == model_name & CRT_methods == "d_CRT_with_r_squared" & error_type %in% error_type_name) %>%
    group_by(cov_conditioned_on, Method) %>%
    summarise(
      upper_error = min(mean(error) + 1.96 * sd(error)/sqrt(simulation.iter), rep(1,length(sparsity_cutoff_m))),
      lower_error = max(mean(error) - 1.96 * sd(error)/sqrt(simulation.iter), rep(0,length(sparsity_cutoff_m))),
      error = mean(error),
      upper_power = min(mean(power) + 1.96 * sd(power)/sqrt(simulation.iter), rep(1,length(sparsity_cutoff_m))),
      lower_power = max(mean(power) - 1.96 * sd(power)/sqrt(simulation.iter), rep(0,length(sparsity_cutoff_m))),
      power = mean(power)
    )
  df.summary = subset(df.summary, cov_conditioned_on %in% (m-sparsity_cutoff_m_values))
  
  if(error_type_name == 'fdr') df.summary = subset(df.summary, Method %in% fdr_method_shortlist)
  if(error_type_name == 'fwer') df.summary = subset(df.summary, Method %in% fwer_method_shortlist)
  
  n_methods = ifelse(error_type_name == 'fdr', length(fdr_method_shortlist), length(fwer_method_shortlist))
  
  p = ggplot(df.summary, aes_string("cov_conditioned_on", tolower(name))) +
    geom_line(aes(group = Method, color = Method, linetype = Method), data = df.summary) +
    geom_errorbar(aes(ymin = get(paste0("lower_", tolower(name))), ymax = get(paste0("upper_", tolower(name))), 
                      color = Method), width = 0, show.legend = FALSE) +
    theme_classic() + 
    labs(x = expression("# Covariates Conditioned on (|" * italic(D)^c * "|)"),
         y = ifelse(name == 'Power', 'Average Power', toupper(error_type_name)),
         color = "Methods ", linetype = "Methods ")+
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text.align = 0) 
  
  if(name == 'ERROR'){
    p=p+coord_cartesian(ylim = c(0, 4*alpha)) + geom_hline(aes(yintercept = 0.1), linetype = "dashed") 
  } 
  if(name == 'Power'){
    p=p+coord_cartesian(ylim = c(0, 0.4))
  } 
  unique_methods <- as.character(unique(df.summary$Method))
  labels_expr <- sapply(unique_methods, convert_to_expression)
  
  p = p + scale_color_manual(name = 'Methods', values = cbPalette,
                             labels = labels_expr)+
    scale_linetype_manual(name = 'Methods', values = c("dotted", "dashed", "longdash", "twodash", "solid"),
                          labels = labels_expr)
  
  p
}
power_plot_fdr = my_plot(name = "Power", model_name = 'dirichlet_multinomial', error_type_name = 'fdr')
power_plot_fwer = my_plot(name = "Power", model_name = 'dirichlet_multinomial', error_type_name = 'fwer')
power_plot = ggarrange(power_plot_fdr, power_plot_fwer, common.legend = TRUE, legend = "right")
power_plot
ggsave(paste0("plots/power_dirichlet_multinomial.pdf"),
       plot=power_plot, device = "pdf", width = 7, height = 3, units = "in")

error_plot_fdr = my_plot(name = "ERROR", model_name = 'dirichlet_multinomial', error_type_name = 'fdr')
error_plot_fwer = my_plot(name = "ERROR", model_name = 'dirichlet_multinomial', error_type_name = 'fwer')
error_plot = ggarrange(error_plot_fdr, error_plot_fwer, common.legend = TRUE, legend = "right")
error_plot
ggsave(paste0("plots/error_dirichlet_multinomial.pdf"),
       plot=error_plot, device = "pdf", width = 7, height = 3, units = "in")

saveRDS(power_plot_fdr, "plots/power_plot_fdr.rds")
saveRDS(power_plot_fwer, "plots/power_plot_fwer.rds")
saveRDS(error_plot_fdr, "plots/error_plot_fdr.rds")
saveRDS(error_plot_fwer, "plots/error_plot_fwer.rds")


