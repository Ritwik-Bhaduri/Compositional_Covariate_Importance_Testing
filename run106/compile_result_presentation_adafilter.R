remove(list = ls())
library(ggplot2); library(stringr)
library(ggpubr)
library(dplyr)
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run103/")

##### Parameters ################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
alpha = 0.1 # target fdr level
CRT.replications = 200 # number of knockoff variables generated
target_simulation.iter = 200*25
simulation.iter = 25 # number of replications of the simulations
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

CRT.method_array = c("d_CRT_with_r_squared"); CRT.method = CRT.method_array[1]
model_array = c('dirichlet', 'normal') # model for generating X
# # fdr_method_shortlist = c("max-BH","adafilter","adafilter_oracle_6","adafilter_oracle_50",
# #                          "DHH_storey_0.3", "DHH_storey_0.5", "DHH_storey_0.7",
# #                          "Simes_oracle_6", "Simes_oracle_7", "Simes_oracle_50", "Simes_oracle_99")
# fdr_method_shortlist = c("max-BH", "univariate-BH",
#                          "adafilter","adafilter_oracle_6",
#                          "Simes_oracle_6", "Simes_oracle_50")
# 
# # fwer_method_shortlist = c("max-Holm", "univariate-Holm", "adaptive-Holm", "adaptive-Holm-aggressive",
# #                           "Simes_holm_oracle_6", "Simes_holm_oracle_50", "Simes_holm_oracle_95", "Simes_holm_oracle_99",
# #                           "adafilter_bonferroni_oracle_6", "adafilter_bonferroni_oracle_50",
# #                           "adafilter_bonferroni_oracle_99")
# 
# fwer_method_shortlist = c("max-Holm", "adaptive-Holm",
#                           "Simes_holm_oracle_6", "Simes_holm_oracle_50",
#                           "adafilter_bonferroni_oracle_6")

# snr_union = c(0.1,0.5,1,2,5,10,25,50,75,100)
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
df <- data.frame(matrix(NA, nrow = 2*length(model_array)*length(snr_array_union)*simulation.iter*length(CRT.method_array), ncol = 7))

file_list = list.files(path="./result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]

result_list = list()
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
  cat(sprintf("k=%d iterations done out of %d\n", k, length(file_list)), file = "progress.txt", append = FALSE)
  result_list[[k]] = rbind(d_fdr,d_fwer)
}
df = do.call(rbind, result_list)
names(df) = names(d_fdr)
df_backup = df

df[df$control_methods=="max-BH", "control_methods"] = "Max-BH"
df[df$control_methods=="adaptive-Holm", "control_methods"] = "Max-Holm"
df[df$control_methods=="univariate-BH", "control_methods"] = "LOO-BH"
df[df$control_methods=="univariate-Holm", "control_methods"] = "LOO-Holm"
df[df$control_methods=="Simes_oracle_11", "control_methods"] = "CB-PCH-oracle"
df[df$control_methods=="Simes_oracle_50", "control_methods"] = "CB-PCH"
df[df$control_methods=="adaptive_Simes_holm_oracle_11", "control_methods"] = "CB-PCH-oracle"
df[df$control_methods=="adaptive_Simes_holm_oracle_50", "control_methods"] = "CB-PCH"
df[df$control_methods=="adafilter_oracle_11", "control_methods"] = "adafilter-BH-10"
df[df$control_methods=="adafilter_oracle_50", "control_methods"] = "adafilter-BH-50"
df[df$control_methods=="adafilter_bonferroni_oracle_11", "control_methods"] = "adafilter-bonferroni-10  "
df[df$control_methods=="adafilter_bonferroni_oracle_50", "control_methods"] = "adafilter-bonferroni-50  "

fdr_method_shortlist = c("adafilter-BH-10", "adafilter-BH-50", "CB-PCH-oracle", "CB-PCH")
fwer_method_shortlist = c("adafilter-bonferroni-10  ","adafilter-bonferroni-50  ", "CB-PCH-oracle", "CB-PCH") 
cbPalette_fdr <- c("#4767CD", "#47B9CD","#FA9400", "#FE0000")
cbPalette_fwer <- c("#4767CD", "#47B9CD","#FA9400", "#FE0000")

model_name = "dirichlet"
snr_array = snr_list[[1]]

my_plot <- function(name = "ERROR", model_name = "dirichlet", error_type_name = 'fdr'){
  df = subset(df, snr %in% snr_array &  model %in% model_name & error_type %in% error_type_name)
  df = subset(df, CRT_methods %in% CRT.method_array)
  
  df.summary <- subset(df, model == model_name & CRT_methods == "d_CRT_with_r_squared" & error_type %in% error_type_name) %>%
    group_by(snr, control_methods) %>%
    summarise(
      upper_error = min(mean(error) + 1.96 * sd(error)/sqrt(target_simulation.iter), rep(1,length(snr))),
      lower_error = max(mean(error) - 1.96 * sd(error)/sqrt(target_simulation.iter), rep(0,length(snr))),
      error = mean(error),
      upper_power = min(mean(power) + 1.96 * sd(power)/sqrt(target_simulation.iter), rep(1,length(snr))),
      lower_power = max(mean(power) - 1.96 * sd(power)/sqrt(target_simulation.iter), rep(0,length(snr))),
      power = mean(power)
    )
  df.summary = subset(df.summary, snr %in% snr_array)
  if(model_name == "normal" && error_type_name == 'fdr'){
    # View(df.summary)
  }
  
  if(error_type_name == 'fdr') df.summary = subset(df.summary, control_methods %in% fdr_method_shortlist)
  if(error_type_name == 'fwer') df.summary = subset(df.summary, control_methods %in% fwer_method_shortlist)
  
  n_methods = ifelse(error_type_name == 'fdr', length(fdr_method_shortlist), length(fwer_method_shortlist))
  
  p = ggplot(df.summary, aes_string("snr", tolower(name))) +
    geom_line(aes(group = control_methods, color = control_methods, linetype = control_methods), data = df.summary) +
    geom_errorbar(aes(ymin = get(paste0("lower_", tolower(name))), ymax = get(paste0("upper_", tolower(name))), 
                      color = control_methods), width = 0, show.legend = FALSE) +
    theme_classic() + 
    labs(x = "SNR", y = ifelse(name == 'Power', 'Power', toupper(error_type_name)), 
         title = ifelse(name == 'Power', 'Power', toupper(error_type_name)),
         color = "Methods ", linetype = "Methods ")+
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  if(name == 'ERROR') p=p+coord_cartesian(ylim = c(0, 6*alpha)) + geom_hline(aes(yintercept = 0.1), linetype = "dashed") 
  if(name == 'Power') p=p+coord_cartesian(ylim = c(0, 1))
  if(error_type_name == 'fdr') p = p + scale_color_manual(values = cbPalette_fdr)
  if(error_type_name == 'fwer') p = p + scale_color_manual(values = cbPalette_fwer)
  
  return(p)
  
}

fdr_plot_1 = my_plot(name = "ERROR", model_name = "dirichlet", error_type_name = 'fdr')
fdr_plot_2 = my_plot(name = "Power", model_name = "dirichlet", error_type_name = 'fdr')
fdr_plot = ggarrange(fdr_plot_1, fdr_plot_2, common.legend = TRUE, legend = "right")
fdr_plot
ggsave("plots/with adafilter/fdr_dirichlet.pdf", plot = fdr_plot, device = "pdf", width = 7, height = 3, units = "in")


fwer_plot_1 = my_plot(name = "ERROR", model_name = "dirichlet", error_type_name = 'fwer')
fwer_plot_2 = my_plot(name = "Power", model_name = "dirichlet", error_type_name = 'fwer')
fwer_plot = ggarrange(fwer_plot_1, fwer_plot_2, common.legend = TRUE, legend = "right")
fwer_plot 
ggsave("plots/with adafilter/fwer_dirichlet.pdf", plot = fwer_plot, device = "pdf", width = 7, height = 3, units = "in")


fdr_plot_1 = my_plot(name = "ERROR", model_name = "normal", error_type_name = 'fdr')
fdr_plot_2 = my_plot(name = "Power", model_name = "normal", error_type_name = 'fdr')
fdr_plot = ggarrange(fdr_plot_1, fdr_plot_2, common.legend = TRUE, legend = "right")
fdr_plot
ggsave("plots/with adafilter/fdr_normal.pdf", plot = fdr_plot, device = "pdf", width = 7, height = 3, units = "in")


fwer_plot_1 = my_plot(name = "ERROR", model_name = "normal", error_type_name = 'fwer')
fwer_plot_2 = my_plot(name = "Power", model_name = "normal", error_type_name = 'fwer')
fwer_plot = ggarrange(fwer_plot_1, fwer_plot_2, common.legend = TRUE, legend = "right")
fwer_plot 
ggsave("plots/with adafilter/fwer_normal.pdf", plot = fwer_plot, device = "pdf", width = 7, height = 3, units = "in")

