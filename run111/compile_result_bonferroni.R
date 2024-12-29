remove(list = ls())
library(ggplot2); library(stringr)
library(ggpubr)
library(dplyr)
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run111/")

##### Parameters ################################### 

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
alpha = 0.05 # target fdr level
CRT.replications = 50000 # number of knockoff variables generated
simulation.iter = 5 # number of replications of the simulations
target_simulation.iter = 500
p_filter_tolerance = 0.9 # 0.9 means if 90% of p vals of a variable is rejected only then it will reject that variable

CRT.method_array = c("d_CRT_with_r_squared")
CRT.method = CRT.method_array[1]
model_array = c('dirichlet', 'normal') # model for generating X
snr_list = list(c(0.01,0.02,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5))

CI_coverage = 90
#################################################### 

snr_array_union = sort(unique(unlist(snr_list)))
result_file_name = function(snr_scale_factor, model,seed){
  paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
         ",alpha=", alpha, ",CRT.replications=", CRT.replications,
         ",simulation.iter=", simulation.iter,
         ",snr=", snr_scale_factor, ",seed=",seed)
}
df <- data.frame(matrix(NA, nrow = 2*length(model_array)*length(snr_array_union)*simulation.iter*length(CRT.method_array), ncol = 7))

file_list = list.files(path="./result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]

result_list = list()
for(k in 1:length(file_list)){
  file = file_list[k]
  result = readRDS(paste0("result files/", file))
  result_fields = c('model', 'snr','PCH_methods','rejection_non_null','rejection_null','m_1','error_type')
  d = setNames(data.frame(matrix(ncol = length(result_fields), nrow = 0)), result_fields)
  for(i in 1:simulation.iter){
    d_temp = data.frame('model' = str_match(file, "model=\\s*(.*?)\\s*,")[2], 
                        'snr' = as.numeric(str_match(file, "snr=\\s*(.*?)\\s*,")[2]),
                        'PCH_methods' = result[[i]][[1]]$pch_methods,
                        'CRT_methods' = CRT.method,
                        'rejection_non_null'  = result[[i]][[1]]$rejection_non_null,
                        'rejection_null'  = result[[i]][[1]]$rejection_null,
                        'm_1' = str_match(file[1], "m_1=\\s*(.*?)\\s*,")[2])
    d <- rbind(d, d_temp)
  }
  cat(paste0("k=", k, "iterations done out of ", length(file_list)), "\r")
  result_list[[k]] = rbind(d)
}
df = do.call(rbind, result_list)
names(df) = names(d)
df_backup = df

df=df_backup
df <- df %>%
  mutate(
    PCH_function = case_when(
      str_detect(PCH_methods, "Simes") ~ "Simes",
      str_detect(PCH_methods, "bonferroni") ~ "Bonferroni",
      str_detect(PCH_methods, "LOO") ~ "LOO",
      TRUE ~ "Other"
    ),
    Method = case_when(
      str_detect(PCH_methods, "Simes") & as.numeric(str_extract(PCH_methods, "\\d+$")) == "99" ~ "BCP(p-1)",
      str_detect(PCH_methods, "Simes") & str_extract(PCH_methods, "\\d+$") == "50" ~ "BCP(p/2)-Simes",
      str_detect(PCH_methods, "Simes") & str_extract(PCH_methods, "\\d+$") == "11" ~ "BCP(s+1)-Simes",
      str_detect(PCH_methods, "bonferroni") & str_extract(PCH_methods, "\\d+$") == "50" ~ "BCP(p/2)-Bonferroni",
      str_detect(PCH_methods, "bonferroni") & str_extract(PCH_methods, "\\d+$") == "11" ~ "BCP(s+1)-Bonferroni",
      # str_detect(PCH_methods, "LOO") ~ "LOO",
      TRUE ~ "Other"
    )
  ) %>%
  filter(PCH_function != "Other", Method != "Other")

df_summary <- df %>%
  mutate(
    # Convert logical to numeric for calculations
    rejection_non_null = as.numeric(rejection_non_null),
    rejection_null = as.numeric(rejection_null)
  ) %>%
  group_by(model, snr, CRT_methods, PCH_function, Method) %>%
  summarize(
    power = mean(rejection_non_null, na.rm = TRUE),
    power_se = sqrt(power * (1 - power) / n()),  # Calculating standard error for power
    type_1_error = mean(rejection_null, na.rm = TRUE),
    type_1_error_se = sqrt(type_1_error * (1 - type_1_error) / n()),  # Standard error for type-1-error
    .groups = 'drop'
  )
df_summary$Method = factor(df_summary$Method, levels = c("BCP(p-1)", 
                                                         "BCP(p/2)-Simes", "BCP(s+1)-Simes", 
                                                         "BCP(p/2)-Bonferroni", "BCP(s+1)-Bonferroni"))

cbPalette <- c("#8A2BE2", "#FA9400", "#FE0000", "#4767CD", "#47B9CD")
snr_array_dirichlet = c(0.01,0.02,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5)

convert_to_expression <- function(label) {
  # Replace 'p_1+1' with 'p[1]+1' etc.
  label <- gsub("p_([0-9]+)", "p[\\1]", label)
  parse(text=label)
}
my_plot <- function(data, model_name = "normal", show_error_bars = FALSE, name = "power"){
  data = subset(data, model == model_name & CRT_methods == "d_CRT_with_r_squared" & snr %in% get(paste0("snr_array_", model_name)))
  
  p = ggplot(data, aes_string("snr", name)) +
    geom_line(aes(color = Method, linetype = Method)) +
    theme_classic() + 
    labs(x = "Coefficient", y = tools::toTitleCase(gsub("_", " ", name)))+
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text.align = 0)
  
  if(show_error_bars==TRUE) p = p+geom_errorbar(aes(ymin = get(name)-2*get(paste0(name,"_se")), ymax = get(name)+2*get(paste0(name,"_se")), color = Method), width = 0, show.legend = FALSE)
  
  if(name == 'type_1_error') p=p+coord_cartesian(ylim = c(0, 6*alpha)) + geom_hline(aes(yintercept = alpha), linetype = "dashed") 
  if(name == 'power') p=p+coord_cartesian(ylim = c(0, 1))
  p = p + scale_color_manual(values = cbPalette, 
                             name = 'Methods',
                             labels = sapply(as.character(unique(df_summary$Method)), convert_to_expression)) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid", "twodash", "longdash"), 
                          name = 'Methods',
                          labels = sapply(as.character(unique(df_summary$Method)), convert_to_expression))
  return(p)
}

dirichlet_plot_1 = my_plot(df_summary, model_name = "dirichlet", show_error_bars = TRUE, name = "type_1_error")
dirichlet_plot_2 = my_plot(df_summary, model_name = "dirichlet", show_error_bars = TRUE, name = "power")
dirichlet_plot = ggarrange(dirichlet_plot_1, dirichlet_plot_2, common.legend = TRUE, legend = "right")
dirichlet_plot
ggsave("plots/dirichlet_single_test_plot.pdf", plot = dirichlet_plot, device = "pdf", width = 7, height = 3, units = "in")

