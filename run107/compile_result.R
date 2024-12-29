remove(list = ls())
library(ggplot2); library(stringr)
library(ggpubr)
library(dplyr)
library(data.table)
library(latex2exp)
library(extrafont)
setwd("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run107")

##### Parameters ###################################

n = 100 # number of rows in data
m = 100 # number of columns in data
m_1 = 10 # proportion of null variables in data
alpha = 0.05 # target fdr level
CRT.replications = 1500 # number of knockoff variables generated
simulation.iter = 10 # number of replications of the simulations
target_simulation.iter = 500

CRT_method_array = c("d_CRT_with_r_squared")
CRT_method = CRT_method_array[1]
model_array = c('dirichlet_multinomial') # model for generating X

snr = 0.5
sparsity_cutoff_m_values = c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)

CI_coverage = 90
####################################################

result_file_name = function(snr_scale_factor, model,seed){
  paste0("result files/Result-", "model=", model, ",n=", n, ",m=", m, ",m_1=", m_1,
         ",alpha=", alpha, ",CRT.replications=", CRT.replications,
         ",simulation.iter=", simulation.iter,
         ",snr=", snr_scale_factor,
         ",p_filter_tolerance=", p_filter_tolerance,",seed=",seed)
}
df <- data.frame(matrix(NA, nrow = 2*length(model_array)*length(sparsity_cutoff_m_values)*simulation.iter*length(CRT_method_array), ncol = 7))

file_list = list.files(path="./result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("Result-", file_list)]
file_list <- file_list[grep(paste0("snr=", snr), file_list)]
file_list <- file_list[grep(paste0("simulation.iter=", simulation.iter,","), file_list)]
file_list <- file_list[grep(paste0("alpha=", alpha), file_list)]

start_time = Sys.time()
result_list = list()
for (k in seq_along(file_list)) {
  file <- file_list[k]
  result <- readRDS(paste0("result files/", file))
  result_fields = c('model', 'sparsity_cutoff_m_value','PCH_methods','rejection_non_null','rejection_null','m_1','CRT_methods')
  num_methods <- length(result[[1]][[1]][[CRT_method]]$pch_methods) # assuming number of pch methods remain same accross simulations
  d = setNames(data.frame(matrix(ncol = length(result_fields), nrow = simulation.iter*length(sparsity_cutoff_m_values)*num_methods)), result_fields)
  
  # assuming most fields remain same accross simulations
  d$model = str_match(file, "model=\\s*(.*?)\\s*,")[2]
  d$sparsity_cutoff_m_value = rep(rep(sparsity_cutoff_m_values, each = num_methods), times = simulation.iter)
  d$PCH_methods = rep(result[[1]][[1]][[CRT_method]]$pch_methods, times = simulation.iter*length(sparsity_cutoff_m_values))
  d$m_1 = str_match(file, "m_1=\\s*(.*?)\\s*,")[2]
  d$CRT_methods = CRT_method
  
  idx1 <- 1
  for (j in 1:simulation.iter) {
    for (i in 1:length(sparsity_cutoff_m_values)) {
      temp_data <- result[[j]][[i]][[CRT_method]]
      idx2 <- idx1 + num_methods - 1
      
      # d[idx1:idx2, "model"] = rep(str_match(file, "model=\\s*(.*?)\\s*,")[2], num_methods)
      # d[idx1:idx2, "sparsity_cutoff_m_value"] = rep(temp_data$sparsity_cutoff_m_value, num_methods)
      # d[idx1:idx2, "PCH_methods"] = temp_data$pch_methods
      # d[idx1:idx2, "CRT_methods"] = rep(CRT_method, num_methods)
      d[idx1:idx2, "rejection_non_null"] = temp_data$rejection_non_null
      d[idx1:idx2, "rejection_null"] = temp_data$rejection_null
      # d[idx1:idx2, "m_1"] = rep(str_match(file, "m_1=\\s*(.*?)\\s*,")[2], num_methods)
      
      idx1 <- idx2 + 1
    }
  }
  cat(paste0("k=", k, " iterations done out of ", length(file_list), "\r"))
  result_list[[k]] <- d
}

df = do.call(rbind, result_list)
print(Sys.time()-start_time)

names(df) = c('model', 'sparsity_cutoff_m_value','PCH_methods','rejection_non_null','rejection_null','m_1','CRT_methods')
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
      str_detect(PCH_methods, "Simes") & as.numeric(str_extract(PCH_methods, "\\d+$")) > 25 ~ "BCP(p-1)",
      str_detect(PCH_methods, "Simes") & str_extract(PCH_methods, "\\d+$") == "25" ~ "BCP(p/4)",
      str_detect(PCH_methods, "Simes") & str_extract(PCH_methods, "\\d+$") == "11" ~ "BCP(s+1)",
      str_detect(PCH_methods, "LOO") ~ "LOO",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Method %in% c("BCP(p-1)","BCP(p/4)","BCP(s+1)"))
df$Method = factor(df$Method, levels = c("BCP(p-1)","BCP(p/4)","BCP(s+1)"))

df_summary <- df %>%
  mutate(
    # Convert logical to numeric for calculations
    rejection_non_null = as.numeric(rejection_non_null),
    rejection_null = as.numeric(rejection_null)
  ) %>%
  group_by(model, sparsity_cutoff_m_value, CRT_methods, PCH_function, Method) %>%
  summarize(
    power = mean(rejection_non_null, na.rm = TRUE),
    power_se = sqrt(power * (1 - power) / n()),  # Calculating standard error for power
    type_1_error = mean(rejection_null, na.rm = TRUE),
    type_1_error_se = sqrt(type_1_error * (1 - type_1_error) / n()),  # Standard error for type-1-error
    .groups = 'drop'
  )

# cbPalette <- c("#FE0000", "#4767CD","#47B9CD","#FA9400")
cbPalette <- c("#FE0000", "#FA9400","#47B9CD")
df_summary$cov_conditioned_on = m - df_summary$sparsity_cutoff_m_value


font_import()  # Only run this once to import all system fonts
loadfonts(device = "pdf")  # Load fonts for PDF output
View(df_summary)
my_plot <- function(data, model_name = "dirichlet_multinomial", show_error_bars = FALSE, name = "power"){
  data = subset(data, model == model_name & CRT_methods == "d_CRT_with_r_squared")
  data$Method <- factor(data$Method, levels = c("BCP(p-1)", "BCP(p/4)", "BCP(s+1)"))
  
  p = ggplot(data, aes_string("cov_conditioned_on", name)) +
    geom_line(aes(color = Method, linetype = Method)) +
    theme_classic() + 
    # labs(x = expression("|" * italic(D) * "|"), y = tools::toTitleCase(gsub("_", " ", name)))+
    labs(x = expression("# Covariates Conditioned on (|" * italic(D)^c * "|)"), 
         y = tools::toTitleCase(gsub("1", "I", gsub("_", " ", name))))+
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text.align = 0)
  
  if(show_error_bars==TRUE) p = p+geom_errorbar(aes(ymin = get(name)-2*get(paste0(name,"_se")), ymax = get(name)+2*get(paste0(name,"_se")), color = Method), width = 0, show.legend = FALSE)
  
  if(name == 'type_1_error') p=p+coord_cartesian(ylim = c(0, 0.4)) + geom_hline(aes(yintercept = 0.05), linetype = "dashed") 
  if(name == 'power') p=p+coord_cartesian(ylim = c(0, 1))
  p = p + scale_color_manual(values = cbPalette, 
                             name = 'Methods',
                             labels = expression("BCP(|" * italic(D) * "|-1)", "BCP(p/4)", "BCP(s + 1)", "LOO")) +
    scale_linetype_manual(values = c("dotted","dashed","solid"), 
                          name = 'Methods',
                          labels = expression("BCP(|" * italic(D) * "|-1)", "BCP(p/4)", "BCP(s + 1)", "LOO"))
  return(p)
}

dirichlet_multinomial_plot_1 = my_plot(df_summary, model_name = "dirichlet_multinomial", show_error_bars = TRUE, name = "type_1_error")
dirichlet_multinomial_plot_2 = my_plot(df_summary, model_name = "dirichlet_multinomial", show_error_bars = TRUE, name = "power")
dirichlet_multinomial_plot = ggarrange(dirichlet_multinomial_plot_1, dirichlet_multinomial_plot_2, common.legend = TRUE, legend = "right")
dirichlet_multinomial_plot
ggsave("plots/dirichlet_multinomial_single_test_plot.pdf", plot = dirichlet_multinomial_plot, device = "pdf", width = 7, height = 3, units = "in")
ggsave("plots/dirichlet_multinomial_single_test_power_plot.pdf", plot = dirichlet_multinomial_plot_2, device = "pdf", width = 3.5, height = 3, units = "in")
ggsave("plots/dirichlet_multinomial_single_test_error_plot.pdf", plot = dirichlet_multinomial_plot_1, device = "pdf", width = 3.5, height = 3, units = "in")

saveRDS(dirichlet_multinomial_plot_1, "plots/dirichlet_multinomial_error.rds")
saveRDS(dirichlet_multinomial_plot_2, "plots/dirichlet_multinomial_power.rds")


####################################################################################
power_single = readRDS("../run107/plots/dirichlet_multinomial_power.rds") + theme(legend.position = "bottom") + guides(color = guide_legend(ncol = 1))
power_fdr = readRDS("../run109/plots/power_plot_fdr.rds") + theme(legend.position = "bottom") + guides(color = guide_legend(ncol = 1))
power_fwer = readRDS("../run109/plots/power_plot_fwer.rds") + theme(legend.position = "bottom") + guides(color = guide_legend(ncol = 1))
power_plot = ggarrange(power_single,power_fwer,power_fdr,ncol = 3, legend = "bottom")
power_plot
ggsave("plots/dirichlet_multinomial_combined_power_plot.pdf", plot = power_plot, device = "pdf", width = 9, height = 5, units = "in")

error_single = readRDS("../run107/plots/dirichlet_multinomial_error.rds") + theme(legend.position = "bottom") + guides(color = guide_legend(ncol = 1))
fdr = readRDS("../run109/plots/error_plot_fdr.rds") + theme(legend.position = "bottom") + guides(color = guide_legend(ncol = 1))
fwer = readRDS("../run109/plots/error_plot_fwer.rds") + theme(legend.position = "bottom") + guides(color = guide_legend(ncol = 1))
error_plot = ggarrange(error_single,fwer,fdr,ncol = 3, legend = "bottom")
error_plot
ggsave("plots/dirichlet_multinomial_combined_error_plot.pdf", plot = error_plot, device = "pdf", width = 9, height = 5, units = "in")

