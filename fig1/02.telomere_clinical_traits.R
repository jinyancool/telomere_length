suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(patchwork)
  library(rms)
  .libPaths(c(.libPaths(), "~/sbin/R/library/4.3.0/"))
  library(plotRCS)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/"
setwd(wdir)

# data
df <- readRDS("doc/patient_data_full.rds")
df_meta <- readRDS("doc/patient_data_meta.rds") %>%
  select(!tab_old)

# continuous type clinical traits
vars <- df_meta %>%
  filter(class == "con") %>%
  pull(var)
dt_con <- lapply(vars, function(var){
  res <- cor.test(df[[var]], df$TelomereLength_qmotif, method = "spearman")
  pval <- res$p.value
  rho <- res$estimate
  return(data.frame(var = var,
                    pval = pval,
                    rho = rho))
}) %>% rbindlist() %>%
  left_join(df_meta)


color_tab <- c("#80B1D3", "#FB8072", "#8DD3C7")
names(color_tab) <- df_meta %>% 
  filter(tab != "Outcome") %>% 
  pull(tab) %>% unique() %>% sort()
p <- ggplot() +
  geom_point(data = dt_con %>% filter(pval >= 0.05), aes(rho, -log10(pval), size = -log10(pval)), color = "grey70") +
  geom_point(data = dt_con %>% filter(pval < 0.05), aes(rho, -log10(pval), color = tab, size = -log10(pval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = color_tab, name = "") +
  labs(x = "Spearman correlation") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("02_telomere_clinical_traits/telomere_vs_clinial_traits_volcano_plot.pdf", width = 7, height = 6)

# categorical type clinical traits with two levels
vars <- df_meta %>%
  filter(class == "bool" | class == "dis") %>%
  pull(var)
dt_bool <- lapply(vars, function(var){
  res <- wilcox.test(df$TelomereLength_qmotif ~ as.factor(df[[var]]))
  pval <- res$p.value
  return(data.frame(var = var,
                    pval = pval))
}) %>% rbindlist() %>%
  left_join(df_meta)

# # categorical type clinical traits - more than two levels
# vars <- df_meta %>%
#   filter(class == "dis") %>%
#   pull(var)
# dt_dis <- lapply(vars, function(var){
#   res <- aov(df$TelomereLength_qmotif ~ as.factor(df[[var]])) %>% summary
#   pval <- unlist(res)[9]
#   return(data.frame(var = var,
#                     pval = pval))
# }) %>% rbindlist() %>%
#   left_join(df_meta)

# all
dt <- rbindlist(list(dt_con, dt_bool), fill = TRUE)
fwrite(dt, "02_telomere_clinical_traits/telomere_vs_clinical_traits.csv")

# # treat all variables as continuous
# vars <- df_meta %>%
#   filter(class %in% c("con", "dis", "bool")) %>%
#   pull(var)
# dt_all <- lapply(vars, function(var){
#   print(var)
#   res <- cor.test(as.numeric(as.factor(df[[var]])), df$TelomereLength_qmotif, method = "spearman")
#   pval <- res$p.value
#   rho <- res$estimate
#   return(data.frame(var = var,
#                     pval = pval,
#                     rho = rho))
# }) %>% rbindlist() %>%
#   left_join(df_meta)
# fwrite(dt_all, "02_telomere_clinical_traits/telomere_vs_clinical_traits.csv")
# 
# color_tab <- RColorBrewer::brewer.pal(8, "Set3")
# names(color_tab) <- dt_all %>% filter(pval < 0.05) %>% pull(tab) %>% unique() %>% sort()
# p <- ggplot() +
#   geom_point(data = dt_all %>% filter(pval >= 0.05), aes(rho, -log10(pval), size = -log10(pval)), color = "grey70") +
#   geom_point(data = dt_all %>% filter(pval < 0.05), aes(rho, -log10(pval), color = tab, size = -log10(pval))) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   scale_color_manual(values = color_tab, name = "") +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank())
# ggsave("02_telomere_clinical_traits/telomere_vs_clinial_traits_volcano_plot.pdf", width = 6, height = 5)

# individual plot
for(i in 1:nrow(dt)){
  var_class <- dt[i, class]
  var_name <- dt[i, var]
  
  if(var_class == "con"){
    p <- ggplot(df, aes(.data[[var_name]], TelomereLength_qmotif)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      stat_cor(method = "spearman") +
      labs(y = "Telomere Length") +
      theme_classic()
    ggsave(paste0("02_telomere_clinical_traits/correlation_plot/correlation_telomere_vs_", var_name, ".pdf"), p, width = 5, height = 5)
  } else if(var_class == "bool" | var_class == "dis"){
    p <- ggboxplot(df, x = var_name, y = "TelomereLength_qmotif", color = var_name, add = "jitter") +
      stat_compare_means() +
      scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
      labs(y = "Telomere Length") +
      guides(color = "none") +
      theme_classic()
    ggsave(paste0("02_telomere_clinical_traits/boxplot/boxplot_telomere_vs_", var_name, ".pdf"), p, width = 4, height = 5)
  } # else if(var_class == "dis"){
  #   n_level <- df[[var_name]] %>% as.factor %>% droplevels() %>% levels() %>% length()
  #   if(n_level > 2){
  #     p <- ggboxplot(df, x = var_name, y = "TelomereLength_qmotif", color = var_name, add = "jitter") +
  #       stat_compare_means(method = "anova") +
  #       scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF")) +
  #       labs(y = "Telomere Length") +
  #       guides(color = "none") +
  #       theme_classic()
  #     ggsave(paste0("02_telomere_clinical_traits/boxplot/boxplot_telomere_vs_", var_name, ".pdf"), p, width = 5, height = 5)
  #   } else {
  #     p <- ggboxplot(df, x = var_name, y = "TelomereLength_qmotif", color = var_name, add = "jitter") +
  #       stat_compare_means() +
  #       scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
  #       labs(y = "Telomere Length") +
  #       guides(color = "none") +
  #       theme_classic()
  #     ggsave(paste0("02_telomere_clinical_traits/boxplot/boxplot_telomere_vs_", var_name, ".pdf"), p, width = 4, height = 5)
  #   }
  # }
}

# ============================= rcs ==========================
# rcs
vars <- df_meta %>%
  filter(class == "con" & tab != "Outcome") %>%
  pull(var)
covars <- c("Age", "Gender", "SmokingHistory", "Immunosuppression", "ChronicLungDiseaseOrAsthma", "SolidTumor", "HM")

ddist <- datadist(df %>% select(c("TelomereLength_qmotif", covars, vars) %>% unique))
options(datadist = "ddist")

dt_rcs <- lapply(vars, function(var_name){
  if(var_name %in% covars){
    covars <- covars[!covars %in% var_name]
  }
  data <- df %>% 
    select(all_of(c("patient_id", "TelomereLength_qmotif", covars, var_name))) %>%
    na.omit()
  colnames(data) <- c("patient_id", "telomere", covars, var_name)
  dt <- tryCatch({
    f_str <- paste("telomere ~ rcs(", var_name, ", 3) + ", 
                   paste(covars, collapse = " + "), sep = "")
    f_obj <- as.formula(f_str)
    fit <- ols(f_obj, data = data)
    res <- anova(fit)
    dt <- data.table(var_name = var_name,
                     pval = res[var_name, "P"],
                     pval_nonlinear = res[" Nonlinear", "P"])
  }, error = function(e){
    dt <- data.table(var_name = var_name,
                     pval = NA,
                     pval_nonlinear = NA)
  })
  return(dt)
}) %>% rbindlist()
fwrite(dt_rcs, "02_telomere_clinical_traits/telomere_vs_clinical_traits_rcs_results.csv")

dt_sig <- dt_rcs %>% filter(pval < 0.05 | pval_nonlinear < 0.05)
for(var_name in dt_sig$var){
  covars <- c("Age", "Gender", "SmokingHistory", "Immunosuppression", "ChronicLungDiseaseOrAsthma", "SolidTumor", "HM")
  if(var_name %in% covars){
    covars <- covars[!covars %in% var_name]
  }
  data <- df %>% 
    select(all_of(c("patient_id", "TelomereLength_qmotif", covars, var_name))) %>%
    na.omit()
  p <- rcsplot(data = data,
               outcome = "TelomereLength_qmotif", 
               exposure = var_name, 
               covariates = covars,
               xlab        = var_name,              
               ylab        = "Difference in Telomere Length (95% CI)", 
               pvalue.digits  = 3,             
               pvalue.position= c(0.02, 0.98), 
               fontsize    = 16,              
               linesize    = 1.5,                        
               ref.value = "median",                      
               knots = c(0.1, 0.5, 0.9),                  
               knots.line  = FALSE,           
               conf.int    = TRUE,             
  ) + theme_classic(base_size = 14) +                         
    theme(                                                    
      axis.title = element_text(size = 14, face = "bold"),   
      axis.text  = element_text(size = 14) ,                 
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
    )
  ggsave(paste0("02_telomere_clinical_traits/rcsplot/rcs_plot_telomere_vs_", var_name, ".pdf"), p, width = 5.5, height = 5)
}

vars <- c("NDVI_average12months", "GDP2020")
for(var_name in vars){
  covars <- c("Age", "Gender", "SmokingHistory", "Immunosuppression", "ChronicLungDiseaseOrAsthma", "SolidTumor", "HM")
  if(var_name %in% covars){
    covars <- covars[!covars %in% var_name]
  }
  data <- df %>% 
    dplyr::select(all_of(c("patient_id", "TelomereLength_qmotif", covars, var_name))) %>%
    na.omit()
  p <- rcsplot(data = data,
               outcome = "TelomereLength_qmotif", 
               exposure = var_name, 
               covariates = covars,
               xlab        = var_name,              
               ylab        = "Difference in Telomere Length (95% CI)", 
               pvalue.digits  = 3,             
               pvalue.position= c(0.02, 0.98), 
               fontsize    = 16,              
               linesize    = 1.5,                        
               ref.value = "median",                      
               knots = c(0.1, 0.5, 0.9),                  
               knots.line  = FALSE,           
               conf.int    = TRUE,             
  ) + theme_classic(base_size = 14) +                         
    theme(                                                    
      axis.title = element_text(size = 14, face = "bold"),   
      axis.text  = element_text(size = 14) ,                 
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
    )
  ggsave(paste0("02_telomere_clinical_traits/rcsplot/rcs_plot_telomere_vs_", var_name, ".pdf"), p, width = 5.5, height = 5)
}
