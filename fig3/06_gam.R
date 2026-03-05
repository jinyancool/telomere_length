suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(mgcv)
  library(clusterProfiler)
})
# =========================== rcs ==============================
wdir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5/"
rds_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5/rds"
plot_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5/gam" %>% checkdir()
setwd(wdir)

# expression mat
mat_vst <- read.csv(glue("{rds_dir}/module_score.csv"), row.names = 1) %>% t() %>% as.data.frame()
genes <- rownames(mat_vst)

# sample info
df_clin_v2 <- readRDS("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/patient_data_full.rds")
outcome <- read.csv("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/patient_outcome_data.csv")
outcome <- subset(outcome, select = -c(VisualAnalogScaleScore, VAP_Presence, VAP_Pathogen,
                                       HAP_Presence, HAP_Date, HAP_Pathogen,
                                       Pneumothorax_Presence, Pneumothorax_Date, NonVentilation_Time))
df_clin <- left_join(df_clin_v2, outcome, by = "patient_id")
df_clin$true_name <- paste0(df_clin$patient_id,"_D1")

intersect_sample <- intersect(df_clin$true_name,names(mat_vst))
df_clin <- df_clin[df_clin$true_name %in% intersect_sample,]
rownames(df_clin) <- NULL
df_clin <- column_to_rownames(df_clin, var = "true_name")

plot.name <- c("Age","Gender","PneumoniaType_CAP","RespiratorySupport_24h","Immunosuppression","MyocardialInfarction",
               "SolidTumor","HM","ModerateToSevereRenalImpairment","MildRenalImpairment","SOFA_24h","TelomereLength_qmotif")
df_clin <- df_clin[, plot.name]
#修改顺序一致
mat_vst <- mat_vst[, intersect_sample]
df_clin <- df_clin[intersect_sample, ]

df_pvals <- data.frame()
mat_pred <- vector()
p_lst <- list()
se_pred <- vector()
mat_vst <- as.matrix(mat_vst)
#mat_vst <- mat_vst %>% t()

for(gene in genes){
  data <- cbind(df_clin, mat_vst[gene, ])
  colnames(data) <- c(colnames(df_clin), gene)
  colnames(data) <- gsub("[ ()-]+", "_", colnames(data))
  colnames(data) <- gsub(",", "", colnames(data))
  gene <- tail(colnames(data), 1)
  pAge_telomereLength <- rcsplot(data = data,
                                 outcome = gene,
                                 exposure = "TelomereLength_qmotif",
                                 covariates = c("Age","Gender","PneumoniaType_CAP","RespiratorySupport_24h","Immunosuppression","MyocardialInfarction",
                                                "SolidTumor","HM","ModerateToSevereRenalImpairment","MildRenalImpairment","SOFA_24h"),
                                 xlab = "Telomere Length",
                                 ylab = gene,
                                 pvalue.digits = 3,
                                 pvalue.position = c(0.05,0.90),
                                 fontsize = 16,
                                 linesize = 1.5,
                                 ref.value = "median",
                                 knots = c(0.1,0.5,0.9),
                                 knots.line = FALSE,
                                 conf.int = TRUE
  )+
    theme_classic(base_size = 14)+
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )+
    ggtitle(glue("Association_between_telomere_length_and_{gene}"))
  
  dd <- datadist(data)
  options(datadist = "dd")
  
  formula <- as.formula(paste0(
    gene, " ~ rcs(TelomereLength_qmotif, 3) + ",  # 这里交换了变量位置
    "Age + Gender + PneumoniaType_CAP + ",
    "RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + SOFA_24h"
  ))
  
  model <- ols(formula, data = data)
  
  # 3. 获取p值
  anova_results <- anova(model)
  
  # 整体p值
  overall_p <- anova_results["TelomereLength_qmotif", "P"]
  # 非线性p值  
  nonlinear_p <- anova_results[" Nonlinear", "P"]
  add_frame <- data.frame(
    overall_p = overall_p,
    nonlinear_p = nonlinear_p
  )
  rownames(add_frame) <- gene
  df_pvals <- rbind(df_pvals,add_frame)
  ggsave(pAge_telomereLength, file = glue("{plot_dir}/{gene}_telomereLength.pdf"), width = 12,height = 8)
}
write.csv(df_pvals, file = glue("{plot_dir}/telomereLength_pval_module.csv"))