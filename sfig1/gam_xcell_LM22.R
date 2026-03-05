
.libPaths(c(.libPaths(), "/cluster/home/yliang_jh/sbin/R/library/4.3.0/"))
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(mgcv)
  library(clusterProfiler)
  library(plotRCS)
  library(rms)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/"
workdir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong"
rds_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/rds"
plot_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/gam" %>% checkdir()
data_dir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/_"

setwd(workdir)

# xcell mat
xcell <- read.csv("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/xcell_res.csv", row.names = 1)

# cybersortx mat
cybersortx <- read.csv("/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/rnaseq/balf/balf_cibersort_lm22_results_mtx_from_log_cpm.csv")
cybersortx$day <- sapply(strsplit(cybersortx$true_name, "_"), tail, 1)
cybersortx$samp_id <- sub("_[^_]+$", "", cybersortx$true_name)


cybersortx_sub <- cybersortx %>%
  group_by(samp_id) %>%
  # 提取D后面的数字并转换为数值
  dplyr::mutate(day_num = as.numeric(gsub("D", "", day))) %>%
  # 筛选出每个患者最小的D值对应的行
  dplyr::filter(day_num == min(day_num)) %>%
  # 移除辅助列
  dplyr::select(-day_num) %>%
  ungroup() %>% 
  as.data.frame()

cybersortx_sub$samp_id <- paste0(cybersortx_sub$samp_id,"_D1")

cybersortx_sub <- cybersortx_sub %>%
  mutate(across(where(is.character), ~ gsub("HZSY", "HZYY", .)))


#cybersortx_sub <- readRDS(file = "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/rds/cybersortx_sub.rds")
# sample info
df_clin <- read_rds(glue("{rds_dir}/df_clin.rds"))
cybersortx_sub <- cybersortx_sub[cybersortx_sub$samp_id %in% rownames(df_clin),]
rownames(cybersortx_sub) <- NULL
cybersortx_sub <- cybersortx_sub %>% column_to_rownames(.,"samp_id")
cybersortx_sub <- cybersortx_sub[, !colnames(cybersortx_sub) %in% c("true_name", "P.value", "Correlation","RMSE","day")]

#统一行名
xcell <- xcell %>% t() %>% as.data.frame()
xcell <- xcell[rownames(df_clin),]
cybersortx_sub <- cybersortx_sub[rownames(df_clin),]


#xcell gam res----------------------------------------------------------------------
genes <- names(xcell)
xcell <- as.data.frame(xcell)
#identical(rownames(xcell), rownames(df_clin))
df_pvals <- data.frame()
mat_pred <- vector()
p_lst <- list()
se_pred <- vector()
mat_vst <- as.matrix(xcell)
#mat_vst <- mat_vst %>% t()
#Epithelial_cells + ImmuneScore

for(gene in genes){
  tryCatch({
  data <- cbind(df_clin, mat_vst[,gene])
  colnames(data) <- c(colnames(df_clin), gene)
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
    ggtitle(glue("Association_between_telomere_length_and_xcell_{gene}"))
  
  dd <- datadist(data)
  options(datadist = "dd")
  
  formula <- as.formula(paste0(
    "TelomereLength_qmotif ~ rcs(", gene, ", 3) + ",
    "Age + Gender + PneumoniaType_CAP + ",
    "RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + SOFA_24h"
  ))
  
  model <- ols(formula, data = data)
  
  # 3. 获取p值
  anova_results <- anova(model)
  
  # 整体p值
  overall_p <- anova_results[gene, "P"]
  # 非线性p值  
  nonlinear_p <- anova_results[" Nonlinear", "P"]
  add_frame <- data.frame(
    overall_p = overall_p,
    nonlinear_p = nonlinear_p
  )
  rownames(add_frame) <- gene
  df_pvals <- rbind(df_pvals,add_frame)
  ggsave(pAge_telomereLength, file = glue("{plot_dir}/xcell/{gene}_telomereLength.pdf"), width = 12,height = 8)
  
  }, error = function(e) {
    # message(paste("跳过基因:", gene, "- 错误:", conditionMessage(e)))
  })
}
write.csv(df_pvals, file = glue("{plot_dir}/xcell/telomereLength_pval_xcell.csv"))



#cybersortx_sub gam res----------------------------------------------------------------------
genes <- names(cybersortx_sub)
cybersortx_sub <- as.data.frame(cybersortx_sub)
#identical(rownames(xcell), rownames(df_clin))
df_pvals <- data.frame()
mat_pred <- vector()
p_lst <- list()
se_pred <- vector()
mat_vst <- as.matrix(cybersortx_sub)
#mat_vst <- mat_vst %>% t()
#Epithelial_cells + ImmuneScore

for(gene in genes){
  tryCatch({
    data <- cbind(df_clin, mat_vst[,gene])
    colnames(data) <- c(colnames(df_clin), gene)
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
      ggtitle(glue("Association_between_telomere_length_and_LM22_{gene}"))
    
    dd <- datadist(data)
    options(datadist = "dd")
    
    formula <- as.formula(paste0(
      "TelomereLength_qmotif ~ rcs(", gene, ", 3) + ",
      "Age + Gender + PneumoniaType_CAP + ",
      "RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + SOFA_24h"
    ))
    
    model <- ols(formula, data = data)
    
    # 3. 获取p值
    anova_results <- anova(model)
    
    # 整体p值
    overall_p <- anova_results[gene, "P"]
    # 非线性p值  
    nonlinear_p <- anova_results[" Nonlinear", "P"]
    add_frame <- data.frame(
      overall_p = overall_p,
      nonlinear_p = nonlinear_p
    )
    rownames(add_frame) <- gene
    df_pvals <- rbind(df_pvals,add_frame)
    ggsave(pAge_telomereLength, file = glue("{plot_dir}/LM22/{gene}_telomereLength.pdf"), width = 12,height = 8)
    
  }, error = function(e) {
    # 出错时简单跳过，不打印信息
    # 如果需要可以取消注释下面的行
    # message(paste("跳过基因:", gene, "- 错误:", conditionMessage(e)))
  })
}
write.csv(df_pvals, file = glue("{plot_dir}/LM22/telomereLength_pval_LM22.csv"))
