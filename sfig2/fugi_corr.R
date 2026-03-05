pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "vroom", "jhtools", "glue", "readxl", "ggsci", "patchwork", "cowplot", 
          "tidyverse", "dplyr", "WGCNA", "pheatmap","survival","openxlsx","caret","compositions")


for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


workdir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva" %>% checkdir() %>% setwd()
rds_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/rds" %>% checkdir()

##fugi
fugi <- read.csv("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/microbiome_relative_cfu_new_all_levels.csv")
fugi <- fugi[fugi$level %in% "Genus",]
rownames(fugi) <- NULL
fugi <- fugi %>% column_to_rownames(.,"taxon") %>% .[,-c(1, 2)]
colnames(fugi) <- gsub("_D\\d+$", "_D1", colnames(fugi))

spinfo <- readRDS("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5/rds/sampleinfo_subsample_with_survival_time.rds")
fugi <- as.data.frame(fugi)
spinfo <- spinfo[spinfo$true_name %in% names(fugi),]

fugi <- fugi[,spinfo$true_name]
zero_ratio <- rowSums(fugi == 0) / ncol(fugi)
fugi <- fugi[zero_ratio <= 0.4, ]

df_clr <- clr(t(fugi))  # clr() 默认对行计算，需先转置
MEs_t <- t(df_clr)   # 转置回原始维度

spinfo1 <- data.frame(
  sample_id = spinfo$true_name,
  OS_month = spinfo$time,
  OS_status = as.numeric(spinfo$status))


spinfo1 <- spinfo1[spinfo1$sample_id %in% names(MEs_t),]
#----------------------------------------------------------------------------------
#  Step 2: Module score ~ OS (coxph, prognosis)
#----------------------------------------------------------------------------------
module.name <- rownames(MEs_t)
plot.info <- NULL
plot.data <- spinfo1
plot.info.os <- NULL
## (1) Module score ~ OS
for (i in 1:length(module.name)){
  # data:  module score -> Low/High level (median cutoff)
  plot.data <- plot.data %>% 
    mutate(ME_score = MEs_t[module.name[i], plot.data$sample_id] %>% as.numeric())
  res.cut <- surv_cutpoint(plot.data, 
                           time = "OS_month", 
                           event = "OS_status",
                           variables = "ME_score",
                           minprop = 0.2)
  plot.data <- plot.data %>% 
    mutate(ME_level = ifelse(ME_score >= res.cut$cutpoint$cutpoint, "High", "Low"),
           ME_level = factor(ME_level, levels = c("Low", "High")))
  # coxph
  fit <- coxph(Surv(OS_month, OS_status) ~ ME_level, data = plot.data)
  fit.info <- summary(fit)
  # fit.p <- survdiff(Surv(OS_month, OS_status) ~ ME_level, data = plot.data)$pvalue
  fit.info.out <- c(module.name[i], "OS", fit.info$conf.int[1, 1], fit.info$coefficients[1,5], "Log-rank test")
  plot.info.os <- rbind(plot.info.os, fit.info.out)
}
plot.info.os <- as.data.frame(plot.info.os)
colnames(plot.info.os) <- c("Module", "Factor", "Type","P", "Test")
# result
plot.info.os <- plot.info.os %>% 
  mutate(Type=as.numeric(Type),P=as.numeric(P),
         Padj= p.adjust(P, method = "fdr")) %>% 
  mutate(Tag=case_when(Type>1&P<0.05 ~ "Unfavorable",
                       Type<1&P<0.05 ~ "Favorable",
                       .default="notSig"))

#----------------------------------------------------------------------------------
#  Step 3: Module score ~ continuous type clinical trait (cor)
#----------------------------------------------------------------------------------
plot.name <- c("BMI","TelomereLength_qmotif","time","NonVentilation_Time","SurvivalTimeWithin28Days","hsCRP","PCT","NDVI_average12months","gb75_cumulative_PM25")

plot.name[! plot.name %in% colnames(spinfo)]
plot.info.cat <- NULL
for (kkk in 1:length(plot.name)) {
  plot.data <- spinfo
  plot.data$Plot <- plot.data[, which(colnames(plot.data) == plot.name[kkk])]
  plot.data <- plot.data[which(!is.na(plot.data$Plot)), ]
  plot.info.cat.sub <- NULL
  for (i in 1:length(module.name)) {
    plot.data$ME_score <- MEs_t[module.name[i], plot.data$true_name] %>% as.numeric()
    sub <- c(module.name[i], plot.name[kkk], cor(plot.data$ME_score, plot.data$Plot, method = "spearman"), 
             cor.test(plot.data$ME_score, plot.data$Plot, method = "spearman")$p.value, "Spearman's correlation")
    plot.info.cat.sub <- rbind(plot.info.cat.sub, sub)
  }
  plot.info.cat.sub <- as.data.frame(plot.info.cat.sub)
  colnames(plot.info.cat.sub) <- c("Module", "Factor", "Type","P", "Test")
  plot.info.cat.sub <- plot.info.cat.sub %>% 
    mutate(Type=as.numeric(Type),
           P=as.numeric(P),
           Padj=p.adjust(P,method="fdr")) 
  plot.info.cat <- rbind(plot.info.cat, plot.info.cat.sub)
}
plot.info.cat <- plot.info.cat %>% 
  mutate(Tag=case_when(Type<0&P<0.05 ~ "Negative",
                       Type>0&P<0.05 ~ "Positive",
                       .default="notSig"))
#----------------------------------------------------------------------------------
#  Step 4: Module score ~ category type clinical trait (Wilcoxon rank-sum test or ANOVA)
#----------------------------------------------------------------------------------
plot.name <- c("SmokingHistory","DrinkingHistory","Gender","PneumoniaType","DeathWithin28DaysAfterEnrollment","TransplantHistory","VAP_Presence","HAP_Presence","Pneumothorax_Presence")
plot.info.con <- NULL
for (kkk in 1:length(plot.name)) {
  plot.data <- spinfo
  plot.data$Plot <- plot.data[, which(colnames(plot.data) == plot.name[kkk])]
  plot.data <- plot.data[which(!is.na(plot.data$Plot)), ]
  plot.info.con.sub <- NULL
  if (length(unique(plot.data$Plot)) == 2) { # 2 types category -> wilcox
    for (i in 1:length(module.name)) {
      plot.data$ME_score <- MEs_t[module.name[i], plot.data$true_name] %>% as.numeric()
      sub <- c(module.name[i], plot.name[kkk], NA, wilcox.test(plot.data$ME_score ~ as.factor(plot.data$Plot) )$p.value, "Wilcoxon rank-sum test")
      plot.info.con.sub <- rbind(plot.info.con.sub, sub)
    } 
  } else { # more than 2 types category -> anova
    for (i in 1:length(module.name)) {
      plot.data$ME_score <- MEs_t[module.name[i], plot.data$true_name] %>% as.numeric()
      sub <- c(module.name[i], plot.name[kkk], NA, unlist(summary(aov(plot.data$ME_score ~ as.factor(plot.data$Plot) )))[9], "ANOVA")
      plot.info.con.sub <- rbind(plot.info.con.sub, sub)
    }
  }
  plot.info.con.sub <- as.data.frame(plot.info.con.sub)
  colnames(plot.info.con.sub) <- c("Module", "Factor", "Type", "P", "Test")
  plot.info.con.sub <- plot.info.con.sub %>% 
    mutate(Type=as.numeric(Type),
           P=as.numeric(P),
           Padj=p.adjust(P,method="fdr")) 
  plot.info.con <- rbind(plot.info.con, plot.info.con.sub)
}
plot.info.con <- plot.info.con %>% mutate(Tag=ifelse(P<0.05,"Sig","notSig"))


#----------------------------------------------------------------------------------
#  Step 6: plot
#----------------------------------------------------------------------------------
plot.info <- rbind(plot.info.os, plot.info.cat, plot.info.con)
plot.info <- plot.info %>% mutate(LogP=-log10(P),
                                  Module=factor(as.character(Module), levels = rev(module.name)),
                                  Factor=factor(as.character(Factor), levels = unique(as.character(Factor))))
p <- ggplot(plot.info, aes(Factor, Module)) +
  geom_point(aes(size = LogP, colour = Tag, fill = Tag), shape = 16) + 
  xlab("") + ylab("")
p <- p + theme_base() + scale_size(range = c(0,10))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))
p <- p + theme(axis.text.y = element_text(size = 20))
p <- p + scale_colour_manual(values = c(Favorable = "#00599F", Unfavorable = "#D01910", Negative = "#009632", Positive = "#8f00b7", Sig = "#ed7a00", Basal = "#ea6c59", Classical = "#1288a5", notSig = "#CCCCCC")) + theme(plot.background = element_blank())

pdf("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig3_wgcna/S3/fugi_vs_clinical_traits.pdf", 
    width = 12, height = 15)   
print(p)
dev.off()
write.xlsx(plot.info, "./module_score_vs_clinical_traits.xlsx")

# -------------------------------------------------------------------------
# Renamed sample info
# -------------------------------------------------------------------------
write.xlsx(spinfo, "./clinical_traits_clean.xlsx")

# -------------------------------------------------------------------------
# Ribosome & Ribosome Kaplan-Meier curves
# -------------------------------------------------------------------------
# (1) Sample information and clinical characteristics
MEs_t <- gsva_score
# Overall survival stratified by Ribosome
plot.cohort <- spinfo
plot.cohort <- plot.cohort %>% 
  mutate(Apoptosis_pro = as.numeric(MEs_t["Apoptosis", spinfo$true_name]))
res.cut <- surv_cutpoint(plot.cohort, 
                         time = "time", 
                         event = "status",
                         variables = "Apoptosis_pro",
                         minprop = 0.2)
plot.cohort <- plot.cohort %>% 
  mutate(Apoptosis_pro_l = ifelse(Apoptosis_pro >= res.cut$cutpoint$cutpoint, "High", "Low"),
         Apoptosis_pro_l = factor(Apoptosis_pro_l, levels = c("Low", "High")))


# coxph
info <- summary(coxph(Surv(time, status) ~ Apoptosis_pro_l, data = plot.cohort))
anno.text <- ""
for (i in 1:nrow(info$conf.int)) {
  anno.text <- paste0(anno.text, "\n", paste0(rownames(info$conf.int)[i], " HR=", round(info$conf.int[i, 1], 3), " CI=", round(info$conf.int[i, 3], 3), "-", round(info$conf.int[i, 4], 3), " P=", signif(info$coefficients[i, 5], 4) ))
}
anno.text <- paste0(anno.text, "\nKaplan-Meier P=", signif(survdiff(Surv(time, status) ~ Apoptosis_pro_l, data = plot.cohort)$pvalue, 4) )
anno.text <- str_replace_all(anno.text, "Apoptosis_pro_l", "")
fit <- survfit(Surv(time, status) ~ Apoptosis_pro_l, data = plot.cohort)
p1 <- ggsurvplot(fit, 
                 data = plot.cohort,
                 xlab = 'Time (Months)',
                 pval = TRUE,
                 risk.table = TRUE, 
                 risk.table.height = 0.28,
                 conf.int.alpha = 0.05,
                 conf.int = TRUE, 
                 palette = c("#00599F","#d80700"),
                 axes.offset = TRUE,
                 break.time.by = 12,  xlim = c(0, 24),
                 title= paste0("OS Apoptosis_pro_l \n", anno.text))


# Overall survival stratified by Ribosome
plot.cohort <- spinfo
plot.cohort <- plot.cohort %>% 
  mutate(Ribosome_pro = as.numeric(MEs_t["Ribosome", spinfo$true_name]))
res.cut <- surv_cutpoint(plot.cohort, 
                         time = "time", 
                         event = "status",
                         variables = "Ribosome_pro",
                         minprop = 0.2)
plot.cohort <- plot.cohort %>% 
  mutate(Ribosome_pro_l = ifelse(Ribosome_pro >= res.cut$cutpoint$cutpoint, "High", "Low"),
         Ribosome_pro_l = factor(Ribosome_pro_l, levels = c("Low", "High")))
# coxph
info <- summary(coxph(Surv(time, status) ~ Ribosome_pro_l, data = plot.cohort))
anno.text <- ""
for (i in 1:nrow(info$conf.int)) {
  anno.text <- paste0(anno.text, "\n", paste0(rownames(info$conf.int)[i], " HR=", round(info$conf.int[i, 1], 3), " CI=", round(info$conf.int[i, 3], 3), "-", round(info$conf.int[i, 4], 3), " P=", signif(info$coefficients[i, 5], 4) ))
}
anno.text <- paste0(anno.text, "\nKaplan-Meier P=", signif(survdiff(Surv(time, status) ~ Ribosome_pro_l, data = plot.cohort)$pvalue, 4) )
anno.text <- str_replace_all(anno.text, "Ribosome_pro_l", "")
fit <- survfit(Surv(time, status) ~ Ribosome_pro_l, data = plot.cohort)
p2 <- ggsurvplot(fit, 
                 data = plot.cohort,
                 xlab = 'Time (Months)',
                 pval = TRUE,
                 risk.table = TRUE, 
                 risk.table.height = 0.28,
                 conf.int.alpha = 0.05,
                 conf.int = TRUE, 
                 palette = c("#00599F","#d80700"),
                 axes.offset = TRUE,
                 break.time.by = 12,  xlim = c(0, 24),
                 title= paste0("OS Ribosome_pro_l \n", anno.text))

## output
p <- arrange_ggsurvplots(list(p1, p2), ncol = 2, nrow = 1, print = FALSE)
ggsave(paste0("./Apoptosis_Ribosome_survival_analysis.pdf"), p, width = 12, height = 8)
