suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(ggplot2)
  library(ggpubr)
  library(ggthemes)
  library(survival)
  library(survminer)
  library(glue)
})

workdir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5" %>% setwd()
rds_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5/rds"
MEs <- read.csv(glue("{rds_dir}/module_score.csv"), row.names = 1)
spinfo <- readRDS(glue("{rds_dir}/sampleinfo_subsample_with_survival_time.rds"))

spinfo <- spinfo[spinfo$true_name %in% rownames(MEs), ]
MEs <- MEs[spinfo$true_name,]
MEs_t <- t(MEs)

spinfo1 <- data.frame(
  sample_id = spinfo$true_name,
  OS_month = spinfo$time,
  OS_status = as.numeric(spinfo$status))

spinfo1 <- spinfo1[spinfo1$sample_id %in% rownames(MEs),]
#----------------------------------------------------------------------------------
#  Step 2: Module score ~ OS (coxph, prognosis)
#----------------------------------------------------------------------------------

module.name <- rownames(MEs_t)

plot.info <- NULL
plot.data <- spinfo1
plot.info.os <- NULL

## (1) Module score ~ OS
for (i in 1:length(module.name)) {
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
plot.name <- classify_info[classify_info$class %in% "con",] %>% pull(var)
plot.name <- c(plot.name,"BMI","TelomereLength_qmotif","time","NonVentilation_Time")

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
plot.name <- classify_info[classify_info$class %in% "bool",] %>% pull(var)
plot.name[! plot.name %in% colnames(spinfo)]
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

pdf("./module_score_vs_clinical_traits.pdf", 
    width = 35, height = 15)   
print(p)
dev.off()
write.xlsx(plot.info, "./module_score_vs_clinical_traits.xlsx")



# -------------------------------------------------------------------------
# Renamed sample info
# -------------------------------------------------------------------------
write.xlsx(spinfo, "./clinical_traits_clean.xlsx")

# -------------------------------------------------------------------------
# MEblack & MEblack Kaplan-Meier curves
# -------------------------------------------------------------------------
# (1) Sample information and clinical characteristics 
spinfo <- readRDS(glue("{rds_dir}/sampleinfo_subsample_with_survival_time.rds"))
# (2) Scores of the modules
MEs <- read.csv(glue("{workdir}/module_score.csv"), row.names = 1)
MEs <- MEs[spinfo$true_name,]
MEs_t <- t(MEs)
# Overall survival stratified by MEblack
plot.cohort <- spinfo
plot.cohort <- plot.cohort %>% 
  mutate(MEblue_pro = as.numeric(MEs_t["MEblue", spinfo$true_name]))
res.cut <- surv_cutpoint(plot.cohort, 
                         time = "time", 
                         event = "status",
                         variables = "MEblue_pro",
                         minprop = 0.2)
plot.cohort <- plot.cohort %>% 
  mutate(MEblue_pro_l = ifelse(MEblue_pro >= res.cut$cutpoint$cutpoint, "High", "Low"),
         MEblue_pro_l = factor(MEblue_pro_l, levels = c("Low", "High")))


# coxph
info <- summary(coxph(Surv(time, status) ~ MEblue_pro_l, data = plot.cohort))
anno.text <- ""
for (i in 1:nrow(info$conf.int)) {
  anno.text <- paste0(anno.text, "\n", paste0(rownames(info$conf.int)[i], " HR=", round(info$conf.int[i, 1], 3), " CI=", round(info$conf.int[i, 3], 3), "-", round(info$conf.int[i, 4], 3), " P=", signif(info$coefficients[i, 5], 4) ))
}
anno.text <- paste0(anno.text, "\nKaplan-Meier P=", signif(survdiff(Surv(time, status) ~ MEblue_pro_l, data = plot.cohort)$pvalue, 4) )
anno.text <- str_replace_all(anno.text, "MEblue_pro_l", "")
fit <- survfit(Surv(time, status) ~ MEblue_pro_l, data = plot.cohort)
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
                 title= paste0("OS MEblue_pro_l \n", anno.text))


# Overall survival stratified by MEblack
plot.cohort <- spinfo
plot.cohort <- plot.cohort %>% 
  mutate(MEblack_pro = as.numeric(MEs_t["MEblack", spinfo$true_name]))
res.cut <- surv_cutpoint(plot.cohort, 
                         time = "time", 
                         event = "status",
                         variables = "MEblack_pro",
                         minprop = 0.2)
plot.cohort <- plot.cohort %>% 
  mutate(MEblack_pro_l = ifelse(MEblack_pro >= res.cut$cutpoint$cutpoint, "High", "Low"),
         MEblack_pro_l = factor(MEblack_pro_l, levels = c("Low", "High")))
# coxph
info <- summary(coxph(Surv(time, status) ~ MEblack_pro_l, data = plot.cohort))
anno.text <- ""
for (i in 1:nrow(info$conf.int)) {
  anno.text <- paste0(anno.text, "\n", paste0(rownames(info$conf.int)[i], " HR=", round(info$conf.int[i, 1], 3), " CI=", round(info$conf.int[i, 3], 3), "-", round(info$conf.int[i, 4], 3), " P=", signif(info$coefficients[i, 5], 4) ))
}
anno.text <- paste0(anno.text, "\nKaplan-Meier P=", signif(survdiff(Surv(time, status) ~ MEblack_pro_l, data = plot.cohort)$pvalue, 4) )
anno.text <- str_replace_all(anno.text, "MEblack_pro_l", "")
fit <- survfit(Surv(time, status) ~ MEblack_pro_l, data = plot.cohort)
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
                 title= paste0("OS MEblack_pro_l \n", anno.text))

## output
p <- arrange_ggsurvplots(list(p1, p2), ncol = 2, nrow = 1, print = FALSE)
ggsave(paste0("./MEblue_MEblack_survival_analysis.pdf"), p, width = 12, height = 8)


#---separate plot




