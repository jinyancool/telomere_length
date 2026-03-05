suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(mgcv)
  library(ComplexHeatmap)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/"
fig_dir <- ""
setwd(wdir)

#Ly summary all
dt_all <- read.csv("/cluster/home/yliang_jh/projects/wgs/icu_telomere/07_multiomics/spearman_correlation_res.csv")

LM22_cor <- read.xlsx("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/pearson/module_trait_results_LM22.xlsx", sheet = 1)
LM22_pval <- read.xlsx("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/pearson/module_trait_results_LM22.xlsx", sheet = 2)

LM22 <- data.frame(
  var = LM22_cor$Module,
  pval = LM22_pval$TelomereLength_qmotif,
  rho = LM22_cor$TelomereLength_qmotif,
  tab = "LM22"
)

dt_all <- dt_all[dt_all$tab != "Transcriptome", ]

# gsva
# expression mat
mat_vst <- readRDS("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/rna_gsva_kegg.rds") %>% as.data.frame()
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
saveRDS(df_clin,"/cluster/home/xyzhang_jh/projects/icu/docs/df_clin.rds")
#修改顺序一致
mat_vst <- mat_vst[, intersect_sample]

gene = rownames(mat_vst)[[1]]

dt_rna <- lapply(rownames(mat_vst), function(gene){
  rna <- mat_vst[gene, ]
  tel <- df_clin$TelomereLength_qmotif[match(names(rna), rownames(df_clin))]
  rna <- as.numeric(rna)
  res <- cor.test(rna, tel, method = "spearman")
  pval <- res$p.value
  rho <- res$estimate
  return(data.frame(var = gene,
                    pval = pval,
                    rho = rho, 
                    tab = "GSVA"))
}) %>% rbindlist()

# save results
dt_all <- rbind(dt_all, dt_rna)
fwrite(dt_all, "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig1/spearman_correlation_res.csv")
saveRDS(mat_vst, "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig1/rna_gsva_kegg.rds")
# volcano plot
pt <- dt_all %>%
  mutate(tab = factor(tab, levels = c("Baseline characteristics", "Disease-related clinical variables", 
                                      "Pre-admission exposure", "Bacteria", "Fungi", "Viruses", "GSVA"))) %>%
  arrange(desc(tab))
color_tab <- c("#FDB462", "#FB8072", "#8DD3C7", "#BEBADA", "#FFFFB3", "#B3DE69", "#80B1D3","#FFED6F")
names(color_tab) <- levels(pt$tab)
p <- ggplot() +
  geom_point(data = pt %>% filter(pval >= 0.05), aes(rho, -log10(pval), size = -log10(pval)), color = "grey70") +
  geom_point(data = pt %>% filter(pval < 0.05), aes(rho, -log10(pval), color = tab, size = -log10(pval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = color_tab, name = "") +
  labs(x = "Spearman correlation") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig1/volcano_plot.pdf", width = 7, height = 6)


##pieplot---
dt_all <- read.csv("/cluster/home/yliang_jh/projects/wgs/icu_telomere/07_multiomics/spearman_correlation_res.csv")
dt_total <- rbind(dt_all, dt_rna)

pt <- dt_total %>%
  mutate(tab = factor(tab, levels = c("Baseline characteristics", "Disease-related clinical variables", 
                                      "Pre-admission exposure", "Bacteria", "Fungi", "Viruses", "GSVA", "Transcriptome"))) %>%
  arrange(desc(tab))

# n

dt_n_dplyr <- pt %>%
  group_by(tab) %>%
  summarise(N = n())

pt <- as.data.table(pt) 
dt_n <- pt[, .N, by = tab]

# pie charts
dt_pie <- pt[, .(pct = sum(pval < 0.05) / .N * 100), by = tab][order(tab), ] %>%
  mutate(pct := signif(pct, 3))

color_tab <- c("#FDB462", "#FB8072", "#8DD3C7", "#BEBADA", "#FFFFB3", "#B3DE69", "#80B1D3","#FFED6F")
names(color_tab) <- levels(pt$tab)

p_lst <- lapply(1:nrow(dt_pie), function(i){
  # data
  pct <- dt_pie[i, pct]
  tab <- dt_pie[i, tab]
  data <- data.frame(percentage = c(pct, 100 - pct), type = c("sig", "non-sig"))
  
  # color
  color_fill <- c(color_tab[tab], "grey85")
  names(color_fill) <- c("sig", "non-sig")
  
  # plot
  ggplot(data, aes(2, percentage, fill = type)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = color_fill) +
    theme_void() +
    annotate("text", x = 0, y = 0, label = paste0(tab, "\n(", pct, "%)")) +
    theme(legend.position = "none")
})
p <- gridExtra::grid.arrange(grobs = p_lst, ncol = 2)
ggsave("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig1/pie_chart_spearman_sig_percentage_v2.pdf", p, width = 5, height = 8)


# ============================= loess ================================
library(tidyverse)
library(broom)
library(caret)

set.seed(2025)

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/"
#data full
df <- readRDS("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/patient_data_full.rds")
df_meta <- readRDS("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/patient_data_meta.rds") %>%
  select(!tab_old)
vars <- df_meta %>%
  filter(tab != "Outcome") %>%
  filter(class == "con") %>%
  pull(var)

#mat_bac-------------------------------------------------------------------------------
mat_bac <- readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/microbiome_relative_cfu.xlsx", sheet = "Bacteria") %>%
  tibble::column_to_rownames("taxon") %>%
  as.matrix()

#mat_fungi------------------------------------------------------------------------------
mat_fungi <- readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/microbiome_relative_cfu.xlsx", sheet = "Fungi") %>%
  tibble::column_to_rownames("taxon") %>%
  as.matrix()

#mat_vir------------------------------------------------------------------------------
mat_vir <- readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/microbiome_relative_cfu.xlsx", sheet = "Viruses") %>%
  tibble::column_to_rownames("taxon") %>%
  as.matrix()
#mat_gsva------------------------------------------------------------------------------
mat_gsva <- readRDS("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig1/rna_gsva_kegg.rds")

# data
df1 <- df %>% select(patient_id, sample_id, TelomereLength_qmotif, all_of(vars))
df2 <- mat_bac %>% t() %>% as.data.table(keep.rownames = "sample_id")
df3 <- mat_fungi %>% t() %>% as.data.table(keep.rownames = "sample_id")
df4 <- mat_vir %>% t() %>% as.data.table(keep.rownames = "sample_id")
df5 <- mat_gsva %>% t() %>% as.data.table(keep.rownames = "sample_id") %>%
  mutate(patient_id = stringr::str_extract(sample_id, "^[^_]+_[^_]+")) %>%
  dplyr::select(!sample_id)
df_all <- left_join(df1, df2) %>%
  left_join(df3) %>%
  left_join(df4) %>%
  left_join(df5)






vars_all <- colnames(df_all)[-c(1:3)]

# cross-validation to find the best span
find_best_span <- function(data, span_grid = seq(0.2, 0.8, by = 0.05), k = 5){
  folds <- sample(rep(1:k, length.out = nrow(data)))
  cv_mse <- sapply(span_grid, function(span) {
    mse <- numeric(k)
    for (i in 1:k) {
      train <- data[folds != i, ]
      test  <- data[folds == i, ]
      
      fit <- loess(var ~ TelomereLength_qmotif,
                   data = train,
                   span = span,
                   degree = 1,
                   control = loess.control(surface = "direct"))
      
      pred <- predict(fit, newdata = test$TelomereLength_qmotif)
      mse[i] <- mean((test$var - pred)^2, na.rm = TRUE)
    }
    return(mean(mse))
  })
  best_span <- span_grid[which.min(cv_mse)]
  return(best_span)
}

# loess pred matrix
tel_grid <- seq(500, 3500, by = 50)
pred_lst <- pbmcapply::pbmclapply(vars_all, function(var_test){
  data <- df_all %>%
    select(TelomereLength_qmotif, var_test) %>%
    na.omit() %>%
    dplyr::rename(var = var_test)
  best_span <- find_best_span(data)
  fit <- loess(var ~ TelomereLength_qmotif, data = data,
               span = best_span, degree = 1,
               control = loess.control(surface = "direct"))
  pred <- predict(fit, newdata = tel_grid)
  return(pred)
}, mc.cores = 20)
mat_pred <- do.call(rbind, pred_lst)
rownames(mat_pred) <- vars_all
colnames(mat_pred) <- tel_grid
write.csv(mat_pred, "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig1/mat_pred_loess.csv")

# heatmap
mat_scale <- apply(mat_pred, 1, scale) %>% t()

tel_show <- seq(500, 3500, by = 500)
col_label <- rep("", length(tel_grid))
col_label[tel_grid %in% tel_show] <- as.character(tel_show)

a = intersect(rownames(mat_scale), names(df1))
b = intersect(rownames(mat_scale), names(df2))
c = intersect(rownames(mat_scale), names(df3))
d = intersect(rownames(mat_scale), names(df4))
e = intersect(rownames(mat_scale), names(df5))

# 2. 创建行注释
label_frame <- data.frame(
  label = c(a,b,c,d,e),
  cluster = c(rep("Baseline",length(a)),
              rep("Bacteria",length(b)),
              rep("Fungi",length(c)),
              rep("Viruses",length(d)),
              rep("Gsva",length(e)))
)
label_frame <- label_frame %>% column_to_rownames(var = "label")
label_frame <- label_frame[rownames(mat_scale),]

row_ha <- rowAnnotation(
  Group = label_frame,
  col = list(Group = setNames(c("#FDB462", "#BEBADA", "#FFFFB3", "#B3DE69", "#80B1D3" )
                              , unique(label_frame))),
  annotation_name_side = "top"
)



ht <- Heatmap(mat_scale, name = "Z-score",
              col = circlize::colorRamp2(c(-2, 0, 2), c("#2C7BB6", "#FFFFB3", "#D7191C")),
              left_annotation = row_ha,
              row_split = label_frame,
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              # show_column_names = FALSE,
              column_labels = col_label,
              use_raster = TRUE)
pdf("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig1/heatmap_loess_prediction_v2.pdf", width = 5, height = 8)
draw(ht)
dev.off()






