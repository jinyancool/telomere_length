pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "vroom", "jhtools", "glue", "readxl", "ggsci", "patchwork", "cowplot", 
          "tidyverse", "dplyr", "WGCNA", "pheatmap","openxlsx","rtracklayer","tibble")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/"
workdir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong"
rds_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/rds"
plot_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/pearson" %>% checkdir()
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


# barplot TNM grade-----------------------------------------------------------------
# 定义起始和结束颜色
color.bin <- c("#00599F", "#d80700")

# 创建渐变调色板函数
pal <- colorRampPalette(color.bin)

# 生成 3 个渐变颜色（包括两端）
colors_3 <- pal(3)

# 连续变量相关性分析
# 保持rowname一致
df_clin <- df_clin[rownames(cybersortx_sub), c("Age","TelomereLength_qmotif")]

moduleTraitCor <- cor(cybersortx_sub, df_clin, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 343)

textMatrix <- paste(signif(moduleTraitCor, 2), "(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")

pdf("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/pearson/corre_heatmap_v2_LM22.pdf",width = 12,height = 12)
par(mar = c(11, 11.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(df_clin),
               yLabels = names(cybersortx_sub),
               ySymbols = names(cybersortx_sub),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.5,0.5),
               main = paste("LM22-TelomereLength_qmotif"))
dev.off()

output_file <- list()
output_file[["corre"]] <- as.data.frame(moduleTraitCor) %>%
  tibble::rownames_to_column("Module")

output_file[["pvalue"]] <- as.data.frame(moduleTraitPvalue) %>%
  tibble::rownames_to_column("Module")
write_xlsx(output_file, glue("{plot_dir}/module_trait_results_LM22.xlsx"))

pdf(glue("{plot_dir}/all_modules_scatterplots_LM22.pdf"), width = 12, height = 8)
# 设置图形布局
n_traits <- ncol(df_clin)
n_modules <- ncol(cybersortx_sub)

# 对每个模块绘制散点图矩阵
for(i in 1:n_modules) {
  module_name <- names(cybersortx_sub)[i]
  module_eigengene <- cybersortx_sub[, i]
  
  # 设置每页显示4个性状
  par(mar = c(3, 3, 3, 1))  # 调整边距
  
  for(j in 1:n_traits) {
    trait_name <- names(df_clin)[j]
    
    # 计算相关性和p值
    trait_cor <- cor(df_clin[, j], module_eigengene, use = "pairwise.complete.obs")
    trait_p <- corPvalueStudent(trait_cor, nrow(df_clin))
    
    plot(df_clin[, j], module_eigengene,
         xlab = trait_name, 
         ylab = paste(module_name, "Eigengene"),
         main = paste(module_name, "vs", trait_name, "\n",
                      "Cor =", signif(trait_cor, 3),
                      "P =", signif(trait_p, 3)),
         pch = 16, col = "blue", cex = 0.8)
    
    # 添加回归线
    if(sum(!is.na(df_clin[, j]) & !is.na(module_eigengene)) > 2) {
      abline(lm(module_eigengene ~ df_clin[, j]), col = "red", lwd = 2)
    }
  }
}

dev.off()

#XCELL--------------------------------------------------------------------------
xcell <- read.csv("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/xcell_res.csv", row.names = 1)
xcell <- t(xcell)
color.bin <- c("#00599F", "#d80700")

# 创建渐变调色板函数
pal <- colorRampPalette(color.bin)

# 生成 3 个渐变颜色（包括两端）
colors_3 <- pal(3)

# 连续变量相关性分析
# 保持rowname一致
df_clin <- df_clin[rownames(xcell), c("Age","TelomereLength_qmotif")]

moduleTraitCor <- cor(xcell, df_clin, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 343)

textMatrix <- paste(signif(moduleTraitCor, 2), "(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")

pdf("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/pearson/corre_heatmap_v2_xcell.pdf",width = 15,height = 15)
par(mar = c(15, 15, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(df_clin),
               yLabels = colnames(xcell),
               ySymbols = colnames(xcell),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.5,0.5),
               main = paste("xcell-TelomereLength_qmotif"))
dev.off()

xcell <- as.data.frame(xcell)
output_file <- list()
output_file[["corre"]] <- as.data.frame(moduleTraitCor) %>%
  tibble::rownames_to_column("Module")

output_file[["pvalue"]] <- as.data.frame(moduleTraitPvalue) %>%
  tibble::rownames_to_column("Module")
write_xlsx(output_file, glue("{plot_dir}/module_trait_results_xcells.xlsx"))

pdf(glue("{plot_dir}/all_modules_scatterplots_xcell.pdf"), width = 12, height = 8)
# 设置图形布局
n_traits <- ncol(df_clin)
n_modules <- ncol(xcell)

# 对每个模块绘制散点图矩阵
for(i in 1:n_modules) {
  module_name <- names(xcell)[i]
  module_eigengene <- xcell[, i]
  
  # 设置每页显示4个性状
  par(mar = c(3, 3, 3, 1))  # 调整边距
  
  for(j in 1:n_traits) {
    trait_name <- names(df_clin)[j]
    
    # 计算相关性和p值
    trait_cor <- cor(df_clin[, j], module_eigengene, use = "pairwise.complete.obs")
    trait_p <- corPvalueStudent(trait_cor, nrow(df_clin))
    
    plot(df_clin[, j], module_eigengene,
         xlab = trait_name, 
         ylab = paste(module_name, "Eigengene"),
         main = paste(module_name, "vs", trait_name, "\n",
                      "Cor =", signif(trait_cor, 3),
                      "P =", signif(trait_p, 3)),
         pch = 16, col = "blue", cex = 0.8)
    
    # 添加回归线
    if(sum(!is.na(df_clin[, j]) & !is.na(module_eigengene)) > 2) {
      abline(lm(module_eigengene ~ df_clin[, j]), col = "red", lwd = 2)
    }
  }
}

dev.off()




