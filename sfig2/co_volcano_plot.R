.libPaths(new = c(.libPaths(), "/cluster/home/danyang_jh/sbin/R/R-4.3.0", "/cluster/home/danyang_jh/sbin/R/R-4.2.1"))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "QFeatures", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "MSnbase", "msqrob2", 
          "enrichplot", "org.Hs.eg.db", "parallel", "jhuanglabRNAseq", "jhuanglabGO", "limma", 
          "readxl", "writexl", "ComplexHeatmap", "circlize", "clusterProfiler", "MsCoreUtils", 
          "missForest", "paletteer", "factoextra", "FactoMineR", "openxlsx", "ggrepel","scales")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
wdir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/"
rds_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/rds"
plot_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/gam" %>% checkdir()
data_dir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/_"

output_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/"

file_path <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/gsva/gam/no_celltype/telomereLength_pval.csv"
data_total <- read.csv(file = file_path, row.names = 1)

# 创建颜色列
data_total$label <- with(data_total, ifelse(overall_p > 0.05 & nonlinear_p > 0.05, "Not-sig",
                                            ifelse(overall_p > 0.05 & nonlinear_p < 0.05, "nonlinear only",
                                                   ifelse(overall_p < 0.05 & nonlinear_p > 0.05, "overall only",
                                                          ifelse(overall_p < 0.05 & nonlinear_p < 0.05, "nonlinear & overall", NA)))))

write.xlsx(data_total, file = glue::glue(output_dir,"volcanoplot.xlsx"), rowNames = TRUE)

data_total <- read.xlsx(glue::glue(output_dir,"volcanoplot.xlsx"), rowNames = TRUE)
# 计算绝对值
data_total$abs_sum = abs(data_total$overall_p + data_total$nonlinear_p)
data_total$abs_diff = abs(data_total$overall_p - data_total$nonlinear_p)
data_total$gene_name <- rownames(data_total)

# 选择 abs(FoldChange1 + FoldChange2) 最小的五个基因
top_sum_genes <- data_total[order(data_total$abs_sum), ][1:10, ] %>% .[.$label != "Not-sig", "gene_name"]
# 选择 abs(FoldChange1 - FoldChange2) 最大的五个基因
top_diff_genes <- data_total[order(data_total$overall_p), ][1:25, ] %>% .[.$label == "overall only", "gene_name"]

nonlinear <- data_total %>% .[.$label == "nonlinear only", "gene_name"] %>% .[[1]]

label_frame <- data_frame(gene_name = c(top_sum_genes, top_diff_genes,nonlinear),
                          tag = c(top_sum_genes, top_diff_genes, nonlinear))

data_total <- left_join(data_total, label_frame, by = "gene_name")

#colorframe
col_frame  <- data_frame(label = c("Not-sig", "nonlinear only", "overall only", "nonlinear & overall"),
                         color = c("#a6a6a6", "#ed7456", "#2AA12B", "#b67fd0"))
data_total <- left_join(data_total, col_frame, by = "label")
data_total$label <- factor(data_total$label, levels = c("Not-sig", "nonlinear only", "overall only", "nonlinear & overall"))




# 绘制火山图
# 
volcano_plot <- ggplot(data_total, aes(x = nonlinear_p, y = overall_p)) +
  geom_point(aes(color = label), alpha = 0.7, size = 3) +  # 根据颜色列着色
  labs(title = "Volcano Plot",
       x = "pval (nonlinear)",
       y = "pval (overall)",
       color = "Significant") +
  
  theme_minimal()+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = c(0.95, 0.05),  # 将图例放置在右下角
    legend.justification = c("right", "bottom"),
    axis.line = element_line(colour = "black"),
    axis.ticks.length = unit(0.25, "cm"), 
    axis.ticks = element_line(color = "black"), 
    axis.text = element_text(size = 12, color = "black")
  )+
  scale_x_continuous(trans = c("log10", "reverse"), labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(trans = c("log10", "reverse"), labels = label_number(accuracy = 0.01)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "black") +  # 添加X轴虚线
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +  # 添加Y轴虚线
  
  geom_label_repel(data = data_total, aes(label = tag),
                   colour = data_total$color,
                   size = 6,                           # 设置标签大小
                   box.padding = unit(0.5, "lines"),   # 设置标签内边距
                   point.padding = unit(0.8, "lines"), # 设置标签与点的距离
                   segment.color = "black",            # 设置标签边界线颜色
                   show.legend = FALSE,                # 不显示图例
                   max.overlaps = 15)+
  
  scale_color_manual(values = c("Not-sig" = "#a6a6a6", "nonlinear only" = "#ed7456", "overall only" = "#2AA12B", "nonlinear & overall" = "#b67fd0"))

ggsave(volcano_plot, file = glue::glue(output_dir,"volcano_plot_v2.pdf"),width = 10, height = 10)










