pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "vroom", "jhtools", "glue", "readxl", "ggsci", "patchwork", "cowplot", 
          "tidyverse", "dplyr", "WGCNA", "pheatmap","openxlsx","rtracklayer","tibble", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

workdir <- "~/projects/icu/analysis/huanglingtong/WGCNA_v5" %>% checkdir()
setwd(workdir)
library(FSA, lib.loc = "/cluster/home/zhangjie_jh/my_Rlib")
library(dunn.test, lib.loc = "/cluster/home/zhangjie_jh/my_Rlib")

rds_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5/rds" %>% checkdir()
plot_dir <- "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/WGCNA_v5" %>% checkdir()
data_dir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/_"

# samp_info
samp_info <- read_rds(glue("{data_dir}/patient_data.rds"))
RNA_count <- read.csv(glue("{data_dir}/rna_matrix_vst.csv", row.names = 1))
RNA_count <- column_to_rownames(RNA_count, var = "X")

telomere_length <- read.csv(glue("{data_dir}/telomere_length_pass_one_per_patient.csv"))
telomere_length <- telomere_length[, c("patient_id","TelomereLength_qmotif")]
samp_info <- left_join(samp_info, telomere_length, by = "patient_id")

extracted_name <- sub("_D1$", "", names(RNA_count))
intersect_sample <- intersect(samp_info$patient_id, extracted_name)
samp_info <- samp_info[samp_info$patient_id %in% intersect_sample, ]
saveRDS(samp_info, file = glue("{rds_dir}/sampleinfo_subsample.rds"))

survival_info <- read.csv(glue("{data_dir}/patient_survival_nonVenti_data.csv"))
classify_info <- readRDS(glue("{data_dir}/patient_data_meta.rds"))
samp_info <- left_join(samp_info,survival_info,by = "patient_id")
samp_info$true_name <- paste0(samp_info$patient_id,"_D1")
saveRDS(samp_info,file = glue("{rds_dir}/sampleinfo_subsample_with_survival_time.rds"))

# high variable genes
vars <- apply(RNA_count, 1, var)
n_hvg <- 6000
genes <- names(sort(vars, decreasing = TRUE))[1:n_hvg]

write.csv(RNA_count[genes,], glue("{rds_dir}/clean_rnaseq_nor_code_gene_hvg6000.csv")) 

rnaData <- read.csv(file = glue::glue("{rds_dir}/clean_rnaseq_nor_code_gene_hvg6000.csv"), row.names = 1)
rnaData_t <- rnaData %>% t()

# don't log
# rnaData_log <- log2(rnaData_t + 1)
rnaData_log <- rnaData_t

sampleTree = hclust(dist(rnaData_log), method = "average")

## 异常tree检测
pdf(glue("{plot_dir}/tree.pdf"), width = 20, height = 9)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

pdf(glue("{plot_dir}/tree_check.pdf"), width = 20, height = 9)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
abline(h=235,col="red")
dev.off()

data = rnaData_log

datTraits <- samp_info[,c("patient_id","TelomereLength_qmotif","Age","BMI","APACHEII_24h","hsCRP","PCT","Ddimer_ugL",
                              "APTT","PT", "WBC","BUN","Creatinine","UricAcid","Albumin","TotalBileAcids","TotalBilirubin","DirectBilirubin",
                              "ALT","AST","GGT","ALP","LDH","SOFA_24h","OI_24h")]
datTraits$patient_id <- paste0(datTraits$patient_id,"_D1")
intersect_sample <- intersect(datTraits$patient_id,rownames(data))
data <- data[rownames(data) %in% intersect_sample,]
datTraits <- datTraits[datTraits$patient_id %in% intersect_sample,]
rownames(datTraits) <- NULL
datTraits <- column_to_rownames(datTraits, var = "patient_id")
write.csv(as.data.frame(data), glue("{rds_dir}/WGCNA_input_data.csv"))

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(data, powerVector=powers,
                        networkType="unsigned", verbose=5)
saveRDS(sft, glue("{rds_dir}/WGCNA_select_softThread.rds"))

pdf(glue("{plot_dir}/select_softThread.pdf"), width = 12, height = 8)
par(mfrow = c(1,2))
# 理论此处应该大于等于0.85！！！！！
cex1 = 0.85
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()

power = 14

# 这里冲突cor，需要修改
cor <- WGCNA::cor
net = blockwiseModules(data, power = power, maxBlockSize = ncol(data),
                       TOMType = "signed", minModuleSize = 20, randomSeed = 54321, 
                       numericLabels = TRUE, # reassignThreshold = 0,
                       pamRespectsDendro = FALSE, saveTOMs = TRUE, nThreads = 64, 
                       saveTOMFileBase = glue("{workdir}/blockwiseTOM"), verbose = 3)
cor<-stats::cor

saveRDS(net, glue("{rds_dir}/WGCNA_net.rds"))

# net <- readRDS(glue("{workdir}/Fig1/WGCNA/WGCNA_net.rds"))
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs_tmp = moduleEigengenes(data, moduleColors)
MEs = MEs_tmp$eigengenes
MEs = orderMEs(MEs)

gene_color_map <- data.frame(gene_name = colnames(data),
                             module_color = glue("ME{MEs_tmp$validColors}"))

write.csv(gene_color_map, glue("{rds_dir}/gene_color_map.csv"))

write.csv(MEs, glue("{rds_dir}/module_score.csv"))

color = moduleColors[net$blockGenes[[1]]]

pdf(glue("{plot_dir}/WGCNA_module.pdf"), width = 8, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cluster Dendrogram",
                    font.main=3,cex.main=2.5,
                    cex.colorLabels=1, font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
dev.off()

# barplot TNM grade-----------------------------------------------------------------
# 定义起始和结束颜色
MEs <- read.csv(glue("{rds_dir}/module_score.csv"), row.names = 1)
color.bin <- c("#00599F", "#d80700")

# 创建渐变调色板函数
pal <- colorRampPalette(color.bin)

# 生成 3 个渐变颜色（包括两端）
colors_3 <- pal(3)

# 连续变量相关性分析
# 保持rowname一致
datTraits <- datTraits[rownames(MEs), ]

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 863)

textMatrix <- paste(signif(moduleTraitCor, 2), "(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")

pdf(glue("{plot_dir}/corre_heatmap_v5.pdf"),width = 18,height = 9)
par(mar = c(8, 16, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.5,0.5),
               main = paste("module-TelomereLength_qmotif"))
dev.off()

output_file <- list()
output_file[["corre"]] <- as.data.frame(moduleTraitCor) %>%
  tibble::rownames_to_column("Module")

output_file[["pvalue"]] <- as.data.frame(moduleTraitPvalue) %>%
  tibble::rownames_to_column("Module")
write_xlsx(output_file, glue("{rds_dir}/module_trait_results.xlsx"))
pdf(paste0(glue("{plot_dir}/all_modules_scatterplots.pdf")), width = 12, height = 8)

# 设置图形布局
n_traits <- ncol(datTraits)
n_modules <- ncol(MEs)

# 对每个模块绘制散点图矩阵
for(i in 1:n_modules) {
  module_name <- names(MEs)[i]
  module_eigengene <- MEs[, i]
  
  # 设置每页显示4个性状
  par(mar = c(3, 3, 3, 1))  # 调整边距
  
  for(j in 1:n_traits) {
    trait_name <- names(datTraits)[j]
    
    # 计算相关性和p值
    trait_cor <- cor(datTraits[, j], module_eigengene, use = "pairwise.complete.obs")
    trait_p <- corPvalueStudent(trait_cor, nrow(datTraits))
    
    plot(datTraits[, j], module_eigengene,
         xlab = trait_name, 
         ylab = paste(module_name, "Eigengene"),
         main = paste(module_name, "vs", trait_name, "\n",
                      "Cor =", signif(trait_cor, 3),
                      "P =", signif(trait_p, 3)),
         pch = 16, col = "blue", cex = 0.8)
    
    # 添加回归线
    if(sum(!is.na(datTraits[, j]) & !is.na(module_eigengene)) > 2) {
      abline(lm(module_eigengene ~ datTraits[, j]), col = "red", lwd = 2)
    }
  }
}

dev.off()



