
library(xCell)
data_dir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/_"
RNA_count <- read.csv(glue("{data_dir}/rna_matrix_vst.csv", row.names = 1))
RNA_count <- column_to_rownames(RNA_count, var = "X")
RNA_count$gene_name <- rownames(RNA_count)

ids <- clusterProfiler::bitr(RNA_count$gene_name,'ENSEMBL','SYMBOL','org.Hs.eg.db')
RNA_count_anno <- merge(RNA_count, ids, by.x='gene_name', by.y='ENSEMBL')
RNA_count_anno <- distinct(RNA_count_anno, SYMBOL, .keep_all = TRUE)
RNA_count_anno <- column_to_rownames(RNA_count_anno, var = "SYMBOL")
RNA_count_anno <- RNA_count_anno[, !names(RNA_count_anno) %in% c("SYMBOL", "gene_name")]

xcell_results <- xCellAnalysis(RNA_count_anno)
write.csv(xcell_results, file = "/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/xcell/xcell_res.csv")

