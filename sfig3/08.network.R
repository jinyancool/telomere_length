suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(MetaNet)
  library(ggpubr)
  library(igraph)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/"
setwd(wdir)

# data
mat_gsva <- readRDS("doc/rna_gsva_kegg.rds")
colnames(mat_gsva) <- stringr::str_extract(colnames(mat_gsva), "^[^_]+_[^_]+")

mat_micro <- read.csv("doc/microbiome_relative_cfu_new_all_levels.csv") %>%
  filter(level == "Species") %>%
  tibble::column_to_rownames("taxon") %>%
  dplyr::select(!level:kingdom) %>%
  as.matrix()
colnames(mat_micro) <- stringr::str_extract(colnames(mat_micro), "^[^_]+_[^_]+")

# match data
sps <- intersect(colnames(mat_gsva), colnames(mat_micro))
mat_gsva <- mat_gsva[, sps]
mat_micro <- mat_micro[, sps]

# CLR of micro
clr <- function(x, pseudocount = 1) {
  x <- x + pseudocount
  logx <- log(x)
  logx - mean(logx)
}
pseudo <- min(mat_micro[mat_micro > 0]) / 2
mat_micro_clr <- apply(mat_micro, 2, clr, pseudocount = pseudo)

# telomere group
dt_tl <- fread("doc/telomere_length_pass_one_per_patient.csv")
df_meta <- data.frame(sample = sps) %>%
  left_join(dt_tl[, .(patient_id, TelomereLength_qmotif)], join_by(sample == patient_id)) %>%
  mutate(group = ifelse(TelomereLength_qmotif > 2000, "High", "Low")) %>%
  select(!TelomereLength_qmotif)

# vertex size
gsva_size <- abs(rowSums(mat_gsva))
df_gsva <- data.frame(name = names(gsva_size), 
                      gsva = gsva_size)
micro_size <- rowSums(mat_micro)
df_micro <- data.frame(name = names(micro_size), 
                      size = micro_size)

# build network using spearman correlation
for(gp in c("High", "Low")){
  sp <- df_meta %>% filter(group == gp) %>% pull(sample)
  corr_lst <- lapply(1:100, function(i){
    sp_sub <- sample(sp, 100)
    mat1 <- mat_gsva[, sp_sub]
    mat2 <- mat_micro_clr[, sp_sub]
    corr <- c_net_calculate(t(mat1), t(mat2), method = "spearman")
  })
  r_lst <- lapply(corr_lst, function(x) x$r)
  r_mat <- Reduce("+", r_lst) / 100
  p_lst <- lapply(corr_lst, function(x) x$p.value)
  p_mat <- Reduce("+", p_lst) / 100
  corr <- list(r = r_mat, p.value = p_mat)
  class(corr) <- "corr"
  
  co_net <- c_net_build(corr, r_threshold = 0.3, p_threshold = 0.05, delete_single = T)
  co_net <- c_net_set(co_net, 
                      data.frame("gsva_score" = abs(rowSums(mat1))),
                      data.frame("abundance" = rowSums(mat2)),
                      vertex_size = c("gsva_score", "abundance"))
  pdf(paste0("08_network/network_telomere_group_", gp, ".pdf"))
  plot(co_net); dev.off()
  pdf(paste0("08_network/network_telomere_group_", gp, "_with_label.pdf"))
  plot(co_net, vertex.label = V(co_net)$label); dev.off()
  saveRDS(co_net, paste0("08_network/network_telomere_group_", gp, ".rds"))
  df_v <- get_v(co_net) %>% 
    left_join(net_par(co_net, mode = "v")$v_index)
  fwrite(df_v, paste0("08_network/network_nodes_attri_telomere_group_", gp, ".csv"))
}

# compare topological indexes
df_topo <- lapply(c("High", "Low"), function(gp){
  sp <- df_meta %>% filter(group == gp) %>% pull(sample)
  
  df_lst <- lapply(1:100, function(i){
    sp_sub <- sample(sp, 100)
    mat1 <- mat_gsva[, sp_sub]
    mat2 <- mat_micro_clr[, sp_sub]
    corr <- c_net_calculate(t(mat1), t(mat2), method = "spearman")
    co_net <- c_net_build(corr, r_threshold = 0.3, p_threshold = 0.05, delete_single = T)
    net_par(co_net)$n_index
  })
  df <- do.call(rbind, df_lst) %>%
    mutate(group = gp) %>%
    relocate(group)
  return(df)
}) %>% do.call(rbind, .)
fwrite(df_topo, "08_network/network_topology_telomere_group.csv")

pt_topo <- df_topo %>%
  dplyr::select(!c(Diameter, Clustering_coefficient)) %>% 
  mutate(group = factor(group, levels = c("Low", "High"))) %>%
  as.data.table() %>%
  melt(id.vars = c("group"), variable.name = "topo", value.name = "value")
p <- ggboxplot(pt_topo, x = "group", y = "value", color = "group", add = "jitter", 
               add.params = list(size = 0.5)) +
  facet_wrap(~ topo, scales = "free", ncol = 4) +
  stat_compare_means(comparisons = list(c("High", "Low")), vjust = -0.2) +
  scale_color_manual(values = c("#80B1D3", "#FB8072")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_classic()
ggsave("08_network/boxplot_network_topology_telomere_group.pdf", p, width = 10, height = 8)

# telomere group vs age
df <- readRDS("doc/patient_data_full.rds")
pt <- df %>% select(patient_id, TelomereLength_qmotif, Age) %>%
  mutate(group = ifelse(TelomereLength_qmotif > 2000, "High", "Low")) %>%
  select(!TelomereLength_qmotif)
p <- ggboxplot(pt, x = "group", y = "Age", color = "group", add = "jitter") +
  stat_compare_means(comparisons = list(c("High", "Low")), method = "t.test") +
  scale_color_manual(values = c("#80B1D3", "#FB8072")) +
  theme_classic()
ggsave("08_network/boxplot_age_vs_telomere_group.pdf", p)

# ========================= age adjusted ============================
library(ppcor)

pcor_calculate_dt <- function(mat1, mat2, z, method = "spearman"){
  dt <- data.table(names1 = rep(colnames(mat1), ncol(mat2)),
                   names2 = rep(colnames(mat2), each = ncol(mat1)))
  dt[, c("r", "p.value") := as.list(pcor.test(mat1[, names1], mat2[, names2], z, method = method))[1:2], by = .I]
  
  mat_r <- dcast(dt, names1 ~ names2, value.var = "r") %>%
    tibble::column_to_rownames("names1") %>% as.matrix()
  mat_p <- dcast(dt, names1 ~ names2, value.var = "p.value") %>%
    tibble::column_to_rownames("names1") %>% as.matrix()
  
  corr <- list(r = mat_r, p.value = mat_p)
  class(corr) <- "corr"
  return(corr)
}
pcor_calculate_mat <- function(mat1, mat2, z, method = "spearman"){
  X <- as.matrix(mat1)
  Y <- as.matrix(mat2)
  Z <- as.matrix(z)
  n <- nrow(Z)
  Zc <- cbind(1, z)
  Pz <- Zc %*% solve(crossprod(Zc)) %*% t(Zc)
  Mz <- diag(n) - Pz
  Xr <- Mz %*% X
  Yr <- Mz %*% Y
  rownames(Xr) <- rownames(mat1)
  rownames(Yr) <- rownames(mat2)
  corr <- c_net_calculate(Xr, Yr, method = "spearman")
  return(corr)
}
pcor_calculate_lm <- function(mat1, mat2, z, method = "spearman"){
  rz1 <- residuals(lm(mat1 ~ z))
  rz2 <- residuals(lm(mat2 ~ z))
  corr <- c_net_calculate(rz1, rz2, method = "spearman")
  return(corr)
}

for(gp in c("High", "Low")){
  sp <- df_meta %>% filter(group == gp) %>% pull(sample)
  corr_lst <- lapply(1:100, function(i){
    sp_sub <- sample(sp, 100)
    mat1 <- mat_gsva[, sp_sub]
    mat2 <- mat_micro_clr[, sp_sub]
    z <- df %>%
      filter(patient_id %in% sp_sub) %>%
      arrange(match(patient_id, sp_sub)) %>%
      pull(Age)
    corr <- pcor_calculate_lm(t(mat1), t(mat2), z, method = "spearman")
    co_net <- c_net_build(corr, r_threshold = 0.3, p_threshold = 0.05, delete_single = T)
    return(c(corr, as.list(net_par(co_net)$n_index)))
  })
  r_lst <- lapply(corr_lst, function(x) x$r)
  r_mat <- Reduce("+", r_lst) / 100
  p_lst <- lapply(corr_lst, function(x) x$p.value)
  p_mat <- Reduce("+", p_lst) / 100
  corr <- list(r = r_mat, p.value = p_mat)
  class(corr) <- "corr"
  
  co_net <- c_net_build(corr, r_threshold = 0.3, p_threshold = 0.05, delete_single = T)
  co_net <- c_net_set(co_net, 
                      data.frame("gsva_score" = abs(rowSums(mat1))),
                      data.frame("abundance" = rowSums(mat2)),
                      vertex_size = c("gsva_score", "abundance"))
  pdf(paste0("08_network/age_adjusted_network_telomere_group_", gp, ".pdf"))
  plot(co_net); dev.off()
  pdf(paste0("08_network/age_adjusted_network_telomere_group_", gp, "_with_label.pdf"))
  plot(co_net, vertex.label = V(co_net)$label); dev.off()
  saveRDS(co_net, paste0("08_network/age_adjusted_network_telomere_group_", gp, ".rds"))
  df_v <- get_v(co_net) %>% 
    left_join(net_par(co_net, mode = "v")$v_index)
  fwrite(df_v, paste0("08_network/age_adjusted_network_nodes_attri_telomere_group_", gp, ".csv"))
  
  df_topo <- lapply(corr_lst, function(x) as.data.frame(x[3:14])) %>% 
    do.call(rbind, .) %>%
    mutate(group = gp) %>%
    relocate(group)
  fwrite(df_topo, paste0("08_network/age_adjusted_network_topology_telomere_group_", gp, ".csv"))
}

pt_topo <- fread("08_network/age_adjusted_network_topology_telomere_group_High.csv") %>%
  rbind(fread("08_network/age_adjusted_network_topology_telomere_group_Low.csv")) %>%
  dplyr::select(!c(Diameter, Clustering_coefficient)) %>% 
  mutate(group = factor(group, levels = c("Low", "High"))) %>%
  melt(id.vars = c("group"), variable.name = "topo", value.name = "value")
p <- ggboxplot(pt_topo, x = "group", y = "value", color = "group", add = "jitter", 
               add.params = list(size = 0.5)) +
  facet_wrap(~ topo, scales = "free", ncol = 4) +
  stat_compare_means(comparisons = list(c("High", "Low")), vjust = -0.2) +
  scale_color_manual(values = c("#80B1D3", "#FB8072")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_classic()
ggsave("08_network/age_adjusted_boxplot_network_topology_telomere_group.pdf", p, width = 10, height = 8)

# =============================== segment ================================
dt_seg <- fread("/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/output/rna/nonlinear_gene_telomere_segment_result.csv")
p <- ggplot(dt_seg, aes(breakpoint)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 20, color = "black", fill = "white") +
  geom_density(color = "red") +
  theme_classic()
ggsave("08_network/nonlinear_gene_telomere_breakpoint_density_plot.pdf", p, width = 4, height = 4)
