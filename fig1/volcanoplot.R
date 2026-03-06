suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(patchwork)
  library(rms)
  .libPaths(c(.libPaths(), "~/sbin/R/library/4.3.0/"))
  library(plotRCS)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/"
setwd(wdir)

# data
df <- readRDS("doc/patient_data_full.rds")
df_meta <- readRDS("doc/patient_data_meta.rds") %>%
  select(!tab_old)

# continuous type clinical traits
vars <- df_meta %>%
  filter(class == "con") %>%
  pull(var)
dt_con <- lapply(vars, function(var){
  res <- cor.test(df[[var]], df$TelomereLength_qmotif, method = "spearman")
  pval <- res$p.value
  rho <- res$estimate
  return(data.frame(var = var,
                    pval = pval,
                    rho = rho))
}) %>% rbindlist() %>%
  left_join(df_meta)


color_tab <- c("#80B1D3", "#FB8072", "#8DD3C7")
names(color_tab) <- df_meta %>% 
  filter(tab != "Outcome") %>% 
  pull(tab) %>% unique() %>% sort()
p <- ggplot() +
  geom_point(data = dt_con %>% filter(pval >= 0.05), aes(rho, -log10(pval), size = -log10(pval)), color = "grey70") +
  geom_point(data = dt_con %>% filter(pval < 0.05), aes(rho, -log10(pval), color = tab, size = -log10(pval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = color_tab, name = "") +
  labs(x = "Spearman correlation") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("02_telomere_clinical_traits/telomere_vs_clinial_traits_volcano_plot.pdf", width = 7, height = 6)