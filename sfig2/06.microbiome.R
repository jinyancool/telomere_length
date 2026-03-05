suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(rms)
  library(ggrepel)
  library(scales)
  .libPaths(c(.libPaths(), "~/sbin/R/library/4.3.0/"))
  library(plotRCS)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/"
setwd(wdir)

# =========================== rcs ==============================
# data
df_pat <- readRDS("doc/patient_data_full.rds")
df_bac <- readxl::read_xlsx("doc/microbiome_relative_cfu.xlsx", sheet = "Bacteria") %>%
  tibble::column_to_rownames("taxon") %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id")
df_fungi <- readxl::read_xlsx("doc/microbiome_relative_cfu.xlsx", sheet = "Fungi") %>%
  tibble::column_to_rownames("taxon") %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id")
df_vir <- readxl::read_xlsx("doc/microbiome_relative_cfu.xlsx", sheet = "Viruses") %>%
  tibble::column_to_rownames("taxon") %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id")
covars <- c("Age", "Gender", "PneumoniaType_CAP", "RespiratorySupport_24h", 
            "Immunosuppression", "MyocardialInfarction", "SolidTumor", "HM", 
            "ModerateToSevereRenalImpairment", "MildRenalImpairment", "SOFA_24h")
df <- df_pat %>%
  dplyr::select(patient_id, sample_id, TelomereLength_qmotif, all_of(covars)) %>%
  left_join(df_bac) %>%
  left_join(df_fungi) %>% 
  left_join(df_vir)

# rcs
ddist <- datadist(df)
options(datadist = "ddist")

vars <- colnames(df)[(length(covars) + 4):ncol(df)]
dt_rcs <- lapply(vars, function(var_name){
  data <- df %>% 
    select(all_of(c("patient_id", "TelomereLength_qmotif", covars, var_name))) %>%
    na.omit()
  dt <- tryCatch({
    f_str <- paste(var_name, " ~ rcs(TelomereLength_qmotif, 3) + ", 
                   paste(covars, collapse = " + "), sep = "")
    f_obj <- as.formula(f_str)
    fit <- ols(f_obj, data = data)
    res <- anova(fit)
    dt <- data.table(var_name = var_name,
                     pval = res["TelomereLength_qmotif", "P"],
                     pval_nonlinear = res[" Nonlinear", "P"])
  }, error = function(e){
    dt <- data.table(var_name = var_name,
                     pval = NA,
                     pval_nonlinear = NA)
  })
  return(dt)
}) %>% rbindlist()
fwrite(dt_rcs, "06_microbiome/microbiome_vs_telomere_rcs_results.csv")

# pval plot
pt_rcs <- dt_rcs %>%
  mutate(sig = ifelse(pval >= 0.05 & pval_nonlinear >= 0.05, "Not-sig",
                        ifelse(pval >= 0.05 & pval_nonlinear < 0.05, "Nonlinear only",
                               ifelse(pval < 0.05 & pval_nonlinear >= 0.05, "Overall only",
                                      ifelse(pval < 0.05 & pval_nonlinear < 0.05, "Nonlinear & overall", NA)))))
color_sig <- c("Not-sig" = "#a6a6a6",
               "Nonlinear only" = "#ed7456",
               "Overall only" = "#2AA12B",
               "Nonlinear & overall" = "#b67fd0")
pt_top <- pt_rcs %>%
  filter(sig %in% c("Overall only", "Nonlinear & overall")) %>%
  arrange(pval) %>%
  slice_head(n = 6)
p <- ggplot(pt_rcs, aes(x = pval_nonlinear, y = pval, color = sig)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(trans = c("log10", "reverse"), labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(trans = c("log10", "reverse"), labels = label_number(accuracy = 0.01)) +
  scale_color_manual(values = color_sig) +
  geom_text_repel(data = pt_top, aes(label = var_name),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3) +
  labs(x = "p-value (non-linear)", y = "p-value (overall)", color = "Significance") +
  theme_classic()
ggsave("06_microbiome/scatter_pval_vs_nonlinear_pval.pdf", p, width = 8, height = 8)

# rcs plot 
dt_sig <- dt_rcs %>% filter(pval < 0.05 | pval_nonlinear < 0.05)
for(var_name in dt_sig$var){
  data <- df %>% 
    select(all_of(c("patient_id", "TelomereLength_qmotif", covars, var_name))) %>%
    na.omit()
  p <- rcsplot(data = data,
               outcome = var_name, 
               exposure = "TelomereLength_qmotif", 
               covariates = covars,
               xlab        = "TelomereLength_qmotif", 
               pvalue.digits  = 3,             
               pvalue.position= c(0.02, 0.98), 
               fontsize    = 16,              
               linesize    = 1.5,                        
               ref.value = "median",                      
               knots = c(0.1, 0.5, 0.9),                  
               knots.line  = FALSE,           
               conf.int    = TRUE,             
  ) + theme_classic(base_size = 14) +                         
    theme(                                                    
      axis.title = element_text(size = 14, face = "bold"),   
      axis.text  = element_text(size = 14) ,                 
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
    )
  ggsave(paste0("06_microbiome/rcsplot/rcs_plot_", var_name, "_vs_telomere.pdf"), p, width = 5.5, height = 5)
}

# microbe composition
mat_bac <- read.csv("07_multiomics/mat_pred_loess.csv", row.names = 1) %>%
  .[colnames(df_bac)[-1], ]
mat_bac[mat_bac < 0] <- 0
sums <- rowSums(mat_bac)
bac10 <- sort(-sums) %>% head(10) %>% names()
dt_bac_cp <- mat_bac[!rownames(mat_bac) %in% bac10, ] %>% colSums() %>%
  rbind(mat_bac[bac10, ]) %>%
  tibble::rownames_to_column("species") %>% as.data.table() %>%
  melt(id.vars = "species", variable.name = "telomere", value.name = "cfu") %>%
  mutate(species = ifelse(species == "1", "Others", species) %>% factor(levels = rev(c(bac10, "Others"))),
         telomere = as.numeric(sub("X", "", telomere))) %>%
  filter(telomere %in% seq(300, 3800, 100))
color_bac <- RColorBrewer::brewer.pal(11, "Set3")
names(color_bac) <- c(bac10[1:8], "Others", bac10[9:10])

p <- ggplot(dt_bac_cp, aes(x = telomere, y = cfu, fill = species)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color_bac) +
  theme_classic()
ggsave("06_microbiome/barplot_bateria_composition_loess_pred.pdf", p, width = 10, height = 5)

# microbe composition - non-linear
bac_sig <- dt_rcs[pval_nonlinear < 0.05, var_name] %>%
  .[. %in% colnames(df_bac)]
mat_bac_sig <- read.csv("07_multiomics/mat_pred_loess.csv", row.names = 1) %>%
  .[bac_sig, ]
mat_bac_sig[mat_bac_sig < 0] <- 0
sums <- rowSums(mat_bac_sig)
bac_sig10 <- sort(-sums) %>% head(10) %>% names()
dt_bac_sig_cp <- mat_bac_sig[!rownames(mat_bac_sig) %in% bac_sig10, ] %>% colSums() %>%
  rbind(mat_bac_sig[bac_sig10, ]) %>%
  tibble::rownames_to_column("species") %>% as.data.table() %>%
  melt(id.vars = "species", variable.name = "telomere", value.name = "cfu") %>%
  mutate(species = ifelse(species == "1", "Others", species) %>% factor(levels = rev(c(bac_sig10, "Others"))),
         telomere = as.numeric(sub("X", "", telomere))) %>%
  filter(telomere %in% seq(300, 3800, 100))
color_bac_sig <- RColorBrewer::brewer.pal(11, "Set3")
names(color_bac_sig) <- c(bac_sig10[1:8], "Others", bac_sig10[9:10])

p <- ggplot(dt_bac_sig_cp, aes(x = telomere, y = cfu, fill = species)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color_bac_sig) +
  theme_classic()
ggsave("06_microbiome/barplot_bateria_composition_loess_pred_nonlinear.pdf", p, width = 10, height = 5)

# =========================== rcs all levels ============================
# data
df_pat <- readRDS("doc/patient_data_full.rds")
df_micro <- read.csv("doc/microbiome_relative_cfu_new_all_levels.csv")
df_level <- df_micro %>% dplyr::select(taxon:kingdom)
df_micro <- df_micro %>%
  dplyr::select(!level:kingdom) %>%
  tibble::column_to_rownames("taxon") %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id")

covars <- c("Age", "Gender", "PneumoniaType_CAP", "RespiratorySupport_24h", 
            "Immunosuppression", "MyocardialInfarction", "SolidTumor", "HM", 
            "ModerateToSevereRenalImpairment", "MildRenalImpairment", "SOFA_24h")
df <- df_pat %>%
  dplyr::select(patient_id, sample_id, TelomereLength_qmotif, all_of(covars)) %>%
  left_join(df_micro)

# rcs
ddist <- datadist(df)
options(datadist = "ddist")

vars <- colnames(df)[(length(covars) + 4):ncol(df)]
dt_rcs <- lapply(vars, function(var_name){
  data <- df %>% 
    dplyr::select(all_of(c("patient_id", "TelomereLength_qmotif", covars, var_name))) %>%
    na.omit()
  dt <- tryCatch({
    f_str <- paste(var_name, " ~ rcs(TelomereLength_qmotif, 3) + ", 
                   paste(covars, collapse = " + "), sep = "")
    f_obj <- as.formula(f_str)
    fit <- ols(f_obj, data = data)
    res <- anova(fit)
    dt <- data.table(var_name = var_name,
                     pval = res["TelomereLength_qmotif", "P"],
                     pval_nonlinear = res[" Nonlinear", "P"])
  }, error = function(e){
    dt <- data.table(var_name = var_name,
                     pval = NA,
                     pval_nonlinear = NA)
  })
  return(dt)
}) %>% rbindlist()
dt_rcs <- left_join(dt_rcs, df_level, join_by(var_name == taxon))
fwrite(dt_rcs, "06_microbiome/microbiome_vs_telomere_rcs_results_all_levels.csv")

# pval plot
pt_rcs <- dt_rcs %>%
  mutate(sig = case_when(pval >= 0.05 & pval_nonlinear >= 0.05 ~ "Not-sig",
                         pval >= 0.05 & pval_nonlinear <  0.05 ~ "Nonlinear only",
                         pval <  0.05 & pval_nonlinear >= 0.05 ~ "Overall only",
                         pval <  0.05 & pval_nonlinear <  0.05 ~ "Nonlinear & overall")) %>%
  mutate(level = ifelse(sig == "Not-sig", "Not-sig", level))

pt_top <- pt_rcs %>%
  filter(sig %in% c("Overall only", "Nonlinear & overall")) %>%
  arrange(pval) %>%
  group_by(sig) %>%
  slice_head(n = 5)

alpha_sig <- c("Nonlinear & overall" = 1,
               "Overall only"        = 0.6,
               "Nonlinear only"      = 0.3,
               "Not-sig"             = 1)
color_level <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#a6a6a6")
names(color_level) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Not-sig")
p <- ggplot(pt_rcs, aes(x = pval_nonlinear, y = pval, color = level, alpha = sig, shape = kingdom)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(trans = c("log10", "reverse"), labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(trans = c("log10", "reverse"), labels = label_number(accuracy = 0.01)) +
  scale_color_manual(values = color_level) +
  scale_alpha_manual(values = alpha_sig, guide = "none") +
  scale_shape_manual(values = c(Bacteria = 16, Fungi = 17, Viruses = 15)) +
  geom_text_repel(data = pt_top, aes(label = var_name),
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.3) +
  labs(x = "p-value (non-linear)", y = "p-value (overall)", color = "Level") +
  theme_classic()
ggsave("06_microbiome/scatter_pval_vs_nonlinear_pval_all_level.pdf", p, width = 8, height = 8)

# rcs plot 
dt_sig <- dt_rcs %>% filter(pval < 0.05 | pval_nonlinear < 0.05)
for(var_name in dt_sig$var){
  data <- df %>% 
    dplyr::select(all_of(c("patient_id", "TelomereLength_qmotif", covars, var_name))) %>%
    na.omit()
  p <- rcsplot(data = data,
               outcome = var_name, 
               exposure = "TelomereLength_qmotif", 
               covariates = covars,
               xlab        = "TelomereLength_qmotif", 
               pvalue.digits  = 3,             
               pvalue.position= c(0.02, 0.98), 
               fontsize    = 16,              
               linesize    = 1.5,                        
               ref.value = "median",                      
               knots = c(0.1, 0.5, 0.9),                  
               knots.line  = FALSE,           
               conf.int    = TRUE,             
  ) + theme_classic(base_size = 14) +                         
    theme(                                                    
      axis.title = element_text(size = 14, face = "bold"),   
      axis.text  = element_text(size = 14) ,                 
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
    )
  ggsave(paste0("06_microbiome/rcsplot_all_level/rcs_plot_", var_name, "_vs_telomere.pdf"), p, width = 5.5, height = 5)
}

# todo: microbe composition - genus


# ========================= mediation by gsva ==============================
library(mediation)
library(mgcv)

# data
mat_gsva <- readRDS("doc/rna_gsva_kegg.rds")
df_gsva <- mat_gsva %>% t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(patient_id = stringr::str_extract(sample_id, "^[^_]+_[^_]+")) %>%
  dplyr::select(!sample_id)
df_medi <- df %>%
  left_join(df_gsva)

covars <- c("Age", "Gender", "PneumoniaType_CAP", "RespiratorySupport_24h", 
            "Immunosuppression", "MyocardialInfarction", "SolidTumor", "HM", 
            "ModerateToSevereRenalImpairment", "MildRenalImpairment", "SOFA_24h")
microbe_name <- "Enterococcus_faecium"

# telomere - gsva - microbe
res_lst <- pbmcapply::pbmclapply(rownames(mat_gsva), function(gsva_name){
  data <- df_medi %>%
    dplyr::select("patient_id", "TelomereLength_qmotif", all_of(c(covars, gsva_name, microbe_name))) %>%
    rename(microbe = microbe_name,
           gsva = gsva_name) %>%
    na.omit()
  
  # mediator model - gam
  fit_med <- gam(gsva ~ s(TelomereLength_qmotif, k = 5) + Age + Gender + PneumoniaType_CAP + 
                   RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + 
                   SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + 
                   SOFA_24h, data = data, method = "REML")
  
  # outcome model - gam
  fit_out <- gam(microbe ~ s(TelomereLength_qmotif, k = 4) + gsva + Age + Gender + PneumoniaType_CAP + 
                   RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + 
                   SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + 
                   SOFA_24h, data = data, method = "REML")
  
  # mediation analysis
  tryCatch({
    model <- mediate(model.m = fit_med, model.y = fit_out, 
                     treat = "TelomereLength_qmotif", mediator = "gsva",
                     boot = TRUE, sims = 1000)
    dt <- data.table(Effect = c("ACME", "ADE", "Total Effect"),
                     Estimate = c(model$d.avg,
                                  model$z.avg,
                                  model$tau.coef),
                     CI_lower = c(model$d.avg.ci[1],
                                  model$z.avg.ci[1],
                                  model$tau.ci[1]),
                     CI_upper = c(model$d.avg.ci[2],
                                  model$z.avg.ci[2],
                                  model$tau.ci[2]),
                     pvalue = c(model$d.avg.p,
                                model$z.avg.p,
                                model$tau.p)) %>%
      mutate(gsva = gsva_name) %>%
      relocate(gsva)
  }, error = function(e){
    dt <- data.table(gsva = gsva_name)
  })
  return(dt)
}, mc.cores = 20) %>% rbindlist(fill = TRUE)
fwrite(dt, paste0("06_microbiome/mediation_gsva_", microbe, "_results.csv"))

