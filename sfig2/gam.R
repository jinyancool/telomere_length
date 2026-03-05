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
df_pat <- readRDS("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/patient_data_full.rds")
covars <- c("Age", "Gender", "PneumoniaType_CAP", "RespiratorySupport_24h", 
            "Immunosuppression", "MyocardialInfarction", "SolidTumor", "HM", 
            "ModerateToSevereRenalImpairment", "MildRenalImpairment", "SOFA_24h")

df_bac <- readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/microbiome_relative_cfu.xlsx", sheet = "Bacteria") %>%
  tibble::column_to_rownames("taxon") %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id")
df_fungi <- readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/microbiome_relative_cfu.xlsx", sheet = "Fungi") %>%
  tibble::column_to_rownames("taxon") %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id")
df_vir <- readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_telomere/doc/microbiome_relative_cfu.xlsx", sheet = "Viruses") %>%
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




# rcs plot 
dt_rcs <- read.csv(file = "/cluster/home/yliang_jh/projects/wgs/icu_telomere/06_microbiome/microbiome_vs_telomere_rcs_results_all_levels.csv")

for(var_name in intersect(unique(dt_rcs$var_name), names(df))){
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
  ggsave(paste0("/cluster/home/xyzhang_jh/projects/icu/analysis/huanglingtong/fig4/", var_name, "_vs_telomere.pdf"), p, width = 5.5, height = 5)
}

