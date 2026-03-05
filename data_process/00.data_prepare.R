library(data.table)
library(magrittr)
library(dplyr)
library(DESeq2)
library(stringr)

# wgs sample info
dt_sp <- readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/doc/BALF-DNA_0517_tidy-data.xlsx")
dt_sp <- dt_sp %>%
  mutate(across(c(patient_id, sample_id), ~ str_replace_all(., "HZSY", "HZYY"))) %>%
  relocate(sample_id)
fwrite(dt_sp, "doc/wgs_sample_info.csv")

# patient data
df <- readRDS("/cluster/home/tanghl_jh/projects/heatmap/docs/huanglingtong/human/patient_info/20251222_patient_metadata_v5.rds") %>%
  dplyr::rename(patient_id = HumanID) %>%
  mutate(BMI = Weight / (Height/100)^2) %>%
  left_join(readxl::read_xlsx("/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/doc/non_ventilation.xlsx"), by = join_by(patient_id == HumanID))
meta <- readRDS("/cluster/home/tanghl_jh/projects/heatmap/docs/huanglingtong/human/patient_info/20251222_var_class_v5.rds") %>%
  rbind(data.frame(var = "BMI", tab = "Baseline characteristics", class = "con", tab_old = "basic")) %>%
  rbind(data.frame(var = "NonVentilation_Time", tab = "Outcome", class = "con", tab_old = "follow"))
saveRDS(meta, "doc/patient_data_meta.rds")

# patient 603
dt_clin <- fread("/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/doc/clinical_data_for_lmm.csv")
dt_clin_s <- dt_clin[1:603, ][, patient_id := str_replace_all(patient_id, "HZSY", "HZYY")]
df_orig <- df %>%
  filter(patient_id %in% dt_clin_s$patient_id)
saveRDS(df_orig, "doc/patient_data_original603.rds")

# patient pneumonia & surv data & age >= 18
table(dt_clin_s$Diagnosis)  # NotPneumonia 21
table(df603$Age < 18) # 1
table(is.na(df603$DeathWithin28DaysAfterEnrollment)) # 5
df1 <- df_orig %>%
  filter(patient_id %in% dt_clin_s[Diagnosis != "NotPneumonia", patient_id]) %>% # 582
  filter(Age >= 18) %>% # 581
  filter(!is.na(DeathWithin28DaysAfterEnrollment)) # 576

# telomere length pass sample
dt_sp <- fread("doc/wgs_sample_info.csv")
dt <- fread("doc/telomere_length_all.csv") %>%
  left_join(dt_sp) %>%
  filter(patient_id %in% df1$patient_id)
dt$patient_id %>% unique %>% length # 448

dt_qc <- fread("doc/qc_summary_with_spinfo_snp_count.csv")
sps <- dt_qc %>%
  filter(pass) %>%
  pull(Sample)
dt_pass <- dt %>%
  filter(wgs_id %in% sps)
dt_pass$patient_id %>% unique %>% length # 403
fwrite(dt_pass, "doc/telomere_length_pass.csv")

# first day telomere length per patient
dt_one <- mutate(dt_pass, day = factor(day, levels = paste0("D", c(1, 4, 7, 14, 21)))) %>%
  arrange(day) %>%
  group_by(patient_id) %>%
  slice_min(day, n = 1) # 403
fwrite(dt_one, "doc/telomere_length_pass_one_per_patient.csv")

# patient data final 403
df2 <- df %>% filter(patient_id %in% dt_one$patient_id) %>%
  left_join(dt_one) %>%
  select(!c(TelomereLength_telseq, day)) %>%
  mutate(Diagnosis = droplevels(Diagnosis))
saveRDS(df2, "doc/patient_data_full.rds")

# survival/outcome data
df_out <- df2 %>% 
  select(c(patient_id, meta %>% filter(tab == "Outcome") %>% pull(var))) %>%
  dplyr::rename(status = DeathWithin28DaysAfterEnrollment,
                time = SurvivalTimeWithin28Days) %>%
  relocate(patient_id, status, time)
fwrite(df_out, "doc/patient_outcome_data.csv")

# rna - coding gene & one sample per patient
mat <- fread("/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/rnaseq/balf/BALFrnacounts_1018.csv") %>%
  tibble::column_to_rownames("V1") %>%
  as.matrix()
colnames(mat) <- colnames(mat) %>%
  str_replace_all(., "HZSY", "HZYY")
dt_sp <- data.table(sample_id = colnames(mat)) %>%
  mutate(patient_id = str_extract(sample_id, "^[^_]+_[^_]+"),
         day = str_extract(sample_id, "D\\d+"))
dt_sp_s <- dt_sp %>%
  filter(day %in% "D1") %>%
  filter(patient_id %in% dt_one$patient_id) # 343

gene_coding <- fread("/cluster/home/yliang_jh/ref/coding_genes_ensembl_id.txt") %>%
  pull("Gene stable ID") %>%
  unique
mat_s <- mat[rownames(mat) %in% gene_coding, dt_sp_s$sample_id]
write.csv(mat_s, "doc/rna_matrix_raw_counts.csv")
fwrite(dt_sp_s, "doc/rna_matrix_sample_info.csv")

# rna - normalize
dds <- DESeqDataSetFromMatrix(countData = mat_s,
                              colData = data.frame(sample_id = colnames(mat_s),
                                                   row.names = colnames(mat_s)),
                              design = ~ 1)
# dds <- dds[rowSums(counts(dds) > 0) > ncol(dds)/5, ]
vsd <- vst(dds, blind = TRUE)
mat_vst <- assay(vsd)
write.csv(mat_vst, "doc/rna_matrix_vst.csv")

# rna - hvg
vars <- apply(mat_vst, 1, var)
n_hvg <- 2000
genes <- names(sort(vars, decreasing = TRUE))[1:n_hvg]
write.csv(mat_vst[genes, ], "doc/rna_matrix_vst_hvg2000.csv")

# rna gsva
gsva <- readRDS("/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/rnaseq/balf/gsva_score_kegg.rds")
colnames(gsva) <- colnames(gsva) %>%
  str_replace_all(., "HZSY", "HZYY")
gsva_s <- gsva[, colnames(mat_s)]
saveRDS(gsva_s, "doc/rna_gsva_kegg.rds")

# microbiom data
dt_one <- fread("doc/telomere_length_pass_one_per_patient.csv")
## bacteria
df1 <- readxl::read_xlsx("/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/metagenomics/rel_cfu/co_occur/cfu_rel_all_days.xlsx", sheet = "Bacteria")
df1_s <- df1 %>% 
  dplyr::select(taxon, any_of(dt_one$sample_id)) %>%
  filter(rowSums(. > 0) > (ncol(.) - 1)/5)
## Fungi
df2 <- readxl::read_xlsx("/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/metagenomics/rel_cfu/co_occur/cfu_rel_all_days.xlsx", sheet = "Fungi")
df2_s <- df2 %>% 
  dplyr::select(taxon, any_of(dt_one$sample_id)) %>%
  filter(rowSums(. > 0) > (ncol(.) - 1)/5)
## Viruses
df3 <- readxl::read_xlsx("/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/metagenomics/rel_cfu/co_occur/cfu_rel_all_days.xlsx", sheet = "Viruses")
df3_s <- df3 %>% 
  dplyr::select(taxon, any_of(dt_one$sample_id)) %>%
  filter(rowSums(. > 0) >= (ncol(.) - 1)/5)

library(openxlsx)
lst <- list(df1_s, df2_s, df3_s)
names(lst) <- c("Bacteria", "Fungi", "Viruses")
jhtools::list2excel(lst, "doc/microbiome_relative_cfu.xlsx")

# microbe all level
f_rlt <- "/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/metagenomics/rel_cfu/co_occur/cfu_rel_all_days.xlsx"
f_tax <- "/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/metagenomics/rel_cfu/co_occur/tax_info_lst.xlsx"

df_lst <- lapply(c("Bacteria", "Fungi", "Viruses"), function(sheet){
  df1 <- readxl::read_xlsx(f_rlt, sheet = sheet)
  df1_tax <- readxl::read_xlsx(f_tax, sheet = sheet)
  
  sps <- colnames(df1)[-1]
  level_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  df_lst <- lapply(level_all, function(level){
    df_level <- df1 %>%
      left_join(df1_tax, join_by(taxon == sp)) %>%
      group_by(.data[[level]]) %>%
      summarise(across(all_of(sps), \(x) sum(x, na.rm = TRUE))) %>%
      ungroup() %>%
      rename(taxon = all_of(level)) %>%
      mutate(level = level)
    return(df_level)
  })
  df_level <- do.call(rbind, df_lst) %>%
    dplyr::select(taxon, level, any_of(dt_one$sample_id)) %>%
    filter(rowSums(. > 0) > (ncol(.) - 2)/5)
  return(df_level)
})  

names(df_lst) <- c("Bacteria", "Fungi", "Viruses")
jhtools::list2excel(df_lst, "doc/microbiome_relative_cfu_all_levels.xlsx")

f_abs <- "/cluster/home/danyang_jh/projects/infectious/analysis/huanglingtong/human/metagenomics/rel_cfu/co_occur/cfu_abs_all_days.xlsx"
df_lst <- lapply(c("Bacteria", "Fungi", "Viruses"), function(sheet){
  df1 <- readxl::read_xlsx(f_abs, sheet = sheet)
  df1_tax <- readxl::read_xlsx(f_tax, sheet = sheet)
  
  sps <- colnames(df1)[-1]
  level_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  df_lst <- lapply(level_all, function(level){
    df_level <- df1 %>%
      left_join(df1_tax, join_by(taxon == sp)) %>%
      group_by(.data[[level]]) %>%
      summarise(across(all_of(sps), \(x) sum(x, na.rm = TRUE))) %>%
      ungroup() %>%
      rename(taxon = all_of(level)) %>%
      mutate(level = level)
    return(df_level)
  })
  df_level <- do.call(rbind, df_lst) %>%
    dplyr::select(taxon, level, any_of(dt_one$sample_id)) %>%
    filter(rowSums(. > 0) > (ncol(.) - 2)/5)
  return(df_level)
})  

names(df_lst) <- c("Bacteria", "Fungi", "Viruses")
jhtools::list2excel(df_lst, "doc/microbiome_absolute_cfu_all_levels.xlsx")  

# recalibrate relative cfu
df_lst <- lapply(c("Bacteria", "Fungi", "Viruses"), function(sheet){
  df1 <- readxl::read_xlsx(f_abs, sheet = sheet)
  df1_tax <- readxl::read_xlsx(f_tax, sheet = sheet)
  
  sps <- colnames(df1)[-1]
  level_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  df_lst <- lapply(level_all, function(level){
    df_level <- df1 %>%
      left_join(df1_tax, join_by(taxon == sp)) %>%
      group_by(.data[[level]]) %>%
      summarise(across(all_of(sps), \(x) sum(x, na.rm = TRUE))) %>%
      ungroup() %>%
      rename(taxon = all_of(level)) %>%
      mutate(level = level)
    return(df_level)
  })
  df_level <- do.call(rbind, df_lst) %>%
    mutate(kingdom = sheet) %>%
    dplyr::select(taxon, level, kingdom, any_of(dt_one$sample_id))
  return(df_level)
})
df_new <- do.call(rbind, df_lst) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
  filter(rowSums(. > 0) > (ncol(.) - 2)/5) %>%
  group_by(level) %>%
  mutate(across(where(is.numeric), ~ . / sum(.))) %>%
  ungroup() %>%
  mutate(taxon = gsub(" ", "_", taxon) %>% gsub("\\[|\\]", "", .))
fwrite(df_new, "doc/microbiome_relative_cfu_new_all_levels.csv")

# extubation
dt <- fread("/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/doc/extubation.tsv")
fwrite(dt, "doc/extubation_data.csv")
