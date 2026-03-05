suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(patchwork)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere/"
setwd(wdir)

# telomere length qmotif vs TelSeq
dt_tel <- fread("doc/telomere_length_pass_one_per_patient.csv")
p <- ggplot(dt_tel, aes(log2(TelomereLength_telseq), log2(TelomereLength_qmotif))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "spearman") +
  labs(x = "TelSeq log2 ratio", y = "qmotif log2 ratio") +
  theme_classic()
ggsave("01_basics/correlation_qmotif_telseq.pdf", p, width = 5, height = 5)

# telomere length histogram
p <- ggplot(dt_tel, aes(TelomereLength_qmotif)) +
  geom_histogram(bins = 20, fill = "white", color = "black") +
  labs(x = "Telomere Length") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(10, 20, 30)) +  
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000))
ggsave("01_basics/telomere_length_hist.pdf", p, width = 5, height = 5)

# patient info
dt_pat <- readRDS("doc/patient_data.rds")
dt_m <- dt_tel %>%
  left_join(dt_pat)

# telomere length histogram with age
p <- ggplot(dt_m, aes(TelomereLength_qmotif, fill = Gender)) +
  geom_histogram(bins = 20, color = "black") +
  scale_fill_manual(values = c(male = "#A6CEE3", female = "#FDBF6F")) +
  labs(x = "Telomere Length") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(10, 20, 30)) +  
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000))
ggsave("01_basics/telomere_length_hist_gender.pdf", p, width = 6, height = 5)

# telomere length histogram with lung disease
p <- ggplot(dt_m, aes(TelomereLength_qmotif, fill = factor(ChronicLungDiseaseOrAsthma))) +
  geom_histogram(bins = 20, color = "black") +
  scale_fill_manual(values = c("#B2DF8A", "#EE6A50FF")) +
  labs(x = "Telomere Length", fill = "Chronic Lung Disease") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(10, 20, 30)) +  
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000))
ggsave("01_basics/telomere_length_hist_lungdisease.pdf", p, width = 6.5, height = 5)

# hospital vs telomere length
p <- ggplot(dt_m, aes(Hospital, TelomereLength_qmotif, color = Hospital)) +
  geom_boxplot(width = 0.7) +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                                "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", 
                                "#7E6148FF", "#B09C85FF", "#FF6B6BFF", "#8F7700FF", 
                                "#7876B1FF", "#1A936FFF", "#114B5FFF")) +
  labs(y = "Telomere Length", x = "") +
  guides(color = "none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +  
  scale_y_continuous(breaks = c(1000, 2000, 3000)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
ggsave("01_basics/telomere_length_hospital.pdf", p, width = 7, height = 5)

# age histogram with gender
p <- ggplot(dt_m, aes(Age, fill = Gender)) +
  geom_histogram(bins = 20, color = "black") +
  scale_fill_manual(values = c(male = "#A6CEE3", female = "#FDBF6F")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("01_basics/age_hist_gender.pdf", p, width = 6, height = 5)

# depth vs telomere
dt_qc <- fread("/cluster/home/yliang_jh/projects/wgs/icu_huanglingtong/output/qc/qc_summary_with_spinfo.csv")
dt_tel <- fread("doc/telomere_length_pass_one_per_patient.csv")
pt <- dt_tel %>%
  left_join(dt_qc, by = join_by(wgs_id == Sample)) %>%
  dplyr::select(wgs_id, TelomereLength_qmotif, mean_coverage)

p <- ggplot(pt, aes(mean_coverage, TelomereLength_qmotif)) +
  geom_point() +
  labs(x = "Sequencing depth", y = "Telomere length") +
  scale_x_continuous(trans = "log10") +
  theme_classic()
ggsave("01_basics/qc_depth_vs_telomere_length.pdf", p, width = 5, height = 5)
