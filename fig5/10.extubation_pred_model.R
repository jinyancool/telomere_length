suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(mgcv)
  library(randomForest)
  library(xgboost)
  library(pROC)
  library(rms)
})

wdir <- "/cluster/home/yliang_jh/projects/wgs/icu_telomere//"
setwd(wdir)

# ======================== Age vs telomere ROC =========================
# data
covars <- c("Gender", "PneumoniaType_CAP", "RespiratorySupport_24h", 
            "Immunosuppression", "MyocardialInfarction", "SolidTumor", "HM", 
            "ModerateToSevereRenalImpairment", "MildRenalImpairment", "SOFA_24h")
df_extub <- fread("doc/extubation_data.csv")
df <- readRDS("doc/patient_data_full.rds") %>%
  left_join(df_extub) %>%
  select(patient_id, Age, TelomereLength_qmotif, extub_1:extub_4, all_of(covars)) %>%
  tibble::column_to_rownames("patient_id")

# loop over 4 extub index
extub <- "extub_1"
for(extub in paste0("extub_", 1:4)){
  data <- df %>% select(Age, TelomereLength_qmotif, all_of(extub), all_of(covars)) %>%
    dplyr::rename(extub_status = extub) %>%
    mutate(extub_status = as.factor(extub_status), 
           Gender = as.factor(Gender)) %>%
    na.omit()
  
  # data
  set.seed(2026)
  train_idx <- sample(1:nrow(data), 0.8 * nrow(data))
  train_data <- data[train_idx, ]
  test_data <- data[-train_idx, ]
  
  # 1. random forest
  rf_telo <- randomForest(extub_status ~ TelomereLength_qmotif + Gender + PneumoniaType_CAP + 
                            RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + 
                            SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + 
                            SOFA_24h, data = train_data)
  rf_age <- randomForest(extub_status ~ Age + Gender + PneumoniaType_CAP + 
                           RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + 
                           SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + 
                           SOFA_24h, data = train_data)
  
  prob_rf_telo <- randomForest:::predict.randomForest(rf_telo, test_data, type = "prob")[,2]
  prob_rf_age <- randomForest:::predict.randomForest(rf_age, test_data, type = "prob")[,2]
  
  roc_rf_telo <- roc(test_data$extub_status, prob_rf_telo)
  roc_rf_age <- roc(test_data$extub_status, prob_rf_age)
  
  rf_roc_test <- roc.test(roc_rf_telo, roc_rf_age, method = "delong")
  
  pdf(paste0("10_extubation/roc_age_vs_telomere_rf_", extub, ".pdf"))
  plot(roc_rf_age, col="#FB8072", main = "Extubation ROC")
  plot(roc_rf_telo, col="#80B1D3", add = TRUE)
  legend("bottomright", legend=c("Age model", "Telomere model"), 
         col=c("#FB8072", "#80B1D3"), lwd=2)
  text(0.9, 1, paste0("P-value = ", signif(rf_roc_test[["p.value"]], 3) ))
  dev.off()
  
  # 2. gam
  gam_age <- gam(extub_status ~ s(Age, k = 5) + Gender + PneumoniaType_CAP + 
                   RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + 
                   SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + 
                   SOFA_24h, data = train_data, family = "binomial")
  gam_telo <- gam(extub_status ~ s(TelomereLength_qmotif, k = 5) + Gender + PneumoniaType_CAP + 
                    RespiratorySupport_24h + Immunosuppression + MyocardialInfarction + 
                    SolidTumor + HM + ModerateToSevereRenalImpairment + MildRenalImpairment + 
                    SOFA_24h, data = train_data, family = "binomial")
  
  prob_gam_age <- mgcv:::predict.gam(gam_age, test_data, type = "response")
  prob_gam_telo <- mgcv:::predict.gam(gam_telo, test_data, type = "response")
  
  roc_gam_age <- roc(test_data$extub_status, prob_gam_age)
  roc_gam_telo <- roc(test_data$extub_status, prob_gam_telo)
  
  rf_gam_test <- roc.test(roc_gam_age, roc_gam_telo, method = "delong")
  
  pdf(paste0("10_extubation/roc_age_vs_telomere_gam_", extub, ".pdf"))
  plot(roc_gam_age, col="#FB8072", main = "Extubation ROC")
  plot(roc_gam_telo, col="#80B1D3", add = TRUE)
  legend("bottomright", legend=c("Age model", "Telomere model"), 
         col=c("#FB8072", "#80B1D3"), lwd=2)
  text(0.9, 1, paste0("P-value = ", signif(rf_gam_test[["p.value"]], 3) ))
  dev.off()
  
  # 3. xgboost
  train_data_xg <- train_data %>%
    mutate(Gender = as.numeric(Gender) - 1,
           extub_status = as.numeric(extub_status) -1)
  test_data_xg <- test_data %>%
    mutate(Gender = as.numeric(Gender) - 1,
           extub_status = as.numeric(extub_status) -1)
  
  train_matrix <- as.matrix(train_data_xg %>% select(!Age))
  test_matrix <- as.matrix(test_data_xg %>% select(!Age))
  xgb_telo <- xgboost(data = train_matrix[, colnames(train_matrix) != "extub_status"], 
                      label = train_matrix[, "extub_status"], 
                      nrounds = 100, objective = "binary:logistic", verbose = 0)
  prob_xgb_telo <- predict(xgb_telo, test_matrix[, colnames(train_matrix) != "extub_status"])
  roc_xgb_telo <- roc(test_matrix[, "extub_status"], prob_xgb_telo)
  
  train_matrix <- as.matrix(train_data_xg %>% select(!TelomereLength_qmotif))
  test_matrix <- as.matrix(test_data_xg %>% select(!TelomereLength_qmotif))
  xgb_age <- xgboost(data = train_matrix[, colnames(train_matrix) != "extub_status"], 
                     label = train_matrix[, "extub_status"], 
                     nrounds = 100, objective = "binary:logistic", verbose = 0)
  prob_xgb_age <- predict(xgb_age, test_matrix[, colnames(train_matrix) != "extub_status"])
  roc_xgb_age <- roc(test_matrix[, "extub_status"], prob_xgb_age)
  
  rf_xgb_test <- roc.test(roc_gam_age, roc_gam_telo, method = "delong")
  
  pdf(paste0("10_extubation/roc_age_vs_telomere_xgb_", extub, ".pdf"))
  plot(roc_xgb_age, col="#FB8072", main = "Extubation ROC")
  plot(roc_xgb_telo, col="#80B1D3", add = TRUE)
  legend("bottomright", legend=c("Age model", "Telomere model"), 
         col=c("#FB8072", "#80B1D3"), lwd=2)
  text(0.9, 1, paste0("P-value = ", signif(rf_xgb_test[["p.value"]], 3) ))
  dev.off()
}

# ===================== multi-omics prediction model ========================


