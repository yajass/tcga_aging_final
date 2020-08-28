library(TCGAutils)

source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
folder <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA_CIBERSORT_ABS_WCM/"
files <- list.files(folder)
cancer_codes <- c("BRCA","THCA","LGG","COAD","OV","UCEC")

resDF <- setNames(data.frame(matrix(ncol = 4, nrow = length(cancer_codes) + 1),
                             row.names = c(cancer_codes,"COMBINED")), c("Beta", "ci_lower", "ci_upper","pval"))

# regression by tumor type
################################################################
for (i in 1:length(cancer_codes)) {
  dat0 <- data.table::fread(paste0(folder, files[stringr::str_detect(files, cancer_codes[i])]))[,-c(24:27)] %>%
    tibble::column_to_rownames(., "TCGA_barcode")
  
  # Preproc
  dat0 <- dat0[(TCGAbiospec(rownames(dat0))$sample_definition == "Primary Solid Tumor" |
                  TCGAbiospec(rownames(dat0))$sample_definition == 
                  "Primary Blood Derived Cancer - Peripheral Blood"), ]
  if (sum(duplicated(TCGAbiospec(rownames(dat0))$submitter_id)) > 0){
    dat0 <- dat0[-c(which(duplicated(TCGAbiospec(rownames(dat0))$submitter_id))), ]
  }
  rownames(dat0) <- substring(rownames(dat0), 1, 12)
  
  # Get scores and merge
  immune_scores <- data.frame(bcr_patient_barcode = rownames(dat0),
                              lymphocyte_score = rowSums(dat0[, c(1:12)]))
  surv1 <- survival_data[, c(2, 3, 36)]
  res <- merge(surv1, immune_scores) %>%
    dplyr::filter(Quartiles %in% c(1, 4)) %>%
    dplyr::mutate(Quartiles = factor(ifelse(Quartiles == 1, "Young", "Old"),
                                     levels = c("Young", "Old")))
  
  # regression + write results
  reg_mod <- glm(Quartiles ~ lymphocyte_score, data = res, family = "binomial")
  resDF$Beta[i] <- summary(reg_mod)$coefficients[2,1]
  resDF$pval[i] <- summary(reg_mod)$coefficients[2,4]
  resDF$ci_lower[i] <- confint(reg_mod)[2,1]
  resDF$ci_upper[i] <- confint(reg_mod)[2,2]
}

# Combined cancer regression
################################################################
for (i in 1:length(cancer_codes)) {
  if (i == 1){
    dat0 <- data.table::fread(paste0(folder, files[stringr::str_detect(files, cancer_codes[i])]))[,-c(24:27)] %>%
      tibble::column_to_rownames(., "TCGA_barcode")
  }
  if (i > 1){
    dat0 <- rbind(dat0, data.table::fread(paste0(folder, files[stringr::str_detect(files, cancer_codes[i])]))[,-c(24:27)] %>%
                    tibble::column_to_rownames(., "TCGA_barcode"))
  }
}
# Preproc
dat0 <- dat0[(TCGAbiospec(rownames(dat0))$sample_definition == "Primary Solid Tumor" |
                TCGAbiospec(rownames(dat0))$sample_definition == 
                "Primary Blood Derived Cancer - Peripheral Blood"), ]
if (sum(duplicated(TCGAbiospec(rownames(dat0))$submitter_id)) > 0){
  dat0 <- dat0[-c(which(duplicated(TCGAbiospec(rownames(dat0))$submitter_id))), ]
}
rownames(dat0) <- substring(rownames(dat0), 1, 12)
# Get scores and merge
immune_scores <- data.frame(bcr_patient_barcode = rownames(dat0),
                            lymphocyte_score = rowSums(dat0[, c(1:12)]))
surv1 <- survival_data[, c(2, 3, 36)]
res <- merge(surv1, immune_scores) %>%
  dplyr::filter(Quartiles %in% c(1, 4)) %>%
  dplyr::mutate(Quartiles = factor(ifelse(Quartiles == 1, "Young", "Old"),
                                   levels = c("Young", "Old")))
# regression + write results
reg_mod <- glm(Quartiles ~ lymphocyte_score, data = res, family = "binomial")
resDF[7, ] <- c(summary(reg_mod)$coefficients[2,1],confint(reg_mod)[2,1],confint(reg_mod)[2,2],summary(reg_mod)$coefficients[2,4])


# Plot
################################################################
resDF <- tibble::rownames_to_column(resDF, "cancer") %>%
  dplyr::arrange(Beta) %>%
  dplyr::mutate(FDR = p.adjust(pval, "fdr"), 
                Age = factor(ifelse(FDR < 0.05 & Beta > 0, "Old",
                                    ifelse(FDR < 0.05 & Beta < 0, "Young", "Neither")),
                             levels = c("Young","Old","Neither")),
                cancer = factor(cancer, levels = rev(c("LGG","COAD","OV","THCA","BRCA","UCEC","COMBINED"))))
ggplot(resDF, aes(x = Beta, xmin = ci_lower, xmax = ci_upper, y = cancer, color = Age)) +
  geom_point(shape = 18, size = 6) + 
  geom_errorbar(width = 0.5, size = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Lymphocyte Score", y = NULL, color = "Association") +
  xlim(c(-30, 30)) +
  scale_color_jco() +
  theme_pubclean(25)
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S2/lymphocyte_score_cibersort_absolute_FDR_05.eps")
