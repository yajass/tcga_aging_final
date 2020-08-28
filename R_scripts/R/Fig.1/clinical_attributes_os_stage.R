library(ggsci)
library(data.table)
library(dplyr)
library(openxlsx)
library(survival)
library(survminer)

survival_data <- read.xlsx("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Survival_Data/TCGA-CDR-SupplementalTableS1.xltx")
survival_data <- survival_data[!is.na(survival_data$age_at_initial_pathologic_diagnosis), ]
survival_data$type[survival_data$type == "READ"] <- "COAD"

# create cancer type based age quartiles
setDT(survival_data)[,Subtype_Quartile := cut(age_at_initial_pathologic_diagnosis, quantile(age_at_initial_pathologic_diagnosis, probs = 0:4/4),
                                              labels = FALSE, include.lowest = TRUE), by = type]
setDT(survival_data)[,Quartiles := cut(age_at_initial_pathologic_diagnosis, quantile(age_at_initial_pathologic_diagnosis, probs = 0:4/4),
                                              labels = FALSE, include.lowest = TRUE), by = type]
survival_data$Subtype_Quartile <- paste0("Q",survival_data$Subtype_Quartile)



# remove NA OS instances/time
survival_data <- survival_data[!is.na(survival_data$OS),]
survival_data <- survival_data[!is.na(survival_data$OS.time),]

# loop to calculate OS survival
cancer_code <- unique(survival_data$type)
results_df <- data.frame(matrix(nrow = length(cancer_code),
                                ncol = 7))
colnames(results_df) <- c("Cancer", "Surv_N", "Surv_events", "hazard_ratio",
                          "lower_ci", "upper_ci", "Surv_pval")

for (i in 1:length(cancer_code)) {
  cancerDF <- survival_data[survival_data$type == cancer_code[i], ]
  
  # cancerDF <- cancerDF[cancerDF$Quartiles == "Q1" | cancerDF$Quartiles == "Q4", ]
  
  surv_result <- summary(coxph(Surv(time= OS.time ,event = OS ) ~ age_at_initial_pathologic_diagnosis, 
                               data = cancerDF))
  
  results_df$Surv_pval[i] <- coef(surv_result)[1,5]
  results_df$hazard_ratio[i] = round(surv_result$coefficients[2],5)
  results_df$Surv_N[i] <- surv_result$n
  results_df$Surv_events[i] <- surv_result$nevent
  results_df$lower_ci[i] <- surv_result$conf.int[3]
  results_df$upper_ci[i] <- surv_result$conf.int[4]
  results_df$Cancer[i] <- cancer_code[i]
  
  # survfit(Surv(time= OS.time,event = OS )~ Quartiles,
  #         data = cancerDF)
  # 
  # ggsurvplot(survfit(Surv(time= OS.time,event = OS )~ Quartiles,
  #                    data = cancerDF),data = cancerDF)
}

# FDR Correct
results_df$Surv_pval <- p.adjust(results_df$Surv_pval, "fdr")
sig_cancers <- results_df$Cancer[results_df$Surv_pval<0.01]
results_df <- mutate(results_df, sig = ifelse(results_df$Surv_pval<0.01, "Sig.", "Not Sig."))
results_df <- results_df %>% arrange(sig, hazard_ratio)
levs <- results_df$Cancer
results_df$Cancer <- factor(results_df$Cancer, levels = levs)
results_df$sig <- ifelse(results_df$sig == "Sig.", "< 0.01", "≥ 0.01")

ggplot(data = results_df, aes(x = Cancer, y = hazard_ratio, ymin = lower_ci, ymax = upper_ci)) + 
  geom_pointrange(aes(col = sig), shape = 18, size = 2) +
  ylim(c(0.7,1.3)) +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, col = sig), width = 0.75, lwd = 1.5) +
  geom_hline(yintercept=1, lty=2, size = 2) +
  scale_color_aaas()+
  labs(y = "Hazard Ratio", x = NULL, title = NULL, color = "FDR") +
  theme_pubr(base_size = 35) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 22),
        legend.direction = "vertical",
        legend.justification = c("left", "top"),
        legend.position = c(.05, .95),
        axis.text.x = element_text(size = 30))
# ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.1/overall_survival_hr.jpg", dpi = 320, height = 11, width = 11)
# write.xlsx(x = results_df,
#            file = "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - Overall Survival.xlsx")


# write table for quartile age limits
# survival_data %>%
#   dplyr::mutate(Age_Group = ifelse(Subtype_Quartile == "Q1", "Young",
#                                    ifelse(Subtype_Quartile == "Q4", "Old", "Middle Aged"))) %>%
#   dplyr::select(type, age_at_initial_pathologic_diagnosis,Age_Group) %>%
#   dplyr::filter(Age_Group %in% c("Young", "Old"),
#                 type %in% as.character(results_df$Cancer[results_df$sig == "< 0.01"])) %>%
#   dplyr::group_by(type, Age_Group) %>%
#   dplyr::summarise(Minimum = min(age_at_initial_pathologic_diagnosis),
#                    Maximum = max(age_at_initial_pathologic_diagnosis),
#                    Mean = mean(age_at_initial_pathologic_diagnosis),
#                    Median = median(age_at_initial_pathologic_diagnosis)) %>%
#   write.xlsx(., file = "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - Age Quartiles.xlsx")


# plot KM for types of interest
plt.surv <- survival_data %>%
  dplyr::filter(type %in% c("BRCA", "COAD", "LGG", "THCA", "OV", "UCEC"),
                Quartiles %in% c(1,4)) %>%
  dplyr::mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))
fit <- survfit(Surv(OS.time, OS) ~ Quartiles + type, data = plt.surv)
ggsurvplot(fit, linetype = rep(c("solid", "dashed"), 6),
           palette = rep(ggsci::pal_aaas()(6), each = 2), legend = "none",
           size = 2.5, font.x =  16, font.y =  16, font.ticks = 12, theme = "transparent")
# ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/graphical.abstract/survival_curves.jpg", bg = "transparent")

# check if age is associated with tumor stage
# need to collate tumor stage columns
surv_data_1 <- survival_data[survival_data$Subtype_Quartile %in% c("Q1", "Q4") &
                               survival_data$type %in% c('BRCA','COAD','THCA','UCEC',
                                                         'OV','LGG'), ]
surv_data_1$ajcc_pathologic_tumor_stage <- ifelse(startsWith(surv_data_1$ajcc_pathologic_tumor_stage, "Stage IV"), "Stage IV",
                                                  ifelse(startsWith(surv_data_1$ajcc_pathologic_tumor_stage, "Stage III"), "Stage III",
                                                         ifelse(startsWith(surv_data_1$ajcc_pathologic_tumor_stage, "Stage II"), "Stage II",
                                                                ifelse(startsWith(surv_data_1$ajcc_pathologic_tumor_stage, "Stage I"), "Stage I",
                                                                       surv_data_1$ajcc_pathologic_tumor_stage))))
surv_data_1$clinical_stage <- ifelse(startsWith(surv_data_1$clinical_stage, "Stage IV"), "Stage IV",
                                     ifelse(startsWith(surv_data_1$clinical_stage, "Stage III"), "Stage III",
                                            ifelse(startsWith(surv_data_1$clinical_stage, "Stage II"), "Stage II",
                                                   ifelse(startsWith(surv_data_1$clinical_stage, "Stage I"), "Stage I",
                                                          surv_data_1$clinical_stage))))
surv_data_1$ajcc_pathologic_tumor_stage[surv_data_1$type %in% c('OV','SARC','LGG','UCEC')] <- 
  surv_data_1$clinical_stage[surv_data_1$type %in% c('OV','SARC','LGG','UCEC')]
surv_data_1$Subtype_Quartile <- factor(surv_data_1$Subtype_Quartile, levels = c("Q4", "Q1"))

for (i in 1:length(unique(surv_data_1$type))) {
  subsetData <- surv_data_1[surv_data_1$type == unique(surv_data_1$type)[i], ]
  
  if (i == 1){
    fisher.res <- subsetData %>%
      dplyr::select(type, Subtype_Quartile, ajcc_pathologic_tumor_stage) %>%
      dplyr::group_by(type) %>%
      {table(.$ajcc_pathologic_tumor_stage, .$Subtype_Quartile)} %>%
      rstatix::row_wise_fisher_test(., p.adjust.method = "fdr", detailed = TRUE) %>%
      dplyr::mutate(type = unique(subsetData$type))
  }
  
  if (i != 1){
    fisher.res <- rbind(fisher.res, subsetData %>%
                          dplyr::select(type, Subtype_Quartile, ajcc_pathologic_tumor_stage) %>%
                          dplyr::group_by(type) %>%
                          {table(.$ajcc_pathologic_tumor_stage, .$Subtype_Quartile)} %>%
                          rstatix::row_wise_fisher_test(., p.adjust.method = "fdr", detailed = TRUE) %>%
                          dplyr::mutate(type = unique(subsetData$type)))
  }
}

fisher.res$estimate <- log10(fisher.res$estimate)
fisher.res$conf.low <- log10(fisher.res$conf.low)
fisher.res$conf.high <- log10(fisher.res$conf.high)
fisher.res$p.adj.signif <- ifelse(fisher.res$p.adj.signif == "ns", "≥ 0.05", "< 0.05")
fisher.res$p.adj.signif <- factor(fisher.res$p.adj.signif, levels = c("< 0.05", "≥ 0.05"))

ggplot(data = fisher.res, aes(x = group, y = estimate, ymin = conf.low, ymax = conf.high)) + 
  geom_pointrange(aes(col = p.adj.signif), shape = 18, size = 2) +
  facet_wrap(vars(type)) +
  ylim(c(-5,5)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, col =  p.adj.signif), width = 0.75, lwd = 1.5) +
  geom_hline(yintercept=0, lty=2, size = 2) +
  scale_color_aaas()+
  labs(y = "log Odds Ratio", x = "", title = "Tumor Stage", color = "FDR") +
  theme_pubclean(base_size = 15) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15))
# ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S1/age_stage_or.jpg", dpi = 320)
rm(list=c("cancer_code","cancerDF", "fisher.res", "i", "levs", "results_df", "subsetData", "surv_data_1", "surv_result"))
