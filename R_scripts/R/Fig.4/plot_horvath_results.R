library(gtools)
library(RColorBrewer)
library(ggpubr)
library(rstatix)

source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
PATH <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/HG 19 Meth/Horvath/Output/"
files <- list.files(PATH)[endsWith(list.files(PATH), ".csv")]
normals <- files[startsWith(files, "Normal")]
files <- files[!startsWith(files, "Normal")]

# get normals
filtered_age <- survival_data
datN <- read.csv(paste0(PATH, normals), stringsAsFactors = F, header = T, sep = ",")
datN$SampleID <- gsub('\\.','-', datN$SampleID)
filtered_age <- filtered_age[filtered_age$bcr_patient_barcode %in% datN$Sample, ]
filtered_age <- filtered_age[!duplicated(filtered_age$bcr_patient_barcode), ]
mergedNormal <- merge(filtered_age, datN, by.x = "bcr_patient_barcode", by.y = "SampleID")
mergedNormal$Source = "Normal"

combined_data <- data.frame()
for (i in 1:length(files)) {
  dat0 <- read.csv(paste0(PATH, files[i]), stringsAsFactors = F, header = T, sep = ",")
  combined_data <- rbind(combined_data, dat0)
}
combined_data$SampleID <- gsub('\\.','-', combined_data$SampleID)
combined_data$SampleID <- substring(combined_data$SampleID,1,12)
combined_data <- rbind(combined_data, datN)
dat0 <- merge(survival_data, combined_data, by.x = "bcr_patient_barcode", by.y = "SampleID")
dat0$Source = "Cancer"
dat0 <- rbind(dat0, mergedNormal)
dat0$Quartiles <- factor(dat0$Quartiles, levels = c(1,2,3,4))

## Plot horvath results
################################################################################################################
# Remove samples with < 0.8 cor with gold standard and Q2-Q3 samples
dat0 <- dat0[dat0$type %in% c("BRCA","COAD","LGG","THCA","UCEC"), ]
df_save <- data.frame(Barcode = dat0$bcr_patient_barcode,
                      Age = dat0$age_at_initial_pathologic_diagnosis,
                      DNAmAge = dat0$DNAmAge,
                      DNAmmeanAbsDifferenceSampleVSgoldstandard = dat0$meanAbsDifferenceSampleVSgoldstandard,
                      DNAmAgeAccDiff = dat0$AgeAccelerationDiff,
                      Status = dat0$Source)
write.xlsx(df_save, "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - DNAm Age.xlsx")
dat0 <- dat0[dat0$corSampleVSgoldstandard >= 0.8, ]

pdf("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.4/Epigenetic_Age.pdf",
    width = 9, height = 9)
ggscatter(
  data = dat0,
  x = "Age",
  y = "DNAmAge",
  alpha = 1/3,
  color = 'Source',
  # color = "dodgerblue2",
  palette = "aaas",
  title = "Epigenetic Age", 
  add = "reg.line",
  fullrange = TRUE
) +
  labs(title = NULL) +
  theme_pubr(base_size = 35) +
  theme(legend.direction = "vertical",
        legend.justification = c("left", "top"),
        legend.position = c(.01, .99)) +
  stat_cor(aes(color = Source), method = "pearson", show.legend = FALSE,
           inherit.aes = TRUE, label.x.npc = 0.6, label.y.npc = 0.95, size = 7)
dev.off()

# Compare groups
dat0$Quartiles <- ifelse(dat0$Quartiles == 1, "Young", 
                         ifelse(dat0$Quartiles == 4,"Old","Middle Aged"))
dat0$Quartiles <- factor(dat0$Quartiles, levels = c("Young", "Middle Aged","Old"))

stat.test <- dat0 %>%
  t_test(AgeAccelerationDiff ~ Quartiles) %>%
  mutate(y.position = 40)
stat.test <- stat.test[-c(1,3), ]

ggline(
  data = dat0,
  x = "Quartiles",
  y = "AgeAccelerationDiff",
  color = "type",
  # color = "Source",
  add = "mean_ci",
  linetype = 2,
  size = 2,
  palette = "jama",
  xlab = FALSE
) +
  theme_pubr(base_size = 35) +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        legend.box.spacing = unit(0, "line")) +
  labs(color = "Cancer", title = NULL, y = "DNAm AgeAccelerationDiff") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", size = 10, bracket.size = 1.5, position = position_dodge(0.8))
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.4/Quartile_TumorType_DNAm_Age_Acceleration.eps", dpi = 320,
       width = 9, height = 9)

stat.test1 <- dat0 %>%
  group_by(Quartiles) %>%
  t_test(AgeAccelerationDiff ~ Source) %>%
  mutate(y.position = c(2,-1,-11))
stat.test1$stars <- stars.pval(stat.test1$p)

dat0$Quartiles <- factor(dat0$Quartiles, levels = c("Young", "Middle Aged", "Old"))
ggline(
  data = dat0,
  x = "Quartiles",
  y = "AgeAccelerationDiff",
  add = "mean_ci",
  color = "Source",
  linetype = 2,
  size = 2,
  xlab = FALSE,
  palette = "aaas"
) +
  theme_pubr(base_size = 35, 
             legend = "none") +
  labs(y = "DNAm AgeAccelerationDiff") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  stat_compare_means(method = "t.test", aes(`cols` = Source), label = "p.signif",
                     label.y = c(2,-1,-11), size = 10, hide.ns = TRUE)
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.4/Quartile_Overall_DNAm_Age_Acceleration.eps", dpi = 320,
       width = 9, height = 9)


# ggline(
#   data = dat0,
#   x = "Quartiles",
#   y = "AgeAccelerationDiff",
#   add = "mean_ci",
#   color = "Source",
#   linetype = 2,
#   xlab = FALSE,
#   palette = "jco",
#   facet.by = "type"
# ) +
#   stat_compare_means(method = "t.test", aes(`cols` = Source), label = "p.signif",
#                      label.y = c(2,-5,-13))
# 
# 
# dat1 <- dat0
# dat1$CASE_ID <- substr(dat1$CASE_ID,1,12)
# dat1 <- dat1[dat1$CASE_ID %in% dat1$CASE_ID[dat1$Source == "Normal"], ]
# ggline(
#   data = dat1,
#   x = "Quartiles",
#   y = "AgeAccelerationDiff",
#   add = "mean_ci",
#   color = "Source",
#   linetype = 2,
#   xlab = FALSE,
#   palette = "jco",
#   facet.by = "type",
#   scales = "free",
#   title = "Paired Age Acceleration"
# ) +
#   stat_compare_means(method = "t.test", comparisons = list(c("Young", "Old")), tip.length = 0.01) +
#   stat_compare_means(method = "t.test", aes(`cols` = Source), label = "p.signif")
# ggsave("./Final/Figures/Figure 4/Paired_AgeAccelDiff.jpg", dpi = 320, width = 13)