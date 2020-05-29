library(TCGAbiolinks)
library(ggpubr)
library(dplyr)

# Data from https://www.nature.com/articles/nature12113
# Integrated genomic characterization of endometrial carcinoma
# Integrative cluster based on CNV, RNASeq, mutations and DNA methylation
data <- readxl::read_xls("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA UCEC/datafile.S1.1.KeyClinicalData.xls")

# Get clinical data and quartile
tumor_clinical <- GDCquery_clinic(project = "TCGA-UCEC",
                                  type = "Clinical")
tumor_clinical$primary_diagnosis <- "UCEC"
colnames(tumor_clinical)[1] <- "bcr_patient_barcode"
tumor_clinical$age_at_diagnosis <- ntile(tumor_clinical$age_at_diagnosis, 4)
tumor_clinical$age_at_diagnosis <- ifelse(tumor_clinical$age_at_diagnosis == 1, "Young",
                                          ifelse(tumor_clinical$age_at_diagnosis == 4, "Old", NA))
tumor_clinical <- tumor_clinical[!is.na(tumor_clinical$age_at_diagnosis), ]

# Integrate datasets
young_samples <- tumor_clinical$bcr_patient_barcode[tumor_clinical$age_at_diagnosis == "Young"]
old_samples <- tumor_clinical$bcr_patient_barcode[tumor_clinical$age_at_diagnosis == "Old"]
data$Quartile <- NA
data[data$bcr_patient_barcode %in% young_samples, ]$Quartile <- "Young"
data[data$bcr_patient_barcode %in% old_samples, ]$Quartile <- "Old"
data <- na.omit(data)
data$IntegrativeCluster

# Dataframe to plot
df <- data.frame(table(data$Quartile, data$IntegrativeCluster))
df$Var1 <- factor(df$Var1, levels = c("Young", "Old"))

# visualize
ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_fill_brewer() +
  theme_pubclean(base_size = 35) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(fill = "UCEC \nIntegrative \nCluster")
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/uecc_integrative_clusters.eps",
       dpi = 320, height = 11, width = 6)
