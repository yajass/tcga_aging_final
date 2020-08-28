library(gtable)
library(grid)
library(gtools)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(TCGAbiolinks)
library(maftools)

source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
drivers <- readxl::read_xlsx("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA Drivers/driver.genes.xlsx",
                             sheet = 2, trim_ws = TRUE, skip = 3)[,1:2]

cancer_code <- c("BRCA","COAD","THCA","UCEC","LGG","OV")

# collect results
################################################################################
# results DF
res <- setNames(data.frame(matrix(ncol = 9)),
                c("Hugo_Symbol", "Young", "Old", "pval",
                  "or", "ci.up", "ci.low", "adjPval", "cancercode"))
# Sample Sumamry DF
ssumary <- setNames(data.frame(matrix(ncol = 3)),  c("Cohort","SampleSize","cancercode"))
# Collet results
collectRes <- list(results = res, SampleSummary = ssumary)


for (i in 1:length(cancer_code)) {
  
  # download MAF
  tumor_maf <- GDCquery_Maf(tumor = cancer_code[i],
                            pipelines = "mutect2")
  
  # Download clinical data
  tumor_clinical <- survival_data %>% dplyr::filter(type == cancer_code[i])
  
  # Quartile age
  colnames(tumor_clinical)[2] <- "Tumor_Sample_Barcode"
  tumor_clinical$Quartiles <- ifelse(tumor_clinical$Quartiles == 1, "Young",
                                     ifelse(tumor_clinical$Quartiles == 4, "Old", NA))
  tumor_clinical <- tumor_clinical[!is.na(tumor_clinical$Quartiles), ]
  
  # Get young and old barcodes
  young_samples <- tumor_clinical$Tumor_Sample_Barcode[tumor_clinical$Quartiles == "Young"]
  old_samples <- tumor_clinical$Tumor_Sample_Barcode[tumor_clinical$Quartiles == "Old"]
  
  # Subset the downloaded maf
  tumor_maf$Tumor_Sample_Barcode <- substr(tumor_maf$Tumor_Sample_Barcode,1,12)
  young_maf <- tumor_maf[tumor_maf$Tumor_Sample_Barcode %in% young_samples, ]
  old_maf <- tumor_maf[tumor_maf$Tumor_Sample_Barcode %in% old_samples, ]
  
  # Setup results
  if (i == 1){
    overall_young_maf <- young_maf
    overall_old_maf <- old_maf
    young_clinical <- tumor_clinical[tumor_clinical$Quartiles == "Young", ]
    old_clinical <- tumor_clinical[tumor_clinical$age_at_diagnosis == "Old", ]
  }
  if (i != 1){
    overall_young_maf <- plyr::rbind.fill(overall_young_maf, young_maf)
    overall_old_maf <- plyr::rbind.fill(overall_old_maf, old_maf)
    young_clinical <- plyr::rbind.fill(young_clinical, tumor_clinical[tumor_clinical$Quartiles == "Young", ])
    old_clinical <- plyr::rbind.fill(old_clinical, tumor_clinical[tumor_clinical$Quartiles == "Old", ])
  }
}

overall_old_maf$Age <- "Old"
overall_young_maf$Age <- "Young"
combined_maf <- rbind(overall_young_maf, overall_old_maf)


xx <- combined_maf %>%
  dplyr::filter(Hugo_Symbol %in% c("TP53","PIK3CA","PTEN","APC","BRAF","ARID1A","IDH1","KMT2C","ATRX","RYR2")) %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::count(Tumor_Sample_Barcode, Age, Hugo_Symbol)

pvals <- numeric(length = length(unique(xx$Hugo_Symbol)))
for (i in 1:length(unique(xx$Hugo_Symbol))) {
  df <- xx %>% dplyr::filter(Hugo_Symbol == unique(xx$Hugo_Symbol)[i])
  tab <- setNames(data.frame(matrix(c(sum(df$n[df$Age == "Young"] == 1),
                                      sum(df$n[df$Age == "Old"] == 1),
                                      sum(df$n[df$Age == "Young"] > 1),
                                      sum(df$n[df$Age == "Old"] > 1)),
                                    nrow = 2, ncol = 2),
                             row.names = c("Young", "Old")),
                  c("Single", "Multi"))
  print(unique(xx$Hugo_Symbol)[i])
  print(tab)
  pvals[i]<- fisher.test(t(tab))$p.value
  names(pvals)[i] <- unique(xx$Hugo_Symbol)[i]
}

p.adjust(pvals, "fdr")
