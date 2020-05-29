library(maftools)
library(plyr)
library(reshape2)
library(dplyr)
library(ggpubr)

# read clinical
clinical <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_clinical_patient copy.txt",
                     header = T, stringsAsFactors = F, sep = "\t")
clinical <- clinical[-c(1,2),]
clinical$Age <- ifelse(clinical$Age.at.Diagnosis <= 49, "Young",
                       ifelse(clinical$Age.at.Diagnosis >= 68, "Old", "Middle.Aged"))

# total maf
metabric_mut <- read.maf(maf = "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_mutations_extended.txt",
                         isTCGA = F, clinicalData = clinical,
                         gisticScoresFile = "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_CNA.txt")

# young vs old
young_samples <- clinical$Tumor_Sample_Barcode[clinical$Age == "Young"]
old_samples <- clinical$Tumor_Sample_Barcode[clinical$Age == "Old"]
young_maf <- subsetMaf(maf = metabric_mut, query = 'Tumor_Sample_Barcode %in% young_samples')
old_maf <- subsetMaf(maf = metabric_mut, query = 'Tumor_Sample_Barcode %in% old_samples')

# Co-oncoplot
gene_list <- c("TP53", "PIK3CA", "KMT2C", "SF3B1")
postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S5/METABRIC_ONCOP.eps",
           width = 10, height = 4, horizontal = F)
coOncoplot(m1 = young_maf, m2 = old_maf, m1Name = "Young", m2Name = "Old",
           genes = gene_list, sepwd_genes1 = 0, sepwd_genes2 = 0,
           titleFontSize = 2, geneNamefont = 1, sepwd_samples1 = 0, sepwd_samples2 = 0)
dev.off()

# lollipop
postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S5/METABRIC_lollipop.eps",
           width = 10, height = 4, horizontal = F)
lollipopPlot2(m1 = young_maf, m2 = old_maf,
              gene = "TP53", AACol1 = "HGVSp", AACol2 = "HGVSp", 
              m1_name = "Young", m2_name = "Old",
              legendTxtSize = 2, pointSize = 2, domainLabelSize = 2, labPosSize = 2)
dev.off()



