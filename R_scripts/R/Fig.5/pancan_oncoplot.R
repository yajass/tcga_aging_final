library(TCGAbiolinks)
library(maftools)
library(dplyr)

# Load data
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
cancer_code <- c("BRCA","COAD","THCA","UCEC","LGG","OV")
# Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
          "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
          "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")
survival_data$Quartiles <- ifelse(survival_data$Quartiles == 1, "Young",
                                  ifelse(survival_data$Quartiles == 4, "Old", NA))
survival_data <- survival_data[!is.na(survival_data$Quartiles), ] %>%
  dplyr::filter(type %in% cancer_code)
colnames(survival_data)[2] <- "Tumor_Sample_Barcode"

# download MAF and combine
for (i in 1:length(cancer_code)) {
  if (i ==1){
    tumor_maf <- GDCquery_Maf(tumor = cancer_code[i],
                              pipelines = "mutect2")
  }
  if (i != 1){
    downloadedMAF <- GDCquery_Maf(tumor = cancer_code[i],
                                  pipelines = "mutect2")
    tumor_maf <- plyr::rbind.fill(tumor_maf, downloadedMAF)
    rm(downloadedMAF)
  }
}

# Create MAF file for samples of interest after removing FLAGS
tumor_maf$Tumor_Sample_Barcode <- substr(tumor_maf$Tumor_Sample_Barcode,1,12)
tumor_maf <- tumor_maf %>% 
  dplyr::filter(!Hugo_Symbol %in% flags) %>% 
  dplyr::filter(Tumor_Sample_Barcode %in% survival_data$Tumor_Sample_Barcode)
colnames(survival_data)[3] <- "Tumor Type"
colnames(survival_data)[36] <- "Age"
MAF <- read.maf(maf = tumor_maf, clinicalData = survival_data)

age_col <- setNames(c(pal_jco()(2)[1], pal_jco()(2)[2]), c("Young", "Old"))
cancer_col <- setNames(RColorBrewer::brewer.pal(n = 7,name = 'Dark2'), cancer_code)


# plot
par(mar=c(1,1,1,1))
# postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/pancan_oncoplot.eps",
#            width = 7, height = 7)
tiff("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/pancan_oncoplot.tiff",
     width = 7, height = 7, units = "in", res = 320)
oncoplot(maf =MAF, 
         top = 10, 
         # bgCol = "white",
         legendFontSize = 1,sepwd_samples = 0,draw_titv = T,logColBar = T,
         sortByAnnotation = T,annotationFontSize = 1, fontSize = 1, titleFontSize = 1,
         annotationColor = list(Age = age_col, `Tumor Type` = cancer_col),
         gene_mar = 7, showTitle = FALSE,
         clinicalFeatures = c("Age","Tumor Type"), # which columns should be plotted
         annotationDat = survival_data) # data frame with metadata
dev.off()
