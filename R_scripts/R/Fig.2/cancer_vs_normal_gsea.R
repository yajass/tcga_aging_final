###########################################################################################################
###########################################################################################################
####                                     Load packages and essentials                                  ####
###########################################################################################################
###########################################################################################################

library(pheatmap)
library(dplyr)
library(tibble)
library(fgsea)
library(edgeR)
library(limma)

pathways.hallmark <- gmtPathways("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GSEA_pathways/h.all.v7.0.symbols.gmt")

# Converts ENSEMBL IDs (row names) to unique gene symbols
getGeneNamesFromENSELBL <- function(countData){
  library(org.Hs.eg.db)
  # rownames(countData) <- gsub("\\..*","",rownames(countData))
  countData[,ncol(countData)+1] <- mapIds(org.Hs.eg.db, keys = gsub("\\..*","",rownames(countData)),
                                          keytype = "ENSEMBL", column="SYMBOL", "first")
  countData <- na.omit(countData)
  countData <- countData[-which(duplicated(countData[,ncol(countData)])),]
  rownames(countData) <- countData[,ncol(countData)]
  countData <- countData[,-ncol(countData)]
  return(countData)
}


###########################################################################################################
###########################################################################################################
####                                        GTEx Breast FGSEA                                          ####
###########################################################################################################
###########################################################################################################

# Load data
###########################################################################################################
meta1 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t")
meta2 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t")
gtex_data <- read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/Tissues/GTEx_Breast - Mammary Tissue.txt")

# Setup for DE
###########################################################################################################
# match column names with meta2
colnames(gtex_data) <- gsub('\\.', '-', colnames(gtex_data))
colnames(gtex_data) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", colnames(gtex_data))

# convert transcript rownames to ensembl
gtex_data <- gtex_data[!duplicated(gsub("\\..*","",rownames(gtex_data))), ]
rownames(gtex_data) <- gsub("\\..*","",rownames(gtex_data))

age_meta <- meta2[meta2$SUBJID %in% colnames(gtex_data), ]
age_meta$AGE <- ifelse(age_meta$AGE == "20-29", "Young",
                       ifelse(age_meta$AGE == "30-39", "Young",
                              ifelse(age_meta$AGE == "40-49", "Young",
                                     ifelse(age_meta$AGE == "50-59", NA, "Old"))))
age_meta <- na.omit(data.frame(SUBJID = age_meta$SUBJID, AGE = age_meta$AGE))
gtex_data <- gtex_data[, colnames(gtex_data) %in% age_meta$SUBJID]

# Verify this is true
all(colnames(gtex_data) == age_meta$SUBJID)

# Limma
###########################################################################################################
age_meta$AGE <- factor(age_meta$AGE, levels = c("Old","Young"))
condition <- c("Old","Young")
comparison <- "Old-Young"
design <- model.matrix(~0+ age_meta$AGE)
colnames(design) <- condition
contmatrix <- makeContrasts(as.character(comparison),levels=design)

gtex_data<-gtex_data[which(rowSums(cpm(gtex_data)> 2) >= 2) ,]
print(dim(gtex_data))
Data.voom <- voom(gtex_data,plot=TRUE)
fit <- lmFit(Data.voom,design)
fit2 <- contrasts.fit(fit,contmatrix)
fit2 <-eBayes(fit2)
tt <- topTable(fit2,n=Inf)
tt <- getGeneNamesFromENSELBL(tt)


# FGSEA
###########################################################################################################
set.seed(2019)
tt <- tt[order(tt$t, decreasing = T),]
ranks <- tt$t
names(ranks) <- rownames(tt)
fgseaRes_GTx <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes_GTx$NES <- ifelse(fgseaRes_GTx$padj > 0.05, 0, fgseaRes_GTx$NES)


# Keep essentials
###########################################################################################################
rm(list=setdiff(ls(), c("fgseaRes_GTx", "getGeneNamesFromENSELBL", "pathways.hallmark")))



###########################################################################################################
###########################################################################################################
####                                        TCGA Normal Breast FGSEA                                   ####
###########################################################################################################
###########################################################################################################

# Download data
###########################################################################################################
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  legacy = FALSE,
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Solid Tissue Normal")
GDCdownload(query)
temp <- GDCprepare(query)
normal_data <- data.frame(assay(temp))
colnames(normal_data) <- gsub('\\.', '-', substr(colnames(normal_data), 1, 12))

# Survival data
###########################################################################################################
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
brca_samples <- survival_data$bcr_patient_barcode[survival_data$type == "BRCA"]

# Select data
###########################################################################################################
survival_data <- survival_data[survival_data$type %in% "BRCA", ]
survival_data <- survival_data[survival_data$Quartiles == 1 | 
                                 survival_data$Quartiles == 4]

normal_data <- normal_data[, colnames(normal_data) %in% survival_data$bcr_patient_barcode]
normal_data <- normal_data[, sort(colnames(normal_data))]
survival_data <- survival_data[survival_data$bcr_patient_barcode %in% colnames(normal_data),]

input <- data.frame(Patient = as.character(survival_data$bcr_patient_barcode),
                    Age = survival_data$Quartiles,
                    row.names = as.character(survival_data$bcr_patient_barcode))
input <- input[sort(rownames(input)), ]
input$Age <- ifelse(input$Age == 1, "Young", "Old")

# Verify this is true
all(input$Patient == colnames(normal_data))


# Limma
###########################################################################################################
input$Age <- factor(input$Age, levels = c("Old","Young"))
condition <- c("Old","Young")
comparison <- "Old-Young"
design <- model.matrix(~0+ input$Age)
colnames(design) <- condition
contmatrix <- makeContrasts(as.character(comparison),levels=design)

normal_data<-normal_data[which(rowSums(cpm(normal_data)> 2) >= 2) ,]
print(dim(normal_data))
Data.voom <- voom(normal_data,plot=TRUE)
fit <- lmFit(Data.voom,design)
fit2 <- contrasts.fit(fit,contmatrix)
fit2 <-eBayes(fit2)
tt <- topTable(fit2,n=Inf)
tt <- getGeneNamesFromENSELBL(tt)


# FGSEA
###########################################################################################################
set.seed(2019)
tt <- tt[order(tt$t, decreasing = T),]
ranks <- tt$t
names(ranks) <- rownames(tt)
fgseaRes_NAT <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes_NAT$NES <- ifelse(fgseaRes_NAT$padj > 0.05, 0, fgseaRes_NAT$NES)


# Keep essentials
###########################################################################################################
rm(list=setdiff(ls(), c("fgseaRes_NAT", "fgseaRes_GTx", "getGeneNamesFromENSELBL", "pathways.hallmark")))



###########################################################################################################
###########################################################################################################
####                                        TCGA Breast Cancer FGSEA                                   ####
###########################################################################################################
###########################################################################################################

# Survival data
###########################################################################################################
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
brca_samples <- survival_data$bcr_patient_barcode[survival_data$type == "BRCA"]

# Load data
###########################################################################################################
load("~/Documents/Elemento/AgeTCGA-master/DATA/NewTCGAcounts/Full Matrix.RData")
filtered_counts <- RNAseq_counts
library(TCGAutils)
filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Solid Tumor" |
                                        TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Blood Derived Cancer - Peripheral Blood")]
biospec_table <- TCGAbiospec(colnames(filtered_counts))
biospec_table <- biospec_table[order(biospec_table$plate,decreasing = TRUE), ]
filtered_counts <- filtered_counts[,colnames(filtered_counts)[order(biospec_table$plate, decreasing = TRUE)]]

if (sum(duplicated(biospec_table$submitter_id)) > 0){
  filtered_counts <- filtered_counts[, !duplicated(biospec_table$submitter_id)]
}

# change column names to match metadata
colnames(filtered_counts) <- substr(colnames(filtered_counts), 1, 12)


# Select data
###########################################################################################################
survival_data <- survival_data[survival_data$bcr_patient_barcode %in% brca_samples, ]
survival_data <- survival_data[survival_data$Quartiles %in% c(1,4), ]
common <- intersect(colnames(filtered_counts), survival_data$bcr_patient_barcode)
filtered_counts <- filtered_counts[,colnames(filtered_counts) %in%  common]
filtered_counts <- filtered_counts[, -which(!colnames(filtered_counts) %in% common)]
survival_data <- survival_data[survival_data$bcr_patient_barcode %in% common, ]

# Set up for limma
###########################################################################################################
filtered_counts <- filtered_counts[, sort(colnames(filtered_counts))]

input <- data.frame(Patient = as.character(survival_data$bcr_patient_barcode),
                    Age = survival_data$Quartiles,
                    row.names = as.character(survival_data$bcr_patient_barcode))
input <- input[sort(rownames(input)), ]
input$Age <- ifelse(input$Age == 1, "Young", "Old")

# Verify this is true
all(input$Patient == colnames(filtered_counts))


# Limma
###########################################################################################################
input$Age <- factor(input$Age, levels = c("Old","Young"))
condition <- c("Old","Young")
comparison <- "Old-Young"
design <- model.matrix(~0+ input$Age)
colnames(design) <- condition
contmatrix <- makeContrasts(as.character(comparison),levels=design)

filtered_counts<-filtered_counts[which(rowSums(cpm(filtered_counts)> 2) >= 2) ,]
print(dim(filtered_counts))
Data.voom <- voom(filtered_counts,plot=TRUE)
fit <- lmFit(Data.voom,design)
fit2 <- contrasts.fit(fit,contmatrix)
fit2 <-eBayes(fit2)
tt <- topTable(fit2,n=Inf)
tt <- getGeneNamesFromENSELBL(tt)


# FGSEA
###########################################################################################################
set.seed(2019)
tt <- tt[order(tt$t, decreasing = T),]
ranks <- tt$t
names(ranks) <- rownames(tt)
fgseaRes_BRCA <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes_BRCA$NES <- ifelse(fgseaRes_BRCA$padj > 0.05, 0, fgseaRes_BRCA$NES)


# Keep essentials
###########################################################################################################
rm(list=setdiff(ls(), c("fgseaRes_BRCA", "fgseaRes_NAT", "fgseaRes_GTx", "getGeneNamesFromENSELBL", "pathways.hallmark")))




###########################################################################################################
###########################################################################################################
####                                            METABRIC FGSEA                                         ####
###########################################################################################################
###########################################################################################################


source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
survival_data <- survival_data[survival_data$type == "BRCA", ]
survival_data <- survival_data[survival_data$Quartiles %in% c(1, 4), ]
tcga_q1_max <- max(survival_data$age_at_initial_pathologic_diagnosis[survival_data$Quartiles == 1])
tcga_q4_min <- min(survival_data$age_at_initial_pathologic_diagnosis[survival_data$Quartiles == 4])

# Convert metabric into young and old
metabric_clinical <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_clinical_patient.txt",
                              sep = "\t", header = T, stringsAsFactors = F)
metabric_clinical <- metabric_clinical[-c(1,2), ]
metabric_clinical$Quartile = ifelse(metabric_clinical$Age.at.Diagnosis <= tcga_q1_max, "Young", 
                                    ifelse(metabric_clinical$Age.at.Diagnosis >= tcga_q4_min, "Old", NA))
metabric_clinical <- metabric_clinical[!is.na(metabric_clinical$Quartile), ]

# Load metabric gene expression matrix
metabric_expression <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_expression.txt", 
                                sep = "\t", header = T, stringsAsFactors = F)
metabric_expression <- metabric_expression[, -2]
metabric_expression <- column_to_rownames(metabric_expression, "Hugo_Symbol")
colnames(metabric_expression) <- gsub('\\.', '-', colnames(metabric_expression))

# identify common cases
common_cases <- intersect(colnames(metabric_expression), metabric_clinical$X.Patient.Identifier)

# Keep common samples
metabric_expression <- metabric_expression[,common_cases]
metabric_clinical <- metabric_clinical[metabric_clinical$X.Patient.Identifier %in% common_cases, ]


# Set up for limma
###########################################################################################################
metabric_expression <- metabric_expression[, sort(colnames(metabric_expression))]

input <- data.frame(Patient = as.character(metabric_clinical$X.Patient.Identifier),
                    Age = metabric_clinical$Quartile,
                    row.names = as.character(metabric_clinical$X.Patient.Identifier))
input <- input[sort(rownames(input)), ]

# Verify this is true
all(input$Patient == colnames(metabric_expression))

# Limma
###########################################################################################################
input$Age <- factor(input$Age, levels = c("Old","Young"))
condition <- c("Old","Young")
comparison <- "Old-Young"
design <- model.matrix(~0+ input$Age)
colnames(design) <- condition
contmatrix <- makeContrasts(as.character(comparison),levels=design)

# filtered_counts<-metabric_expression[which(rowSums(cpm(metabric_expression)> 2) >= 2) ,]
print(dim(metabric_expression))
# Data.voom <- voom(metabric_expression,plot=TRUE)
# fit <- lmFit(Data.voom,design)
fit <- lmFit(metabric_expression,design)
fit2 <- contrasts.fit(fit,contmatrix)
fit2 <-eBayes(fit2)
tt <- topTable(fit2,n=Inf)
# tt <- getGeneNamesFromENSELBL(tt)


# FGSEA
###########################################################################################################
set.seed(2019)
tt <- tt[order(tt$t, decreasing = T),]
ranks <- tt$t
names(ranks) <- rownames(tt)
fgseaRes_METABRIC <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes_METABRIC$NES <- ifelse(fgseaRes_METABRIC$padj > 0.05, 0, fgseaRes_METABRIC$NES)


# Keep essentials
###########################################################################################################
rm(list=setdiff(ls(), c("fgseaRes_METABRIC", "fgseaRes_BRCA", "fgseaRes_NAT", "fgseaRes_GTx", "getGeneNamesFromENSELBL", "pathways.hallmark")))



###########################################################################################################
###########################################################################################################
####                                        TCGA Thyroid Cancer FGSEA                                  ####
###########################################################################################################
###########################################################################################################

# Survival data
###########################################################################################################
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
thca_samples <- survival_data$bcr_patient_barcode[survival_data$type == "THCA"]

# Load data
###########################################################################################################
load("~/Documents/Elemento/AgeTCGA-master/DATA/NewTCGAcounts/Full Matrix.RData")
filtered_counts <- RNAseq_counts
filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Solid Tumor" |
                                        TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Blood Derived Cancer - Peripheral Blood")]
biospec_table <- TCGAbiospec(colnames(filtered_counts))
biospec_table <- biospec_table[order(biospec_table$plate,decreasing = TRUE), ]
filtered_counts <- filtered_counts[,colnames(filtered_counts)[order(biospec_table$plate, decreasing = TRUE)]]

if (sum(duplicated(biospec_table$submitter_id)) > 0){
  filtered_counts <- filtered_counts[, !duplicated(biospec_table$submitter_id)]
}

# change column names to match metadata
colnames(filtered_counts) <- substr(colnames(filtered_counts), 1, 12)


# Select data
###########################################################################################################
survival_data <- survival_data[survival_data$bcr_patient_barcode %in% thca_samples, ]
survival_data <- survival_data[survival_data$Quartiles %in% c(1,4), ]
common <- intersect(colnames(filtered_counts), survival_data$bcr_patient_barcode)
filtered_counts <- filtered_counts[,colnames(filtered_counts) %in%  common]
# filtered_counts <- filtered_counts[, -which(!colnames(filtered_counts) %in% common)]
survival_data <- survival_data[survival_data$bcr_patient_barcode %in% common, ]

# Set up for limma
###########################################################################################################
filtered_counts <- filtered_counts[, sort(colnames(filtered_counts))]

input <- data.frame(Patient = as.character(survival_data$bcr_patient_barcode),
                    Age = survival_data$Quartiles,
                    row.names = as.character(survival_data$bcr_patient_barcode))
input <- input[sort(rownames(input)), ]
input$Age <- ifelse(input$Age == 1, "Young", "Old")

# Verify this is true
all(input$Patient == colnames(filtered_counts))


# Limma
###########################################################################################################
input$Age <- factor(input$Age, levels = c("Old","Young"))
condition <- c("Old","Young")
comparison <- "Old-Young"
design <- model.matrix(~0+ input$Age)
colnames(design) <- condition
contmatrix <- makeContrasts(as.character(comparison),levels=design)

filtered_counts<-filtered_counts[which(rowSums(cpm(filtered_counts)> 2) >= 2) ,]
print(dim(filtered_counts))
Data.voom <- voom(filtered_counts,plot=TRUE)
fit <- lmFit(Data.voom,design)
fit2 <- contrasts.fit(fit,contmatrix)
fit2 <-eBayes(fit2)
tt <- topTable(fit2,n=Inf)
tt <- getGeneNamesFromENSELBL(tt)


# FGSEA
###########################################################################################################
set.seed(2019)
tt <- tt[order(tt$t, decreasing = T),]
ranks <- tt$t
names(ranks) <- rownames(tt)
fgseaRes_THCA <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes_THCA$NES <- ifelse(fgseaRes_THCA$padj > 0.05, 0, fgseaRes_THCA$NES)


# Keep essentials
###########################################################################################################
rm(list=setdiff(ls(), c("fgseaRes_METABRIC", "fgseaRes_BRCA", "fgseaRes_NAT", 
                        "fgseaRes_GTx", "getGeneNamesFromENSELBL", "pathways.hallmark", "fgseaRes_THCA")))


###########################################################################################################
###########################################################################################################
####                                        GTEx Thyroid FGSEA                                         ####
###########################################################################################################
###########################################################################################################

# Load data
###########################################################################################################
meta1 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t")
meta2 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t")
gtex_data <- read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/Tissues/GTEx_Thyroid.txt")

# Setup for DE
###########################################################################################################
# match column names with meta2
colnames(gtex_data) <- gsub('\\.', '-', colnames(gtex_data))
colnames(gtex_data) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", colnames(gtex_data))

# convert transcript rownames to ensembl
gtex_data <- gtex_data[!duplicated(gsub("\\..*","",rownames(gtex_data))), ]
rownames(gtex_data) <- gsub("\\..*","",rownames(gtex_data))

age_meta <- meta2[meta2$SUBJID %in% colnames(gtex_data), ]
age_meta$AGE <- ifelse(age_meta$AGE == "20-29", "Young",
                       ifelse(age_meta$AGE == "30-39", "Young",
                              ifelse(age_meta$AGE == "40-49", "Young",
                                     ifelse(age_meta$AGE == "50-59", NA, "Old"))))
age_meta <- na.omit(data.frame(SUBJID = age_meta$SUBJID, AGE = age_meta$AGE))
gtex_data <- gtex_data[, colnames(gtex_data) %in% age_meta$SUBJID]

# Verify this is true
all(colnames(gtex_data) == age_meta$SUBJID)

# Limma
###########################################################################################################
age_meta$AGE <- factor(age_meta$AGE, levels = c("Old","Young"))
condition <- c("Old","Young")
comparison <- "Old-Young"
design <- model.matrix(~0+ age_meta$AGE)
colnames(design) <- condition
contmatrix <- makeContrasts(as.character(comparison),levels=design)

gtex_data<-gtex_data[which(rowSums(cpm(gtex_data)> 2) >= 2) ,]
print(dim(gtex_data))
Data.voom <- voom(gtex_data,plot=TRUE)
fit <- lmFit(Data.voom,design)
fit2 <- contrasts.fit(fit,contmatrix)
fit2 <-eBayes(fit2)
tt <- topTable(fit2,n=Inf)
tt <- getGeneNamesFromENSELBL(tt)


# FGSEA
###########################################################################################################
set.seed(2019)
tt <- tt[order(tt$t, decreasing = T),]
ranks <- tt$t
names(ranks) <- rownames(tt)
fgseaRes_Thy_GTx <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes_Thy_GTx$NES <- ifelse(fgseaRes_Thy_GTx$padj > 0.05, 0, fgseaRes_Thy_GTx$NES)


# Keep essentials
###########################################################################################################
rm(list=setdiff(ls(), c("fgseaRes_METABRIC", "fgseaRes_BRCA", "fgseaRes_NAT", 
                        "fgseaRes_GTx", "getGeneNamesFromENSELBL", "pathways.hallmark",
                        "fgseaRes_THCA","fgseaRes_Thy_GTx")))



###########################################################################################################
###########################################################################################################
####                                           Visualize                                               ####
###########################################################################################################
###########################################################################################################

df <- data.frame(GTEx.Breast = fgseaRes_GTx$NES,
                 BRCA.NAT = fgseaRes_NAT$NES,
                 GTEx.Thyroid = fgseaRes_Thy_GTx$NES,
                 TCGA.BRCA = fgseaRes_BRCA$NES,
                 METABRIC = fgseaRes_METABRIC$NES,
                 TCGA.THCA = fgseaRes_THCA$NES,
                 row.names = fgseaRes_BRCA$pathway)
rownames(df) <- substring(rownames(df), 10)
rownames(df) <- gsub('\\_', ' ', rownames(df))

# select rows
rowsToKeep <- rowSums(df != 0) > 3
df <- df[rowsToKeep, ]

# annotation
column_anno <- data.frame(Source = c(rep("Healthy",3), rep("Tumor", 3)),
                          row.names = colnames(df))

postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.2/ALL_cancer_vs_normal_pathways.eps")
pheatmap(mat = df,border_color = "white", color = pal_gsea("default")(12), cellheight = 15, cellwidth = 15,
         annotation_col = column_anno, annotation_names_col = F,  fontsize = 15, 
         treeheight_row = 15, treeheight_col = 15, legend = F, cutree_cols = 2,
         annotation_colors = list(Source = c(`Tumor` = "black", `Healthy` = "gray")))
dev.off()


df <- df[-c(2,3,5,7,8,11,13,14,15,16,17,19,20,21,22,25), ]
postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.2/SELECTED_cancer_vs_normal_pathways.eps")
pheatmap(mat = df,border_color = "white", color = pal_gsea("default")(12), cellheight = 20, cellwidth = 20,
         annotation_col = column_anno, annotation_names_col = F,  fontsize = 20, 
         treeheight_row = 15, treeheight_col = 15, legend = F, cutree_cols = 2,
         annotation_colors = list(Source = c(`Tumor` = "black", `Healthy` = "gray")))
dev.off()
