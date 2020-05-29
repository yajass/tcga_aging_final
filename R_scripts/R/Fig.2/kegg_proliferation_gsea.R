library(ggsci)
library(dplyr)
library(ggcorrplot)
library(RColorBrewer)
library(corrplot)
library(GeneOverlap)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(TCGAbiolinks)
library(tibble)
library(SummarizedExperiment)
library(wateRmelon)
library(limma)
library(edgeR)
library(ggpubr)
library(TCGAbiolinks)
library(TCGAutils)
library(fgsea)

# Data
########################################################################################################################
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
survival_data <- survival_data[survival_data$Quartiles %in% c(1,4), ]
survival_data$Quartiles <- ifelse(survival_data$Quartiles == 1, "Young", "Old")
pathways <- gmtPathways("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GSEA_pathways/KEGG-Cancer_Pathways_Proliferation")

# Setup functions
########################################################################################################################
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

differentialEpression <- function(){
  # to collect results
  results_rnaseq_genes <- setNames(replicate(length(cancer_codes),data.frame()),cancer_codes)
  
  for (i in 1:length(cancer_codes)){
    # Load age data
    survival_data <- read.table("/Users/Yajas/Downloads/survival_data.txt", 
                                sep = "\t", header = T, row.names = 1)
    survival_data <- survival_data[survival_data$Quartiles %in% c(1,4), ]
    survival_data$Quartiles <- ifelse(survival_data$Quartiles == 1, "Young", "Old")
    # Select data
    filtered_counts <- RNAseq_counts
    filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Solid Tumor" |
                                            TCGAbiospec(colnames(filtered_counts))$sample_definition == 
                                            "Primary Blood Derived Cancer - Peripheral Blood")]
    ## remove duplicates if they exist
    if (sum(duplicated(TCGAbiospec(colnames(filtered_counts))$submitter_id)) > 0){
      filtered_counts <- filtered_counts[,-c(which(duplicated(TCGAbiospec(colnames(filtered_counts))$submitter_id)))]
    }
    ## making sure both datasets have the same case id format
    colnames(filtered_counts) <- substr(colnames(filtered_counts), start = 1, stop = 12)
    ## remove CASE ID factors from age metadata
    survival_data$bcr_patient_barcode <- as.character(survival_data$bcr_patient_barcode)
    ## selecting case IDs common to both datasets and discarding the rest
    common_cases <- intersect(colnames(filtered_counts), survival_data$bcr_patient_barcode)
    filtered_age <- survival_data[survival_data$bcr_patient_barcode %in% common_cases,]
    filtered_counts <- filtered_counts[,colnames(filtered_counts) %in% common_cases]
    ## sorting case IDs
    filtered_age <- filtered_age[order(filtered_age$bcr_patient_barcode),]
    filtered_counts <- filtered_counts[,order(colnames(filtered_counts))]
    
    # Limma
    age_subset <- filtered_age[which(filtered_age$type == cancer_codes[i]),  c("bcr_patient_barcode","Quartiles")]
    rna_subset <- filtered_counts[,which(colnames(filtered_counts) %in% age_subset$bcr_patient_barcode)]
    input <- as.data.frame(age_subset[,c("bcr_patient_barcode","Quartiles")])
    input$Quartiles <- factor(input$Quartiles, levels = c("Old","Young"))
    design <- model.matrix(~0+ input$Quartiles)
    colnames(design) <- condition
    contmatrix <- makeContrasts(as.character(comparison),levels=design)
    Count_input <- rna_subset
    Count_input<-Count_input[which(rowSums(cpm(Count_input)> 2) >= 2) ,]
    print(dim(Count_input))
    Data.voom <- voom(Count_input,plot=FALSE)
    Voom_output <- Data.voom
    fit <- lmFit(Data.voom,design)
    fit2 <- contrasts.fit(fit,contmatrix)
    fit2 <-eBayes(fit2)
    tt <- topTable(fit2,n=Inf)
    tt <- getGeneNamesFromENSELBL(tt)
    
    # Run fgsea
    set.seed(2019)
    tt <- tt[order(tt$t, decreasing = T),]
    ranks_TCGA <- tt$t
    names(ranks_TCGA) <- rownames(tt)
    results_rnaseq_genes[[i]] <- fgsea(pathways=pathways, stats=ranks_TCGA, nperm=1000)
  }
  return(results_rnaseq_genes)
}


# Load data
########################################################################################################################
# comparisions
condition <- c("Old","Young")
comparison <- "Old-Young"
# For expression
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/NewTCGAcounts/Full Matrix.RData")
cancer_codes <- c("BRCA","COAD","LGG","OV","THCA","UCEC")


# Run differntial analyses
########################################################################################################################
results_gsea <- differentialEpression()

# Combine and plot
########################################################################################################################
nes <- data.frame(lapply(results_gsea, function (x) x$NES),
                  row.names = gsub('\\_',' ', substring(results_gsea$BRCA$pathway, 6)),
                  stringsAsFactors = F)
fdr <- data.frame(lapply(results_gsea, function (x) x$padj),
                  row.names = gsub('\\_',' ', substring(results_gsea$BRCA$pathway, 6)),
                  stringsAsFactors = F)
nes[fdr >= 0.05] <- 0
fdr[fdr >= 0.05] <- 1
fdr <- fdr[abs(rowSums(nes)) != 0, ]
nes <- nes[abs(rowSums(nes)) != 0, ]

postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S2/kegg_pathway_gsea_proliferation.eps")
pheatmap(nes, border_color = "white", color = pal_gsea()(10),
         cellheight = 15, cellwidth = 15, treeheight_row = 15, treeheight_col = 15,
         fontsize_number = 15, fontsize = 15)
dev.off()
