# run on R 3.6.3
# Load everything
################################################################################################################
# packages
library(stringr)
library(TCGAbiolinks)
library(TCGAutils)
library(SummarizedExperiment)
library(dplyr)
library(ggpubr)
library(tibble)
library(tidyr)
library(org.Hs.eg.db)
library(biomaRt)
library(data.table)
library(CePa)
library(GSVA)

# scripts
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")

# functions
getGeneNamesFromENSELBL <- function(countData){
  rownames(countData) <- gsub("\\..*","",rownames(countData))
  countData[,ncol(countData)+1] <- mapIds(org.Hs.eg.db, keys = gsub("\\..*","",rownames(countData)),
                                          keytype = "ENSEMBL", column="SYMBOL", "first")
  countData <- na.omit(countData)
  countData <- countData[-which(duplicated(countData[,ncol(countData)])),]
  rownames(countData) <- countData[,ncol(countData)]
  countData <- countData[,-ncol(countData)]
  return(countData)
}

# gene signatures
gene_signatures <- GSA::GSA.read.gmt("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GSEA_pathways/senesccence_gene_sets.gmt")
gene_signatures <- setNames(gene_signatures$genesets,
                            gene_signatures$geneset.names)

# The signatures from msigdb include:
# [1] "COURTOIS SENESCENCE TRIGGERS"                                            
# [2] "GO AGING"                                                                
# [3] "GO CELL AGING"                                                           
# [4] "DEMAGALHAES AGING UP"                                                    
# [5] "FRIDMAN SENESCENCE UP"                                                   
# [6] "GO MULTI-\nCELLULAR ORGANISM AGING"                                      
# [7] "GO STRESS INDUCED PREMATURE SENESCENCE"                                  
# [8] "REACTOME CELLULAR SENESCENCE"                                            
# [9] "REACTOME SENESCENCE ASSOCIATED SECRETORY PHENOTYPE SASP"                 
# [10] "REACTOME OXIDATIVE STRESS INDUCED SENESCENCE"                            
# [11] "REACTOME ONCOGENE INDUCED SENESCENCE"                                    
# [12] "REACTOME FORMATION OF SENESCENCE ASSOCIATED HETERO-\nCHROMATIN FOCI SAHF"
# [13] "REACTOME DNA DAMAGE TELOMERE STRESS INDUCED SENESCENCE"  

# TCGA TPM Data
folder <- "/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/TPM/TCGA/"
list.files(folder)

# Get ssgsea
for (i in 1:length(list.files(folder))) {
  dat0 <- read.table(paste0(folder, list.files(folder)[i]), sep = "\t", stringsAsFactors = F, header = T)
  colnames(dat0) <- gsub('\\.','-',colnames(dat0))
  
  ## select smples and change column names
  dat0 <- dat0[ ,(TCGAbiospec(colnames(dat0))$sample_definition == "Primary Solid Tumor" |
                    TCGAbiospec(colnames(dat0))$sample_definition == "Primary Blood Derived Cancer - Peripheral Blood")]
  # Remove duplicates if they exist
  biospec_table <- TCGAbiospec(colnames(dat0))
  biospec_table <- biospec_table[order(biospec_table$plate,decreasing = TRUE), ]
  dat0 <- dat0[,colnames(dat0)[order(biospec_table$plate, decreasing = TRUE)]]
  if (sum(duplicated(biospec_table$submitter_id)) > 0){
    dat0 <- dat0[, !duplicated(biospec_table$submitter_id)]
  }
  
  # change column names to match metadata
  colnames(dat0) <- substr(colnames(dat0), 1, 12)
  dat0 <- getGeneNamesFromENSELBL(dat0[, colnames(dat0) %in% survival_data$bcr_patient_barcode])
  
  if (i == 1){
    set.seed(2019)
    enrich_result <- data.frame(gsva(as.matrix(dat0), gene_signatures, "Gaussian", "ssgsea"),
                                stringsAsFactors = FALSE)
  }
  if (i > 1){
    set.seed(2019)
    enrich_result <- cbind(enrich_result,
                           data.frame(gsva(as.matrix(dat0), gene_signatures, "Gaussian", "ssgsea"),
                                      stringsAsFactors = FALSE))
  }
}
enrich_result <- data.frame(enrich_result, stringsAsFactors = F)
colnames(enrich_result) <- gsub('\\.', '-', colnames(enrich_result))
write.table(enrich_result, "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Aging Enrichment/tcga_tpm_ssgsea_enrichment.txt",
            sep = "\t", quote = F)
