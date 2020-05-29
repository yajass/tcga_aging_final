library(stringr)
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

gtex_differentialEpression <- function(){
  # to collect results
  results_rnaseq_genes <- setNames(replicate(length(files),data.frame()),
                                   unlist(purrr::map(str_split(files, "#"), 1)))
  
  for (i in 1:length(files)){
    # Select data
    filtered_counts <- data.table::fread(paste0(folder, files[i])) %>%
      tibble::column_to_rownames("V1") %>%
      dplyr::select(sort(colnames(.)))
    # Load age data
    age_subset <- read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                             stringsAsFactors = F, header = T, sep = "\t") %>%
      dplyr::mutate(AGE = ifelse(AGE == "20-29", "Young",
                                 ifelse(AGE == "30-39", "Young",
                                        ifelse(AGE == "40-49", "Young",
                                               ifelse(AGE == "50-59", NA, "Old"))))) %>%
      dplyr::filter(SUBJID %in% colnames(filtered_counts)) %>%
      dplyr::select(SUBJID, AGE) %>%
      tidyr::drop_na() %>%
      dplyr::arrange(SUBJID)
    
    filtered_counts <- filtered_counts[,colnames(filtered_counts) %in% age_subset$SUBJID]
    
    print(ifelse(all(age_subset$SUBJID == colnames(filtered_counts)), paste0(files[i]," Clear"), paste0(files[i], " does NOT match")))
    # Limma
    age_subset$AGE <- factor(age_subset$AGE, levels = c("Old","Young"))
    design <- model.matrix(~0+ age_subset$AGE)
    colnames(design) <- condition
    contmatrix <- makeContrasts(as.character(comparison),levels=design)
    Count_input <- filtered_counts
    Count_input <- Count_input[which(rowSums(cpm(Count_input)> 2) >= 2) ,]
    print(dim(Count_input))
    Data.voom <- voom(Count_input,plot=FALSE)
    fit <- lmFit(Data.voom,design)
    fit2 <- contrasts.fit(fit,contmatrix)
    fit2 <-eBayes(fit2)
    tt <- topTable(fit2,n=Inf)
    tt <- getGeneNamesFromENSELBL(tt)
    tt <- rownames_to_column(tt, "Gene")
    results_rnaseq_genes[[i]] <- tt
  }
  return(results_rnaseq_genes)
}

tcga_differentialEpression <- function(){
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
    tt <- rownames_to_column(tt, "Gene")
    results_rnaseq_genes[[i]] <- tt
  }
  return(results_rnaseq_genes)
}


# Load data
########################################################################################################################
# comparisions
condition <- c("Old","Young")
comparison <- "Old-Young"
# For GTEX
folder <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/TissueBySubtype/"
files <- list.files(folder)
# For TCGA
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/NewTCGAcounts/Full Matrix.RData")
cancer_codes <- c("BRCA","COAD","THCA","LGG","OV","UCEC")

# Run differntial analyses
########################################################################################################################
tcga_deg <- tcga_differentialEpression(); rm(RNAseq_counts)
gtex_deg <- gtex_differentialEpression()

# Reorder results and keep what's useful
########################################################################################################################
tcga_deg <- tcga_deg[c(1,2,5,3,6)]


# Create Lists and Overlap
########################################################################################################################
gene_lists <- setNames(rep(list(setNames(rep(list(NA), 4), c("GTEx_Old", "GTEx_Young", "TCGA_Old", "TCGA_Young"))), 5), names(tcga_deg))

for (i in 1:length(gtex_deg)) {
  gene_lists[[i]][[1]] <- gtex_deg[[i]] %>%
    dplyr::mutate(Regulation = ifelse(logFC > 0, "Old", "Young")) %>%
    dplyr::filter(adj.P.Val < 0.05,
                  Regulation =="Old") %>%
    .$Gene
  gene_lists[[i]][[2]] <- gtex_deg[[i]] %>%
    dplyr::mutate(Regulation = ifelse(logFC > 0, "Old", "Young")) %>%
    dplyr::filter(adj.P.Val < 0.05,
                  Regulation =="Young") %>%
    .$Gene
  gene_lists[[i]][[3]] <- tcga_deg[[i]] %>%
    dplyr::mutate(Regulation = ifelse(logFC > 0, "Old", "Young")) %>%
    dplyr::filter(adj.P.Val < 0.05,
                  Regulation =="Old") %>%
    .$Gene
  gene_lists[[i]][[4]] <- tcga_deg[[i]] %>%
    dplyr::mutate(Regulation = ifelse(logFC > 0, "Old", "Young")) %>%
    dplyr::filter(adj.P.Val < 0.05,
                  Regulation =="Young") %>%
    .$Gene
}

names(gene_lists) <- c("Breast", "Colon", "Ovary", "Thyroid", "Uterus")
myCol <- brewer.pal(4, "Pastel2")
for (i in 1:length(names(gene_lists))) {
  VennDiagram::venn.diagram(x = gene_lists[[i]],
                            paste0("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S2/venn_",names(gene_lists)[i],".tiff"),
                            imagetype="tiff" ,
                            height = 1200 , 
                            width = 1200 , 
                            resolution = 300,
                            compression = "lzw",
                            fill = myCol,
                            cex = .6,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 0.6,
                            cat.fontface = "bold",
                            cat.default.pos = "outer",
                            cat.fontfamily = "sans",
                            main = names(gene_lists)[i],
                            main.cex = 1.2,
                            main.fontface = "bold",
                            main.fontfamily = "sans")
}
