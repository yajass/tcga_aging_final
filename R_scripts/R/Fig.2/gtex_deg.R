# Load everything
################################################################################################################
# packages
library(dplyr)
library(purr)
library(tibble)
library(tidyr)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(biomaRt)
library(CePa)
# data
gtex_counts <- read.gct("/athena/khuranalab/scratch/yas4002/Aging/gtex_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
gtex_meta1 <- read.delim("/athena/khuranalab/scratch/yas4002/Aging/gtex_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                         stringsAsFactors = FALSE, header = TRUE, sep = "\t")
gtex_meta1 <- gtex_meta1[gtex_meta1$SMAFRZE == "RNASEQ", ]
gtex_meta2 <- read.delim("/athena/khuranalab/scratch/yas4002/Aging/gtex_data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                         stringsAsFactors = FALSE, header = TRUE, sep = "\t")

tissues_of_interest <- c("Breast - Mammary Tissue", "Thyroid", "Colon - Sigmoid", "Ovary", "Uterus")


# setup functoins
################################################################################################################
# Converts ENSEMBL IDs (row names) to unique gene symbols
getGeneNamesFromENSELBL <- function(countData){
  # rownames(countData) <- gsub("\\..*","",rownames(countData))
  countData[,ncol(countData)+1] <- mapIds(org.Hs.eg.db, keys = gsub("\\..*","",rownames(countData)),
                                          keytype = "ENSEMBL", column="SYMBOL", "first")
  countData <- na.omit(countData)
  countData <- countData[-which(duplicated(countData[,ncol(countData)])),]
  rownames(countData) <- countData[,ncol(countData)]
  countData <- countData[,-ncol(countData)]
  return(countData)
}

# annotate age groups
annotateAgeGroups <- function(input){
  input <- data.frame(SUBJID = input$SUBJID,
                      AGE = input$AGE, stringsAsFactors = FALSE)
  input$AGE <- ifelse(input$AGE == "20-29", "Young",
                      ifelse(input$AGE == "30-39", "Young",
                             ifelse(input$AGE == "40-49", "Young",
                                    ifelse(input$AGE == "50-59", NA, "Old"))))
  input <- na.omit(input)
  return(input)
}

# differential expression
runDEG <- function(lstData){
  countMatrix <- lstData[[1]]
  input <- lstData[[2]]
  
  input$AGE <- factor(input$AGE, levels = c("Old","Young"))
  design <- model.matrix(~0+ input$AGE)
  colnames(design) <- condition
  contmatrix <- makeContrasts(as.character(comparison),levels=design)
  
  countMatrix <- countMatrix[which(rowSums(cpm(countMatrix)> 2) >= 2) ,]
  print(dim(countMatrix))
  Data.voom <- voom(countMatrix, plot=FALSE)
  fit <- lmFit(Data.voom,design)
  fit2 <- contrasts.fit(fit,contmatrix)
  fit2 <-eBayes(fit2)
  tt <- topTable(fit2, n=Inf)
  tt <- getGeneNamesFromENSELBL(tt)
  
  return(tt)
}


# select samples for each tissue type
selectData <- function(counts = gtex_counts, meta1 = gtex_meta1, meta2 = gtex_meta2){
  sub.gtex_counts <- data.frame(counts)
  tissue <- tissues_of_interest[i]
  
  # convert transcript rownames to ensembl
  sub.gtex_counts <- sub.gtex_counts[!duplicated(gsub("\\..*","",rownames(sub.gtex_counts))), ]
  rownames(sub.gtex_counts) <- gsub("\\..*","",rownames(sub.gtex_counts))
  
  # match count colnames with meta1
  colnames(sub.gtex_counts) <- gsub('\\.', '-', colnames(sub.gtex_counts))
  
  # select samples for tissue of interest --> subset counts
  sub.meta1 <- meta1[meta1$SMTSD == tissue, ]
  sub.gtex_counts <- sub.gtex_counts[ ,colnames(sub.gtex_counts) %in% sub.meta1$SAMPID]
  
  # match count colnames with meta2 --> subset meta2
  colnames(sub.gtex_counts) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", colnames(sub.gtex_counts))
  sub.meta2 <- meta2[meta2$SUBJID %in% colnames(sub.gtex_counts), ]
  sub.gtex_counts <- sub.gtex_counts[, sub.meta2$SUBJID]
  
  # order the data
  sub.gtex_counts <- sub.gtex_counts[, sort(colnames(sub.gtex_counts))]
  sub.meta2 <- sub.meta2[order(sub.meta2$SUBJID), ]
  
  # check if everything is fine
  if (all(colnames(sub.gtex_counts) == sub.meta2$SUBJID)) {
    print(paste0(tissues_of_interest[i], " matches up"))
  }
  if (!all(colnames(sub.gtex_counts) == sub.meta2$SUBJID)) {
    print(paste0(tissues_of_interest[i], " does NOT match up"))
  }
  
  return(list(sub.gtex_counts, sub.meta2))
}


# setup comparisions and results
################################################################################################################
condition <- c("Old","Young")
comparison <- "Old-Young"
gtex_meta2 <- annotateAgeGroups(gtex_meta2)
de_results <- setNames(c(rep(list(NA), length(tissues_of_interest))), make.names(tissues_of_interest, unique = TRUE))
gtex_counts <- as.data.frame(gtex_counts)

# run DE
################################################################################################################
for (i in 1:length(tissues_of_interest)) {
  data_for_DE <- selectData(counts = gtex_counts,
                            meta1 = gtex_meta1,
                            meta2 = gtex_meta2)
  # get age distribution
  ages <- data_for_DE[[2]]
  if (i == 1){
    age_groups <- data.frame(table(ages$AGE))
    colnames(age_groups) <- c("Age", make.names(tissues_of_interest[i]))
  }
  if (i > 1){
    age_groups <- data.frame(cbind(age_groups), matrix(data.frame(table(ages$AGE))[,2], ncol = 1))
    colnames(age_groups)[i+1] <- make.names(tissues_of_interest[i])
  }
  
  tt <- runDEG(data_for_DE)
  de_results[[i]] <- tt
}


# remove everything but the output
rm(list = setdiff(ls(), c("de_results", "age_groups")))
saveRDS(de_results, "/athena/khuranalab/scratch/yas4002/Aging/GTEx_aging_DEG_all.RDS")