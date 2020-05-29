library(TCGAutils)
# Load data
########################################################################################################################
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/General Functions/immunophenotypeAnalysisFunctions.R")
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/NewTCGAcounts/Full Matrix.RData")

# Load functions
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

# Subset data
########################################################################################################################
gene.interest <- c("MKI67","PCNA","MYBL2","PLK1","CCND1","TOP2A","FOXM1")
tcga_counts <- getGeneNamesFromENSELBL(RNAseq_counts)[gene.interest, ]
colnames(tcga_counts) <- gsub('\\.','-',colnames(tcga_counts))

tcga_counts <- tcga_counts[ ,(TCGAbiospec(colnames(tcga_counts))$sample_definition == "Primary Solid Tumor" |
                                TCGAbiospec(colnames(tcga_counts))$sample_definition == 
                                "Primary Blood Derived Cancer - Peripheral Blood")]
## remove duplicates if they exist
if (sum(duplicated(TCGAbiospec(colnames(tcga_counts))$submitter_id)) > 0){
  tcga_counts <- tcga_counts[,-c(which(duplicated(TCGAbiospec(colnames(tcga_counts))$submitter_id)))]
}
## making sure both datasets have the same case id format
colnames(tcga_counts) <- substr(colnames(tcga_counts), start = 1, stop = 12)

# Set parameters
########################################################################################################################
tumor_types <- c('BRCA','UCEC','THCA','COAD','LGG','OV')
stat.type <- "spearman"
comp.type <- "quartiles"

# get correlation matrix
########################################################################################################################
for (i in 1:length(tumor_types)) {
  if (i == 1){
    # Temp for labelling rows
    temp_res <- deconvCorrelateStat(deconv_matrix = tcga_counts,
                                    comparision = comp.type,
                                    surv_data = survival_data,
                                    statistic = stat.type)
    
    # Initialize results dataframe
    cor_res_stat <- data.frame(matrix(ncol = 6, nrow = length(temp_res)), row.names = names(temp_res))
    colnames(cor_res_stat)[i] <- tumor_types[i]
    cor_res_p <- cor_res_stat
    colnames(cor_res_p)[i] <- tumor_types[i]
    
    # Write results
    cor_res_stat[,i] <- sapply(temp_res, "[[", "estimate")
    cor_res_p[,i] <- sapply(temp_res, "[[", "p.value")
  }
  
  if (i > 1){
    colnames(cor_res_stat)[i] <- tumor_types[i]
    colnames(cor_res_p)[i] <- tumor_types[i]
    
    temp_res <- deconvCorrelateStat(deconv_matrix = tcga_counts,
                                    comparision = comp.type,
                                    surv_data = survival_data,
                                    statistic = stat.type)
    
    # Write results
    cor_res_stat[,i] <- sapply(temp_res, "[[", "estimate")
    cor_res_p[,i] <- sapply(temp_res, "[[", "p.value")
  }
}

# Process and plot
########################################################################################################################
cor_res_fdr <- apply(cor_res_p, 2, p.adjust, "fdr")
postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S2/proliferation_aggressiveness_correlation.eps")
corrplot::corrplot(t(as.matrix(cor_res_stat)), p.mat = t(cor_res_fdr),
                   "color", is.corr = F,
                   insig = "label_sig", sig.level = c(.001, .01, .05),
                   pch.cex = 2, pch.col = "white", order = "hclust",
                   cl.lim = c(-0.35,0.35), tl.cex = 1.5,
                   col = colorRampPalette(c("blue","white","red"))(200), tl.col = "black", mar = c(0,0,0,2))
p <- par('usr')
text(p[2], mean(p[3:4])-0.55, labels = 'Spearman', xpd = NA, srt = -90, cex = 2)
dev.off()
