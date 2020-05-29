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

# Data
########################################################################################################################
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
survival_data <- survival_data[survival_data$Quartiles %in% c(1,4), ]
survival_data$Quartiles <- ifelse(survival_data$Quartiles == 1, "Young", "Old")

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

differentialMethylation <- function(){
  # create dataframes to store results
  results_meth_genes <- setNames(replicate(length(tumor_type),data.frame()),tumor_type)
  
  for (i in 1:length(files)) {
    beta_values <- read.table(paste0(folder, files[i]), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
    beta_values <- beta_values[-1, ]
    rnames <- rownames(beta_values)
    beta_values <- as.data.frame(sapply(beta_values, as.numeric),
                                 row.names = rnames)
    
    m_values <- Beta2M(beta_values)
    colnames(m_values) <- substr(colnames(m_values), 1, 12)
    colnames(m_values) <- gsub('\\.', '-', colnames(m_values))
    
    # find common cases
    input <- survival_data[survival_data$type == tumor_type[i], ]
    input <- data.frame(Subject = input$bcr_patient_barcode,
                        Age = input$Quartiles)
    common <- intersect(input$Subject, colnames(m_values))
    
    # keep common cases and reorder df
    input <- input[input$Subject %in% common, ]
    m_values <- m_values[, common]
    input <- input[order(input$Subject), ]
    m_values <- m_values[, order(colnames(m_values))]
    input$Subject <- as.character(input$Subject)
    
    # check
    print(ifelse(all(input$Subject == colnames(m_values)), paste0(tumor_type[i], " clear"), paste0(tumor_type[i], " ERROR")))
    
    # limma
    input$Age <- factor(input$Age, levels = c("Old","Young"))
    design <- model.matrix(~0+ input$Age)
    colnames(design) <- condition
    contmatrix <- makeContrasts(as.character(comparison),levels=design)
    fit <- lmFit(m_values,design)
    fit2 <- contrasts.fit(fit,contmatrix)
    fit2 <-eBayes(fit2)
    tt_meth <- topTable(fit2,n=Inf)
    tt_meth <- rownames_to_column(tt_meth, "Gene")
    
    # write results
    results_meth_genes[[i]] <- tt_meth
  }
  return(results_meth_genes)
}

differentialEpression <- function(){
  # to collect results
  results_rnaseq_genes <- setNames(replicate(length(cancer_codes),data.frame()),cancer_codes)
  
  for (i in 1:length(cancer_codes)){
    # Load age data
    # survival_data <- read.table("/Users/Yajas/Downloads/survival_data.txt", 
    #                             sep = "\t", header = T, row.names = 1)
    # survival_data1 <- survival_data[survival_data$Quartiles %in% c(1,4), ]
    # survival_data1$Quartiles <- ifelse(survival_data1$Quartiles == 1, "Young", "Old")
    survival_data1 <- survival_data
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
    survival_data1$bcr_patient_barcode <- as.character(survival_data1$bcr_patient_barcode)
    ## selecting case IDs common to both datasets and discarding the rest
    common_cases <- intersect(colnames(filtered_counts), survival_data1$bcr_patient_barcode)
    filtered_age <- survival_data1[survival_data1$bcr_patient_barcode %in% common_cases,]
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
# For methylation
folder <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_gene_beta/"
files <- list.files(folder)
tumor_type <- gsub("\\..*","", files)
# For expression
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/NewTCGAcounts/Full Matrix.RData")
cancer_codes <- tumor_type


# Run differntial analyses
########################################################################################################################
results_meth_genes <- differentialMethylation()
results_rnaseq_genes <- differentialEpression()

# Number of up/down genes
########################################################################################################################
gene_count_rna <- bind_rows(results_rnaseq_genes, .id = "Cancer")
gene_count_rna %>%
  dplyr::mutate(Regulation = ifelse(logFC > 0, "Old", "Young")) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::group_by(Cancer) %>%
  dplyr::count(Regulation) %>%
  tidyr::spread(Regulation, n)


# Gene overlap - RNASeq Pairwise
########################################################################################################################
up_genes_rna <- lapply(results_rnaseq_genes, function(u)u$Gene[u$adj.P.Val < 0.05 & u$logFC > 0])
down_genes_rna <- lapply(results_rnaseq_genes, function(u)u$Gene[u$adj.P.Val < 0.05 & u$logFC < 0])

# find genes common - referenced in text
names(table(plyr::ldply(down_genes_rna, cbind)[,2]))[table(plyr::ldply(down_genes_rna, cbind)[,2]) > 2]
names(table(plyr::ldply(up_genes_rna, cbind)[,2]))[table(plyr::ldply(up_genes_rna, cbind)[,2]) > 2]

# Upregulated genes
or_data_up <- log10(getMatrix(newGOM(up_genes_rna), name = "odds.ratio")+1)
pval_data_up_before_adjust <- getMatrix(newGOM(up_genes_rna),name = "p")
pval_data_up <- matrix(p.adjust(pval_data_up_before_adjust, "fdr"), nrow = nrow(pval_data_up_before_adjust), ncol = ncol(pval_data_up_before_adjust))
dimnames(pval_data_up) <- dimnames(pval_data_up_before_adjust)

# Downregulated genes
or_data_down <- log10(getMatrix(newGOM(down_genes_rna), name = "odds.ratio")+1)
pval_data_down_before_adjust <- getMatrix(newGOM(down_genes_rna),name = "p")
pval_data_down <- matrix(p.adjust(pval_data_down_before_adjust, "fdr"), nrow = nrow(pval_data_down_before_adjust), ncol = ncol(pval_data_down_before_adjust))
dimnames(pval_data_down) <- dimnames(pval_data_up_before_adjust)

postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S2/DEG_DEG_Overlap.eps")
par(mfrow=c(1,2))
corrplot(or_data_down, is.corr = F, method = "color", p.mat = pval_data_down, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = 1.5, pch.col = "white", tl.col = "black", type = "upper",
         tl.cex = 1.5, outline = "white", col = brewer.pal(9, "Blues"), mar = c(0,0,0,2))
mtext(text = "Young", side = 3, line = -3, cex = 2)
corrplot(or_data_up, is.corr = F, method = "color", p.mat = pval_data_up, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = 1.5, pch.col = "white", tl.col = "black", type = "upper",
         tl.cex = 1.5, outline = "white", col = brewer.pal(9, "Reds"), mar = c(0,0,0,2))
mtext(text = "Old", side = 3, line = -3, cex = 2)
dev.off()

# Gene overlap - RNASeq and DNAm
########################################################################################################################
up_genes_rna <- lapply(results_rnaseq_genes, function(u)u$Gene[u$adj.P.Val < 0.05 & u$logFC > 0])
down_genes_rna <- lapply(results_rnaseq_genes, function(u)u$Gene[u$adj.P.Val < 0.05 & u$logFC < 0])
hyper_meth_genes <- lapply(results_meth_genes, function(u)u$Gene[u$adj.P.Val < 0.05 & u$logFC > 0])
hypo_meth_genes <- lapply(results_meth_genes, function(u)u$Gene[u$adj.P.Val < 0.05 & u$logFC < 0])

names(up_genes_rna) <- paste0(names(up_genes_rna), " mRNA")
names(down_genes_rna) <- paste0(names(down_genes_rna), " mRNA")
names(hyper_meth_genes) <- paste0(names(hyper_meth_genes), " DNAm")
names(hypo_meth_genes) <- paste0(names(hypo_meth_genes), " DNAm")

gom.object_up <- newGOM(lapply(up_genes_rna, unlist),lapply(hypo_meth_genes, unlist))
or_data_up <- log10(getMatrix(gom.object_up, name = "odds.ratio")+1)
pval_data_up_before_adjust <- getMatrix(gom.object_up,name = "p")
pval_data_up <- matrix(p.adjust(pval_data_up_before_adjust, "fdr"), nrow = nrow(pval_data_up_before_adjust), ncol = ncol(pval_data_up_before_adjust))
dimnames(pval_data_up) <- dimnames(pval_data_up_before_adjust)

gom.object_down <- newGOM(lapply(down_genes_rna, unlist),lapply(hyper_meth_genes, unlist))
or_data_down <- log10(getMatrix(gom.object_down, name = "odds.ratio")+1)
pval_data_down_before_adjust <- getMatrix(gom.object_down,name = "p")
pval_data_down <- matrix(p.adjust(pval_data_down_before_adjust, "fdr"), nrow = nrow(pval_data_down_before_adjust), ncol = ncol(pval_data_down_before_adjust))
dimnames(pval_data_down) <- dimnames(pval_data_up_before_adjust)

newDF <- data.frame(OR = c(diag(or_data_down), diag(or_data_up)),
                    FDR = c(p.adjust(diag(pval_data_down_before_adjust), "fdr"),
                            p.adjust(diag(pval_data_up_before_adjust), "fdr")),
                    Tumor = rep(gsub('.{5}$', '', colnames(or_data_down)), 2),
                    Type = c(rep("Young", 6), rep("Old", 6))) %>%
  dplyr::mutate(Sig = ifelse(FDR < 0.05, "< 0.05", "≥ 0.05"),
                Type = factor(Type, levels = c("Young", "Old")))

ggbarplot(data = newDF, x = 'Tumor', y = 'OR', fill = 'Sig',
          palette = 'aaas') +
  facet_wrap(vars(Type), nrow = 2) +
  labs(x = NULL, y = expression(log[10]~(OR+1)), fill = "FDR") +
  theme_pubr(base_size = 30) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.2/DEG_DMR_Overlap.jpeg",
       dpi = 320, width = 8, height = 12)
dev.off()

# Plot differential methylation
########################################################################################################################
# mutate results
results_meth_genes <- bind_rows(results_meth_genes, .id = "CancerCode")
results_meth_genes <- results_meth_genes %>%
  filter(adj.P.Val < 0.05) %>%
  # filter(abs(logFC) > 0.5) %>%
  group_by(CancerCode) %>%
  mutate(Status = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated")) %>%
  mutate(Case = paste0(CancerCode,"-",Status)) %>%
  mutate(Count = 1)
results_rnaseq_genes <- bind_rows(results_rnaseq_genes, .id = "CancerCode")
results_rnaseq_genes <- results_rnaseq_genes %>%
  filter(adj.P.Val < 0.05) %>%
  # filter(abs(logFC) > 0.5) %>%
  group_by(CancerCode) %>%
  mutate(Status = ifelse(logFC > 0, "Upregulated", "Downregulated")) %>%
  mutate(Case = paste0(CancerCode,"-",Status)) %>%
  mutate(Count = 1)

# plot meth results
df1 <- as.data.frame(table(results_meth_genes$CancerCode))
df1 <- df1[order(df1$Freq), ]
plotOrder <- factor(df1$Var1, levels = df1$Var1)
plot1 <- ggbarplot(data = df1, x = 'Var1', y = 'Freq', order = plotOrder, fill = "black") +
  labs(x = NULL, y = "DMG") +
  theme_pubclean(base_size = 35) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank())

df2 <- as.data.frame(table(results_meth_genes$CancerCode, results_meth_genes$Status))
df2$Var2 <- factor(df2$Var2, levels = c("Hypomethylated","Hypermethylated"))
plot2 <- ggbarplot(data = df2, x = 'Var1', y = 'Freq', fill = 'Var2', position = position_fill(),
                   order = plotOrder) +
  labs(x = NULL, y = "Fraction of DMG", fill = "Methylation") +
  scale_fill_jco(labels = c("Young","Old")) + 
  theme_pubclean(base_size = 35) +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank())

ggarrange(plot1,plot2, ncol = 1, align = "h", common.legend = T, heights = c(1.3,2))
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S2/differential_meth_gene_level.eps", dpi = 320)


# clusterProfiler Pathway Analysis
########################################################################################################################
# combine meth and RNASeq results
colnames(results_meth_genes)[3:ncol(results_meth_genes)] <- paste0(colnames(results_meth_genes)[3:ncol(results_meth_genes)], "_DNAm")
combined_results <- results_meth_genes %>% right_join(results_rnaseq_genes, by = c("CancerCode", "Gene"))
combined_results <- na.omit(data.frame(CancerCode = combined_results$CancerCode,
                                       Gene = combined_results$Gene, 
                                       DNAm = combined_results$Status_DNAm,
                                       Expression = combined_results$Status))
dna.rna.interaction <- combined_results %>%
  dplyr::group_by(CancerCode) %>%
  dplyr::count(DNAm, Expression)
openxlsx::write.xlsx(dna.rna.interaction,
                     "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - TCGA_DNAm-RNASeq_intersection_ALL.xlsx")

combined_results$Intersection <- paste0(combined_results$DNAm,"-",combined_results$Expression)
accepted <- c("Hypermethylated-Downregulated", "Hypomethylated-Upregulated")
combined_results <- combined_results[combined_results$Intersection %in% accepted, ]
dna.rna.interaction <- combined_results %>%
  dplyr::group_by(CancerCode) %>%
  dplyr::count(DNAm, Expression)
openxlsx::write.xlsx(dna.rna.interaction,
                     "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - TCGA_DNAm-RNASeq_intersection_FILTERED.xlsx")

# Genes up/downregulated across cancers
combined_results %>% dplyr::filter(Intersection == "Hypermethylated-Downregulated") %>% dplyr::count(Gene)
combined_results %>% dplyr::filter(Intersection == "Hypomethylated-Upregulated") %>% dplyr::count(Gene)

## Cluster Profiler
combined_results$CancerCode.Intersection <- paste0(combined_results$CancerCode, "-", combined_results$Intersection)
combined_results$CancerCode.Intersection <- as.character(combined_results$CancerCode.Intersection)
keys <- unique(combined_results$CancerCode.Intersection)

# convert to entrez
mart <- bitr(combined_results$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Gene", "ENTREZ")
combined_results <- merge(combined_results, mart, "Gene")

gene_list <- setNames(rep(list(NA), length(keys)), keys)
for (i in 1:length(keys)) {
  gene_list[[i]] <- as.character(combined_results$ENTREZ[combined_results$CancerCode.Intersection == keys[i]])
}
names(gene_list) <- gsub('\\-',' ', names(gene_list))

# Overall young vs old clusters
combined_results$Age <- ifelse(combined_results$Expression == "Upregulated", "Old", "Young")
formula_res <- compareCluster(ENTREZ~Age, data=combined_results, fun="enrichPathway")
dotplot(formula_res) + 
  scale_x_discrete(labels = wrap_format(10)) +
  scale_y_discrete(labels = wrap_format(40)) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(5,15)) +
  theme_pubclean(base_size = 25) +
  theme(legend.position = "right", 
        axis.ticks = element_blank())
# ggsave("/Users/Yajas/Downloads/RNA-DNAm_reactome_pathway_clusters.jpg", dpi = 320, width = 11, height = 8)

# Young vs old by tumor type
formula_res <- compareCluster(ENTREZ~CancerCode+Age, data=combined_results, fun="enrichPathway")
openxlsx::write.xlsx(formula_res@compareClusterResult,
                     "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - DEG_DMG_intersection_clusterprofiler.xlsx")
# filter results --> plot
tryRows <- c(1,6,12,18,25,26,29,37,47,57,60,61,65,66,68,69,71,80,87,89,92,95,98,101,102)
dat0 <- as.data.frame(formula_res)
dat0$Description <- factor(dat0$Description, levels = unique(dat0$Description))
dat0$GeneRatio <- round(as.numeric(sub("\\/.*", "", dat0$GeneRatio))/as.numeric(sub('.*/', '', dat0$GeneRatio)), 2)
dat0$Cluster <- gsub('\\.', ' ', dat0$Cluster)
dat0$Cluster <- factor(dat0$Cluster, levels = unique(dat0$Cluster))

ggplot(data = dat0[tryRows, ], aes(x = Cluster, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_size_continuous(range = c(5,15)) +
  scale_x_discrete(labels = wrap_format(4))+
  scale_y_discrete(labels = wrap_format(45)) +
  theme_pubclean(base_size = 25) +
  theme(legend.position = "right") +
  ylab(NULL) + xlab(NULL) + scale_color_viridis_c() + labs(color = "FDR")
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.2/RNA-DNAm_reactome_pathway_clusters.eps",
       dpi = 320, width = 16, height = 12)
