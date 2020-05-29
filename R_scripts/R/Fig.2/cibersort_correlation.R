library(gtools)
library(reshape2)
library(openxlsx)
library(TCGAutils)
library(matrixStats)
library(tidyverse)

source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/General Functions/immunophenotypeAnalysisFunctions.R")
tcga_signatures <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Immune Landscape of Cancer/Scores_160_Signatures.tsv",
                            header = T, stringsAsFactors = F, sep = "\t")
tcga_scores <- tcga_signatures[,3:ncol(tcga_signatures)]

sig_cancers <- c("BRCA", "COAD", "LGG", "OV", "THCA", "UCEC")

index_Cibersort <- which(tcga_signatures$Source == "Cibersort")
tcga_signatures <- tcga_signatures[index_Cibersort,]
tcga_scores <- tcga_scores[index_Cibersort,]
colnames(tcga_scores) <- gsub('\\.', '-', colnames(tcga_scores))

tcga_scores <- tcga_scores[ ,(TCGAbiospec(colnames(tcga_scores))$sample_definition == "Primary Solid Tumor" |
                                TCGAbiospec(colnames(tcga_scores))$sample_definition == "Primary Blood Derived Cancer - Peripheral Blood")]

# Remove duplicates if they exist
biospec_table <- TCGAbiospec(colnames(tcga_scores))
biospec_table <- biospec_table[order(biospec_table$plate,decreasing = TRUE), ]
tcga_scores <- tcga_scores[,colnames(tcga_scores)[order(biospec_table$plate, decreasing = TRUE)]]

if (sum(duplicated(biospec_table$submitter_id)) > 0){
  tcga_scores <- tcga_scores[, !duplicated(biospec_table$submitter_id)]
}

# change column names to match metadata
colnames(tcga_scores) <- substr(colnames(tcga_scores), 1, 12)
survival_data <- survival_data %>%
  dplyr::filter(Quartiles %in% c(1, 4),
                type %in% sig_cancers)

common_cases <- intersect(survival_data$bcr_patient_barcode, colnames(tcga_scores))
tcga_scores <- data.frame(tcga_scores[, common_cases])
survival_data <- survival_data[survival_data$bcr_patient_barcode %in% common_cases, ]

# Cibersort mean values in young vs old
################################################################################################

survival_data1 <- as.data.frame(survival_data)
tcga_scores1 <- as.data.frame(tcga_scores)

colnames(tcga_scores1) <- gsub('\\.','-',colnames(tcga_scores1))
tcga_scores1 <- tcga_scores1[, order(colnames(tcga_scores1))]
survival_data1 <- survival_data1 %>%
  dplyr::arrange(bcr_patient_barcode)

all(survival_data1$bcr_patient_barcode == colnames(tcga_scores1))


tumor_types <- c('BRCA','UCEC','OV','THCA','LGG','COAD')

# setup comparision - correlation
stat.type <- "pearson"
comp.type <- "continuous"

# get correlation matrix
for (i in 1:length(tumor_types)) {
  if (i == 1){
    # Temp for labelling rows
    temp_res <- deconvCorrelateStat(deconv_matrix = tcga_scores1,
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
    
    temp_res <- deconvCorrelateStat(deconv_matrix = tcga_scores1,
                                    comparision = comp.type,
                                    surv_data = survival_data,
                                    statistic = stat.type)
    
    # Write results
    cor_res_stat[,i] <- sapply(temp_res, "[[", "estimate")
    cor_res_p[,i] <- sapply(temp_res, "[[", "p.value")
  }
}

# adjust p for cell types of interest
# fdr for tumor types
cor_res_fdr <- apply(cor_res_p, 2, p.adjust, "fdr")

res <- cbind(melt(as.data.frame(cor_res_stat)), melt(cor_res_fdr))
colnames(res)[2] <- "estimate"
colnames(res)[5] <- "fdr"

res$Var1 <- as.character(res$Var1)
res$Var1 <- gsub('\\.',' ', res$Var1)
res$neglogfdr <- -log10(res$fdr)
res$Var1 <- factor(res$Var1)
res$Var2 <- factor(res$Var2)
res$direction <- factor(ifelse(res$estimate > 0, "Old", "Young"), levels = c("Young", "Old"))
res$absest <- abs(res$estimate)
res$combined <- paste0(res$Var1, " - ",res$Var2)
res <- na.omit(res)
res <- res %>% dplyr::arrange(estimate)
res$combined <- factor(res$combined, levels = res$combined)

ggbarplot(res[res$fdr < 0.05, ], x = 'combined', y = 'estimate', fill = 'direction', palette = "jco", position = position_dodge()) +
  scale_x_discrete(labels = function (x) stringr::str_wrap(x, width = 20)) +
  ylim(-0.31, 0.31) +
  scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.1),
                     labels = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3)) +
  labs(x = NULL, y = "Pearson (R)",
       fill = "Group") +
  coord_flip() +
  theme_bw(base_size = 25) +
  guides(fill = guide_legend(override.aes = list(size=15))) +
  theme(legend.position = c(0.75, 0.2),
        legend.title = element_blank(),
        panel.grid = element_blank())
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.2/CIBERSORT_correlation.eps",
       height = 12)
