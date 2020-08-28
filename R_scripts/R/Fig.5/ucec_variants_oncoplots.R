library(gtable)
library(grid)
library(gtools)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(TCGAbiolinks)
library(maftools)

source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
drivers <- readxl::read_xlsx("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA Drivers/driver.genes.xlsx",
                             sheet = 2, trim_ws = TRUE, skip = 3)[,1:2]
survival_data <- survival_data %>%
  dplyr::filter(type == "UCEC") %>%
  dplyr::filter(Quartiles %in% c(1,4)) %>%
  dplyr::mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))
colnames(survival_data)[2] <- "Tumor_Sample_Barcode"

tumor_maf <- GDCquery_Maf(tumor = "UCEC",
                          pipelines = "mutect2")
tumor_maf$Tumor_Sample_Barcode <- substr(tumor_maf$Tumor_Sample_Barcode,1,12)


# Get young and old barcodes
young_samples <- survival_data$Tumor_Sample_Barcode[survival_data$Quartiles == "Young"]
old_samples <- survival_data$Tumor_Sample_Barcode[survival_data$Quartiles == "Old"]

# clinical subsets
young_clinical <- survival_data %>% filter(Quartiles == "Young")
old_clinical <- survival_data %>% filter(Quartiles == "Old")

# maf subsets
young_maf <- tumor_maf %>% filter(Tumor_Sample_Barcode %in% young_samples)
old_maf <- tumor_maf %>% filter(Tumor_Sample_Barcode %in% old_samples)

# compile maf
young_maf <- read.maf(maf = young_maf, clinicalData = young_clinical)
old_maf <- read.maf(maf = old_maf, clinicalData = old_clinical)

# compare mafs --> forest plot
compare_maf <- mafCompare(m1 = young_maf, m2 = old_maf, m1Name = "Young", m2Name = "Old", minMut = 5)

# co oncoplot
gns <- c("PTEN", "PIK3CA", "TP53", "CTCF", "KRAS")
postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/ucec_cooncoplot.eps",
           width = 6, height = 3)
coOncoplot(m1 = young_maf, m2 = old_maf, genes = gns, m1Name = "Young", m2Name = "Old",sepwd_genes2 = 0,
           sepwd_genes1 = 0)
dev.off()

# co lollipop
postscript("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/ucec_coolollipoplot.eps",
           width = 10, height = 4)
lollipopPlot2(m1 = young_maf, m2 = old_maf, gene = "PTEN", m1_name = "Young", m2_name = "Old",
              AACol1 = "HGVSp", AACol2 = "HGVSp", legendTxtSize = 2, pointSize = 2, domainLabelSize = 1.5)
dev.off()


# some pathways
pdf("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S5/ucec_mut_young_RAS.pdf", height = 9.5)
maftools::PlotOncogenicPathways(young_maf, "RTK-RAS")
mtext(text = "Young", adj = 0.8, padj = 1, outer = TRUE, line = -1, cex = 1.5)
dev.off()
pdf("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S5/ucec_mut_old_RAS.pdf", height = 9.5)
maftools::PlotOncogenicPathways(old_maf, "RTK-RAS")
mtext(text = "Old", adj = 0.8, padj = 1, outer = TRUE, line = -1, cex = 1.5)
dev.off()
pdf("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S5/ucec_mut_young_PI3K.pdf", height = 9.5)
maftools::PlotOncogenicPathways(young_maf, "PI3K")
mtext(text = "Young", adj = 0.8, padj = 1, outer = TRUE, line = -1, cex = 1.5)
dev.off()
pdf("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S5/ucec_mut_old_PI3K.pdf", height = 9.5)
maftools::PlotOncogenicPathways(old_maf, "PI3K")
mtext(text = "Old", adj = 0.8, padj = 1, outer = TRUE, line = -1, cex = 1.5)
dev.off()


# compare boxplots
young <- data.frame(young_maf@variant.classification.summary)
young$Age <- "Young"
old <- data.frame(old_maf@variant.classification.summary)
old$Age <- "Old"
combined <- rbind(young, old)
colnames(combined) <- gsub('\\_',' ',colnames(combined))
combined <- combined[, -c(1,11)]
combined <- melt(combined)
combined$Age <- factor(combined$Age, levels = c("Young", "Old"))
compareTest <- combined %>%
  dplyr::group_by(variable) %>%
  rstatix::wilcox_test(data = ., formula = value ~ Age)
compareTest$p.adj <- p.adjust(compareTest$p, method = "fdr")
compareTest$fdrstar <- stars.pval(compareTest$p.adj)

ggboxplot(combined, x = 'Age', y = 'value', fill = 'Age', palette = 'jco', ylab = "Frequency",
          title = NULL) + 
  facet_wrap(~variable, labeller = label_wrap_gen(width=15), scales = "free_y") + scale_y_log10() + 
  stat_pvalue_manual(data = compareTest, y.position = c(2.5,3.2,1.5,3,1.5),
                     label = "fdrstar", hide.ns = T, size = 15, vjust = 1, bracket.size = 0) +
  theme_pubr(base_size = 25) + theme(axis.title.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     legend.position = "none")
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/ucec_variants.eps")
