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

cancer_code <- c("BRCA","COAD","THCA","UCEC","LGG","OV")

# collect results
################################################################################
# results DF
res <- setNames(data.frame(matrix(ncol = 9)),
                c("Hugo_Symbol", "Young", "Old", "pval",
                  "or", "ci.up", "ci.low", "adjPval", "cancercode"))
# Sample Sumamry DF
ssumary <- setNames(data.frame(matrix(ncol = 3)),  c("Cohort","SampleSize","cancercode"))
# Collet results
collectRes <- list(results = res, SampleSummary = ssumary)


for (i in 1:length(cancer_code)) {
  
  # download MAF
  tumor_maf <- GDCquery_Maf(tumor = cancer_code[i],
                            pipelines = "mutect2")
  
  # Download clinical data
  tumor_clinical <- survival_data %>% dplyr::filter(type == cancer_code[i])
  
  # Quartile age
  colnames(tumor_clinical)[2] <- "Tumor_Sample_Barcode"
  tumor_clinical$Quartiles <- ifelse(tumor_clinical$Quartiles == 1, "Young",
                                     ifelse(tumor_clinical$Quartiles == 4, "Old", NA))
  tumor_clinical <- tumor_clinical[!is.na(tumor_clinical$Quartiles), ]
  
  # Get young and old barcodes
  young_samples <- tumor_clinical$Tumor_Sample_Barcode[tumor_clinical$Quartiles == "Young"]
  old_samples <- tumor_clinical$Tumor_Sample_Barcode[tumor_clinical$Quartiles == "Old"]
  
  # Subset the downloaded maf
  tumor_maf$Tumor_Sample_Barcode <- substr(tumor_maf$Tumor_Sample_Barcode,1,12)
  young_maf <- tumor_maf[tumor_maf$Tumor_Sample_Barcode %in% young_samples, ]
  old_maf <- tumor_maf[tumor_maf$Tumor_Sample_Barcode %in% old_samples, ]
  
  # Setup results
  if (i == 1){
    overall_young_maf <- young_maf
    overall_old_maf <- old_maf
    young_clinical <- tumor_clinical[tumor_clinical$Quartiles == "Young", ]
    old_clinical <- tumor_clinical[tumor_clinical$age_at_diagnosis == "Old", ]
  }
  if (i != 1){
    overall_young_maf <- plyr::rbind.fill(overall_young_maf, young_maf)
    overall_old_maf <- plyr::rbind.fill(overall_old_maf, old_maf)
    young_clinical <- plyr::rbind.fill(young_clinical, tumor_clinical[tumor_clinical$Quartiles == "Young", ])
    old_clinical <- plyr::rbind.fill(old_clinical, tumor_clinical[tumor_clinical$Quartiles == "Old", ])
  }
  
  # Make MAFs for young/old for current tumor type
  young_maf_compare <- read.maf(maf = young_maf, clinicalData = tumor_clinical[tumor_clinical$Quartiles == "Young", ])
  old_maf_compare <- read.maf(maf = old_maf, clinicalData = tumor_clinical[tumor_clinical$Quartiles == "Old", ])
  
  # Go to nexr loop iteration if theres a fisher test error
  possibleError <- tryCatch(
    young_vs_old <- mafCompare(m1 = young_maf_compare, m2 = old_maf_compare, m1Name = 'Young', m2Name = 'Old', minMut = 5),
    error=function(e) e
  )
  if(inherits(possibleError, "error")) next
  
  # Considering only genes which are mutated in at-least in 5
  # samples in one of the cohort to avoid bias due to genes mutated in single sample.
  young_vs_old <- mafCompare(m1 = young_maf_compare, m2 = old_maf_compare, 
                             m1Name = 'Young', m2Name = 'Old', minMut = 5)
  young_vs_old$results <- data.frame(young_vs_old$results)
  young_vs_old$SampleSummary <- data.frame(young_vs_old$SampleSummary)
  
  young_vs_old$results[,9] <- cancer_code[i]
  colnames(young_vs_old$results)[9] <- "cancercode"
  young_vs_old$SampleSummary[,3] <- toupper(cancer_code[i])
  colnames(young_vs_old$SampleSummary)[3] <- "cancercode"
  
  # Aggregate results
  collectRes$results <- rbind(collectRes$results, young_vs_old$results)
  collectRes$SampleSummary <- rbind(collectRes$SampleSummary, young_vs_old$SampleSummary)
}
collectRes$results <- collectRes$results[-1, ]
collectRes$SampleSummary <- collectRes$SampleSummary[-1, ]

write.xlsx(collectRes$results[collectRes$results$adjPval < 0.05, ],
           "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - Young vs Old Mutations Fisher.xlsx")

# Significant results by age and tumor type
na.omit(collectRes$results) %>% dplyr::group_by(cancercode) %>% dplyr::summarize(total = length(adjPval), 
                                                                                 SigYoung = sum(adjPval < 0.05 & or > 1),
                                                                                 SigOld = sum(adjPval < 0.05 & or < 1))

setDT(collectRes$results)
setDT(collectRes$SampleSummary)

# Get driver genes
driver.gene.compare <- collectRes
pancan.drivers <- drivers$Gene[drivers$Cancer == "PANCAN"]
driver.gene.compare$results$Hugo_Symbol <- ifelse(driver.gene.compare$results$Hugo_Symbol %in% pancan.drivers,
                                                  paste0(driver.gene.compare$results$Hugo_Symbol, "-", "PANCAN"),
                                                  driver.gene.compare$results$Hugo_Symbol)
pancan.indices <- which(stringr::str_detect(driver.gene.compare$results$Hugo_Symbol, "-"))
driver.gene.compare$results$Hugo_Symbol[-pancan.indices] <- paste0(driver.gene.compare$results$Hugo_Symbol[-pancan.indices],
                                                                   "-",
                                                                   driver.gene.compare$results$cancercode[-pancan.indices])
drivers$Combined <- paste0(drivers$Gene, "-", drivers$Cancer)
driver.gene.compare$results <- merge(driver.gene.compare$results,
                                     drivers, 
                                     by.x = "Hugo_Symbol",
                                     by.y = "Combined")
pancan.indices <- which(stringr::str_detect(driver.gene.compare$results$Hugo_Symbol, "-"))
driver.gene.compare$results$Hugo_Symbol <- stringr::str_remove(driver.gene.compare$results$Hugo_Symbol, "-PANCAN")
driver.gene.compare$results$Hugo_Symbol[pancan.indices] <- paste0(driver.gene.compare$results$Hugo_Symbol[pancan.indices],
                                                                  "-",
                                                                  driver.gene.compare$results$cancercode[pancan.indices])
drivers.to.plot <- driver.gene.compare$results       
drivers.to.plot[,5:7] <- -log10(drivers.to.plot[,5:7])
drivers.to.plot$Enrichment <- ifelse(drivers.to.plot$or < 0, "Young", "Old")
drivers.to.plot$Enrichment <- factor(drivers.to.plot$Enrichment, levels = c("Young", "Old"))
drivers.to.plot$cancercode <- factor(drivers.to.plot$cancercode, levels = unique(drivers.to.plot$cancercode))
# drivers.to.plot$Hugo_Symbol <- gsub("\\-.*","",drivers.to.plot$Hugo_Symbol)
drivers.to.plot <- drivers.to.plot[order(drivers.to.plot$or), ]
drivers.to.plot$Hugo_Symbol <- factor(drivers.to.plot$Hugo_Symbol, levels = unique(drivers.to.plot$Hugo_Symbol))


ax.label <- expression(log[10]~Odds~Ratio)

plt <- ggplot(data = drivers.to.plot[drivers.to.plot$adjPval < 0.005,],
              aes(x = Hugo_Symbol, y = or, ymin = ci.low, ymax = ci.up)) + 
  geom_pointrange(aes(color = Enrichment), shape = 18, size = 2, position = position_dodge(width = 1)) +
  geom_errorbar(aes(color = Enrichment, ymin=ci.low, ymax=ci.up),width=0.5, position = position_dodge(width = 1)) +
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5,1.5)) +
  geom_hline(yintercept=0, lty=2) +
  coord_flip() +
  labs(color = "Enrichment", y = ax.label, x = NULL) +
  scale_color_jco() + theme_pubclean(base_size = 35) +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "none")

p.df <- data.frame(Gene = drivers.to.plot[drivers.to.plot$adjPval < 0.005,]$Hugo_Symbol,
                   pval = drivers.to.plot[drivers.to.plot$adjPval < 0.005,]$pval,
                   p.stars =stars.pval(drivers.to.plot[drivers.to.plot$adjPval < 0.005,]$pval))

drivers.to.plot$p.stars <- stars.pval(drivers.to.plot$pval)

# tables
tab_base <- ggplot(drivers.to.plot[drivers.to.plot$adjPval < 0.005,], aes(y=Hugo_Symbol)) +
  ylab(NULL) + xlab("  ") +
  theme(plot.title = element_text(hjust = 0.5, size=28), ## centering title on text
        axis.text.x=element_text(color="white"), ## need text to be printed so it stays aligned with figure but white so it's invisible
        axis.line=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
tab1 <- tab_base + 
  geom_text(aes(x=1, label=Young), size = 10) + 
  ggtitle("Young")
tab2 <- tab_base + 
  geom_text(aes(x=1, label=Old), size = 10) + 
  ggtitle("Old")
tab3 <- tab_base + 
  geom_text(aes(x=1, label=p.stars), size = 10) + 
  ggtitle("p")


ggarrange(plt, tab1,tab2,tab3, widths = c(15,1,1,1), nrow = 1, align = "h")
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/pancan_forest_plot.eps",
       dpi = 320, width = 15)


drivers.to.plot$Hugo_Symbol <- as.character(drivers.to.plot$Hugo_Symbol)
ucec_drivers <- drivers.to.plot %>% 
  filter(adjPval < 0.05) %>% 
  filter(Enrichment == "Young") %>%
  pull(Hugo_Symbol)
ucec_drivers <- ucec_drivers[endsWith(ucec_drivers, "-UCEC")]
ucec_drivers <- sub("-UCEC", "", ucec_drivers)
ucec_drivers <- sub("-UCEC", "", ucec_drivers)

# cluster profiler
library(clusterProfiler)
library(ReactomePA)
mart <- bitr(ucec_drivers, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Hugo_Symbol", "ENTREZ")
ucec_drivers <- mart$ENTREZ

results <- enrichPathway(gene = ucec_drivers, organism = "human")

jpeg("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/ucec_mutations_young_reactome.jpeg",
     height = 7, width = 11, res = 320, units = "in")
emapplot(results, label = "none") + scale_color_viridis_c() + labs(color = "FDR")
dev.off()

write.xlsx(results@result[results@result$p.adjust < 0.05, ],
           "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - UCEC young patient driver mutation reactome pathway analysis.xlsx")

# repeat for old
ucec_drivers <- drivers.to.plot %>% 
  filter(adjPval < 0.05) %>% 
  filter(Enrichment == "Old") %>%
  pull(Hugo_Symbol)
ucec_drivers <- ucec_drivers[endsWith(ucec_drivers, "-UCEC")]
ucec_drivers <- sub("-UCEC", "", ucec_drivers)
ucec_drivers <- sub("-UCEC", "", ucec_drivers)

mart <- bitr(ucec_drivers, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Hugo_Symbol", "ENTREZ")
ucec_drivers <- mart$ENTREZ

results <- enrichPathway(gene = ucec_drivers, organism = "human")

jpeg("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S5/ucec_mutations_old_reactome.jpeg",
     height = 10, width = 11, res = 320, units = "in")
emapplot(results, label = "none") + scale_color_viridis_c() + labs(color = "FDR")
dev.off()

write.xlsx(results@result[results@result$p.adjust < 0.05, ],
           "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - UCEC old patient driver mutation reactome pathway analysis.xlsx")
