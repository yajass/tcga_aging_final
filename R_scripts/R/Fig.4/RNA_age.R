# Load everything
################################################################
# load rna age prection model
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/Forked_Packages/RNAAgeCalc-master/R/predict_age.R")
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/Forked_Packages/RNAAgeCalc-master/R/makeplot.R")
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Forked_Packages/RNAAgeCalc-master/R/sysdata.rda")

# load packages
library(impute)
library(rstatix)
library(ggpubr)
library(dplyr)
library(stringr)

# parametes
sig = "Dev"

# load tcga metadata
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")

# Schultz batch correct data
folder <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Schultz/"
tcgaindex <- grepl(pattern = "tcga", list.files(folder))

# Setup data frame for results
################################################################
res <- data.frame(Barcode = character(),
                  RNAAge = integer(),
                  ChronAge = integer(),
                  AgeAccelResid = integer(),
                  Tissue = factor(),
                  Status = factor(),
                  Dataset = factor())

# impute tcga samples
################################################################
for (i in list.files(folder)[tcgaindex]) {
  
  # extract info
  tissue <- str_split(i, fixed("-"))[[1]][1]
  tumor.type <- str_split(i, fixed("-"))[[1]][2]
  status <- ifelse(grepl(pattern = "tcga-t", i), "Tumor", "Normal")
  
  # load and preproc data
  dat0 <- read.table(paste0(folder,i), stringsAsFactors = FALSE, header = T, row.names = 1)[,-1]
  colnames(dat0) <- gsub('\\.','-',substring(colnames(dat0),1,12))
  dat0 <- dat0[,!duplicated(colnames(dat0))]
  
  # tcga metadata
  metadata <- data.frame(UID = survival_data$bcr_patient_barcode[survival_data$bcr_patient_barcode %in% colnames(dat0)],
                         chonage = survival_data$age_at_initial_pathologic_diagnosis[survival_data$bcr_patient_barcode %in% colnames(dat0)],
                         stringsAsFactors = FALSE)
  metadata <- metadata[order(metadata$UID), ]
  
  # make sure tcga fpkm is fine
  dat0 <- dat0[, metadata$UID]
  dat0 <- dat0[, sort(colnames(dat0))]
  
  # check if everthything matches up
  if (!all(metadata$UID == colnames(dat0))) print(paste0("Error: ",i, " does not match up"))
  
  # get age
  res0 <- predict_age(exprdata = dat0, tissue = tissue, exprtype = "FPKM", idtype = "SYMBOL", 
                      signature = sig, stype = "caucasian", chronage = metadata)
  res0 <- tibble::rownames_to_column(res0, "Barcode")
  res0$Barcode <- gsub('\\.','-',res0$Barcode)
  # res0 <- res0[,-3]
  res0$Tissue <- tissue
  res0$Status <- status
  res0$Dataset <- "TCGA"
  
  # merge with existing results
  res <- data.frame(plyr::rbind.fill(res, res0), stringsAsFactors = FALSE)
  
  print(paste0("Finished: ", i))
}

# save tcga samples
res2 <- res
write.xlsx(res, "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - TCGA RNA Age.xlsx")

# impiute GTEx
################################################################
gtexindex <- grepl(pattern = "gtex", list.files(folder))
meta2 <- read.csv("/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t")
for (i in list.files(folder)[gtexindex]) {
  
  # extract info
  tissue <- str_split(i, fixed("-"))[[1]][1]
  tumor.type <- str_split(i, fixed("-"))[[1]][2]
  status <- "GTEx"
  
  # load and preproc data
  dat0 <- read.table(paste0(folder,i), stringsAsFactors = FALSE, header = T, row.names = 1)[,-1]
  colnames(dat0) <- gsub("^([^-]*-[^-]*)-.*$",
                         "\\1",
                         gsub('\\.', '-', colnames(dat0)))
  dat0 <- dat0[,!duplicated(colnames(dat0))]
  
  # tcga metadata
  metadata <- data.frame(Barcode = meta2$SUBJID[meta2$SUBJID %in% colnames(dat0)],
                         chonage = meta2$AGE[meta2$SUBJID %in% colnames(dat0)],
                         stringsAsFactors = FALSE)
  metadata <- metadata[order(metadata$Barcode), ]
  
  # make sure tcga fpkm is fine
  dat0 <- dat0[, metadata$Barcode]
  dat0 <- dat0[, sort(colnames(dat0))]
  
  # check if everthything matches up
  if (!all(metadata$Barcode == colnames(dat0))) print(paste0("Error: ",i, " does not match up"))
  
  # get age
  res0 <- predict_age(exprdata = dat0, tissue = tissue, exprtype = "FPKM", idtype = "SYMBOL", 
                      signature = sig, stype = "caucasian")
  res0 <- tibble::rownames_to_column(res0, "Barcode")
  res0$Barcode <- gsub('\\.','-',res0$Barcode)
  res0$ChronAge <- metadata$chonage
  
  res0$Tissue <- tissue
  res0$Status <- status
  res0$Dataset <- "GTEx"
  
  # merge with existing results
  res <- data.frame(plyr::rbind.fill(res, res0), stringsAsFactors = FALSE)
  
  print(paste0("Finished: ", i))
}

write.xlsx(res, "/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Tables/Supplementary Table SXX - GTEx RNA Age.xlsx")

# Discretize results
################################################################
res$ChronAge[res$Dataset == "TCGA"] <- cut(as.numeric(res$ChronAge[res$Dataset == "TCGA"]),
                                           breaks = c(10,20,30,40,50,60,70,80,90,100),
                                           labels = F)
res$ChronAge[res$Dataset == "TCGA"] <- ifelse(res$ChronAge[res$Dataset == "TCGA"] == 1, "10-19",
                                              ifelse(res$ChronAge[res$Dataset == "TCGA"] == 2, "20-29",
                                                     ifelse(res$ChronAge[res$Dataset == "TCGA"] == 3, "30-39",
                                                            ifelse(res$ChronAge[res$Dataset == "TCGA"] == 4, "40-49",
                                                                   ifelse(res$ChronAge[res$Dataset == "TCGA"] == 5, "50-59",
                                                                          ifelse(res$ChronAge[res$Dataset == "TCGA"] == 6, "60-69",
                                                                                 ifelse(res$ChronAge[res$Dataset == "TCGA"] == 7, "70-79",
                                                                                        ifelse(res$ChronAge[res$Dataset == "TCGA"] == 8, "80-89",
                                                                                               "90-99"))))))))
res$ChronAge <- factor(res$ChronAge, levels = sort(unique(res$ChronAge)))

ggscatter(data = res,
          x = "ChronAge", y = "RNAAge", title = sig) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_grid(vars(Tissue), vars(Status), scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################################################################
ggscatter(data = res2,
          x = "ChronAge", y = "RNAAge", title = "RNA Age") +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_grid(vars(Tissue), vars(Status), scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

res2$AgeAccelerationDiff <- res2$RNAAge - res2$ChronAge
metadata <- data.frame(Barcode = survival_data$bcr_patient_barcode,
                       Quartile = survival_data$Quartiles,
                       stringsAsFactors = FALSE)
res2 <- merge(res2, metadata, "Barcode")
res2$Quartile <- ifelse(res2$Quartile == 1, "Young",
                        ifelse(res2$Quartile == 4, "Old", "Middle Aged"))
res2$Quartile <- factor(res2$Quartile, levels = c("Young", "Middle Aged", "Old"))
paired_barcodes <- res2 %>%
  dplyr::count(Barcode) %>%
  dplyr::filter(n > 1) %>%
  dplyr::select(Barcode)



# residual normality
ggqqplot(data = res2, x = 'AgeAccelResid') +
  facet_grid(Quartile ~ Status + Tissue)
res2 %>%
  filter(Tissue != "uterus" & Status != "Normal") %>%
  group_by(Quartile,Status, Tissue) %>%
  shapiro_test(AgeAccelResid) %>%
  rstatix::adjust_pvalue(output.col = "p.adj", method = "fdr")

# plots
signif <- res2 %>% 
  wilcox_test(AgeAccelerationDiff ~ Quartile, comparisons = list(c("Young", "Old")),p.adjust.method = "fdr")
signif2 <- res2[res2$Barcode %in% paired_barcodes$Barcode, ] %>% 
  group_by(Quartile) %>%
  wilcox_test(AgeAccelerationDiff ~ Status, comparisons = list(c("Tumor", "Normal")), paired = TRUE) %>%
  rstatix::adjust_pvalue(method = "fdr") %>% 
  rstatix::add_significance("p.adj")
res2$Status <- factor(res2$Status, levels = c("Tumor", "Normal"))
ggline(
  data = res2,
  x = "Quartile",
  y = "AgeAccelerationDiff",
  ylab = "RNA AgeAccelerationDiff",
  add = "mean_ci",
  color = "Status",
  linetype = 2,
  size = 2,
  palette = "aaas",
)  + 
  scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) +
  labs(x = NULL) +
  stat_pvalue_manual(data = signif, y.position = 17, size = 15, bracket.size = 2, vjust = 0.6) +
  stat_pvalue_manual(data = signif2, y.position = c(10,-5,-20), size = 5, x = "Quartile") +
  theme_pubr(base_size = 35)
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.4/AgeAccDiff(Unpaired Samples).eps",
       dpi = 320)



ggpaired(data = res2[res2$Barcode %in% paired_barcodes$Barcode & res2$Quartile != "Middle Aged", ], 
         x = 'Status', y = 'AgeAccelerationDiff', fill = 'Quartile', id = 'Barcode',
         line.color = "gray80", palette = 'jco', facet.by = 'Quartile') +
  labs(x = NULL, y = "RNA AgeAccelerationDiff", fill = "Age") +
  stat_compare_means(method = "wilcox", paired = TRUE, label = "p.signif",label.x.npc = 0.5, label.y.npc = 0.7, hide.ns = TRUE, size = 15) +
  theme_pubr(base_size = 35) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.4/PairedAgeAccDiff(ggpaired).eps",
       dpi = 320)


tcga <- ggscatter(data = res2,
                  x = "ChronAge", y = "RNAAge") +
  geom_smooth(method = "lm") +
  labs(x = NULL, y = "RNA Age") +
  stat_cor(method = "pearson", size = 7) +
  facet_grid(vars(Tissue), vars(Status), scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_pubr(base_size = 35) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gtex <- ggscatter(data = res[res$Dataset == "GTEx", ],
                  x = "ChronAge", y = "RNAAge") +
  labs(x = NULL, y = NULL) +
  facet_grid(vars(Tissue), vars(Status), scales = "free") +
  theme_pubr(base_size = 35) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

jpeg("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.S4/CorrelationScatterPlot.jpeg",
     width = 18, height = 10, units = "in", res = 320)
annotate_figure(ggarrange(tcga, gtex, nrow = 1, widths = c(2,1), align = "h"),
                bottom = text_grob("Chronological Age", size = 35, vjust = -0.25))
dev.off()
