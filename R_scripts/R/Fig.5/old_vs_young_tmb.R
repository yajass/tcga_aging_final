library(gtools)
library(maftools)
library(dplyr)
library(rstatix)

tcga.cohort = system.file('extdata', 'tcga_cohort.txt.gz', package = 'maftools')
tcga.cohort = data.table::fread(cmd = paste('zcat <', tcga.cohort), sep = '\t', stringsAsFactors = FALSE)
tcga.cohort$total <- tcga.cohort$total/50

source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
colnames(tcga.cohort)[2] <- "bcr_patient_barcode"
combinedDF <- merge(x = tcga.cohort, y = survival_data, by = 'bcr_patient_barcode')

combinedDF <- combinedDF[combinedDF$cohort %in% c("BRCA","UCEC","COAD","LGG","OV","THCA"), ]
combinedDF <- combinedDF[combinedDF$Quartiles %in% c(1,4), ]
combinedDF$Quartiles <- ifelse(combinedDF$Quartiles == 1, "Young", "Old")
combinedDF$Quartiles <- factor(combinedDF$Quartiles, levels = c("Young", "Old"))
combinedDF$cohort <- factor(combinedDF$cohort, levels = rev(c("THCA", "LGG", "BRCA", "OV", "UCEC", "COAD")))


stat.test <- combinedDF %>%
  group_by(cohort) %>%
  wilcox_test(total ~ Quartiles) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(y.position = 3)
stat.test


p <- ggboxplot(data = combinedDF, x = 'Quartiles', y = 'total', fill = 'Quartiles', palette = "jco") +
  facet_wrap(~ cohort, nrow = 1) +
  scale_y_log10(name = "TMB",
                breaks = c(0.1,1,10,100,1000),
                labels = c(expression(10^1),expression(10^2),expression(10^3),
                           expression(10^4),expression(10^5)),
                limits = c(0.01,1000)) +
  labs(fill = 'Age') +
  theme_pubr(base_size =  35) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "lines")
  )


stat.test$rounded <-  round(stat.test$p.adj, 10)
stat.test$fdrstar <- stars.pval(p.value = stat.test$rounded)
p + stat_pvalue_manual(stat.test, label = "fdrstar", hide.ns = T, label.size = 10)
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.5/old_vs_young_tmb.eps",
       dpi = 320, width = 18)
