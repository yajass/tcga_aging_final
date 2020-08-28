

source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/General Functions/immunophenotypeAnalysisFunctions.R")
enrich_result <- read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Aging Enrichment/tcga_tpm_ssgsea_enrichment.txt",
                            header = T, stringsAsFactors = F)
colnames(enrich_result) <- gsub('\\.','-',colnames(enrich_result))
# Try t-test/wilcox
################################################################################################################################
# setup comparision 
stat.type <- "t.test"
tumor_types <- c('BRCA','UCEC','OV','THCA','LGG','COAD')

# get stat matrix
# function does NOT inherently adjust p values
for (i in 1:length(tumor_types)) {
  if (i == 1){
    # Temp for labelling rows
    temp_res <- deconvCompareStat(deconv_matrix = enrich_result,
                                  surv_data = survival_data,
                                  statistic = stat.type) %>%
      dplyr::mutate(Type = tumor_types[i])
  }
  
  if (i > 1){
    temp_res <- rbind(temp_res, deconvCompareStat(deconv_matrix = enrich_result,
                                                  surv_data = survival_data,
                                                  statistic = stat.type) %>%
                        dplyr::mutate(Type = tumor_types[i]))
  }
}


temp_res <- temp_res %>%
  dplyr::group_by(Type) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  # dplyr::filter(p.adj < 0.05) %>%
  dplyr::mutate(neglogfdr = -log10(p.adj))

temp_res$Type <- factor(temp_res$Type)
temp_res$variable <- factor(temp_res$variable)
if (stat.type == "t.test"){
  temp_res$Higher_Med <- ifelse(temp_res$statistic > 0, "Old", "Young")
}
temp_res$Higher_Med <- factor(temp_res$Higher_Med, levels = c("Young", "Old"))
temp_res$variable <- gsub('\\_',' ', temp_res$variable)
temp_res$variable <- ifelse(temp_res$variable == "GO MULTICELLULAR ORGANISM AGING",
                            "GO MULTI-\nCELLULAR ORGANISM AGING", temp_res$variable)
temp_res$variable <- ifelse(temp_res$variable == "REACTOME FORMATION OF SENESCENCE ASSOCIATED HETEROCHROMATIN FOCI SAHF",
                            "REACTOME FORMATION OF SENESCENCE ASSOCIATED HETERO-\nCHROMATIN FOCI SAHF", temp_res$variable)

ggscatter(data = temp_res[temp_res$variable != "COURTOIS SENESCENCE TRIGGERS", ], x = 'variable', y = 'neglogfdr',
          shape = 'Higher_Med', color = 'Type',
          size = 'neglogfdr', palette = "jama", stroke = 5) +
  geom_hline(yintercept = -log10(0.05), lty = 2, lwd = 1) +
  scale_size(range = c(1, 20), guide = FALSE) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x = NULL, y = expression(-log[10]~FDR~Young/Old), shape = "Enrichment", color = "Cancer Type") +
  # coord_flip() +
  theme_pubclean(44) +
  theme(legend.position = "top",
        legend.direction = "horizontal") +
  guides(color = guide_legend(override.aes = list(size = 20),
                              nrow = 1, by.row = TRUE),
         shape = guide_legend(override.aes = list(size = 20)))
ggsave("/Users/Yajas/Documents/Elemento/tcga_aging_final/Results/Figures/Fig.2/Senescence_ssGSEA_TPM_TCGA_new.eps", dpi = 320,
       height = 14.5, width = 44)

