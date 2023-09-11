#               +--------------------------------+
#               |    analysis IGFBP_6 and PE     | 
#               +--------------------------------+
#read the Insulin-like growth factor-binding protein 6 GWAS into R
library(vroom)
IGFBP6_gwas <- vroom("F:\\Aging_gwas\\IGFBP_6\\IGFBP_6.csv",
                  col_names = TRUE)

colnames(IGFBP6_gwas)
IGFBP6_gwas$phenotype <- "IGFBP_6"

library(TwoSampleMR)
library(MRInstruments)

IGFBP6_gwas$neg_p <- -IGFBP6_gwas$LP
IGFBP6_gwas$p_value <- 10^IGFBP6_gwas$neg_p

##ÌáÈ¡p£¼5e-08µÄsnp
IGFBP6_p <- subset(IGFBP6_gwas, p_value < 5e-08)
