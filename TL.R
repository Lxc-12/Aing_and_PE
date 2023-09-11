#               +---------------------------+
#               |     analysis TL and PE    | 
#               +---------------------------+
## import telomere length GWAS data into R.
library(vroom)
TL_gwas <- vroom("F:\\Aging_gwas\\Telomere length\\TL.csv",
                col_names = TRUE)

colnames(TL_gwas)
TL_gwas$phenotype <- "telomere length"

library(TwoSampleMR)
library(MRInstruments)

TL_gwas$neg_p <- -TL_gwas$LP
TL_gwas$p_value <- 10^TL_gwas$neg_p

## extract SNPs with p£¼5e-08.
TL_p <- subset(TL_gwas, p_value < 5e-08)

## format the data
TL_format <- format_data(TL_p, type = "exposure",
                         snps = NULL,
                         header = TRUE,
                         phenotype_col = "phenotype",
                         snp_col = "ID",
                         beta_col = "ES",
                         se_col = "SE",
                         eaf_col = "AF",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF",
                         pval_col = "p_value",
                         units_col = "units",
                         ncase_col = "ncase",
                         ncontrol_col = "ncontrol",
                         samplesize_col = "N",
                         gene_col = "gene",
                         id_col = "id",
                         min_pval = 0,
                         z_col = "z",
                         info_col = "info",
                         chr_col = "CHROM",
                         pos_col = "POS",
                         log_pval = FALSE)

## clump the data
exposure_data <- clump_data(TL_format, clump_kb = 10000, clump_r2 = 0.001)

## extract those SNPs included in exposure_data on PE GWAS data.
outcome_data <- read_outcome_data(snps = exposure_data$SNP,
                                 filename = "F:\\PE_gwas\\PE_gwas.txt",
                                 sep = "\t",
                                 snp_col = "SNP",
                                 beta_col = "Effect",
                                 se_col = "StdErr",
                                 effect_allele_col = "Allele1",
                                 other_allele_col = "Allele2",
                                 eaf_col = "Freq1",
                                 pval_col = "P_value",
                                 units_col = "Units",
                                 gene_col = "Gene",
                                 samplesize_col = "n")

## to harmonize data
dat <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
## MR analysis
mr(dat)
generate_odds_ratios(mr_res = mr(dat))
## forest plot
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
## scatter plot
mr_scatter_plot(mr_results = mr(dat),dat)
## heterogeneity test
mr_heterogeneity(dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
## pleiotropy test
mr_pleiotropy_test(dat)
## sensitivity analysis
mr_leaveoneout(dat)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))

## MR-PRESSO analysis
sumarydat <- subset(dat, mr_keep=='TRUE')
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE, data = sumarydat, NbDistribution = 3000,  
          SignifThreshold = 0.05)