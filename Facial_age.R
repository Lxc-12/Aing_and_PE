#               +---------------------------------+
#               |  analysis Facial Aging and PE   | 
#               +---------------------------------+
## import the Facial Aging GWAS data into R.
library(vroom)
FA_gwas <- vroom("F:\\Aging_gwas\\Facial Aging\\Facial_Aging.csv",
                     col_names = TRUE)

colnames(FA_gwas)
FA_gwas$phenotype <- "Facial_Aging"

library(TwoSampleMR)
library(MRInstruments)

## Convert -log10 p-value to P-value
FA_gwas$neg_p <- -FA_gwas$LP
FA_gwas$p_value <- 10^FA_gwas$neg_p

## extract SNPs with p��5e-08.
FA_p <- subset(FA_gwas, p_value < 5e-08)

## format the data.
FA_format <- format_data(FA_p, type = "exposure",
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
                             min_pval = 1e-500,
                             z_col = "z",
                             info_col = "INFO",
                             chr_col = "CHROM",
                             pos_col = "POS",
                             log_pval = FALSE)

## clump the data
exposure_data <- clump_data(FA_format, clump_kb = 10000, clump_r2 = 0.001)

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
sumarydat <- subset(dat, mr_keep == 'TRUE')
library(MRPRESSO)
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE, data = dat, NbDistribution = 3000,  
          SignifThreshold = 0.05)