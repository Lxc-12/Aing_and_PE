#        +---------------------------------------+
#        |     Epigenetic Age Acceleration       |
#        +---------------------------------------+
## import the EEAA(Extrinsic Epigenetic Age Acceleration) GWAS data into R.
library(vroom)
EEAA_gwas <- vroom('F:\\Aging_gwas\\EEA\\EEAA_metaPlusGS.csv',
                   col_names = TRUE)

colnames(EEAA_gwas)
EEAA_gwas$phenotype <- "EEAA"

library(TwoSampleMR)
library(MRInstruments)

## extract SNPs with p£¼5e-08.
EEAA_p <- subset(EEAA_gwas, metaPval<5e-08)

## format the data.
EEAA_format <- format_data(EEAA_p, type = "exposure",
                             snps = NULL,
                             header = TRUE,
                             phenotype_col = "phenotype",
                             snp_col = "SNP",
                             beta_col = "meta_b",
                             se_col = "meta_se",
                             eaf_col = "metaFreq",
                             effect_allele_col = "metaA1",
                             other_allele_col = "metaA2",
                             pval_col = "metaPval",
                             units_col = "units",
                             ncase_col = "ncase",
                             ncontrol_col = "ncontrol",
                             samplesize_col = "N",
                             gene_col = "gene",
                             id_col = "id",
                             min_pval = 1e-500,
                             z_col = "z",
                             info_col = "info",
                             chr_col = "Chr",
                             pos_col = "bp",
                             log_pval = FALSE)
## clump the data
exposure_data1 <- clump_data(EEAA_format, clump_kb = 10000,clump_r2 = 0.001)

#---------------------------------------------------------------------

## import the IEAA(Intrinsic Epigenetic Age Acceleration) GWAS data into R.
library(vroom)
IEAA_gwas <- vroom('F:\\Aging_gwas\\EEA\\IEAA_metaPlusGS.csv',
                   col_names = TRUE)

colnames(IEAA_gwas)
IEAA_gwas$phenotype <- "IEAA"

library(TwoSampleMR)
library(MRInstruments)

## extract SNPs with p£¼5e-08.
IEAA_p <- subset(IEAA_gwas, metaPval<5e-08)

## format the data.
IEAA_format <- format_data(IEAA_p, type = "exposure",
                             snps = NULL,
                             header = TRUE,
                             phenotype_col = "phenotype",
                             snp_col = "SNP",
                             beta_col = "meta_b",
                             se_col = "meta_se",
                             eaf_col = "metaFreq",
                             effect_allele_col = "metaA1",
                             other_allele_col = "metaA2",
                             pval_col = "metaPval",
                             units_col = "units",
                             ncase_col = "ncase",
                             ncontrol_col = "ncontrol",
                             samplesize_col = "N",
                             gene_col = "gene",
                             id_col = "id",
                             min_pval = 1e-500,
                             z_col = "z",
                             info_col = "info",
                             chr_col = "Chr",
                             pos_col = "bp",
                             log_pval = FALSE)
## clump the data
exposure_data2 <- clump_data(IEAA_format, clump_kb = 10000,clump_r2 = 0.001)

## set the EEAA exposure data and IEAA exposure data.
exposure_data <- rbind(exposure_data1, exposure_data2)

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
## sctter plot
mr_scatter_plot(mr_results = mr(dat),dat)
## heterogeneity test
mr_heterogeneity(dat)
run_mr_presso(dat,NbDistribution = 1000)
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
          DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
          SignifThreshold = 0.05)