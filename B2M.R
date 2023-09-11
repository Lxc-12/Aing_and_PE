#               +---------------------------+
#               |     analysis B2M and PE   | 
#               +---------------------------+
## import the beta-2-microglobulin GWAS data into R.
library(vroom)
B2M_gwas <- vroom("F:\\Aging_gwas\\B2M\\B2M.csv",
                  col_names = TRUE)

colnames(B2M_gwas)
B2M_gwas$phenotype <- "beta-2-microglobulin"

library(TwoSampleMR)
library(MRInstruments)

## Convert -log10 p-value to P-value
B2M_gwas$neg_p <- -B2M_gwas$LP
B2M_gwas$p_value <- 10^B2M_gwas$neg_p

## extract SNPs with p£¼5e-08.
B2M_p <- subset(B2M_gwas, p_value < 5e-08)

## format the data.
B2M_format <- format_data(B2M_p, type = "exposure",
                         snps = NULL,
                         header = TRUE,
                         phenotype_col = "phenotype",
                         snp_col = "ID",
                         beta_col = "ES",
                         se_col = "SE",
                         eaf_col = "",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF",
                         pval_col = "p_value",
                         units_col = "units",
                         ncase_col = "ncase",
                         ncontrol_col = "ncontrol",
                         samplesize_col = "SS",
                         gene_col = "gene",
                         id_col = "id",
                         min_pval = 0,
                         z_col = "z",
                         info_col = "info",
                         chr_col = "CHROM",
                         pos_col = "POS",
                         log_pval = FALSE)

## clump the data
exposure_data <- clump_data(B2M_format, clump_kb = 10000, clump_r2 = 0.001)

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


## Only one SNP was conducted MR analysis, Not applicable for MRPRESSO analysis.