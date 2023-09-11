#               +-------------------------------------+
#               |    analysis Pleiotrophin and PE     | 
#               +-------------------------------------+
#read the Pleiotrophin GWAS into R
library(vroom)
Pleiotrophin_gwas <- vroom("F:\\Aging_gwas\\Pleiotrophin\\Pleiotrophin.csv",
                      col_names = TRUE)

colnames(Pleiotrophin_gwas)
Pleiotrophin_gwas$phenotype <- "Pleiotrophin"

library(TwoSampleMR)
library(MRInstruments)

Pleiotrophin_gwas$neg_p <- -Pleiotrophin_gwas$LP
Pleiotrophin_gwas$p_value <- 10^Pleiotrophin_gwas$neg_p

##提取p＜5e-08的snp
Pleiotrophin_p <- subset(Pleiotrophin_gwas, p_value < 5e-08)

##数据转换为后缀名为.exposure的标准化数据
Pleiotrophin_format <- format_data(Pleiotrophin_p, type = "exposure",
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
                              info_col = "info",
                              chr_col = "CHROM",
                              pos_col = "POS",
                              log_pval = FALSE)

##去除连锁不平衡
exposure_data <- clump_data(Pleiotrophin_format, clump_kb = 10000,clump_r2 = 0.001)


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

##harmonise数据
dat <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
##进行MR分析
mr(dat)
generate_odds_ratios(mr_res = mr(dat))
##森林图
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
##散点图
mr_scatter_plot(mr_results = mr(dat),dat)
##异质性
mr_heterogeneity(dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
##多效性
mr_pleiotropy_test(dat)
##敏感性分析
mr_leaveoneout(dat)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))


## MR-PRESSO analysis
# Not applicable for MR PRESSO