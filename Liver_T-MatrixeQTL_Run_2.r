
rm(list=ls());options(stringsAsFactors=FALSE)
library(MatrixEQTL)

# datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/2_Liver_Genes/1_PheGen/0_Data"
datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/2_Liver_Trpts/1_PheGen/0_Data"
# datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/2_Liver_Genes/1_PheGen/0_Data"
# datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/2_Liver_Trpts/1_PheGen/0_Data"

# outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/2_Liver_Genes"; setwd(outdir)
outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/2_Liver_Trpts"; setwd(outdir)
# outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/2_Liver_Genes"; setwd(outdir)
# outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/2_Liver_Trpts"; setwd(outdir)

SNP_file_name = paste0(outdir, "/output_SNP.txt");
expression_file_name = paste0(outdir, "/output_Phe.txt");
covariates_file_name = paste0(outdir, "/output_Cov.txt");

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 10000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) { cvrt$LoadFile(covariates_file_name); }

errorCovariance = as.matrix(read.table(paste0(datdir, "/Data_clean_kinship_EMMAX.hIBS.txt"), header=F)) # positive definite #
# errorCovariance = as.matrix(read.table(paste0(datdir, "/Data_clean_kinship_GCTA.Additive.txt"), header=F)) # not positive definite #
# errorCovariance = as.matrix(read.table(paste0(datdir, "/Data_clean_kinship_GCTA.Dominant.txt"), header=F)) # positive definite #
# errorCovariance = as.matrix(read.table(paste0(datdir, "/Data_clean_kinship_GEMMA.cXX.txt"), header=F)) # not positive definite #

pvOutputThreshold = 0.05/nrow(snps)
# pvOutputThreshold = 0.05/(nrow(snps)*nrow(gene))

me.anova = Matrix_eQTL_engine(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name = "results_ANOVA_eQTLs.txt",
pvOutputThreshold = pvOutputThreshold,
useModel = modelANOVA, 
errorCovariance = errorCovariance, 
verbose = TRUE,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

