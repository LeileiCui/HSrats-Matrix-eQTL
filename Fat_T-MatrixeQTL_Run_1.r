
rm(list=ls());options(stringsAsFactors=FALSE)

library(data.table); library(parallel); library(bigmemory)
read <- function(...) as.data.frame(fread(header=F,colClasses="double",...))

num_cov = 41

# datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/1_Fat_Genes/1_PheGen/0_Data"
datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/1_Fat_Trpts/1_PheGen/0_Data"
# datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/2_Liver_Genes/1_PheGen/0_Data"
# datdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/2_Liver_Trpts/1_PheGen/0_Data"

# outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/1_Fat_Genes"; setwd(outdir)
outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/1_Fat_Trpts"; setwd(outdir)
# outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/2_Liver_Genes"; setwd(outdir)
# outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu_MatrixeQTL/2_Liver_Trpts"; setwd(outdir)

### [1] Generate the SNP file ###
system(paste0("plink --bfile ", datdir, "/Data_clean --recode12 --out tmp_clean"))
Chip_ped = read.big.matrix("tmp_clean.ped", sep = " ", type = "integer")
Chip_fam = read(paste0(datdir, "/Data_clean.fam")); Chip_map = read("tmp_clean.map")

allele_to_genotype_012 <- function(i){
	x_sum = rowSums(ped_allele[,(2*i-1):(2*i)])
	x_result = x_sum - 2; x_result[x_result < 0] = 0 # allele code: 0/1/2 -> genotype code: 0/1/3/4 -> genotype code minus 2: -2/-1/1/2 #
	return(x_result)
}

ped_allele = Chip_ped[,c(-1:-6)]
ped_geno_lma = do.call('cbind', mclapply(1:((dim(ped_allele)[2])/2), allele_to_genotype_012, mc.cores=1))
rownames(ped_geno_lma) = paste0(Chip_fam[,1], "_", Chip_fam[,2])
all_SNP = cbind(c("id",paste0(Chip_map[,1], "_", Chip_map[,4])), rbind(rownames(ped_geno_lma), t(ped_geno_lma)))

### [2] Generate the Phe file ###
Chip_phe = read.table(paste0(datdir, "/../2_PheClean/Phe_Clean.txt"),header=T)
all_Phe = t(Chip_phe)[(num_cov+2):ncol(Chip_phe),]
all_Phe = rbind(Chip_phe[,1], all_Phe)
all_Phe = cbind(c("id", rownames(all_Phe)[-1]), all_Phe)
all_Phe[-1,1] = colnames(Chip_phe)[(num_cov+2):ncol(Chip_phe)]

### [3] Generate the Cov file ###
all_Cov = t(Chip_phe)[2:(num_cov+1),]
all_Cov = rbind(Chip_phe[,1], all_Cov)
all_Cov = cbind(c("id", rownames(all_Cov)[-1]), all_Cov)
all_Cov[-1,1] = colnames(Chip_phe)[2:(num_cov+1)]

write.table(all_Phe, file="output_Phe.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(all_Cov, file="output_Cov.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(all_SNP, file="output_SNP.txt", row.names=F, col.names=F, quote=F, sep="\t")

system("rm tmp*")

### [4] Prepare the kinship matrix ###

Kinship_type = "EMMAX"
if(Kinship_type == "EMMAX"){
	setwd(datdir)
	if(!file.exists("tmp_emmax.tped")){system("plink --noweb --silent --bfile Data_clean --recode12 --output-missing-genotype 0 --transpose --out tmp_emmax")}
	system("emmax-kin -v -h -s -d 10 tmp_emmax") # *.hIBS, IBS matrix  #
	system("mv tmp_emmax.hIBS.kinf Data_clean_kinship_EMMAX.hIBS.txt")
	#system("emmax-kin -v -h -d 10 tmp_emmax") # *.BN, BN matrix #
	#system("mv tmp_emmax.hBN.kinf Data_clean_kinship_EMMAX.hBN.txt")
	system("rm tmp_*")
}

