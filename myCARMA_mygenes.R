pkgs = c("data.table", "magrittr", "dplyr", "devtools","R.utlis")
pkgs.na = pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(pkgs.na) > 0) {
  install.packages(pkgs.na)
}
if (!"CARMA" %in% installed.packages()[, "Package"]) {
  devtools::install_github("igarcia17/myCARMA")
}

library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(myCARMA)

trial <- 'ANK3_10000'
path_files <- paste0("/home/ines/CARMA_tryouts/myCARMA_genes/", trial, '/')

#My genes: ANK3
inputsumstats<- paste0(path_files, trial, '.sumstats')
sumstat0 <- fread(file = inputsumstats, sep="\t", header=T, check.names = F,
                  data.table = F, stringsAsFactors = F)

inputLD <- paste0(path_files, trial, '_LD.ld')
ld0 <- fread(file = inputLD,
            sep = "\t", header = F, check.names = F, data.table = F,
            stringsAsFactors = F)

#Load list of snps in LD panel
input_var_inLD <- paste0(path_files, trial, '_SNPs4LD.txt')
var_inLD <- fread(file= input_var_inLD, sep="\t", header=F)

#Filter sumstats file so it only has variants that are present in the LD panel
sumstat <- sumstat0[sumstat0$SNP %in% var_inLD$V1,]

#Make a submatrix of LD with only the SNPs on sumstat file
#Make subset
vector <- var_inLD$V1 %in% sumstat$SNP
ld <- ld0[vector, vector]

#Calculate Z score of summary statistics
#https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html
sumstat$Z <- (log(sumstat$OR))/sumstat$SE

#Preparation of lists
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1 #TODO: mirar como usar GATC-COJO para captar numero de señales independientes

#Run CARMA
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                     outlier.switch=T) #outlier switch is TRUE if external LD panel
#Get PIP and credible sets
sumstat.result <- sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}

#Save results
outputsumstats <- paste0(path_files, "res_noannot.txt")
fwrite(x = sumstat.result,
       file = outputsumstats,
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T)

#Check credible sets
CARMA.results[[1]]$`Credible set`[[2]]

#Check credible models
CARMA.results[[1]]$`Credible model`[[3]]


###TODO
##Run CARMA with annotations
inputannot <- paste0(path_files, "Sample_data/sumstats_chr1_200937832_201937832_annotations.txt.gz")
annot=fread(file = inputannot,
            sep="\t", header = T, check.names = F, data.table = F,
            stringsAsFactors = F)
#z list and ld list stay the same

annot.list<-list()
annot.list[[1]]<-annot
#Run CARMA with annotations
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list=annot.list,
                     outlier.switch=F)
#Get PIP and credible sets
sumstat.result <- sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}
#Save
outputsumstats_withannot <- paste0(path_files, "reswithannot.txt")
fwrite(x = sumstat.result,
       file = outputsumstats_withannot,
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")
#Get credible sets
CARMA.results[[1]]$`Credible set`[[2]]
#Get credible models
CARMA.results[[1]]$`Credible model`[[3]]







