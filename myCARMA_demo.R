path_files <- "/home/ines/CARMA_tryouts"

pkgs = c("data.table", "magrittr", "dplyr", "devtools","R.utlis")
pkgs.na = pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(pkgs.na) > 0) {
  install.packages(pkgs.na)
}
if (!"CARMA" %in% installed.packages()[, "Package"]) {
  devtools::install_github("ZikunY/CARMA")
}

library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(myCARMA)

inputsumstats <- paste0(path_files, "Sample_data/sumstats_chr1_200937832_201937832.txt.gz")
sumstat<- fread(file = inputsumstats,
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)
inputLD <- paste0(path_files, "Sample_data/sumstats_chr1_200937832_201937832_ld.txt.gz")
ld <- fread(file = inputLD,
           sep = "\t", header = F, check.names = F, data.table = F,
           stringsAsFactors = F)

print(head(sumstat))
print(head(ld))

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

#Run CARMA
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                     outlier.switch=F)
#Get PIP and credible sets
sumstat.result <- sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}

#Save results
outputsumstats <- paste0(path_files, "Sample_data/prueba1_sumstats_chr1_200937832_201937832_carma.txt.gz")
fwrite(x = sumstat.result,
       file = outputsumstats,
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")

#Check credible sets
CARMA.results[[1]]$`Credible set`[[2]]

#Check credible models
CARMA.results[[1]]$`Credible model`[[3]]

##Run CARMA with annotations
inputannot <- paste0(path_files, "Sample_data/sumstats_chr1_200937832_201937832_annotations.txt.gz")
annot=fread(file = inputannot,
            sep="\t", header = T, check.names = F, data.table = F,
            stringsAsFactors = F)
#z list and ld list stay the same

#z.list<-list()
#ld.list<-list()
#lambda.list<-list()
#z.list[[1]]<-sumstat$Z
#ld.list[[1]]<-as.matrix(ld)
#lambda.list[[1]]<-1

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
outputsumstats_withannot <- paste0(path_files, "Sample_data/prueba1_sumstats_chr1_200937832_201937832_carma_annot.txt.gz")
fwrite(x = sumstat.result,
       file = outputsumstats_withannot,
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")
#Get credible sets
CARMA.results[[1]]$`Credible set`[[2]]
#Get credible models
CARMA.results[[1]]$`Credible model`[[3]]

