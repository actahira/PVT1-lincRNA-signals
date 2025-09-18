#Input and output
##################

output_prefix <- "glist"

counts <- read.delim(file="AR/8.DESeq/manualcounts_AR_consensus.tsv")
reference_counts <- read.delim(file="8.greenlist_norm/all_samples_counts.tsv",header=T)
metadata <- "diffex_meta.tsv"

test_variables <- c("Condition")
control_values <- c("Control")
# Define column(s) in the metadata table to be used as
# test variable(s) of interest, as well as defining the
# values for each to be used as control (in the same order,
# one value per variable)

covariates <- c("Replicate")
# Define column(s) in the metadata table to be used as
# confounding variable(s), ie. factors that are not of 
# interest but might be affecting the observed results 



####################################
####################################
#Loading working packages
############################
library(sva)
library(DESeq2)
library(dplyr)
library(tidyr)

#Formatting metadata, preparing count table and design
##################
metadata <- read.delim(file = metadata)
meta <- metadata[metadata$Factor == "AR",]
meta <- data.frame(meta,row.names = "SampleID")
meta$Replicate <- as.factor(meta$Replicate)
filt <- subset(reference_counts[,1:8])
counts_with_glist <- rbind(filt,counts)

dds <- DESeqDataSetFromMatrix(countData = counts_with_glist,
                              colData = meta,
                              design = ~ Replicate + Condition)

dds$Condition <- relevel(dds$Condition, ref = control_values[1])

#Artifact-normalized test
##################
dds <- estimateSizeFactors(dds, controlGenes=1:nrow(filt))
sf = sizeFactors(dds)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Replicate + Condition)
dds$Condition <- relevel(dds$Condition, ref = control_values[1])
sizeFactors(dds) <- sf


#Filtering low counts
##################
keep <- rowSums(counts(dds) >= 10) >= 1
dds <- dds[keep,]



#Generating DE lists from tests
##################

DEtest = DESeq(dds)
padj = 'fdr'
res_co.v.Horm = results(DEtest, contrast = c("Condition","hormone_only","Control"),cooksCutoff=F,alpha=0.05,pAdjustMethod=padj)
res_co.v.KD = results(DEtest, contrast = c("Condition","knockdown_only","Control"),cooksCutoff=F,alpha=0.05,pAdjustMethod=padj)
res_Horm.v.comb = results(DEtest, contrast = c("Condition","Combination","hormone_only"),cooksCutoff=F,alpha=0.05,pAdjustMethod=padj)

#Saving/reporting results
##################
save = T
savefolder = "AR/8.DESeq/AR"
if (save==F) {

  print("Co vs Horm")
  summary(res_co.v.Horm)
  print(paste("Control v. Hormone,",nrow(res_co.v.Horm[!is.na(res_co.v.Horm$padj),][res_co.v.Horm$padj[!is.na(res_co.v.Horm$padj)]<=0.05,]),"DE peaks found",sep=' '))
  print("========================")
    
  print("Co vs KD")
  summary(res_co.v.KD)
  print(paste("Control v. KD,",nrow(res_co.v.KD[!is.na(res_co.v.KD$padj),][res_co.v.KD$padj[!is.na(res_co.v.KD$padj)]<=0.05,]),"DE peaks found",sep=' '))
  print("========================")
    
  print("Horm vs Comb")
  summary(res_Horm.v.comb)
  print(paste("Hormone v. Combination,",nrow(res_Horm.v.comb[!is.na(res_Horm.v.comb$padj),][res_Horm.v.comb$padj[!is.na(res_Horm.v.comb$padj)]<=0.05,]),"DE peaks found",sep=' '))
  print("========================")
    
} else {
  # Co vs horm
  write.table(data.frame(res_co.v.Horm),
              row.names=T,col.names=F,quote=F, sep="\t",na="1",
              file=paste(savefolder,output_prefix,"co.vs.horm_allpeaks.tsv",sep='_'))
  write.table(data.frame(res_co.v.Horm[!is.na(res_co.v.Horm$padj),][res_co.v.Horm$padj[!is.na(res_co.v.Horm$padj)]<=0.05,]),
              row.names=T,col.names=F,quote=F, sep="\t",na="1",
              file=paste(savefolder,output_prefix,"co.vs.horm_DEpeaks.tsv",sep='_'))
  ###
  # Co vs KD
  write.table(data.frame(res_co.v.KD),
              row.names=T,col.names=F,quote=F, sep="\t",na="1",
              file=paste(savefolder,output_prefix,"co.vs.kd_allpeaks.tsv",sep='_'))
  write.table(data.frame(res_co.v.KD[!is.na(res_co.v.KD$padj),][res_co.v.KD$padj[!is.na(res_co.v.KD$padj)]<=0.05,]),
              row.names=T,col.names=F,quote=F, sep="\t",na="1",
              file=paste(savefolder,output_prefix,"co.vs.kd_DEpeaks.tsv",sep='_'))
  ###
  # Horm vs Comb
  write.table(data.frame(res_Horm.v.comb),
              row.names=T,col.names=F,quote=F, sep="\t",na="1",
              file=paste(savefolder,output_prefix,"horm.vs.comb_allpeaks.tsv",sep='_'))
  write.table(data.frame(res_Horm.v.comb[!is.na(res_Horm.v.comb$padj),][res_Horm.v.comb$padj[!is.na(res_Horm.v.comb$padj)]<=0.05,]),
              row.names=T,col.names=F,quote=F, sep="\t",na="1",
              file=paste(savefolder,output_prefix,"horm.vs.comb_DEpeaks.tsv",sep='_'))
}







