#Input and output
##################

output_prefix <- "3.WGCNA/allsamps"

meta <- "3.WGCNA/allsamps_meta.tsv"

counts <- "3.WGCNA/allsamps_counts.tsv"



####################################
####################################
# Loading working packages
############################
library(WGCNA)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(edgeR)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 20)

# Formatting count data and normalizing
##################

meta <- read.delim(file = meta)
counts <- read.delim(file = counts)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~Condition) #Change for meta-variables of interest

#Normalize
dds = DESeq(dds)
rcount = assay(rlog(dds,blind=F))
#Transpose table, because WGCNA needs you to
tcount = t(rcount)
#Check for outliers, genes/samples with too many zeroes, etc
check = goodSamplesGenes(tcount)
if (check$allOK==T) {
  print('All genes and samples looking good!')
} else {
  if (sum(!gsg$goodGenes)>0)
    print(paste("Bad genes:", paste(names(tcount)[!check$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Bad samples:", paste(rownames(tcount)[!check$goodSamples], collapse = ", ")));
}


# Picking soft thresholding power
##################

#Scale-free topology approach; copy paste from the tutorial
#####
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(tcount, powerVector = powers)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#addfiglab("A")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");
abline(h=500,col="red")
#addfiglab("B")


####

# Manually picking the thresholding power (important)
##################

pwr = 12


# One-step network construction and module detection
##################

temp_cor <- cor
cor <- WGCNA::cor
net = blockwiseModules(tcount, power = pwr,
                       TOMType = "signed", networkType = 'signed',
                       maxBlockSize = 25000,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       deepSplit = 2,
                       numericLabels = TRUE, #pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = paste(output_prefix,'signed',sep='_'))
cor <- temp_cor

mergedColors = labels2colors(net$colors)
# Plot ugly dendrogram
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#Save modules as separate variable
modules = data.frame(gene_id=names(net$colors),
                     m_number = net$colors,
                     colors=labels2colors(net$colors))
write.table(modules,file=paste(output_prefix,'signed_modules.tsv',sep='_'),sep='\t',quote=F,row.names =F,col.names = F)

# Get and plot eigengenes
##################
MEs = moduleEigengenes(tcount, mergedColors)$eigengenes
MEs <- orderMEs(MEs)
module_order = names(MEs) %>% gsub("ME","", .)
MEs$treatment = row.names(MEs)
mME = MEs %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order),
    treatment = factor(treatment,levels = sample_order)
  )

ggplot(mME, aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(low = "blue",high = "red", mid = "white",
                       midpoint = 0,limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

table(as.factor(modules$colors))

# Export network as edge/vertices file if needed
##################
# Takes a VERY long time, output is VERY big
edge_list = TOMsimilarityFromExpr(tcount,power=pwr)
row.names(edge_list) = colnames(tcount)
colnames(edge_list) = colnames(tcount)
edge_list = data.frame(edge_list) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(module1 = MEs[gene1,]$colors,module2 = MEs[gene2,]$colors)

write.table(edge_list,file=paste(output_prefix,'edge_list.tsv',sep='_'),sep='\t',quote=F)



# Module preservation
################################

# Load query data
reference_counts <- read.delim(file="/work/verjo103/group/prostata/processamento/cutnrun/8.greenlist_norm/all_samples_counts.tsv",header=T)

AR_counts = read.delim(file='3.WGCNA/2.CnR_data/AR_promoters_manualcounts.tsv',header=T,row.names = 1)
ezh2_counts = read.delim(file='3.WGCNA/2.CnR_data/ezh2_promoters_manualcounts.tsv',header=T,row.names = 1)
k27ac_counts = read.delim(file='3.WGCNA/2.CnR_data/h3k27ac_promoters_manualcounts.tsv',header=T,row.names = 1)
k27me3_counts = read.delim(file='3.WGCNA/2.CnR_data/h3k27me3_promoters_manualcounts.tsv',header=T,row.names = 1)
k4me1_counts = read.delim(file='3.WGCNA/2.CnR_data/h3k4me1_promoters_manualcounts.tsv',header=T,row.names = 1)
k4me3_counts = read.delim(file='3.WGCNA/2.CnR_data/h3k4me3_promoters_manualcounts.tsv',header=T,row.names = 1)
cneg_counts = read.delim(file='3.WGCNA/2.CnR_data/cneg_promoters_manualcounts.tsv',header=T,row.names = 1)
input_counts = read.delim(file='3.WGCNA/2.CnR_data/input_promoters_manualcounts.tsv',header=T,row.names = 1)
mainrna_counts = read.delim(file='3.WGCNA/mainsamps_counts.tsv')
main_meta = read.delim(file='3.WGCNA/mainsamps_meta.tsv')
cnr_meta = read.delim(file='3.WGCNA/2.CnR_data/cnr_meta.tsv')

# Normalize and format query counts
ddds <- DESeqDataSetFromMatrix(countData = AR_counts[rowSums(AR_counts) >= 1,],
                               colData = cnr_meta,
                               design = ~Condition)
sizeFactors(ddds) = estimateSizeFactorsForMatrix(reference_counts[,1:8])
ddds = estimateDispersions(ddds)
AR_tcount = t(counts(ddds,normalized=T))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = ezh2_counts[rowSums(ezh2_counts) >= 1,],
                               colData = cnr_meta,
                               design = ~Condition)
sizeFactors(ddds) = estimateSizeFactorsForMatrix(reference_counts[,17:24])
ddds = estimateDispersions(ddds)
ezh2_tcount = t(counts(ddds,normalized=T))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = k27ac_counts[rowSums(k27ac_counts) >= 1,],
                               colData = cnr_meta,
                               design = ~Condition)
sizeFactors(ddds) = estimateSizeFactorsForMatrix(reference_counts[,25:32])
ddds = estimateDispersions(ddds)
k27ac_tcount = t(counts(ddds,normalized=T))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = k27me3_counts[rowSums(k27me3_counts) >= 1,],
                               colData = cnr_meta,
                               design = ~Condition)
sizeFactors(ddds) = estimateSizeFactorsForMatrix(reference_counts[,33:40])
ddds = estimateDispersions(ddds)
k27me3_tcount = t(counts(ddds,normalized=T))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = k4me1_counts[rowSums(k4me1_counts) >= 1,],
                               colData = cnr_meta[c(1:4,6:8),],
                               design = ~Condition)
sizeFactors(ddds) = estimateSizeFactorsForMatrix(reference_counts[,41:47])
ddds = estimateDispersions(ddds)
k4me1_tcount = t(counts(ddds,normalized=T))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = k4me3_counts[rowSums(k4me3_counts) >= 1,],
                               colData = cnr_meta,
                               design = ~Condition)
sizeFactors(ddds) = estimateSizeFactorsForMatrix(reference_counts[,48:55])
ddds = estimateDispersions(ddds)
k4me3_tcount = t(counts(ddds,normalized=T))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = cneg_counts[rowSums(cneg_counts) >= 1,],
                               colData = cnr_meta,
                               design = ~Condition)
sizeFactors(ddds) = estimateSizeFactorsForMatrix(reference_counts[,9:16])
ddds = estimateDispersions(ddds)
cneg_tcount = t(counts(ddds,normalized=T))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = input_counts[rowSums(input_counts >= 10) >= 1,],
                               colData = cnr_meta,
                               design = ~Condition)
ddds = DESeq(ddds)
input_tcount = t(assay(rlog(ddds,blind=F)))
rm(ddds)
#
ddds <- DESeqDataSetFromMatrix(countData = mainrna_counts,
                               colData = main_meta,
                               design = ~Condition)
ddds = DESeq(ddds)
mainrna_tcount = t(assay(rlog(ddds,blind=F)))
rm(ddds)
#####
data.frame(dim(AR_tcount),dim(cneg_tcount),dim(ezh2_tcount),dim(k27ac_tcount),dim(k27me3_tcount),dim(k4me1_tcount),dim(k4me3_tcount))
# Build multiExp and run; subtract IgG counts from other ab's counts
multiExpr = multiData(Ref = tcount, AR = AR_tcount-cneg_tcount, ezh2 = ezh2_tcount-cneg_tcount, k27ac = k27ac_tcount-cneg_tcount, k27me3 = k27me3_tcount-cneg_tcount, k4me1 = k4me1_tcount-cneg_tcount[rownames(k4me1_tcount),], k4me3 = k4me3_tcount-cneg_tcount, cneg = cneg_tcount)
colorList = list(Ref = set_names(modules$colors,modules$gene_id))

start=Sys.time()
start
cor <- WGCNA::cor
set.seed(1)
cnr_net = modulePreservation(multiExpr,colorList,
                             dataIsExpr=T,networkType = 'signed',corFnc='bicor',
                             referenceNetworks=1,testNetworks=c(2:8),
                             nPermutations = 10000,calculateQvalue=F,
                             randomSeed=1,quickCor=0,checkData=F,
                             maxModuleSize=4000,maxGoldModuleSize=1000,
                             verbose=5,parallelCalculation=T)
cor <- temp_cor
Sys.time()-start
save.image(file='3.WGCNA/preservation_minus.norm_bicor.RData')



###########################################################
##################### Significance analysis


library(ggpubr)


df=NULL
df=as.data.frame(df)
for (j in 1:7){
  for (i in 1:33){
    vec=NULL
    for (k in 1:10000){
      vec=append(vec,cnr_net[["permutationDetails"]][["permutedStatistics"]][[1]][[5]][["regStats"]][i,c(8,9,11,12,15,16,18)[j],k]-cnr_net[["permutationDetails"]][["permutedStatistics"]][[1]][[8]][["regStats"]][i,c(8,9,11,12,15,16,18)[j],])  
      df[i,j]=as.data.frame(table(vec>=cnr_net$preservation$observed$ref.Ref$inColumnsAlsoPresentIn.k27me3[i,c(5,6,7,8,10,11,13)[j]]-cnr_net$preservation$observed$ref.Ref$inColumnsAlsoPresentIn.cneg[i,c(5,6,7,8,10,11,13)[j]])[1]/10000)
    }}}
df=1-df
rownames(df)=rownames(cnr_net[["permutationDetails"]][["permutedStatistics"]][[1]][[5]][["regStats"]][,c(8,9,11,12,15,16,18),1])
colnames(df)=colnames(cnr_net[["permutationDetails"]][["permutedStatistics"]][[1]][[5]][["regStats"]][,c(8,9,11,12,15,16,18),1])
df$veredict=rowSums(df<=0.05)









