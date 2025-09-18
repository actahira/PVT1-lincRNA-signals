#Input and output
##################

output_prefix <- "wgcna_tpm"

meta_in <- "Documents/Mestrado/Resultados/TCGA/combined_meta.tsv"

counts_in <- "combined_tpm.tsv"



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
library(mgcv)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 40)

# Formatting count data and normalizing
##################


meta <- read.delim(file = meta_in,header=T,row.names = NULL,check.names=F,stringsAsFactors = T,na.strings = "N/A")
rownames(meta)=meta$file_name
counts <- read.delim(file = counts_in,header=T,row.names=1,check.names=F)
CVs = apply(counts,1,FUN=function(x) sd(x)/mean(x))
tcount=t(counts[CVs>0.5,])


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
powers = c(c(1:10), seq(from = 12, to=30, by=2))
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
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
####

# Manually picking the thresholding power (important)
##################

pwr = 4


# One-step network construction and module detection
##################

temp_cor <- cor
cor <- WGCNA::cor
net = blockwiseModules(tcount, power = pwr,
                       TOMType = "signed", networkType = 'signed',
                       maxBlockSize = 250000,quickCor=1,
                       #reassignThreshold = 0, 
                       #mergeCutHeight = 0.1,
                       #deepSplit = 2,
                       numericLabels = TRUE, #pamRespectsDendro = FALSE,
                       saveTOMs = F,nThreads=40)
# saveTOMFileBase = paste(output_prefix,'signed',sep='_'))
cor <- temp_cor

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#Save modules as separate variable
modules = data.frame(gene_id=names(net$colors),
                     m_number = net$colors,
                     colors=labels2colors(net$colors))
#write.table(modules,file=paste(output_prefix,'signed_modules.tsv',sep='_'),sep='\t',quote=F,row.names =F,col.names = F)

# Get and plot eigengenes
##################
MEs = moduleEigengenes(tcount, mergedColors)$eigengenes
MEs <- orderMEs(MEs)
module_order = names(MEs) %>% gsub("ME","", .)
smplcor=cor(t(MEs[,1:ncol(MEs)-1]))
hc = hclust(as.dist((1-smplcor)/2))
sample_order=rownames(MEs[hc$order,])
mod_zscore <- as.data.frame(t(apply(MEs, 1, function(row) (row - mean(row)) / sd(row))))
mod_zscore$treatment=row.names(mod_zscore)
MEs$treatment = row.names(MEs)
mME = mod_zscore %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order),
    treatment = factor(treatment, levels = sample_order)
  )

ggplot(mME, aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(low = "blue",high = "red", mid = "white",
                       #midpoint = 0,limit = c(-1,1)
  ) +
  scale_x_discrete(labels=NULL)+
  #theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-sample Relationships", y = "Modules", fill="Z-score",x='Samples')

table(as.factor(modules$colors))

clustered_pvt1_exp=t(counts)[sample_order,'ENSG00000249859.12']
pvt1_expr_df=data.frame(Exp=clustered_pvt1_exp,Order=1:length(clustered_pvt1_exp))
ggplot(pvt1_expr_df, aes(x=Order,y=Exp))+
  geom_line()+
  theme_bw()+
  scale_y_continuous(trans='log2',name = 'TPM')+
  geom_smooth(method='gam')


# Discrete variable testing on norm counts
##################

variable_of_interest <- "gleason_score"

design <- model.matrix(~ 0 + as.factor(meta[[variable_of_interest]]))
colnames(design) <- levels(as.factor(meta[[variable_of_interest]]))
colnames(design) = make.names(colnames(design))

# Fit the linear model
fit <- lmFit(heatmap_data[,!is.na(meta[[variable_of_interest]])], design)

# Specify contrasts (e.g., pairwise comparisons)
contrast_matrix <- makeContrasts(contrasts = paste(colnames(design), collapse = "-"), levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, number = Inf, sort.by = "P")


# Continuous variable testing on norm counts
##################

variable_of_interest <- "psa_most_recent_results"
gene_results <- list()

# Loop through each gene and fit a GAM
for (gene in colnames(tcount[!is.na(meta[[variable_of_interest]]),])) {
  gam_fit <- gam(tcount[,gene] ~ s(meta$psa_most_recent_results,k=length(unique(meta[[variable_of_interest]]))-1), method = "REML")
  summary_fit <- summary(gam_fit)
  p_value <- summary_fit$s.pv[1]  # Extract p-value for the smooth term
  
  # Save the results
  gene_results[[gene]] <- data.frame(
    Gene = gene,
    EDF = summary_fit$s.table[1, "edf"],  # Effective degrees of freedom
    F_stat = summary_fit$s.table[1, "F"],  # F-statistic
    P_value = p_value
  )
}

# Combine results into a data frame
results <- do.call(rbind, gene_results)
results$Adj.P = p.adjust(results$P_value,method = 'fdr')


variable_of_interest <- "brief_overall_T"
gene_results <- list()

# Loop through each gene and fit a GAM
for (gene in colnames(tcount[!is.na(meta[[variable_of_interest]]),])) {
  glm_fit <- glm(tcount[,gene] ~ as.factor(meta[[variable_of_interest]]))
  summary_fit <- summary(aov(glm_fit))

  # Save the results
  gene_results[[gene]] <- data.frame(
    Gene = gene,
    F_stat = summary_fit[[1]][['F value']][1],  # F-statistic
    P_value = summary_fit[[1]][['Pr(>F)']][1]
  )
}

# Combine results into a data frame
results <- do.call(rbind, gene_results)
results$Adj.P = p.adjust(results$P_value,method = 'fdr')

