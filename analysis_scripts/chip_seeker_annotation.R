

input <- "peakseeker_consensus_sheet.tsv"



library(ChIPseeker)
library(clusterProfiler) 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mirbase.db)
library(Homo.sapiens)
options(mc.cores=32)


hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene


meta <- read.delim(file = input)
chip_factors <- unique(meta$Factor)

for (i in chip_factors){
  assign(paste(i,"peaks",sep="_"),
         readPeakFile(paste(meta[meta$Factor == i,]$DiffBind.Folder,
                            meta[meta$Factor == i,]$Peakset,sep="")))
}

for (i in chip_factors){
  assign(paste(i,"annotation",sep="_"),
         annotatePeak(eval(parse(text=paste(i,"peaks",sep="_"))), 
                      annoDb = "org.Hs.eg.db", TxDb = hg38,overlap='all',
                      tssRegion = c(-5000, 5000), flankDistance = 5000))
}




write.table(as.data.frame(AR_annotation@anno)[, c(1,2,3,6,7,15,8,9,10,17,12,13,14,16,18)],
            file = "AR/9.peak_annotation/AR_consensus.custombed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = c("chrom","start","end","peak_name","annotation","distanceToTSS","genechr","genestart","geneend","symbol","geneStrand","entrezID","transcriptID","ensembleID","gene_name"))

write.table(as.data.frame(IgG_annotation@anno)[, c(1,2,3,6,7,15,8,9,10,17,12,13,14,16,18)],
            file = "IgG/9.peak_annotation/IgG_consensus.custombed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = c("chrom","start","end","peak_name","annotation","distanceToTSS","genechr","genestart","geneend","symbol","geneStrand","entrezID","transcriptID","ensembleID","gene_name"))

write.table(as.data.frame(EZH2_annotation@anno)[, c(1,2,3,6,7,15,8,9,10,17,12,13,14,16,18)],
            file = "ezh2/9.peak_annotation/ezh2_consensus.custombed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = c("chrom","start","end","peak_name","annotation","distanceToTSS","genechr","genestart","geneend","symbol","geneStrand","entrezID","transcriptID","ensembleID","gene_name"))

write.table(as.data.frame(H3K27ac_annotation@anno)[, c(1,2,3,6,7,15,8,9,10,17,12,13,14,16,18)],
            file = "h3k27ac/9.peak_annotation/h3k27ac_consensus.custombed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = c("chrom","start","end","peak_name","annotation","distanceToTSS","genechr","genestart","geneend","symbol","geneStrand","entrezID","transcriptID","ensembleID","gene_name"))

write.table(as.data.frame(H3K27me3_annotation@anno)[, c(1,2,3,6,7,15,8,9,10,17,12,13,14,16,18)],
            file = "h3k27me3/9.peak_annotation/h3k27me3_consensus.custombed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = c("chrom","start","end","peak_name","annotation","distanceToTSS","genechr","genestart","geneend","symbol","geneStrand","entrezID","transcriptID","ensembleID","gene_name"))

write.table(as.data.frame(H3K4me1_annotation@anno)[, c(1,2,3,6,7,15,8,9,10,17,12,13,14,16,18)],
            file = "h3k4me1/9.peak_annotation/h3k4me1_consensus.custombed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = c("chrom","start","end","peak_name","annotation","distanceToTSS","genechr","genestart","geneend","symbol","geneStrand","entrezID","transcriptID","ensembleID","gene_name"))

write.table(as.data.frame(H3K4me3_annotation@anno)[, c(1,2,3,6,7,15,8,9,10,17,12,13,14,16,18)],
            file = "h3k4me3/9.peak_annotation/h3k4me3_consensus.custombed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = c("chrom","start","end","peak_name","annotation","distanceToTSS","genechr","genestart","geneend","symbol","geneStrand","entrezID","transcriptID","ensembleID","gene_name"))




write.table(AR_annotation@annoStat,
            file = "AR/9.peak_annotation/summary_consensus.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(IgG_annotation@annoStat,
            file = "co-/9.peak_annotation/summary_consensus.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(EZH2_annotation@annoStat,
            file = "ezh2/9.peak_annotation/summary_consensus.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(H3K27ac_annotation@annoStat,
            file = "h3k27ac/9.peak_annotation/summary_consensus.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(H3K27me3_annotation@annoStat,
            file = "h3k27me3/9.peak_annotation/summary_consensus.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)
            
write.table(H3K4me1_annotation@annoStat,
            file = "h3k4me1/9.peak_annotation/summary_consensus.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(H3K4me3_annotation@annoStat,
            file = "h3k4me3/9.peak_annotation/summary_consensus.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)
