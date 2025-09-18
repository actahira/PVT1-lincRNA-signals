#Input and output
##################

output_prefix <- "4.stat/fdr+sva"

input <- "4.stat/combined_metadata_table.tsv"
# Use table of metadata as input, including "Counts_file"
# column with respective file paths

test_variables <- c("Condition")
control_values <- c("NTC_EtOH")
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
library(edgeR)
library(dplyr)

#Formatting metadata, preparing count table and design
##################
meta <- read.delim(file = input)

design <- model.matrix(~0+Condition, data=meta)
colnames(design) <- c("A","B","C","D")

data <- readDGE(files = c(meta$Counts_file), columns = c(1,7), 
                group = design, 
                header = F, skip = 2)



#Filtering low expression RNAs
##################
keep <- filterByExpr(data,design =
                        design, min.prop=0.5)
data <- data[keep, , keep.lib.sizes=F]

#Calculate surrogate variables
##################

fullmodel = model.matrix(~Condition+Replicate, data=meta)
nullmodel = model.matrix(~Replicate, data=meta)
sva_obj <- svaseq(dat=data$counts, mod=fullmodel, mod0=nullmodel)
design <- cbind(design,sva_obj$sv)
colnames(design) <- c("A","B","C","D","sv1","sv2")


#Normalize and calculate dispersion
##################
data <- calcNormFactors(data)
data <- estimateDisp(data, design, robust = T)

#DE statistical test
##################
fit = glmQLFit(data, design, robust = T)

my.contrasts <- makeContrasts(co.vs.horm=B-A, co.vs.kd=C-A, horm.vs.comb=D-B, levels=design)
tests = c()
comparisons = colnames(my.contrasts)
for (i in colnames(my.contrasts)){
  assign(i,glmQLFTest(fit, contrast=my.contrasts[,i]))
  assign(paste(i,"_table",sep=""),eval(parse(text=paste(i,"$table",sep=""))))
  tests <- append(tests,paste(i,"_table",sep=""))
}

#p-value correction
##################
for (i in tests){
  assign(i, mutate(eval(parse(text = i)),
                   p_Adjust = p.adjust(eval(parse(text = paste(i,"$PValue",sep=""))), method = "bonferroni"),
                   FDR = p.adjust(eval(parse(text = paste(i,"$PValue",sep=""))), method = "fdr")))
}
rm(i)


#Save output and graphs
##################
fc_cutoff = 0
pval_cutoff = 0.01
for (i in tests){
  write.table(eval(parse(text = i)),quote=F, sep="\t", col.names=NA,
              file=paste(output_prefix,comparisons[match(i,tests)],"allgenes.tsv",sep="_"))

   assign(paste("DEgenes_",comparisons[match(i,tests)],sep=""),
         filter(eval(parse(text = i)),abs(logFC)>=fc_cutoff & FDR<=pval_cutoff))
  write.table(filter(eval(parse(text = i)),abs(logFC)>=fc_cutoff & FDR<=pval_cutoff),
              quote=F, sep="\t", col.names=NA,
              file=paste(output_prefix,comparisons[match(i,tests)],"DEgenes.tsv",sep="_"))
  print(paste("Filtering",comparisons[match(i,tests)],
              "by FDR,",
              nrow(eval(parse(text=paste("DEgenes_",comparisons[match(i,tests)],sep="")))),"DE genes found",
              sep=" "))

}
rm(i)

####################################

