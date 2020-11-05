# LOAD packages and libraries  =========================================================================

# ... from CRAN
list.of.packages <- c("R.cache","xlsx", "devtools", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')

# ... from Bioconductor
list.of.bioconductor.packages <- c("QDNAseq", "QDNAseq.hg19", "CGHregions", "CGHcall", "GeneBreak")
new.bioconductor.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.bioconductor.packages)) {
  #source("http://bioconductor.org/biocLite.R") now is BiocManager which is installed with CRAN
  #biocLite(new.bioconductor.packages)
  BiocManager::install(new.bioconductor.packages)
}



# NEED to EDIT: Set your working directory!

#setwd("/home/christian/surfdrive/Rosa_van_der_Kaaij/Rosa-Willem-combined/30kb-bins_called")
#setwd("/DATA/c.rausch/Willem2019/testData")
setwd("/Users/l.gonzalez/testing_QDNAseq_Flow/mixdata/genebreak/")
library(QDNAseq)
#copyNumbersCalled=readRDS("allmarkdup.dewave-copyNumbersCalled-30kb-bins.rds")
#copyNumbersCalled=readRDS("h38-copyNumbersCalled-30kb-bins.rds")
#Load CopyNumberCalled rds
copyNumbersCalled=readRDS("../called30kb_dewaving/mixdata30_dew.dewave-copyNumbersCalled-30kb-bins.rds")


library(CGHregions)

cghcall <- makeCgh(copyNumbersCalled)

cghcall_df <- data.frame(featureNames(cghcall), chromosomes(cghcall), bpstart(cghcall), bpend(cghcall), segmented(cghcall))
names(cghcall_df)[1] <- paste("Name")
names(cghcall_df)[2] <- paste("Chromosome")
names(cghcall_df)[3] <- paste("Start")
names(cghcall_df)[4] <- paste("End")
samples=colnames(cghcall_df)[-c(1:4)]
samples_shorter=gsub(".markdup$", "", samples, perl = T)
samples=samples_shorter
colnames(cghcall_df)[-c(1:4)]=samples
# ----------------------------------
# Run GeneBreak
# CHECK if NEEDED:  restrict samples to samples which are not junk
# restrict samples to samples which are not junk
#library(readxl)
#allButJunk_sequenced <- read_excel("/media/sf_surfdrive/tcga_esophagus/allButJunk-sequenced.xlsx", 
#                                   sheet = "allButJunk-sequenced")
#sequencedAndNotJunk = allButJunk_sequenced$sequenced_Rnaming[allButJunk_sequenced$sequencedAndNotJunk == "Yes"]
#samples_shorter=gsub(".markdup$", "", samples, perl = T)

#cghcall=cghcall[,samples_shorter %in% sequencedAndNotJunk]
# End of 'CHECK if NEEDED'

library(GeneBreak)
library(CGHcall)
breakpoints <- getBreakpoints( data = cghcall )

breakpointsFiltered <- bpFilter( breakpoints, filter = "CNA-ass" )

data( "ens.gene.ann.hg19") 
breakpointsAnnotated <- addGeneAnnotation( breakpointsFiltered, ens.gene.ann.hg19 )

#Tho show the names of associated features of e.g.  the ”PCMTD2” gene, give:
featuresPerGene ( breakpointsAnnotated , geneName = "PCMTD2" )

geneFeatures <- geneInfo( breakpointsAnnotated )
head( geneFeatures[ , c("Gene", "Chromosome", "Start", "End", "featureTotal", "featureNames", "remarks") ] )
breakpointGenes <- bpGenes( breakpointsAnnotated )

result_BreakpointGenes <- geneInfo ( breakpointGenes )
head( result_BreakpointGenes[ which ( result_BreakpointGenes$sampleCount > 0 ) , c( "Gene", "Chromosome", "Start", "End", "featureTotal", "nrOfBreakLocations", "sampleCount", "sampleNamesWithBreakpoints") ] )

# Detection of recurrent breakpoint genes
breakpointStatistics <- bpStats( breakpointGenes, level = "gene", method = "Gilbert" )

head( recurrentGenes( breakpointStatistics ) )
recur=recurrentGenes( breakpointStatistics )

breakpointStatisticsInFeatures <- bpStats(breakpointStatistics, level = "feature", method = "BH" )

bpPlot( breakpointStatisticsInFeatures, fdr.threshold = 0.1, plot.chr=20 )
library(xlsx)
write.xlsx(recur, "recur.xlsx")
# END Genebreak

# ----------------------------------
# start converter for RUBIC format

names(cghcall_df)[1] <- paste("Name")
names(cghcall_df)[2] <- paste("Chromosome")
names(cghcall_df)[3] <- paste("Start")
names(cghcall_df)[4] <- paste("End")

markers=data.frame(cghcall_df$Name, cghcall_df$Chromosome, (cghcall_df$End+cghcall_df$Start-1)/2)
colnames(markers)=c("Name", "Chromosome","Position")


write.table(markers, "30kb_markers_h19.tsv", sep = "\t", col.names = T, row.names = F, append = F, quote = F) 
#write.table(markers, "30kb_markers_h38.tsv", sep = "\t", col.names = T, row.names = F, append = F, quote = F) 



# CHECK if NEEDED:  restrict samples to samples which are not junk
# restrict samples to samples which are not junk
#library(readxl)
#allButJunk_sequenced <- read_excel("/media/sf_surfdrive/tcga_esophagus/allButJunk-sequenced.xlsx", 
 #                                  sheet = "allButJunk-sequenced")
#sequencedAndNotJunk = allButJunk_sequenced$sequenced_Rnaming[allButJunk_sequenced$sequencedAndNotJunk == "Yes"]
#samples_shorter=gsub(".markdup$", "", samples, perl = T)

#samples_notJunk=samples[samples_shorter %in% sequencedAndNotJunk]
#samples=samples_notJunk

# END: CHECK if NEEDED:  restrict samples to samples which are not junk

# start for loop over samples here

shift <- function(x, n){
  c(NA, x[-length(x)])
}


onedfsmall_colnames = paste("Sample", "Chromosome", "Start", "End", "MarkerCount", "LogRatio",sep="\t")
write.table(onedfsmall_colnames, "30kb_segcna.tsv", sep="\t", quote=F, col.names = F, row.names = F, append = F)
#write.table(onedfsmall_colnames, "30kb_segcna_h38.tsv", sep="\t", quote=F, col.names = F, row.names = F)


for (sample in samples) {
  # test: sample=samples[1]
  onedf=cghcall_df[,c("Name", "Chromosome", "Start", "End",sample)]
  
  onedf$Chromosome=factor(onedf$Chromosome)
  binSize=onedf$End[1]-onedf$Start[1]+1
  
  # in dataframe onedf: aggregate Start positions by Chromosome, using the min.
  chrStart=aggregate(Start ~ Chromosome, data=onedf, FUN=min)
  chrEnd=aggregate(End ~ Chromosome, data=onedf, FUN=max)
  
  chrStartEnd=data.frame(chrStart$Chromosome,chrStart$Start,chrEnd$End)
  colnames(chrStartEnd)=c("Chromosome", "Start", "End")
  
  
  colnames(onedf)[5]="logvalue"
  
 
  onedf$logvalueshift=onedf$logvalue
  onedf$logvalueshift<- shift(onedf$logvalueshift, 1)
  #rows = apply(onedf[, 5:6], 1, function(i) any(i[-1] != i[1])) # for more comparisons
  rows = apply(onedf[, 5:6], 1, function(i) (i[1] != i[-1]))
  rows[is.na(rows)] = TRUE
  onedfsmall = onedf[rows,]
  newend=NA
  maxrow=dim(onedfsmall)[1]
  
  for (row in c(1:maxrow)) {
    if(row<maxrow & onedfsmall$Chromosome[row]==onedfsmall$Chromosome[row+1]) {
      if(onedfsmall$Chromosome[row]==onedfsmall$Chromosome[row+1]) {
        onedfsmall$newend[row]=onedfsmall$Start[row+1]-1
      }
    }
    else { # so the Chromosome of next row is different
      onedfsmall$newend[row]=chrStartEnd$End[chrStartEnd$Chromosome==onedfsmall$Chromosome[row]]
    }
  }
  
  onedfsmall$Name=sample
  # before setting sample names remove undesired suffix
  sample_shorter=gsub(".markdup$", "", sample, perl = T)
  sample=sample_shorter
  onedfsmall$Name=sample
  onedfsmall$End <- NULL
  onedfsmall$logvalueshift <- NULL
  onedfsmall$numberOfBins = (onedfsmall$newend-onedfsmall$Start+1)/binSize
  
  # End for loop around all samples
  # colnames(onedfsmall) is: "Name"         "Chromosome"   "Start"        "logvalue"     "newend"       "numberOfBins"
  # colnames(onedfsmall) is: "Sample"         "Chromosome"   "Start"        "LogRatio"     "End"       "MarkerCount"
                            
  # not needed: colnames have been written to output file already colnames(onedfsmall)= c("Sample", "Chromosome", "Start", "LogRatio", "End", "MarkerCount")
  write.table(onedfsmall[,c(1,2,3,5,6,4)], "30kb_segcna.tsv", sep = "\t", col.names = F, row.names = F, append = T, quote = F)
 # write.table(onedfsmall[,c(1,2,3,5,6,4)], "30kb_segcna_h38.tsv", sep = "\t", col.names = F, row.names = F, append = T, quote = F)
}


write.table(samples, "30kb_mixsamples.tsv", sep = "\t", col.names = F, row.names = F, append = F, quote = F)
#write.table(samples, "30kb_samples_h38.tsv", sep = "\t", col.names = F, row.names = F, append = F, quote = F) 
getwd()
