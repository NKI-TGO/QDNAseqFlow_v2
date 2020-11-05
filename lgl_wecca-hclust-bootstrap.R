# This program will use called copy number data, do 100 bootstrapping rounds and does a clustering for each of the 100 bootstrapping data sets.
# As clustering methods it uses...


#setwd("/DATA/c.rausch/bea_psc_ibd/comparison_sporadicIntervalCancer_IBDcancer_TumorIBDpsc_excludeSamplesJan2019/IBD-PSC_IBD-cancers/QDNAseq_ibd-psc_ibd-cancers/clustering/") # edit this working directory
setwd("/Users/l.gonzalez/testing_QDNAseq_Flow/mixdata/clustering/Nodewaved/") # edit this working directory


# Function 
WECCA.heatmapChocolateRedWhiteBlueDarkblue <-
  function (cghdata.regioned, dendrogram,...){
    ######################################################################################################
    # heatmap with hard calls and chrosome indicator 
    ######################################################################################################
    
    # calculated number of classes used
    nclass <- length(unique(as.numeric(cghdata.regioned$hardcalls))) #28/7/2011: adjustment MvdW
    
    # Preparation of the plotting of the heatmap
    # Generate alternating colors for chromosomes.
    chr.color <- rep("blue", dim(cghdata.regioned$hardcalls)[1])
    ids <- ((cghdata.regioned$ann[, 1]%%2) == 0)
    chr.color[ids] <- c("yellow")
    
    # Generate labels for begin points of chromosomes.    
    Y <- rep(FALSE, dim(cghdata.regioned$hardcalls)[1])
    for (i in 2:(dim(cghdata.regioned$ann)[1])) {
      if ((cghdata.regioned$ann[i - 1, 1] != cghdata.regioned$ann[i, 1])) {
        Y[i] <- TRUE
      }
    }
    Y[1] <- TRUE
    begin.chr <- rep("", dim(cghdata.regioned$ann)[1])
    begin.chr[Y] <- cghdata.regioned$ann[Y, 1]
    
    # The heatmap is plotted.    
    color.coding <- c("green", "red", "white", "blue", "darkblue")[1:nclass]
    heatmap(cghdata.regioned$hardcalls, Colv = as.dendrogram(dendrogram), 
            Rowv = NA, col = color.coding, labRow = begin.chr,RowSideColors = chr.color, 
            scale = "none",...) #28/7/2011: adjustment MvdW, added ... to allow for using heatmap options in the mother function
  }
### END Function definition.


#setwd("/DATA/c.rausch/Copynumbers2019/clustering/")


#Load the devtools package.
#install.packages("glue")
#install.packages("usethis")
#install.packages("devtools")


#In most cases, you just use install_github("author/package"). For example, with my R/broman package, 
#which exists at github.com/kbroman/broman, you’d type
library(devtools)
library("glue")
#Run this part just one time
#install_github("tgac-vumc/WECCA")

#There’s some extra fanciness that you need to do if the version you want sits on some branch of the repository, 
#or if the package is in a subdirectory of the main repository.
#For example, Bill Engels has an R package HWxtest, but the package actually sits in the pkg subdirectory. 
#To install his package with install_github(), you’d need to do:
  
#install_github("tgac-vumc/WECCA", subdir="pkg")
  
  
library(WECCA)






# generate object of cghCall-class
#data(WiltingCalled)
# make region data (soft and hard calls)
#WiltingRegioned <- regioning(WiltingCalled)
library(QDNAseq)

#Called rds
#CopynumbersCalled=read.table("all.dewave-copyNumbersCalled-30kb-bins.tab", header = TRUE)
CopynumbersCalled=makeCgh(readRDS("/Users/l.gonzalez/testing_QDNAseq_Flow/mixdata/called30kb_Nodw/30Kb_Nodw-copyNumbersCalled-30kb-bins.rds")) # QDNAseqCopyNumbers object is loaded from RDS file and converted to CGHcall object

#Regions_01percerntSmoothing.tab
cghregions <- CGHregions(CopynumbersCalled, averror=0.01) # = 1 % error rate
cghregions_df <- data.frame( chromosomes(cghregions), bpstart(cghregions), bpend(cghregions), nclone(cghregions), avedist(cghregions), regions(cghregions) )
write.table(cghregions_df, "/Users/l.gonzalez/testing_QDNAseq_Flow/mixdata/called30kb_Nodw/30Kb_Nodw-30kb-bins.point01percentsmoothing.CGHregions.tab", quote=F, sep="\t")

CopynumbersRegioned <- regioning(CopynumbersCalled)
# clustering with hard.calls
#dendrogramhc <- WECCAhc(CopynumbersRegioned)
# clustering with soft.calls


dendrogramsc <-  WECCAsc(CopynumbersRegioned, dist.measure="ordinal", weight.type="all.equal", linkage="ward.D2")
#dendrogramsc <- WECCAsc(CopynumbersRegioned, dist.measure="ordinal", weight.type="all.equal", linkage="ward.D2")
#dendrogramsc <- WECCAsc(CopynumbersRegioned, dist.measure="ordinal", weight.type="all.equal", linkage="ward.D2")
#dendrogramsc <- WECCAsc(CopynumbersRegioned, dist.measure="KLdiv", weight.type="all.equal", linkage="ward")
#dendrogramsc <- WECCAsc(CopynumbersRegioned, dist.measure="KLdiv", weight.type="heterogeneity", linkage="ward")

mat=CopynumbersRegioned$softcalls
CopynumbersRegionedBootstrapped=CopynumbersRegioned

install.packages("ape")
library(ape)
#njtrees=list() 
for(i in 1:100) {
  #d<-dist(t(mat[sample(nrow(mat), nrow(mat), replace=TRUE), ]),method="euclidean")
  softcalls=mat[sample(nrow(mat), nrow(mat), replace=TRUE), ]
  CopynumbersRegionedBootstrapped$softcalls=softcalls
  dendrogramscBootstrapped <- WECCAsc(CopynumbersRegionedBootstrapped, dist.measure="ordinal", weight.type="all.equal", linkage="ward.D2")
  dendrogramscBootstrapped$labels=CopynumbersCalled$name
  tree<-dendrogramscBootstrapped
  phy=as.phylo(tree)
  write.tree(phy, file = paste0("bootstraptree_wecca_",i,".phy"), append = FALSE, digits = 10, tree.names = FALSE)
}




#plot(dendrogramscBootstrapped)

#WECCA.heatmap(CopynumbersRegionedBootstrapped, dendrogramscBootstrapped, margins =c(17,11))

#plot(dendrogramsc)
# dendro_sc <- WECCAsc(cghreg, dist.measure="KLdiv", weight.type="all.equal", linkage="ward") # default setting
# dist.measure: The distance measure to be used. This is either ‘"KLdiv"’ (the symmetric Kullback-Leibler divergence) or ‘"ordinal"’ (distance between the cumulative call probability distributions).

# dendro_sc <- WECCAsc(CopynumbersRegioned, dist.measure="ordinal", weight.type="all.equal", linkage="ward") # as used by Josien


#If the input is not an object of class cghCall it should be either a dataframe or a tabseparated textfile
#(textfiles must contain a header). The first three columns should contain the name, chromosome and
#position in bp for each array target respectively. The chromosome and position column must contain
#numbers only. Following these is a column with log2 ratios for each of your samples. If the input
#type is a textfile, missing values should be represented as ’NA’ or an empty field.

#pdf("Copynumbers-WECCA-hard-calls-heatmap.pdf", paper = "a4r", width = 8, height = 11)
#WECCA.heatmap(CopynumbersRegioned, dendrogramhc, margins =c(17,11))



pdf("Copynumbers-WECCA-soft-calls-heatmap_mixdata_Nodewaved.pdf", paper = "a4", width = 5, height = 50)
#a=WECCA.heatmap(CopynumbersRegioned, dendrogramsc, margins =c(17,11))
WECCA.heatmapChocolateRedWhiteBlueDarkblue(CopynumbersRegioned, dendrogramsc, margins =c(17,11))
dev.off()

# specify the number of clusters to be extracted from the dendrogram
#nclusters <- 3
#table.clusters.samples <- sample.cluster.table(CopynumbersRegioned, dendrogramsc, nclusters)
#View(table.clusters.samples)


#---- Bootstrapping with hclust (hierarchical clustering)


#CopyNumberCalled tab file
#calls=read.csv(file="../ibd-psc_ibd-cancers.dewave/30kb-bins/ibd-psc_ibd-cancers.dewave-copyNumbersCalled-30kb-bins.tab", row.names = 1, sep="\t")  ## adjust!!
calls=read.csv(file="/Users/l.gonzalez/testing_QDNAseq_Flow/mixdata/called30kb_Nodw/30Kb_Nodw-copyNumbersCalled-30kb-bins.tab", row.names = 1, sep="\t")  ## adjust!!

calls$Start=calls$Start-1
# remove .markdup from the sample names
colnames(calls)=c(colnames(calls)[1:3],gsub('.{8}$', '', colnames(calls)[4:length(colnames(calls))]))
colnames(calls)[1]=paste0("#",colnames(calls)[1])
write.table(calls, "calls.bed", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
library(readr)
calls <- read_delim("calls.bed", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)

library(gplots)
# previously used: rownames(genes_intersect_calls)=genes_intersect_calls$X4
mat = as.matrix(calls[,c(4:dim(calls)[2])])
#uniq=unique(mat)

#pdf("clustering_from_scratch_WGcalls.pdf")
#heatmap.2(mat, dendrogram=NULL, Colv = NA, Rowv = NA, distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2") , cexRow = 0.1, cexCol=0.1, trace = "none")
#dev.off()

#install.packages("ape")
library(ape)
#njtrees=list() 
for(i in 1:100) {
  d<-dist(t(mat[sample(nrow(mat), nrow(mat), replace=TRUE), ]),method="euclidean")
  tree<-hclust(d,method="ward.D2")
  #njtree <-nj(d)
  #njtrees[[i]]<-njtree
  #phy=as.phylo(njtree)
  phy=as.phylo(tree)
  write.tree(phy, file = paste0("bootstraptree_hierarch_allcalls_",i,".phy"), append = FALSE, digits = 10, tree.names = FALSE)
  
}

