# FUNCTIONS ============================================================================================

# Function to create a frequency plot
frequencyPlotPerInputClass <- function(qdnaseqObject_copyNumbers, clinData, sampleGroupsColumn, selectedClass, filenameOfFrequencyplot) {
  mysamples <- qdnaseqObject_copyNumbers@phenoData@data$name
  # Every element of 'mysamples' is tried to be matched with the pattern(s).
  # first argument of grepl is the pattern and is supposed to be the filenames or unique substrings of the filenames.
  # mymatches is a TRUE/FALSE list of the samples in the copyNumbers-object.
  # the patters can be (unique!) substrings of the filenames in 'mysamples'.
  # Note that each substring needs to be unique and therefore match at max one sequence in the QDNAseq object
  mymatches <- grepl(paste(clinData$Filename[clinData[sampleGroupsColumn]==selectedClass],collapse="|"), mysamples)
 # browser()
  png(file=filenameOfFrequencyplot, width = 1280, height = 1024)
  frequencyPlot(qdnaseqObject_copyNumbers[,mymatches])
  # for testing: frequencyPlot(qdnaseqObject_copyNumbers)
  dev.off()
  filenameOfFrequencyplotTXT = sub("png", "txt", filenameOfFrequencyplot, perl=TRUE)
  write.table(mysamples[mymatches], filenameOfFrequencyplotTXT, quote=F, sep="\t", col.names = FALSE, row.names = FALSE)
}
# Function stopQuietly: create stop function which exits with a blank message (still, an error is generated!)
# TODO: use another way to terminate the program without generating an error
# source: https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()
# End of FUNCTIONS ============================================================================================

# Start info on command line options ==========================================================================

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

#print("Command line arguments set: ", args)
cat("Command line arguments set: ")
print(args)

isRStudio <- Sys.getenv("RSTUDIO") == "1" # running this script from within RStudio is equivalent to using option --interactive
if(isRStudio) {args[1] <- "--interactive"}
print(isRStudio)

if((length(args)==0) || args[1] == "--help" || args[1] == "-h" || !(isRStudio || length(args)==5)) {
  cat(paste("Usage:", 
            "* To show this help: Rscript QDNAseq_FrequencyCGHregionsCGHtestCMDL.R --help OR -h OR without option",
            "* There is no interactive mode for this program version",

            "* Non-interactive mode: Rscript QDNAseq_FrequencyCGHregionsCGHtestCMDL.R requiredArguments" ,
            "   Required arguments for command line input:",
            "   arg1: '/path/to/your/CopyNumbersCalled-30kb-bins.rds copy numbers called file. E.g. CopyNumbersCalled-30kb-bins.rds",
            "   arg2: '/path/to/OutputDir' must exist",
            "   arg3: '/path/to/your/clinicalData.xlsx'",
            "   arg4 'clinicalDataTab' is the tab that will be used on clinicalData.xlsx. It must contain a column named 'Filename'.",
            "                             --> all filenames in the Excel must be present in the CopyNumbersCalled-30kb-bins.rds",
            "   arg5: 'groupName' must be a column containing the group IDs",
            "",
            sep="\n"))
  
  stopQuietly()
}

# End info on command line options ==========================================================================
# LOAD packages and libraries  =========================================================================

# ... from CRAN
list.of.packages <- c("R.cache", "devtools", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')

# ... from Bioconductor
list.of.bioconductor.packages <- c("QDNAseq", "QDNAseq.hg19", "CGHregions", "openxlsx")
new.bioconductor.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.bioconductor.packages)) {
  #source("http://bioconductor.org/biocLite.R") now is BiocManager which is installed with CRAN
  #biocLite(new.bioconductor.packages)
  BiocManager::install(new.bioconductor.packages)
}


# ... from remote file to be downloaded from URL (here: Mark vd Wiels web site: http://www.few.vu.nl/~mavdwiel/CGHtest.html)
if(!"CGHtest" %in% installed.packages()) {
  install.packages('http://www.few.vu.nl/%7Emavdwiel/CGHtest/CGHtest_1.1.tar.gz',repos = NULL)
}
#if(!"CGHtestpar" %in% installed.packages()) {
#  install.packages('http://www.few.vu.nl/~mavdwiel/CGHtest/CGHtestpar_0.0.tar.gz',repos = NULL)
#}

# ... from local file previously downloaded (from Mark vd Wiels web site)
#install.packages('CGHtest_1.1.tar.gz', lib='/home/christian/Rmylocallibs/',repos = NULL)
#library('CGHtest', lib='/home/christian/Rmylocallibs/')


library(CGHtest)
#library(CGHtestpar)
library(QDNAseq)
library(QDNAseq.hg19)
library(R.cache, quietly = TRUE)

library(CGHregions)
library(openxlsx)
# end of LOAD packages and libraries  =========================================================================

# arg1 '/path/to/your/CopyNumbersCalled-30kb-bins.rds copy numbers called file. E.g. CopyNumbersCalled-30kb-bins.rds",
# arg2 '/path/to/OutputDir' must exist",
# arg3 '/path/to/your/clinicalData.xlsx'",
# arg4 'clinicalDataTab' is the tab that will be used on clinicalData.xlsx. It must contain a column named 'Filename'.",
# arg5 'groupName' must be a column containing the group IDs",

args1 <- args[1]
args2 <- args[2]
args3 <- args[3]
args4 <- args[4]
args5 <- args[5]

#args1 <- "/DATA/c.rausch/Willem2019/testData/qdnaseq.dewave-copyNumbersCalled-30kb-bins.rds" # copy numbers RDS file
#args2 <- "/DATA/c.rausch/Willem2019/testData/"                            # name of the xlsx file with clin data
#args3 <- "/DATA/c.rausch/Willem2019/testData/Kopie_van_databaseR+J+cluster-junk_Willem2019.xlsx"
#args4 <- "Database_FIRST" # name of the tab with clinical data
#args5 <- "TRG_1+2_vs_3+4+5"
#args1 <- "/DATA/share/hugo/stat-compar/ovarium.dewave-copyNumbersCalled-30kb-bins.rds"
#args2 <- "/DATA/share/hugo/stat-compar/statClin/"
#args3 <- "/DATA/share/hugo/stat-compar/CFMPB_297_CNVSeq_recodeLabels.xlsx"
#args4 <- "CNVSeq_labels"
#args5 <- "age_cat_recode"


filenameOfCopyNumbersCalled <- args1
outputdir                   <- args2
clinDataTableXLSXName       <- args3 
clinDataTab                 <- args4
sampleGroupsColumn          <- args5


CopyNumbersCalled = readRDS(file=filenameOfCopyNumbersCalled)
clinData <- read.xlsx(xlsxFile = clinDataTableXLSXName, sheet =clinDataTab)

samples <- CopyNumbersCalled@phenoData@data$name


clinData$Filename=gsub(".bam","",clinData$Filename)

clinDataColnames <- colnames(clinData)


#sampleSelection <- read.csv(header = FALSE, file = filenameOfSelection)
selectionClasses = levels(factor(clinData[,clinDataColnames %in% sampleGroupsColumn]))

print(selectionClasses)
typeof(selectionClasses)
length(selectionClasses)


notIgnore = selectionClasses != "ignore"
notIgnore
notUnknownCat = selectionClasses != "unknownCat"
notUnknownCat
remainder = notIgnore & notUnknownCat
remainder
reducedClasses = selectionClasses[remainder]
reducedClasses
selectionClasses=reducedClasses
#stopQuietly()

# do frequency plot for each of the classes
# run cghregions for each of pair of classes
# for each pair of classes run CGHtest for gains and losses



for (selectedClass in selectionClasses) {
  # example: selectedClass = "small_bowel"
  #selectedClass = "65andolder"
  # example: filenameOfCopyNumbersCalled = qdnaseq.dewave-copyNumbersCalled-30kb-bins.rds
  filenameOfFrequencyplot = gsub("copyNumbersCalled-.+kb-bins.rds",  paste0(sampleGroupsColumn, "_class_",selectedClass, ".frequency-plot.png"), filenameOfCopyNumbersCalled, perl=TRUE)
  frequencyPlotPerInputClass(CopyNumbersCalled, clinData, sampleGroupsColumn, selectedClass, filenameOfFrequencyplot)
}


combis = combn(selectionClasses, 2)

for (i in (1:dim(combis)[2])) {
# for testing
#i=1
  print(combis[1:2,i])
  
  # don't do the following because the files of the two classes could be interleaved
  #matches <- grepl(paste(sampleSelection$V1[sampleSelection$V2==combis[1,i]|sampleSelection$V2==combis[2,i]],collapse="|"), samples)
  #cghcall <- makeCgh(CopyNumbersCalled[,matches])
  #cghregions <- CGHregions(cghcall, averror=0.01)
  # browser()
 
  
  matches1 <- grepl(paste(clinData$Filename[clinData[sampleGroupsColumn]==combis[1,i]],collapse="|"), samples)
  matches2 <- grepl(paste(clinData$Filename[clinData[sampleGroupsColumn]==combis[2,i]],collapse="|"), samples)
  cghcall1 <- makeCgh(CopyNumbersCalled[,matches1])
  cghcall2 <- makeCgh(CopyNumbersCalled[,matches2])
  cghcallcombined <-  combine(cghcall1, cghcall2)
  numberOfSamples1 <- dim(cghcall1)[2]
  numberOfSamples2 <- dim(cghcall2)[2]
  cghregions <- CGHregions(cghcallcombined, averror=0.01)
  
  # for testing:
  #cghcall14 <- makeCgh(CopyNumbersCalled[,1:4])
  #cghregions <- CGHregions(cghcall14, averror=0.01)
  datainfo <- data.frame(chromosomes(cghregions), bpstart(cghregions), bpend(cghregions), nclone(cghregions), avedist(cghregions))
  datacgh <- regions(cghregions)
  gfr <- groupfreq(datacgh,group=c(numberOfSamples1,numberOfSamples2),groupnames=c(combis[1,i],combis[2,i]),af=0.1)
  filename_prefix_CGHtestOutput = gsub("_class_.*frequency-plot.png",  paste0("_",combis[1,i], "_vs_",combis[2,i], ".CGHtest"), filenameOfFrequencyplot, perl=TRUE)
  teststats = c("Wilcoxon", "Chi-square")
  # for testing:
  teststats = "Chi-square"
 # j = 1
  for (j in c(-1, 1)) {
    if (j == -1) gainORloss = "loss"
    if (j == 1) gainORloss = "gain"
    
    for (k in (1:length(teststats))) {
      pvs <- pvalstest(datacgh,datainfo, teststats[k] ,group=c(numberOfSamples1,numberOfSamples2),groupnames=c(combis[1,i], combis[2,i]), af=0.1, niter=20000, lgonly = j)
      fdrs <- fdrperm(pvs,mtdirection="stepup")
      filename_CGHtestOutput = paste(filename_prefix_CGHtestOutput, "_", teststats[k], "_", gainORloss, ".tab", sep="")
      write.table(fdrs, file=filename_CGHtestOutput, sep="\t", row.names = FALSE)
    }
  }
  filenameOfCGHregionsOutput = sub("CGHtest", "CGHregions.tab", filename_prefix_CGHtestOutput, perl=TRUE)
  df <- data.frame( chromosomes(cghregions), bpstart(cghregions), bpend(cghregions), nclone(cghregions), avedist(cghregions), regions(cghregions) )
  write.table(df, filenameOfCGHregionsOutput, quote=F, sep="\t")
  
}

