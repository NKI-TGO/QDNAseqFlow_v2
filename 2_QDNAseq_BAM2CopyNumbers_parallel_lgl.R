# FUNCTIONS ============================================================================================
# Function stopQuietly: create stop function which exits with a blank message (still, an error is generated!)
# TODO: use another way to terminate the program without generating an error
# source: https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()

# End of FUNCTIONS =====================================================================================

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

#print("Command line arguments set: ", args)
cat("Command line arguments set: ")
print(args)

isRStudio <- Sys.getenv("RSTUDIO") == "1" # running this script from within RStudio is equivalent to using option --interactive
if(isRStudio) {args[1] <- "--interactive"}
print(isRStudio)

if((length(args)==0) || args[1] == "--help" || args[1] == "-h" || !(isRStudio || args[1] == "--interactive" || length(args)==4)) {
  cat(paste("Usage:", 
            "* To show this help: Rscript QDNAseq_BAM2CopyNumbers.R --help OR -h OR without option",
            "* Interactive mode with GUIs: Rscript QDNAseq_BAM2CopyNumbers.R --interactive",
            "   Requirements",
            "   R must be able to open graphical windows (X11), and TCL & TK must be installed",
            
            "* Non-interactive mode: Rscript QDNAseq_BAM2CopyNumbers.R /path/to/BAMdirectory /path/to/OutputDir binSize projectname" ,
            "   Requirements for non-interactive mode (command line input):",
            "   '/path/to/BAMdirectory' must exist",
            "   '/path/to/OutputDir' must exist",
            "   'binSize' must be one of 1, 5, 10, 15, 30, 50, 100, 500 or 1000 kbs",
            "   'projectname' must not contain spaces or special characters",
            "",
            sep="\n"))
  
stopQuietly()
}


# elegant way to install R.cache package if not yet installed, according
# http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them

list.of.packages <- c("R.cache", "future") # include "BiocManager" if we upgrade R version
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')

list.of.bioconductor.packages <- c("QDNAseq", "QDNAseq.hg19") 
new.bioconductor.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.bioconductor.packages)) {
  source("http://bioconductor.org/install.R")
  BiocManager::install("new.bioconductor.packages")
}

# ... from Bioconductor R 3.6 or later

#list.of.bioconductor.packages <- c("QDNAseq", "QDNAseq.hg19", "CGHregions", "Biobase", "limma", "impute")
#new.bioconductor.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[,"Package"])]
#if(length(new.bioconductor.packages)) {
#source("http://bioconductor.org/biocLite.R")
#biocLite(new.bioconductor.packages)
#  BiocManager::install(new.bioconductor.packages)
#}

library(QDNAseq)
library(QDNAseq.hg19)
library(R.cache, quietly = TRUE)
library("future")
plan(multiprocess)
options(future.globals.maxSize= 1000000000)



if(isRStudio || args[1]=="--interactive") {
  # TODO: Need to check if graphical environment is available and tcltk can be used, if not exit with error message.
  list.of.packages <- c("tcltk", "gWidgetstcltk", "gWidgets") # need to add packages for the GUI
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')
  library(tcltk)
  library(gWidgets)
  options(guiToolkit="tcltk") 
  pathToBAMs = tk_choose.dir(getwd(), "Please select the directory (=folder) which contains your bam files (select with double click and OK)")
  pathToOutputDir = tk_choose.dir(getwd(), "Please select the directory (=folder) in which your output directory will be created (select with double click and OK)")
  binSize <- as.integer(ginput("Enter a desired binsize. It can be: 1, 5, 10, 15, 30, 50, 100, 500 or 1000 kbp", title="Choose binsize", icon ="question"))
  projectname <- ginput("Please enter a project-name or working-name for this run. It will be used as part of all output files", title="Choose project-name", icon ="question")
  args=c(pathToBAMs, pathToOutputDir, binSize, projectname)
} else if(length(args)==4) {
  pathToBAMs <- args[1]
  pathToOutputDir <- args[2]
  binSize <- as.integer(args[3])
  projectname <- args[4] 
} else {
  exit("Invalid command line arguments. Check the possible arguments by running the program with option --help")}


binSizes <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000);

rm(args)

homeDir <- path.expand("~")
qdnaseqDir <- "cached-QDNAseq-binAnnotations"

dir.create(file.path(homeDir, qdnaseqDir), showWarnings = FALSE)

binAnnotationFile = file.path(homeDir, qdnaseqDir, paste("bin",binSize,"kb.rds",sep=""))
if(file.exists(binAnnotationFile)){
  bins <- readRDS(binAnnotationFile)
} else {
  bins <- getBinAnnotations(binSize)
  saveRDS(bins,binAnnotationFile)
}

readCounts <- binReadCounts(bins, path=pathToBAMs, cache=FALSE, chunkSize="longestChromosome")
readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)

# correct for GC content and mappability (applies data previously calculated using estimateCorrection)
copyNumbers <- correctBins(readCountsFiltered)

binSizeDir <- paste(pathToOutputDir, "/",projectname,"/", projectname, "-binSize-", binSize, "/", sep="")
dir.create(binSizeDir, showWarnings = FALSE, recursive = TRUE)
outputfilenameWithoutEnding= paste(binSizeDir,"/", projectname, "-copyNumbers-",binSize,"kb-bins", sep="")
saveRDS(copyNumbers, file = paste(outputfilenameWithoutEnding,".rds", sep=""))

# export as txt file if needed.
#exportBins(copyNumbers, file = paste(binSizeDir,"/", projectname, "-copyNumbers-",binSize,"kb.rds", sep=""))

exportBins(copyNumbers, file = paste(outputfilenameWithoutEnding,".tsv", sep=""), format="tsv", type="copynumber")
system(paste("bzip2 ", outputfilenameWithoutEnding,".tsv", sep=""))
exportBins(copyNumbers, file = paste(outputfilenameWithoutEnding,".igv", sep=""), format="igv", type="copynumber")
system(paste("bzip2 ", outputfilenameWithoutEnding,".igv", sep=""))
