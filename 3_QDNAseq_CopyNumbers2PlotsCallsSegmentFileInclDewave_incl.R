# FUNCTIONS ============================================================================================
# Function sampleNames: return the sample names from a QDNAseq object
sampleNames <- function(largeQDNAseqCopyNumbers){
  sampleNames <- {largeQDNAseqCopyNumbers@phenoData}@data$name
  return(sampleNames)
}

# Function pathToScriptWithScriptname
pathToScriptWithScriptname <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
# End of FUNCTIONS =====================================================================================


# create stop function which exits with a blank message (still, an error is generated!)
# TODO: use another way to terminate the program without generating an error
# source: https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()


#options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(trailingOnly = TRUE)

#print("Command line arguments set: ", args)
cat("Command line arguments set: ")
print(args)

isRStudio <- Sys.getenv("RSTUDIO") == "1" # running this script from within RStudio is equivalent to using option --interactive
if(isRStudio) {args[1] <- "--interactive"}
print(isRStudio)


if((length(args)==0) || args[1] == "--help" || args[1] == "-h" || !(isRStudio || args[1] == "--interactive" || (4 <= length(args)) && (length(args) <= 7))) {
  cat(paste("Usage:", 
            "* To show this help: Rscript QDNAseq_CopyNumbers2PlotsCallsSegmentFileInclDewave.R --help OR -h OR without option",
            "* Interactive mode with GUIs: Rscript QDNAseq_CopyNumbers2PlotsCallsSegmentFileInclDewave.R --interactive",
            "   Requirements",
            "   R must be able to open graphical windows (X11), and TCL & TK must be installed",
            
            "* Non-interactive mode: Rscript QDNAseq_CopyNumbers2PlotsCallsSegmentFileInclDewave.R followed by command line inputs:",
            "  in case of 4 arguments: /path/to/copyNumbers-xxkb-bins.rds /path/to/OutputDir projectname dewave ('yes' or for dewaving, any other character or word for no dewaving)",
            "  in case of 5 arguments: first-4-arguments /path/to/inclusionlist.xls",
            "  in case of 6 arguments: first-4-arguments undo.SD alpha",
            "  in case of 7 arguments: first-4-arguments /path/to/inclusionlist.xls undo.SD alpha.",
            "Note:",
            "  copy number file must look like this: *copyNumbers-xxkb-bins.rds",
            "  inclusionlist.xls must have in the first column a list of unique sample names",
            "",
            sep="\n"))
  
  stopQuietly()
}






# load packages and libraries ...


# elegant way to install R.cache package if not yet installed, according
# http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
# ... from CRAN
list.of.packages <- c("devtools", "R.cache", "MASS", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')

# ... from Bioconductor
list.of.bioconductor.packages <- c("QDNAseq", "QDNAseq.hg19", "CGHregions", "Biobase", "limma", "impute")
new.bioconductor.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.bioconductor.packages)) {
  #source("http://bioconductor.org/biocLite.R")
  #biocLite(new.bioconductor.packages)
  BiocManager::install(new.bioconductor.packages)
}



library(QDNAseq)
library(QDNAseq.hg19)
library(R.cache, quietly = TRUE)


# 4 arguments: classic
# 5 arguments: classic plus inclusion list

# if interactive:
if(isRStudio || args[1]=="--interactive") {
  # TODO: Need to check if graphical environment is available and tcltk can be used, if not exit with error message.
  list.of.packages <- c("tcltk", "gWidgetstcltk", "gWidgets", "gtools", "MASS", "stringr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')
  library(tcltk)
  library(gWidgets)
  options(guiToolkit="tcltk") 
  copynumberfile = tk_choose.files(caption = "Please select your copyNumbers.rds file. It must look like this: *copyNumbers-xxkb-bins.rds")
  pathToOutputDir = tk_choose.dir(getwd(), "Please select the directory (=folder) in which your output directory will be created (select with double click and OK)")
  projectname <- ginput("Please enter a project-name or working-name for this run. It will be used as part of all output files", title="Choose project-name", icon ="question")
  args=c(copynumberfile, pathToOutputDir, projectname)
  dewaving <- ginput("'Dewaving' is currently an experimental and unpublished method in QDNAseq, based on Mark van de Wiel's NoWaves algorithm. If you would like this method to be applied to your data, enter 'yes'. Otherwise enter 'no'.", title="Dewaving", icon ="question")
  inclusionFileQuestion <- ginput("Do you want to continue your analysis with all samples (enter 'all') or do you have an inclusion list (enter 'i')")
  if (inclusionFileQuestion == "i") {
    inclusionFile <- tk_choose.files(caption = "Please select an inclusion file. This should be in XLS or XLSX format and the first column should contain a list of samples you want to keep")
    args=c(copynumberfile, pathToOutputDir, projectname, dewaving, inclusionFile)
    } else {args=c(copynumberfile, pathToOutputDir, projectname, dewaving)}
}



undoSD=1 # default in QDNAseq segmentation
alpha=1e-10 # default in QDNAseq segmentation


# 4 arguments: classic
# 5 arguments: classic plus inclusion list
# 6 arguments: classic plus undoSD and alpha
# 7 arguments: classic plus undoSD and alpha and inclusion list

# 4 arguments: classic
# /path/to/copyNumbers-xxkb-bins.rds /path/to/OutputDir projectname dewaving ('yes' or 'no' for dewaving or not)",
if(length(args)==4) {
  copynumberfile  <- args[1]
  pathToOutputDir <- args[2]
  projectname     <- args[3]
  dewaving     <- args[4]
} else if(length(args)==5) {
  copynumberfile  <- args[1]
  pathToOutputDir <- args[2]
  projectname     <- args[3]
  dewaving     <- args[4]
  inclusionFile   <- args[5]
} else if(length(args)==6) {
  copynumberfile  <- args[1]
  pathToOutputDir <- args[2]
  projectname     <- args[3]
  dewaving     <- args[4]
  undoSD          <- args[5]
  alpha           <- args[6]
} else if(length(args)==7) {
  copynumberfile  <- args[1]
  pathToOutputDir <- args[2]
  projectname     <- args[3]
  dewaving     <- args[4]
  undoSD          <- args[5]
  alpha           <- args[6]
  inclusionFile   <- args[7]
}



if(!file.exists(copynumberfile)|| !file.exists(pathToOutputDir))
{
  stop("Please check: copynumberfile must exist and its filname be of format *copyNumbers-xxkb-bins.rds, /path/to/OutputDir (must exist)")
}

binSize=as.integer(gsub(".*copyNumbers-|kb-bins.rds", "", copynumberfile))
copyNumbers = readRDS(copynumberfile)

if(exists("inclusionFile")) {
  if (!file.exists(inclusionFile)) { exit(paste0("The inclusion file ",inclusionFile," you've provided does not exist"))}
  else {
    list.of.packages <- c("readxl")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')
    library(readxl)
    inclusionlist <- read_excel(inclusionFile, col_names = TRUE, col_types = "text")
    copyNumbersInclusionList =  copyNumbers[,copyNumbers@phenoData@data$name %in% unlist(inclusionlist[,1])]
    copyNumbers = copyNumbersInclusionList
  }
}
library(Biobase)
copyNumbersNormalized <- normalizeBins(copyNumbers) # normalize bins; default is median normalization
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized) # smooth outliers based on previous normalization

if(dewaving == 'yes') {

  # ... from github: e.g.: github.com/tgac-vumc/QDNAseq.dev/tree/dewave
  if(!"QDNAseq.dev" %in% installed.packages()) {
    library(devtools)  
    install_github("tgac-vumc/QDNAseq.dev", ref = "dewave")
  }
  
  projectname = paste(projectname, "dewave",sep=".")
  library(QDNAseq.dev)
  library(Biobase)
  if(binSize == 15) {NormalCalibrationSet_ofcurrentBinSize = NormalCalibrationSet_15kb}
  if(binSize == 30) {NormalCalibrationSet_ofcurrentBinSize = NormalCalibrationSet_30kb}
  if(binSize == 100) {NormalCalibrationSet_ofcurrentBinSize = NormalCalibrationSet_100kb}
  if(binSize == 1000) {NormalCalibrationSet_ofcurrentBinSize = NormalCalibrationSet_1000kb}
  if(!binSize %in% c(15, 30, 100, 1000)) {stop("Your binSize is not one of 15kb, 30kb, 100kb or 1000kb. Dewaving is currently not possible for other bin sizes.")}
  
  copyNumbersDewaved = dewaveBins(copyNumbersSmooth, NormalCalibrationSet_ofcurrentBinSize)
  copyNumbersSmooth = copyNumbersDewaved
}



dir.create(file.path(pathToOutputDir, projectname, paste0(binSize,"kb-bins")), showWarnings = FALSE, recursive=TRUE)
smoothPlotDir <- paste0(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/SmoothPlot/")
dir.create(smoothPlotDir, showWarnings = FALSE, recursive=TRUE)
for (i in 1:ncol(copyNumbersSmooth)) {
  png.name <- paste(smoothPlotDir, sampleNames(copyNumbersSmooth)[i], "_", projectname,"_SmoothPlot_",binSize,"kb.png", sep="")
  png(png.name, width = 1280, height = 1024)
  plot(copyNumbersSmooth[,i])
  dev.off()
}

# TODO!! adjust output filename to input bam, and binSize

# EXAMPLES:
# exportBins(copyNumbersSmooth, file="LGG150.txt")
# exportBins(copyNumbersCalled, file="Calls-50kb.txt")
# exportBins(copyNumbersSmooth, file="LGG150.igv", format="igv")
# exportBins(copyNumbersSmooth, file="LGG150.bed", format="bed")
# copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
# copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="log2", undo.SD=1)

# experimental:
 copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt", undo.SD=undoSD, alpha=alpha, segmentStatistic ="seg.median")

# original: 
#copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt", undo.SD=undoSD, alpha=alpha)




copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

sgmentedPlotDir <- paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/SegmentedPlot/", sep="")
dir.create(sgmentedPlotDir, showWarnings = FALSE, recursive=TRUE)
for (i in 1:ncol(copyNumbersSegmented)) {
  png.name <- paste(sgmentedPlotDir, sampleNames(copyNumbersSegmented)[i], "_", projectname,"_SegmentedPlot_",binSize,"kb.png", sep="")
  png(png.name, width = 1280, height = 1024)
  plot(copyNumbersSegmented[,i])
  dev.off()
}


copyNumbersCalled <- callBins(copyNumbersSegmented)


calledPlotDir <- paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/CalledPlot/", sep="")
dir.create(calledPlotDir, showWarnings = FALSE, recursive=TRUE)
for (i in 1:ncol(copyNumbersCalled)) {
  png.name <- paste(calledPlotDir, sampleNames(copyNumbersCalled)[i], "_", projectname,"_CalledPlot_",binSize,"kb.png", sep="")
  png(png.name, width = 1280, height = 1024)
  plot(copyNumbersCalled[,i])
  dev.off()
}

saveRDS(copyNumbersCalled, file = paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/",projectname,"-copyNumbersCalled-",binSize,"kb-bins.rds", sep=""))
#saveRDS(copyNumbersCalled, file = "copyNumbersCalled-30kb-65bpReads.rds")
# copyNumbersCalled


library(CGHregions)


cghcall <- makeCgh(copyNumbersCalled)


#cghcall_df <- data.frame(featureNames(cghcall), chromosomes(cghcall), bpstart(cghcall), bpend(cghcall), calls(cghcall), copynumber(cghcall), segmented(cghcall))
cghcall_df <- data.frame(featureNames(cghcall), chromosomes(cghcall), bpstart(cghcall), bpend(cghcall), calls(cghcall))
names(cghcall_df)[1] <- paste("Name")
names(cghcall_df)[2] <- paste("Chromosome")
names(cghcall_df)[3] <- paste("Start")
names(cghcall_df)[4] <- paste("End")

# export the calls
write.table(cghcall_df, paste(pathToOutputDir, "/", projectname , "/", binSize, "kb-bins/",projectname,"-copyNumbersCalled-",binSize,"kb-bins.tab", sep=""),col.names=TRUE, row.names=FALSE,quote=F,sep="\t")
# export the the normalized log2 readcounts of the segments
exportBins(copyNumbersSegmented, file = paste(pathToOutputDir, "/", projectname , "/", binSize, "kb-bins/",projectname,"-copyNumbersSegmented-",binSize,"kb-bins.tab", sep=""), format="tsv", type="segments")


frequencyPlot(copyNumbersCalled)
pdf(file=paste(pathToOutputDir, "/", projectname , "/", binSize, "kb-bins/", projectname, "-frequencyPlot-", binSize, "kb-bins.pdf", sep=""))
frequencyPlot(copyNumbersCalled)
dev.off()

# calculate statistics matrix:
OutputDirStats <- paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/stats/", sep="")
dir.create(OutputDirStats, showWarnings = FALSE, recursive=TRUE)

# get the path where this script is in:
currentpathToScriptWithScriptname = pathToScriptWithScriptname()
library(stringr)
#currentpathToScript = str_extract(string = currentpathToScriptWithScriptname , "/.*/")
currentpathToScript = str_extract(string = currentpathToScriptWithScriptname , "/.*/|.*\\\\") # /Unix/ OR C:\Windows\ style of path

source(paste0(currentpathToScript, "/QDNAseq-observedVariance.R")) # this code will take a copyNumbers-object as input, calculate segments, var_expect, var_observed, diffvar, total_reads and return them wrapped as dataframe statsDF

Outlier_removal_in_R_using_IQR_rule <- dget(paste0(currentpathToScript, "/Outlier_removal_in_R_using_IQR_rule.R"))

xlsOutputFile <- paste(projectname, "-statistics-", binSize, "kb-bins.xlsx", sep="")
Outlier_removal_in_R_using_IQR_rule(OutputDirStats, xlsOutputFile, statsDF)



# CGHregions: calculation of regions

filenameOfCGHregionsOutput <- paste0("/", projectname, "-", binSize, "kb-bins.point01percentsmoothing.CGHregions.tab")



# copyNumbersCalled = readRDS(filenameOfCopyNumbersCalledRDSfile) # needed for importing copy numbers here. 
# cghcall <- makeCgh(copyNumbersCalled) # this was done already above
#cghregions_normalRegioning1percent <- CGHregions(cghcall, averror=0.01) # 1% error rate is considered "normal" according to CGHregions paper, when intending to compare groups
#cghregions_lenientRegioningPoint01percent <- CGHregions(cghcall, averror=0.0001) # 0.01 % error rate will result in extremely lenient smoothing
#cghregions_severeRegioning2point5percent <- CGHregions(cghcall, averror=0.025) # 2.5 % error rate will result in extremely dramatic smoothing, when intending to compare groups with <= 10 members
cghregions <- CGHregions(cghcall, averror=0.0001) # =0.01 % error rate will result in extremely lenient smoothing

cghregions_df <- data.frame( chromosomes(cghregions), bpstart(cghregions), bpend(cghregions), nclone(cghregions), avedist(cghregions), regions(cghregions) )
write.table(cghregions_df, paste0(pathToOutputDir, "/", projectname , "/", binSize, "kb-bins/",filenameOfCGHregionsOutput), quote=F, sep="\t")

CGHregionsDF2stats <- dget(paste0(currentpathToScript, "/CGHregionsDF2stats.R"))
CGHregionsDF2stats(cghregions_df, binSize, paste0(OutputDirStats,filenameOfCGHregionsOutput))

