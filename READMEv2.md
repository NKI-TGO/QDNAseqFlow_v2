# QDNAseqFlow 
## TGO group, NKI

### Introduction

QDNAseq is an R package, here we put together one previous step which include bwa alignment (aln or mem), sort files and mark duplicates

The main QDNAseqFlow include two steps: 
- Bam to CopyNumber
- CopyNumber to Called Plots

From there we can use the output of the called plots for: 

- CGHTest to compare two groups

- Clustering usingn WECCA
Hierarchical clustering (WECCA & hclust)
using called copy number data of regions (WECCA) and all bins (hclust). 100 bootstrapping rounds .
Calculation of Strict and Majority consensus clusterings (trees) in Dendroscope

- Recurrent broken genes using GeneBreak
detection of recurrent DNA copy number aberration-associated chromosomal breakpoints within genes

__This information is deprecated but the repository still exist and you can clone or check the information as it is for that version of the pipeline__
--Please have a look at [QDNAseqFlow-Abstract.pdf](https://github.com/NKI-Pathology/QDNAseqFlow/blob/master/QDNAseqFlow-Abstract.pdf) and [QDNAseqFlow\_poster\_ISMB.pdf](https://github.com/NKI-Pathology/QDNAseqFlow/blob/master/QDNAseqFlow_poster_ISMB.pdf) for an introduction.
To better understand how QDNAseq, the main component of this pipeline, works: [Christian-Rausch-QDNAseq-talk\_BioconductorDec2015.pdf](https://github.com/NKI-Pathology/QDNAseqFlow/blob/master/Christian-Rausch-QDNAseq-talk_BioconductorDec2015.pdf)--

This version 2017 is from 
Go to https://github.com/NKI-Pathology/QDNAseqFlow and click button "Clone or Download" --> Download zip.
Unzip the zip file in your Home directory, 'My Documents' or where you are allowed to install programs.

### Installation
New version is from 2020
Go to https://gitlab.rhpc.nki.nl/l.gonzalez/qdnaseq_workflow and click the button "Clone or Download" --> Download zip.
Unzip the zip file in your Home directory, 'My Documents' or where you are allowed to install programs.

**Required R-version etc.**

In the previous version:
R: We have made most tests with R 3.4.1. QDNAseqFlow might work with older R versions but the risk is that some required libraries might not be available.
_The step 2 (Bam to CopyNumber) is still running in the R3.4 version (server)_

R: Update 2020, the new modifications allow to use the latest R version, using BiocManager to install the neccesary packages. Be waare that you need permision to install packages

Java: You need to hava java installed. Make sure it is in your path. On Windows, this is done like here: https://confluence.atlassian.com/doc/setting-the-java_home-variable-in-windows-8895.html

*Note: In the server rjava libraries are missing and I do not have admin rights, Waiting for Rubayte


### Usage

This workflow consists of 3 main R-programs, that need to be run consecutively:

1. QDNAseq\_BAM2CopyNumbers.R
2. QDNAseq\_CopyNumbers2PlotsCallsSegmentFileInclDewave.R
3. QDNAseq\_FrequencyCGHregionsCGHtest.R

Update, the 2019 versions of those steps are not interactluve. Steps 1 and 2, they need to be run from the terminal in the server, and step 3 can be run in a local computer. 
*Note: In the server rjava libraries are missing and I do not have admin rights, Waiting for Rubayte


For the new scrips there is one way to run them

1. From Rstudio: open script file and execute with CTRL-ALT-R

NOt working at this time yet

3. From the command line:

bash script
1- Alignment: bash path/to/script/1_bwaAlnv2_lgl.sh
Need to make the sam index just right after

R script
2- Bam to CopyNumbers
	/path/to/your/r-installation/Rscript QDNAseq\_BAM2CopyNumbers.R with no options or with --help to get an overview of available command line options.

3- CopyNumber to Plots (Called)
    /path/to/your/r-installation/Rscript QDNAseq\_BAM2CopyNumbers.R with no options or with --help to get an overview of available command line options.

CGHTest:
    /path/to/your/r-installation/Rscript QDNAseq\QDNAseq_FrequencyCGHregionsCGHtest_logrank_CMDL.R

Clustering
Open Rstudio and run the script  wecca-hclust-bootstrap.R from there from there


GeneBreak
Open Rstudio and run the script run_GeneBreak-convert2Rubic.R from there


4. Input files
Alignment
You need fastq(fq).gz files

Bam to CopyNumber
You will need the bam and the bai files created in the previous step

CopyNumber to Plots (called)
You will need the rds file generated in the previous step

CGHTest
You will need the rds file with the call generated in the previous step
Excell file with the sample names 

Clustering

GeneBreak


5. How to run 
Alignment 
You need to run the script from the folder/directory where your fastq files are store

Bam to CopyNumber
you need to add in the CLI (space separated)

Program(bin) = DATA/share/pipelines/apps/R-3.4.4/bin/Rscript
Script(path)= /DATA/share/pipelines/PIPELINES/QDNAseqFlow/QDNAseq_BAM2CopyNumbers_parallel.R
sample_bamfiles = Bam/ (path to bam files)# if you are running the analysis form the same directory where the bam are use a '.'
proejectName = project_name 
binSize = 30, 15, 100 (choose one)
batchName = name
Example:
"/DATA/share/pipelines/apps/R-3.4.4/bin/Rscript /DATA/share/pipelines/PIPELINES/QDNAseqFlow/QDNAseq_BAM2CopyNumbers_parallel.R Bam/ IntEnd_batch1 30 batch1"


CopyNumber to Plots (called), Local machine(laptop)
you need to add in the CLI (space separated)

Program = /usr/local/bin/Rscript
Script = /Users/l.gonzalez/PIPELINES_gpMeijer/QDNAseqFlow/QDNAseqFlow-versionMay2019/QDNAseq_CopyNumbers2PlotsCallsSegmentFileInclDewave_incl.R
Inputfile = batch1/batch1-binSize-30/batch1-copyNumbers-30kb-bins.rds
proejectName = out30kb/ # you need to create this folder in advance
batchName = bt1-30kb #different from the previous step
dewaving = yes  #yes/no (options to run) only for 15, 30, 100 bin size
Example "/usr/local/bin/Rscript /Users/l.gonzalez/PIPELINES_gpMeijer/QDNAseqFlow/QDNAseqFlow-versionMay2019/QDNAseq_CopyNumbers2PlotsCallsSegmentFileInclDewave_incl.R  batch1/batch1-binSize-30/batch1-copyNumbers-30kb-bins.rds out30kb/ bt1-30kb yes"



CGHTest  Local machine(laptop)

Program = /usr/local/bin/Rscript
Script =
 /Users/l.gonzalez/PIPELINES_gpMeijer/QDNAseqFlow/QDNAseqFlow-versionMay2019/QDNAseq_FrequencyCGHregionsCGHtestCMDL_X2_IgnoreUnknownCat.R
Input_file = out3step/test2.dewave/30kb-bins/test2.dewave-copyNumbersCalled-30kb-bins.rds
Output_dir = outTest/ 
Excel_path(samples) = test2.dewave_stats_30kb-bins.xlsx 
excel_Tap(name) = stats
excell_column(name) = comparison_CvsN+A

Example:
/usr/local/bin/Rscript /Users/l.gonzalez/PIPELINES_gpMeijer/QDNAseqFlow/QDNAseqFlow-versionMay2019/QDNAseq_FrequencyCGHregionsCGHtestCMDL_X2_IgnoreUnknownCat.R out3step/test2.dewave/30kb-bins/test2.dewave-copyNumbersCalled-30kb-bins.rds outTest/ test2.dewave_stats_30kb-bins.xlsx stats comparison_CvsN+A

Clustering,  Local machine(laptop)
Run in Rstudio
Add the path the the rds calls file 
Add a path for the working directory, where the files will be save 

GeneBreak,  Local machine(laptop)
Run in Rstudio
Add or chage the paths for the:
rds calls file
soomthpercertage tab file
Call file (tab)




6. Out Files
From the alignment, You will get:
sorted bams
mark duplicate bam and bai (after the sam indexing)


From Bam to CopyNumber, you will get 3 files:
projectdir-copyNumbers-30kb-bins.igv.bz2 (raw copy numbers in igv format)
projectdir-copyNumbers-30kb-bins.tsv.bz2 (raw copy numbers as table)
projectdir-copyNumbers-30kb-bins.rds (R data structure with read counts, needed for next 


From the CopyNumber to Plots(called), you will get: This is the *main* output from the pipeline

Directories(dir)
SmoothPlot (dir):  copy number plots
SegmentedPlot: segmented plots
CalledPlot (dir): segmented+called plots
stats(dir): folder with various stats, most info per chromosome arm, gains or losses only. Aberrations are counted or the percentage of the chromosome arm that has gains or losses is given.
stats/projectname-statistics-30kb-bins.xlsx: stats of all plots are listed. Coloring (red or blue) or values if they are outliers according to the 1.5 IQR rule. Empirical criterion to exclude samples with 30kb bins is diff(var) > 0.05.

Files
projectnam30kbbins.point01percentsmoothing.CGHregions.tab 
#output file of CGHregions with ‘information loss = 0.01%
projectname-copyNumbersSegmented-30kb-bins.tab 
#log2 segmented read counts
projectname-copyNumbersCalled-30kb-bins.tab 
#the calls
projectname-copyNumbersCalled-30kb-bins.rds 
#R-data structure, needed in next steps
projectname-frequencyPlot-30kb-bins.pdf 
#frequency plot

CGHTest
Frequency plots for the two populations that you want to compare (pdf format)
Stats from the comparison of the groups


Clustering
heatmap (pdf format)
Bed file 

GeneBreak
3 tsv files to use in RUBIC
excel file with the list of most significatly broken, not by chance, genes 
Right now gene break is running gene centered




