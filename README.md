# QDNAseqFlow 
## TGO group, NKI

### Introduction

QDNAseq is an R package for shaloow Whole Genome Sequency data (sWGS). In this workflow we put together one previous step which include bwa alignment aln, sort files and mark duplicates, then indexing the bam files to handle them to QDNAseq.
This version is an updated version from the previos created by Christian Rhaus. This version has some improvement; mapping and QDNAseq analysis are based in GRCh38 (hg38) genomic assembly, the main workflow is being reduce to two steps (mapping and calling)

The main QDNAseqFlow include two steps: 
- Mapping (sorting and marlduplicates)
- Bam to Called Plots

From there we can use the output of the called plots for: 

- CGHTest to compare two groups

- Clustering usingn WECCA
Hierarchical clustering (WECCA & hclust)
using called copy number data of regions (WECCA) and all bins (hclust). 100 bootstrapping rounds .
Calculation of Strict and Majority consensus clusterings (trees) in Dendroscope

- Recurrent broken genes using GeneBreak
detection of recurrent DNA copy number aberration-associated chromosomal breakpoints within genes

### Installation
New version is from 2020
Go to https://github.com/NKI-TGO/QDNAseqFlow_v2 and click the button "Clone or Download" --> Download zip.
Unzip the zip file in your Home directory, 'My Documents' or where you are allowed to install programs.

**Required R-version etc.**

In the update version you arw allow to use the latest R version, using BiocManager to install the neccesary packages from Bioconductor or Cran. Be waare that you need permision to install packages

Java: You need to hava java installed. Make sure it is in your path. On Windows, this is done like here: https://confluence.atlassian.com/doc/setting-the-java_home-variable-in-windows-8895.html


### Usage

This workflow consists of 2 main scripts, that need to be run consecutively: The workflow is form end to end but you can also run the steps in separate way

1. Aligment, sorting and mark duplicates 
2. QDNAseq_hg38_Bam2CpNCalls.Rmd

Those stepts are not interctive, yo need to take care of few features

- path to the hg38 reference, if you do not have any we recommend to use: [Homo_sapiens.GRCh38.dna.primary_assembly](https://www.ensembl.org/Homo_sapiens/Info/Index)
[direct downloa](ftp://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
- path to your workinf directory, where you want to save the outputs
- Path to the rds files with the hg38 bin annotations 
- Indicate the path to your bam/bai files
- select bin size (15, 30, 100, 1000)
- dewaving not or yes and  if yes, select the bin size accordingly to the the bin size selected
- path to the extra scripts to generate the statsQC, 



#### Only one way to run the scripts for the moment, From the command line:

__Bash__
1. Alignment: bash path/to/script/1_bwaAlnv2_lgl.sh
Need to make the bam index just right after

__R scripts__

2. QDNAseq_hg38_Bam2CpNCalls.Rmd
This is a R markdown script, the best way to run it, and modify is in Rstudio. Once you modifly the path and select you bin size you can use the option Knit to run the entire code and generate an html(pdf etc) to keep a record of the settings that you ran for each study. 

CGHTest:
    /path/to/your/r-installation/Rscript QDNAseq\QDNAseq_FrequencyCGHregionsCGHtest_logrank_CMDL.R

Clustering
Open Rstudio and run the script  wecca-hclust-bootstrap.R from there from there


GeneBreak
Open Rstudio and run the script run_GeneBreak-convert2Rubic.Rmd from there


#### Input files
Alignment
You need fastq(fq).gz files

Bam to CopyNumber to Plots (called)
You will need the bam and the bai files created in the previous step

CGHTest
You will need the rds file calls generated in the previous step
Excell file with the sample names 

Clustering
You will need the rds file calls generated in the previous step

GeneBreak
You will need the rds file calls generated in the previous step


#### How to run 
___Alignment___
You need to run the script from the folder/directory where your fastq files are store

___Bam to CopyNumber Called Plots ___
you need to open and edit the script and run from Rstudio or CLI

Script(path, working directory)= /path/to/the/directory/to/save/outputs
sample_bamfiles = Bam/ (path to bam files)
proejectName = project_name 
binSize = 30, 15, 100, 1000 (choose one)
dewaving = #yes/no if no: comment out dewaving option, if yes: select bin size (15, 30, 100, 1000)


___CGHTest___  Local machine(laptop)

Program = /usr/local/bin/Rscript
Script =
 /Users/l.gonzalez/PIPELINES_gpMeijer/QDNAseqFlow/QDNAseqFlow-versionMay2019/QDNAseq_FrequencyCGHregionsCGHtestCMDL_X2_IgnoreUnknownCat.R
Input_file = copyNumbersCalled-30kb-bins.rds
Output_dir = outTest/ 
Excel_path(samples) = test2_30kb-bins.xlsx 
excel_Tap(name) = stats
excell_column(name) = comparison_CvsN+A

Example:
/usr/local/bin/Rscript /path/QDNAseq_FrequencyCGHregionsCGHtestCMDL_X2_IgnoreUnknownCat.R opath/test2.dewave-copyNumbersCalled-30kb-bins.rds outTest/ test2_30kb-bins.xlsx stats comparison_CvsN+A

___Clustering___,  Local machine(laptop)
Run in Rstudio
Add the path the the rds calls file 
Add a path for the working directory, where the files will be save 

___GeneBreak___,  Local machine(laptop)
Run in Rstudio
Add or chage the paths for the:
rds calls file
soomthpercertage tab file
Call file (tab)



#### Output Files
___From the alignment, You will get:___
sorted bams
mark duplicate bam and bai (after the bam indexing)


___From Bam to CopyNumber CopyNumber to Plots(called), you will get:___

Directories(dir)
_ CalledPlot (dir): segmented+called plots
_ stats(dir): folder with various stats, most info per chromosome arm, gains or losses only. Aberrations are counted or the percentage of the chromosome arm that has gains or losses is given.
stats/projectname-statistics-30kb-bins.xlsx: stats of all plots are listed. Coloring (red or blue) or values if they are outliers according to the 1.5 IQR rule. Empirical criterion to exclude samples with 30kb bins is diff(var) > 0.05.


Files
hg38_raw-30kb.rds 
hg38_copyNumbers-30kb.rds
hg38_copyNumbers-30kb.igv
hg38_segmented-30kb.rds
hg38_segments-30kb.igv
hg38_Original_kb.rds
hg38_calls-30kb.igv
hg38_called-30kb.rds
#frequency plot
projectname-frequencyPlot-30kb-bins.pdf 


___CGHTest___
Frequency plots for the two populations that you want to compare (pdf format)
Stats from the comparison of the groups


___Clustering___
heatmap (pdf format)
Bed file 

___GeneBreak___
3 tsv files to use in RUBIC
excel file with the list of most significatly broken, not by chance, genes 
Right now gene break is running gene centered



## Previous Versions of QDNAseqFlow

__This information is deprecated but the repository still exist and you can clone or check the information as it is for that version of the pipeline__
--Please have a look at [QDNAseqFlow-Abstract.pdf](https://github.com/NKI-Pathology/QDNAseqFlow/blob/master/QDNAseqFlow-Abstract.pdf) and [QDNAseqFlow\_poster\_ISMB.pdf](https://github.com/NKI-Pathology/QDNAseqFlow/blob/master/QDNAseqFlow_poster_ISMB.pdf) for an introduction.
To better understand how QDNAseq, the main component of this pipeline, works: [Christian-Rausch-QDNAseq-talk\_BioconductorDec2015.pdf](https://github.com/NKI-Pathology/QDNAseqFlow/blob/master/Christian-Rausch-QDNAseq-talk_BioconductorDec2015.pdf)--

This version 2017 is from 
Go to https://github.com/NKI-Pathology/QDNAseqFlow and click button "Clone or Download" --> Download zip.
Unzip the zip file in your Home directory, 'My Documents' or where you are allowed to install programs.


In the previous version:
R: We have made most tests with R 3.4.1. QDNAseqFlow might work with older R versions but the risk is that some required libraries might not be available.
_The step 2 (Bam to CopyNumber) is still running in the R3.4 version (server)_
