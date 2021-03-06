
# ViruSeurat
[![GitHub stars](https://img.shields.io/github/stars/cherryamme/ViruSeurat?color=red&logo=Adafruit)](https://github.com/cherryamme/ViruSeurat/stargazers)  ![GitHub repo size](https://img.shields.io/github/repo-size/cherryamme/ViruSeurat?color=yellow&label=Project%20Size&logo=Blueprint)   ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/cherryamme/ViruSeurat?logo=R)   ![GitHub last commit](https://img.shields.io/github/last-commit/cherryamme/ViruSeurat) 
<!-- badges: start -->

![Jesse](./dev_fig/ViruSeurat_logo.png)

<!-- badges: end -->
**Viruseurat** is an R software used to check the transposon expression and virus invasion in single cell data. 

<!-- vscode-markdown-toc -->
* 1. [Introduction](#Introduction)
* 2. [Installation](#Installation)
* 3. [Example](#Example)
	* 3.1. [How to set your own data](#Howtosetyourowndata)
	* 3.2. [How to Use test dat](#HowtoUsetestdat)
	* 3.3. [Main Command](#MainCommand)
* 4. [Release Notes](#ReleaseNotes)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

##  1. <a name='Introduction'></a>Introduction

It is based on the R package-Seurat, and uses virus expression data and transposon expression data. Virus expression data is generated using STAR comparison virus database, and transposon expression data is generated using scTE.

ViruSeurat solves discovering relationship between transposon and Virus in singlecell-level, providing a new view to check virus distribution in human cell. Maybe useful for bioinformatic researcher . 


It can add some tags that comtribute to build some group among these.

It is my first R-package, so it may has lots of bugs that I haven't meet. If you find them , please contact me with email cherryamme@163.com improving this software.

Hope this package will help you.

##  2. <a name='Installation'></a>Installation

You can install the development version of ViruSeurat from [GitHub](https://github.com/) with:

``` r
 if (!require("devtools")){
  install.package("devtools")}
devtools::install_github("cherryamme/ViruSeurat")
```

##  3. <a name='Example'></a>Example
###  3.1. <a name='Howtosetyourowndata'></a>How to set your own data
You must set the working directory in where contains your metadata(sample_info.tsv) transposon/Virus expression tsv/csv data/ folder.

**metadata** must save as sample_info.tsv, split by `\t`,it has two column cointain \<name\>`\t`\<metainfo\>.

**transposon data** must save as CRR703.csv CRR704.csv, a format by \<name\>.csv split by `,`,generated by `scTE`.

**Virus data** must save as CRR703_Exp.tsv CRR704.csv, a format by \<name\>_Exp.tsv split by `\t`,generated by `Viru-track`.

directory is like behind style
```txt
working directory:

.
????????? sample_info.tsv
????????? data
    ????????? CRR703.csv
    ????????? CRR703_Exp.tsv
    ????????? CRR704.csv
    ????????? CRR704_Exp.tsv

```


###  3.2. <a name='HowtoUsetestdat'></a>How to Use test dat

I write a Test data in package source directory, you can simply type `system.file(package = "ViruSeurat")` to get the original directory path. Then you will see a Testdata folder that contains Testdata and metadata.

If you don't want to use original file to build Test data. You can Input Testdata easily by command `ViruSeurat::Input_Testdata()`. It will automatically generate a Test example data.

There are some basic examples which shows you how to solve a common problem:

###  3.3. <a name='MainCommand'></a>Main Command

``` r
library(ViruSeurat)

#Input your own data by below two command, you need to check the directory is correct.
metadata <- Input_meta()
sclist <- Input_Sclist(meta = metadata, inputfiledir = "data/", cellmin = 3, fmin = 20)
vclist <- Input_Vclist(x = viruscount, inputfiledir = "data/", viralname = "NC_003977")

#Check your VirusCounts distribution
Plot_VirusCount(x = vclist, meta = metadata)

#Check your VirusCounts in your sclist or not
CompareVirusCB(vclist, sct = sclist) #check a list
CompareVirusCB(vclist[[1]], sct = sclist) #check a single vector in your data

#Input Virus Counts to Seurat-list
svclist <- VirusinSclist(x = vclist, sct = sclist, meta = metadata)

#Normalize and Integrate Svclist, generate a scRNA : an integrated Seurat object
scRNA <- Normalize_IntegrateSclist(x = svclist, fmethod = "vst", features = 2000)

#Plot_scatter in integrated Seurat object, find corresponding between virus and cells
Plot_Scatter(x = scRNA)

#Calculate Cluster information
ccRNA <- Cluster_ScRNA(x = scRNA, res = c(seq(0.1, 1.6, 0.2)), npcs = 30, dims = 1:30)

#If you are not satisfied about your cluster resolution. Do behind
#ReCalculate Cluster information
ccRNA <- Recluster_CcRNA(x = ccRNA, res = c(seq(0.1, 1.6, 0.2)))

#select a resolution and plot simple information
Plot_Cluster(x = ccRNA, res = 0.5)


```


##  4. <a name='ReleaseNotes'></a>Release Notes
you can check [Release Notes](./News.md).

