% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Input_data.R
\name{Input_Sclist}
\alias{Input_Sclist}
\title{Title Input transposon data as Seurat-list Format}
\usage{
Input_Sclist(meta = metadata, inputfiledir = "data/", cellmin = 3, fmin = 20)
}
\arguments{
\item{meta}{metadata that needs input}

\item{inputfiledir}{transposon data dir}

\item{cellmin}{the min cell you want to keep for}

\item{fmin}{the min feature you want to keep for}
}
\value{
Seurat-list
}
\description{
defaultly,it will automatic detect Data in data/ directory.
The metadata type is like "CRR703  H03P" generate from sample_info.tsv. This tsv-file must be in your workdirectory.
cellmin is the min cell you want to keep. below this value, the observation will be dropped.
fmin is the min feature you want to keep .below this value, the variable will be dropped.
sclist<-Input_Sclist() #it will automatic detect Data in data/ directory.
vclist<-Input_Sclist(data=metadata,cellmin =3,fmin=20) #you can use your own metadata dataset and specific cellmin or fmin.
}
