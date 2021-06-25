#!/usr/bin/Rscript --vanilla

args<-commandArgs(trailingOnly=TRUE)
library(EasyQC)
EasyQC(args[1])
q()

