#!/usr/bin/env Rscript

#install.packages("jpeg")
library(jpeg)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
 stop("
       usage:   Rscript --vanilla Dotplot_n_cov.R <dotplot> <covplot> <picname.jpeg>
       dotplot: e.g. S1phunph_vs_Scer_noMt.ANI.jpeg
       covplot: e.g. CPOPY0_AllSacch2_coverage_per1kb_ScSk.jpeg
       picname: e.g. S1phunph_CPOPY0_dotplot_n_cov.jpeg",
      call.=FALSE)
}

cmd = paste(commandArgs(), collapse = " ")
write(cmd, stderr())

path1 <- ""
path2 <- ""
path3 <- ""
  
file1 <- paste(path1, as.character(args[1]), sep="") # "../Data/Assembly/fltrdSortedPhased/S1phunph_vs_Scer_noMt.ANI.jpeg"
file2 <- paste(path2, as.character(args[2]), sep="") # "../Data/AllSaccMapping/CPOPY0_AllSacch2_coverage_per1kb_ScSk.jpeg"
file3 <- paste(path3, as.character(args[3]), sep="") # "Dotplot_n_cov.jpeg"

covstr  <- strsplit(as.character(args[2]), "_")
covtype <- covstr[[1]][1]
print(covtype)

dotplot <- readJPEG(file1, TRUE)
covplot <- readJPEG(file2, TRUE)

jpeg(file3, width=297, height=210, units="mm", quality=80, res=600);
par(mar=c(0,0,0,0))

plot(x=c(0,6), y=c(0,2), xlim=c(0,6), ylim=c(0,2), xlab="", ylab="", axes=F, pch=3, col="white")

rasterImage(covplot, 6, 0, 10.6, 0.55, angle=90)
rasterImage(dotplot, 0, 0.04, 5, 2.01)

text(x=5.45, y=2, labels=covtype, adj=0.5)

dev.off()