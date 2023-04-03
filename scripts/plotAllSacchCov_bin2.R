#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
 stop("\n
  usage:    Rscript --vanilla plotAllSacchCov_bin2.R <tag> <path> <winlen>
  example:  Rscript --vanilla plotAllSacchCov_bin2.R AllSacch2 CPOPY0_AllSacch2 1kb\n\n", 
  call.=FALSE)
}
cmd = paste(commandArgs(), collapse = " ")
write(cmd, stderr())

#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("DNAcopy")
library(DNAcopy)

tag = as.character(args[1])
prefix = as.character(args[2])
winlen = as.character(args[3])
suffix = paste("_coverage_per",as.character(args[3]),sep="")

plot_title = paste(prefix, "_coverage_per", winlen, "_ScSe", ".jpeg", sep="")
plot_mtit  = paste(prefix, "_coverage_per", winlen, "_ScSe", sep="")

tag = paste(prefix,prefix,sep = "/")
prefix = paste(prefix,"/",sep = "")

lngstLen = 12500000
logORbin = "bin"
mYlim <- 4

infile <- paste(tag,".vcf",sep = "")
cntgLines <- grep('##contig',readLines(infile))
tab <- read.table(infile, comment.char = "~", sep = ",",
                  skip = cntgLines[1] - 1, nrows = length(cntgLines))
Loci <- gsub(pattern = "##contig=<ID=", replacement = "", x = tab$V1)
Loci_Len <- gsub(pattern = "length=", replacement = "", x = tab$V2)
Loci_Len <- gsub(pattern = ">", replacement = "", x = Loci_Len)
lenTab <- data.frame(Loci, Loci_Len)

Loci_IDs <- seq(1,length(cntgLines))
Loci_Labs <- paste(format(as.numeric(Loci_Len) / 1000000, digits = 1, nsmall = 2),"Mb",sep="")

strains = strsplit(Loci,"_")
strain.df <- data.frame(matrix(unlist(strains), nrow=length(strains), byrow=TRUE))
strains = as.vector(levels(factor(strain.df$X1)))

# force order
strains = c(strains[2],strains[3],strains[7],strains[1],strains[5],strains[4],strains[6],strains[8])

strain_coverage_tab <- paste(tag, suffix, sep="")
full_cov_table <- read.table(strain_coverage_tab, sep = "\t", header = T, stringsAsFactors = F,
                             colClasses = c("character", "integer", "numeric", "numeric", "numeric"))
full_cov_table$Observed <- full_cov_table$Observed + 1
full_cov_table$Division <- full_cov_table$Observed / full_cov_table$average
full_cov_table$Log2division <- log2(full_cov_table$Division)

jpeg(file=plot_title, width=333, height=100, units="mm", quality=80, res=600);
par(mfrow=c(2,1), mar=c(2,1,1,0), oma=c(0,0,2,2))

#for (j in 1:length(strains)){
for (j in 1:2){
 #j=2
 int_cov_table <- subset(full_cov_table, grepl(strains[j], Loc))
 # force Scer len as reference
 #int_lenTab <- subset(lenTab, grepl(strains[1], Loci))
 int_lenTab <- subset(lenTab, grepl(strains[j], Loci))
 # replace ID
 int_lenTab$Loci <- gsub("Scer", strains[j], int_lenTab$Loci)
 int_lenTab$cumLen <- cumsum(int_lenTab$Loci_Len)
 int_lenTab$midCumLen <- int_lenTab$cumLen - (as.numeric(int_lenTab$Loci_Len)/rep(2,times=nrow(int_lenTab)))
 
 for (i_cnt in 1:nrow(int_cov_table)){
   for (j_cnt in 2:length(int_lenTab$Loci)){
     k_cnt=j_cnt-1;
     if (grepl(pattern = int_lenTab$Loci[j_cnt], x = int_cov_table[i_cnt,"Loc"], fixed = T) == T){
       int_cov_table[i_cnt,"Pos"] = int_cov_table[i_cnt,"Pos"] + int_lenTab[k_cnt,"cumLen"]
     }
   }
 }
 
 if (logORbin == "log") {
   CNA.object <- CNA(int_cov_table$Log2division, int_cov_table$Loc, int_cov_table$Pos, 
                     sampleid = strains[j], data.type = "logratio", presorted = T)
   smoothed.CNA.object <- smooth.CNA(CNA.object)
 }
 if (logORbin == "bin") {
   CNA.object <- CNA(int_cov_table$Log2division, int_cov_table$Loc, int_cov_table$Pos, 
                     sampleid = strains[j], data.type = "binary", presorted = T)
   smoothed.CNA.object <- CNA.object
 }
 segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
 cov_table <- data.frame(segment.smoothed.CNA.object$data, stringsAsFactors = F, check.names = F)
 cov_table_bup <- cov_table
 #reverse values
 cov_table[,strains[j]] <- ifelse(cov_table[,strains[j]] < 0,
                                  gsub("-","",cov_table[,strains[j]],fixed=TRUE),
                                  paste0("-",cov_table[,strains[j]]))
 cov_table[,strains[j]] <- as.numeric(cov_table[,strains[j]])

 seg_cov_table <- segment.smoothed.CNA.object$output
 seg_cov_table_bup <- seg_cov_table
 #reverse values
 seg_cov_table[,"seg.mean"] <- ifelse(seg_cov_table[,"seg.mean"] < 0,
                                      gsub("-","",seg_cov_table[,"seg.mean"],fixed=TRUE),
                                      paste0("-",seg_cov_table[,"seg.mean"]))
 seg_cov_table[,"seg.mean"] <- as.numeric(seg_cov_table[,"seg.mean"])
 
 plot(cov_table$maploc, cov_table[,strains[j]], type = "p", pch = 20, cex = 0.5, col = "gray50",
      #reverse yaxis
      ylim=c(-mYlim,0), xlim=c(0,lngstLen), axes=F, main="", ylab = "", xlab = "")

 axis(side=1, at=c(0,tail(int_lenTab$cumLen, n=1)), labels=FALSE, padj = -0.75)
 axis(side=1, at=int_lenTab$cumLen, labels=FALSE, padj = -0.75)
 if (j==2){
  axis(side=1, at=int_lenTab$midCumLen, labels=FALSE, padj=-0.75, tick=FALSE)
  text(x=int_lenTab$midCumLen, y=-5, labels=seq(1,nrow(int_lenTab)), srt=-90, xpd=TRUE)
 }
 #reverse y axis
 axis(side=2, line=-2, at=seq(-mYlim,0), labels=FALSE, las=2)
 text(x=-300000, y=seq(-mYlim,0), labels=c("4","3","2","1","0"), srt=-90)

 text(x=-600000, y=-2, labels = strains[j], srt=-90, xpd=TRUE)
 if (j == 1) { title(main = plot_mtit, outer = T) }
 
 for (k in 1:nrow(seg_cov_table)){
  lines(x=c(seg_cov_table[k,"loc.start"],seg_cov_table[k,"loc.end"]),
        y=c(seg_cov_table[k,"seg.mean"],seg_cov_table[k,"seg.mean"]), col="red", lwd = 3)
 }
 
 abline(v=int_lenTab$cumLen, col = "gray", lty=1, lwd = 0.5)
 lines(x=c(0,lngstLen), y=c(0,0), col = "gray", lty=1, lwd = 0.5)
 lines(x=c(0,lngstLen), y=c(-1,-1), col = "gray", lty=1, lwd = 0.5)
 lines(x=c(0,lngstLen), y=c(-2,-2), col = "gray", lty=1, lwd = 0.5)
 lines(x=c(0,lngstLen), y=c(-3,-3), col = "gray", lty=1, lwd = 0.5)
 
 if (j==3){
   featable <- read.table("../../references/genomes/Scer_noMt.feat.coords", sep = " ", header = T, stringsAsFactors = F,
                          colClasses = c("numeric", "character", "numeric", "numeric", "numeric"))
   featable$cumLen <- cumsum(featable$Len)
   tmpVec <- c(0,featable[,"cumLen"])
   featable$corCENcoords <- tmpVec[-length(tmpVec)] + featable$CENcoords
   abline(v=featable$corCENcoords, col="gray", lty=5)
 }
 
}
dev.off()
