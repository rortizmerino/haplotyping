#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
 stop("\n
 usage:    Rscript --vanilla plotVCF.SNVs_n_cov_bin_32_i.R <tag> <path> <length>
 example:  Rscript --vanilla plotVCF.SNVs_n_cov_bin_32_i.R Scen_Seub CPOPY0_Scen_Seub 1kb
 optional: <first scaff> <last scaff> <lngstlen> <log|bin>
 example:  Rscript --vanilla plotVCF.SNVs_n_cov_bin.R S1_np_ph_v0 CPOPY0_S1_np_ph_v0 1kb 1 100 1600000\n\n", 
 call.=FALSE)
}

#install.packages("tidyr")
library(dplyr)
library(tidyr)
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("DNAcopy")
library(DNAcopy)

tag = as.character(args[1])
prefix = as.character(args[2])
winlen = as.character(args[3])
suffix = paste("_coverage_per",as.character(args[3]),sep="")

ori_tag = tag
tag = paste(prefix,prefix,sep = "/")
prefix = paste(prefix,"/",sep = "")

frst = as.numeric(args[4])
last = as.numeric(args[5])
lngstLen = as.numeric(args[6])
logORbin = "bin"
#logORbin = toString(args[7])
write(logORbin, stderr())
cmd = paste(commandArgs(), collapse = " ")
write(cmd, stderr())

# define ylims in covplots
#if (!is.null(logORbin) || logORbin == "log") {
#  logORbin = "log"; mYlim <- 2
#} else if (logORbin == "bin") { logORbin = "bin"; mYlim <- 4 }
#msg = paste("Cov plot version chosen:", logORbin, sep = " ")
#write(msg, stderr())
mYlim <- 2
#mYlim <- 4

infile <- paste(tag,".SNVs.vcf",sep = "")
cntgLines <- grep('##contig',readLines(infile))
tab <- read.table(infile, comment.char = "~", sep = ",",
                  skip = cntgLines[1] - 1, nrows = length(cntgLines))
Loci <- gsub(pattern = "##contig=<ID=", replacement = "", x = tab$V1)
Loci_Len <- gsub(pattern = "length=", replacement = "", x = tab$V2)
Loci_Len <- gsub(pattern = ">", replacement = "", x = Loci_Len)
Loci_IDs <- seq(1,length(cntgLines))
Loci_Labs <- paste(format(as.numeric(Loci_Len) / 1000000, digits = 1, nsmall = 2),"Mb",sep="")

if (!is.na(frst) && !is.na(last)) {
} else { frst = 1; last = length(Loci) }
diff_check = last - frst
msg = paste(c(last," - ",frst," = ",diff_check), collapse = " ")
write(msg, stderr())
if( diff_check > 100){ 
  write("won't plot more than 100 sequences, choose range accordingly", stderr()); stop() 
}

# might help select nice longest length to use for all SNV and covplots
#max(as.numeric(Loci_Len))
#roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
#  if(length(x) != 1) stop("'x' must be of length 1")
#  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
#}
#lngstLen <- roundUpNice(max(as.numeric(Loci_Len)))
#max(sort(Loci_Len[5]))

if (is.na(lngstLen)) { lngstLen <- 1600000 }
msg = paste("Longest sequence lenght seen:", max(sort(Loci_Len[5])), sep = " ")
write(msg, stderr())
msg = paste("Longest sequence lenght chosen:", lngstLen, sep = " ")
write(msg, stderr())

infile <- paste(tag,".SNVs.txt",sep = "")
full_table <- read.table(infile, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
full_table$REF <- NULL; full_table$ALT <- NULL;  
colnames(full_table)[1] <- "Loc"; colnames(full_table)[2] <- "Pos"

strains <- colnames(full_table[c(3:ncol(full_table))])
full_table <- full_table[rowSums(is.na(full_table)) != length(strains),]
strain_Labs <- strains

# force nicer labels; lefthand side of histograms
nice_lbls = c("Scen I", "Scen II", "Scen III", "Scen IV", "Scen V", "Scen VI", "Scen VII", "Scen VIII",
              "Scen IX", "Scen X", "Scen XI", "Scen XII", "Scen XIII", "Scen XIV", "Scen XV", "Scen XVI",
              "Seub I", "Seub II", "Seub III", "Seub IV", "Seub V", "Seub VI", "Seub VII", "Seub VIII",
              "Seub IX", "Seub X", "Seub XI", "Seub XII", "Seub XIII", "Seub XIV", "Seub XV", "Seub XVI")
# force order to show interleaved chromosome order
strain_order = c(1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31,16,32)

for (j in 1:length(strains)){
 #j = 1
 plt_rng <- paste(frst,last,sep = "-")
 ext <- paste("SNVs_n_cov_",logORbin,".i.jpeg",sep="")
 #ext <- "SNVs_n_cov_log.jpeg" 
 plot_title <- paste(strains[j], plt_rng, ext, sep=".")
 
 write(plot_title, stderr())
 
 #jpeg(file=plot_title, width=333, height=999, units="mm", quality=80, res=600);
 jpeg(file=plot_title, width=333, height=1000, units="mm", quality=80, res=600);
 #par(mfrow=c(length(Loci),3), mar=c(0.5,2,1.5,1), oma=c(5,3.5,2,1))
 par(mfrow=c(32,3), mar=c(0.5,2,1.5,1), oma=c(5,9,2,1))

 #strain_coverage_tab <- paste(prefix, strains[j], suffix, sep="")
 strain_coverage_tab <- paste(prefix, ori_tag, suffix, sep="")
 full_cov_table <- read.table(strain_coverage_tab, sep = "\t", header = T, stringsAsFactors = F,
                          colClasses = c("character", "integer", "numeric", "numeric", "numeric"))
 full_cov_table$Observed <- full_cov_table$Observed + 1
 full_cov_table$Division <- full_cov_table$Observed / full_cov_table$average
 full_cov_table$Log2division <- log2(full_cov_table$Division)
 
 if (logORbin == "log") {
  CNA.object <- CNA(full_cov_table$Log2division, full_cov_table$Loc, full_cov_table$Pos, 
                    sampleid = strains[j], data.type = "logratio", presorted = T)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
 }
 if (logORbin == "bin") {
  CNA.object <- CNA(full_cov_table$Log2division, full_cov_table$Loc, full_cov_table$Pos, 
                    sampleid = strains[j], data.type = "binary", presorted = T)
  smoothed.CNA.object <- CNA.object
 }
 segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
 cov_table <- data.frame(segment.smoothed.CNA.object$data, stringsAsFactors = F, check.names = F)
 
 if (winlen == "1kb") { cov_table$maploc_s <- cov_table$maploc - 1000 }
 else { cov_table$maploc_s <- cov_table$maploc - 10000 }
 
 seg_cov_table <- segment.smoothed.CNA.object$output

 colnames(cov_table)[3] <- gsub("\\.", "-",colnames(cov_table)[3])
 seg_cov_table$ID <- gsub("\\.", "-", seg_cov_table$ID)
 
 #if (exists(frst) && exists(last)) {}
 #else {frst = 1; last = length(Loci)}
 
 for (k in frst:last){
 #for (k in length(strain_order)){
  i = strain_order[k]
  #write(Loci[i], stderr())
  strt <- 0; stop <- Loci_Len[i];

  int_table <- subset(full_table, grepl(Loci[i], Loc))
  int_table <- int_table[ which( int_table$Pos >= strt & int_table$Pos < stop), ]
  sub_table <- data.frame(int_table[,c("Pos",strains[j])], stringsAsFactors =F, check.names = F)
  colnames(sub_table) <- c("Pos", "strain")

  if (length(sub_table$Pos) > 1){
   row.names(sub_table) <- c(1:length(sub_table$Pos))

   sub_table <- sub_table %>% separate(strain, c("A", "B"), ",")
   sub_table <- data.frame(lapply(sub_table, function(x) {gsub("NA", NA, x)}), stringsAsFactors =F)
   sub_table0 <- data.frame(c(as.character(sub_table$Pos), as.character(sub_table$Pos)), stringsAsFactors =F)
   sub_table <- cbind(sub_table0, 
                      data.frame(c(as.character(sub_table$A), as.character(sub_table$B)), stringsAsFactors =F))

   colnames(sub_table) <- c("Pos", strains[j])
   sub_table <- sub_table[order(as.numeric(sub_table$Pos)),]
   sub_table <- sub_table[complete.cases(sub_table), ]
   sub_table <- sapply( sub_table, as.numeric )
   sub_table <- data.frame(sub_table, stringsAsFactors =F, check.names = F)
  }
  
# HISTOGRAMS #
  if (exists("sub_table") == TRUE && dim(sub_table)[1] > 1){
   #write("H sub_table check", stderr())
   H<-hist(sub_table[,strains[j]], breaks = 33,  include.lowest = T, plot = F)#, freq = T)
   sub_table2 <- sub_table[ which( sub_table[,strains[j]] < 0.85), ]
  
   mtit<-paste("TAF < 0.85:",as.character(length(which( sub_table[,strains[j]] < 0.85))),"      ",
               "TAF >= 0.85:",as.character(length(which( sub_table[,strains[j]] >= 0.85))),sep=" ")
   VarsNo <- length(which( sub_table[,strains[j]] < 0.85))
  }
  if (exists("sub_table2") == TRUE && dim(sub_table2)[1] != 0 ){
   #write("H2 sub_table2 check", stderr())
   H2<-hist(sub_table2[,strains[j]], breaks = 33, include.lowest = T,
            ylab = "", xlab = "", 
            xlim=c(0.15,0.85), 
            main = "", axes = F,  freq = T, border = "gray50",
            col=ifelse((VarsNo < (length(sub_table[,strains[j]]) * 0.1)), "white", "gray50") )
            #col=ifelse((VarsNo < (length(sub_table[,strains[j]]) * 0.1) | VarsNo < 1000), "white", "gray50") )
   yaxlim <- max(H2$counts)
   abline(v=0.75, col = "DarkGreen", lty=2, lwd = 3); abline(v=0.25, col = "DarkGreen", lty=2, lwd = 3)
   abline(v=0.66, col = "blue", lty=2, lwd = 3); abline(v=0.33, col = "blue", lty=2, lwd = 3)
   abline(v=0.5, col = "purple", lty=2, lwd = 3)
  }
  else{ 
   yaxlim <- 1
   plot(1, 1, type = "p", pch = 20, cex = 0.5, main = "", 
        ylim=c(0,1), ylab = "", xlim=c(0.15,0.85), xlab = "", axes = FALSE)
  }

  mtext(text = nice_lbls[i], side = 2, las=1, line = 1, cex = 1.75)
  axis(side=2, at=c(0,yaxlim), las=1)
  axis(side=1, at=c(0.15, 0.25, 0.33, 0.5, 0.66, 0.75, 0.85), labels=rep("",7))
  if (i == length(Loci) ) {axis(side=1, at=c(0.15, 0.25, 0.33, 0.5, 0.66, 0.75, 0.85), las=2, cex.axis=2.25)}
  
  if (exists("H") == TRUE) { rm(H) }
  if (exists("H2") == TRUE) { rm(H2) }
  if (exists("sub_table2") == TRUE) { rm(sub_table2) }
# HISTOGRAMS #

# SNVs #  
  if (exists("sub_table") == TRUE && dim(sub_table)[1] > 1){
   plot(sub_table$Pos, sub_table[,strains[j]], type = "p", pch = 20, cex = 0.5, xlab = "", 
        main = "", col = "gray50",
        ylim=c(0,1), ylab = "", xlim=c(0,lngstLen), axes = FALSE)
#  if (i == 1) {
#   mtext(text = strains[j], side = 3, adj = 0.5, line = 0.75, cex = 2)
#   abline(v=c(276488,281259), col = "green", lty=3) #MatA
#   abline(v=c(18669,23440), col = "green", lty=3) #HML
#   abline(v=c(920174,921173), col = "red", lty=2, lwd = 3) #CEN
#  }
   #axis(side=2, at=c(0, 0.33, 0.66, 1), las=1)
   #axis(side=4, at=c(0.25, 0.5, 0.75), las=1, line = -1)
   #axis(side=1, at=seq(0, lngstLen, 250000), labels = F)
  }
  #else{ write("SNV sub_table check", stderr()); write(Loci[i], stderr()); plot.new()}
  else{ plot(lngstLen, 1, type = "p", pch = 20, cex = 0.5, xlab = "", 
             main = "", col = "gray50",
             ylim=c(0,1), ylab = "", xlim=c(0,lngstLen), axes = FALSE) }
  axis(side=2, at=c(0, 0.33, 0.66, 1), las=1)
  axis(side=4, at=c(0.25, 0.5, 0.75), las=1, line = -1)
  axis(side=1, at=seq(0, lngstLen, 250000), labels = F)
# SNVs #

# COVERAGE #  
  par(mar=c(0.5,3.5,1.5,1))
  sub_cov_table <- subset(cov_table, grepl(Loci[i], chrom))
  seg_sub_table <- subset(seg_cov_table, grepl(Loci[i], chrom))
  plot(sub_cov_table$maploc, sub_cov_table[,strains[j]], type = "p", pch = 20, cex = 0.5, col = "gray50",
       ylim=c(-mYlim,mYlim), xlim=c(0,lngstLen), axes=F, main="", ylab = "", xlab = "")
  axis(side=1, at=c(0,Loci_Len[i]), labels = c("",Loci_Labs[i]), padj = -0.75)
  axis(side=2, line=-0.5, at=c(-mYlim,-1.5,-1,-0.5,0,0.5,1,1.5,mYlim), 
       labels=c(paste("-",mYlim,sep=""),-1.5,-1,-0.5,0,"+0.5","+1","+1.5",paste("+",mYlim,sep="")), 
       las = 1)
  
  for (k in 1:nrow(sub_cov_table)){
    lines(x=c(sub_cov_table[k,"maploc_s"],sub_cov_table[k,"maploc"]),
          y=c(sub_cov_table[k,3],sub_cov_table[k,3]), col="gray50")
  }
  lines(x=c(0,Loci_Len[i]),y=c(0,0), col="blue", lty=1) 
  for (k in 1:nrow(seg_sub_table)){
   lines(x=c(seg_sub_table[k,"loc.start"],seg_sub_table[k,"loc.end"]),
         y=c(seg_sub_table[k,"seg.mean"],seg_sub_table[k,"seg.mean"]), col="red", lwd = 3)
  }
  lines(x=c(0,Loci_Len[i]), y=c(0.5,0.5), lty=2, lwd = 0.5)
  lines(x=c(0,Loci_Len[i]), y=c(1,1), lty=2, lwd = 0.5)
  lines(x=c(0,Loci_Len[i]), y=c(1.5,1.5), lty=2, lwd = 0.5)
  lines(x=c(0,Loci_Len[i]), y=c(2,2), lty=2, lwd = 0.5)
  lines(x=c(0,Loci_Len[i]), y=-c(0.5,0.5), lty=2, lwd = 0.5)
  lines(x=c(0,Loci_Len[i]), y=-c(1,1), lty=2, lwd = 0.5)
  lines(x=c(0,Loci_Len[i]), y=-c(1.5,1.5), lty=2, lwd = 0.5)
  lines(x=c(0,Loci_Len[i]), y=-c(2,2), lty=2, lwd = 0.5)
  
#  if (i == 1) {
#   abline(v=c(276488,281259), col = "green", lty=3) #MatA
#   abline(v=c(18669,23440), col = "green", lty=3) #HML
#   abline(v=c(920174,921173), col = "red", lty=2, lwd = 3) #CEN
#  }
  
  #axis(side=1, line=-0.5, at=0, las = 1, padj = 1, cex.axis=2.25)
  axis(side=1, line=-0.5, at=0, las = 1, padj = 0)
  par(mar=c(0.5,2,1.5,1))
# COVERAGE #
  
 }
 #dev.off()
}
dev.off()
