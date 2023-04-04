#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("\n
  usage:    Rscript --vanilla ANIhist.R <infile> <outfileprefix>
  example:  Rscript --vanilla ANIhist.R assembly/LM07/default/LM07.r2cat.canu.ids_n_stats.txt assembly/LM07/default/LM07\n\n", 
       call.=FALSE)
}
cmd = paste(commandArgs(), collapse = " ")
write(cmd, stderr())

infile  = as.character(args[1])
outfile = as.character(args[2])

fulltable <- read.table(infile, sep = " ", header = F, stringsAsFactors = F)

divergence0 <-  100-fulltable$V7

figname=paste(outfile,"fastANIhist.jpeg",sep=".")
jpeg(file=figname, width=300, height=300, units="mm", quality=80, res=600, pointsize = 20);

H1 <- hist(divergence0, breaks = 100,  include.lowest = T, freq = T, xlab = "% divergence",
           xlim = c(0,25),
           #ylim = c(0,1500), 
           main = gsub(".jpeg","",figname), 
           axes = F)

axis(side = 1, at = seq(0,25,5))
axis(side = 1, at = seq(0,25,1), lwd.ticks=0.5, labels=F)
axis(side = 2)

#strng <- length(divergence0[divergence0 < 0.3])
#strng <- paste("% divergence < 0.3 = ", strng, sep = "")
#mtext(strng, side = 3, adj = 1, cex=0.66)

dev.off()