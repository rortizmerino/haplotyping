#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("\n
  usage:    Rscript --vanilla plot_dS.R <infile> <outfileprefix>
  example:  Rscript --vanilla plot_dS.R Scer_vs_LM07gt100kb.canu.full_dS\n\n",
       call.=FALSE)
}
cmd = paste(commandArgs(), collapse = " ")
write(cmd, stderr())

filename = as.character(args[1])

table <- read.table(paste(filename,"txt",sep = "."), sep = "\t", header = F, stringsAsFactors = F)
df = table[- grep("^Y", table$V2),]
df = df[grep("^Y", df$V1),]
df = df[- grep("^Seho_", df$V2),]

figname=paste(filename,"jpeg",sep = ".")
jpeg(file=figname, width=300, height=300, units="mm", quality=80, res=600, pointsize = 20);

H_breaks <- hist(df$V3, breaks = 100,  include.lowest = T, plot = F)$breaks
H_colors <- rep("#785EF0", length(H_breaks))              # ScerA purple
H_colors[H_breaks >= 0.25 & H_breaks < 2.75] <- "#DC267F" # ScerB pink
H_colors[H_breaks >= 2.75] <- "#FE6100"                    # Skud orange

H1 <- hist(df$V3, breaks = H_breaks,  include.lowest = T, freq = T, xlab = "dS",
           xlim = c(0,4),
           ylim = c(0,1250), 
           col  = H_colors,
           axes = F,
           main = gsub(".jpeg","",figname))
abline(v=0.25, lty = 2)
abline(v=2.75, lty = 2)

axis(side = 1, at = c(0,0.25,2.75,4), labels = c("", "0.25", "2.75", ""))
axis(side = 2, at = c(0,1250))

dev.off()

