#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
 stop("usage: Rscript --vanilla DotPlot_ANI.R <query> <ref> <picname.jpeg>", call.=FALSE)
}

file1 = paste(as.character(args[1]),".coords",sep="")
file2 = paste(as.character(args[2]),".coords",sep="")
file3 = paste(as.character(args[2]),".comparison.coords.txt",sep="")
file4 = as.character(args[3])
file5 = paste("../../../references/genomes/",as.character(args[2]),".feat.coords",sep="")
file6 = paste(as.character(args[1]),".ids_n_stats.txt",sep="")

tit = paste(args[1], "vs", ylab=args[2], sep=" ") 

# Rscript --vanilla DotPlot_ANI.R S1phunph Scer_noMt S1phunph_vs_Scer_noMt.jpeg
#file1 = paste(as.character("S1phunph"),".coords",sep="")
#file2 = paste(as.character("Scer_noMt"),".coords",sep="")
#file3 = paste(as.character("Scer_noMt"),".comparison.coords.txt",sep="")
#file4 = as.character("S1phunph_vs_Scer_noMt.ANI.jpeg")
#file5 = paste(as.character("Scer_noMt"),".feat.coords",sep="")
#file6 = paste(as.character("S1phunph"),".ids.txt",sep="")
#
#tit = paste("S1phunph", "vs", ylab="Scer_noMt", sep=" ") 

cmd = paste(file1,file2,file3,file4,file5,sep = " ")
write(cmd, stderr())

# Ref, X
table0 <- read.table(file1, sep = " ", header = T, stringsAsFactors = F,
                     colClasses = c("character", "character", "numeric", "numeric",
                                    "numeric", "character", "character"))
table0$cumLen <- cumsum(table0$Len)

# Query, Y
table1 <- read.table(file2, sep = " ", header = T, stringsAsFactors = F,
                     colClasses = c("character", "character", "numeric", "numeric",
                                    "numeric", "character", "character"))
table1$cumLen <- cumsum(table1$Len)

alntable <- read.table(file3, sep = " ", header = T, stringsAsFactors = F,
                       colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                      "character", "character", "character", "character"))

featable <- read.table(file5, sep = " ", header = T, stringsAsFactors = F,
                       colClasses = c("numeric", "character", "numeric", "numeric", "numeric"))
featable$cumLen <- cumsum(featable$Len)
tmpVec <- c(0,featable[,"cumLen"])
featable$corCENcoords <- tmpVec[-length(tmpVec)] + featable$CENcoords

ANItable <- read.table(file6, sep = " ", header = F, stringsAsFactors = F,
                       colClasses = c("character", "character", "character", "character", "numeric", 
                                      "numeric", "numeric"))
ANItable <- ANItable[,c(1,7)]
colnames(ANItable) <- c("G2", "ANI")
#testtable <- merge(alntable, ANItable, by = "G2", sort = F)
alntable <- merge(alntable, ANItable, by = "G2", sort = F)

for (i in 1:nrow(alntable)){
 for (j in 2:length(table1$ID)){
  k=j-1;
  if (grepl(pattern = table1$Order[j], x = alntable[i,"O1"], fixed = T) == T){
   alntable[i,"S1"] = alntable[i,"S1"] + table1[k,"cumLen"]
   alntable[i,"E1"] = alntable[i,"E1"] + table1[k,"cumLen"]
  }
 }
}

for (i in 1:nrow(alntable)){
 for (j in 2:length(table0$ID)){
  k=j-1;
  if (grepl(pattern = table0$Order[j], x = alntable[i,"O2"], fixed = T) == T){
   alntable[i,"S2"] = alntable[i,"S2"] + table0[k,"cumLen"]
   alntable[i,"E2"] = alntable[i,"E2"] + table0[k,"cumLen"]
  }
 }
}

table0$midCumLen <- table0$cumLen - (as.numeric(table0$Len)/rep(2,times=nrow(table0)))
table1$midCumLen <- table1$cumLen - (as.numeric(table1$Len)/rep(2,times=nrow(table1)))

#jpeg(file=file4, width=333, height=333, units="mm", quality=80, res=600)
jpeg(file=file4, width=230, height=152, units="mm", quality=80, res=600)
par(mar=c(2,2.5,1.5,1), cex=0.5)
plot(table0$cumLen[length(table0$cumLen)], table1$cumLen[length(table1$cumLen)], 
     type="p", main=tit, col="gray", axes = FALSE, xaxs = "i", yaxs="i",
     xlim=c(0,table0$cumLen[length(table0$cumLen)]), 
     ylim=c(0,table1$cumLen[length(table1$cumLen)]), xlab="", ylab="")

abline(v=c(0,table0$cumLen), col="gray")
abline(h=c(0,table1$cumLen), col="gray")

# centromeres
abline(h=featable$corCENcoords, col="gray", lty=5)

axis(side=1, at=c(0,table0$cumLen), labels=FALSE, pos=0)
axis(side=1, at=table0$midCumLen, labels=seq(1,nrow(table0)), tick=FALSE, pos=0)
axis(side=2, at=c(0,table1$cumLen), labels=FALSE, pos=0)
axis(side=2, at=table1$midCumLen, labels=seq(1,nrow(table1)), tick=FALSE, las = 1, pos=0)

#axis_offset = -250000
axis_offset = 0
for (i in 1:nrow(table0)){
 if (grepl(pattern = ",", x = table0[i,"Nblock_coords"]) == T){
  all_N_coords <- strsplit(table0[i,"Nblock_coords"], ",")
  for (j in 1:length(all_N_coords[[1]])){
   N_coords <- strsplit(all_N_coords[[1]][j], "-")
   if (i==1) {
    lines(x=c(N_coords[[1]][1],N_coords[[1]][2]), y=c(axis_offset,axis_offset), col="red", type="p", pch="|")
    lines(x=c(N_coords[[1]][1],N_coords[[1]][2]), y=c(axis_offset,axis_offset), col="red", type="l")
   }
   if (i>1) {
    k=i-1
    strt=as.numeric(N_coords[[1]][1])+table0$cumLen[k]
    stp=as.numeric(N_coords[[1]][2])+table0$cumLen[k]
    lines(x=c(strt,stp), y=c(axis_offset,axis_offset), col="red", type="p", pch="|")
    lines(x=c(strt,stp), y=c(axis_offset,axis_offset), col="red", type="l")
   }
  }
 }
}
for (i in 1:nrow(alntable)){

# lines(x=c(alntable[i,"S2"],alntable[i,"E2"]), y=c(alntable[i,"S1"],alntable[i,"E1"]))
# if (alntable[i,"CQ"] > 25){  
  if (alntable[i,"ANI"] >= 98 ) {
   lines(x=c(alntable[i,"S2"],alntable[i,"E2"]),
         y=c(alntable[i,"S1"],alntable[i,"E1"]), col="#785EF0" ) # ScerA purple
  }
  else if (alntable[i,"ANI"] < 98 && alntable[i,"ANI"] > 82) {
   lines(x=c(alntable[i,"S2"],alntable[i,"E2"]),
         y=c(alntable[i,"S1"],alntable[i,"E1"]), col="#DC267F" ) # ScerB pink
  }
  else {
   lines(x=c(alntable[i,"S2"],alntable[i,"E2"]),
         y=c(alntable[i,"S1"],alntable[i,"E1"]), col="#FE6100" ) # Skud orange
  }
# } 
  
}
dev.off()

#hist(alntable$CQ, breaks = 100,  include.lowest = F, freq = T, xlim = c(25,100), ylim = c(0,10))
