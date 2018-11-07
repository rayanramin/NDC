## IMPORTING ALL THE DEPTH FILES ( from $ samtools depth -aa control.alignment.bam > control.depth ) ##
Control <- read.table("control.depth")
Treatment <- read.table("treatment.depth")
mav = function(x,n=100){filter(x,rep(1/n,n),sides=2)}
mC <- mean(Control[,3])
mT <- mean(Treatment[,3])
#calculating NDC :
Treatment[,4] <-  ((mav(Treatment[,3]) / mT ) - (mav(Control[,3]) / mC ))
#calculateing the 5 sigma threshold
FSD <- sd(T[,4], na.rm=TRUE)*5
names(T) <- c("Chr", "POS", "COV" ,"NDC")

## NDC PLOTS
options(bitmapType='cairo')

png(filename = "NDC.png" ,height = 6, width = 24, units = "in", res = 400)
par( mai = c(1,1,0.7,0.2))
xl <- pretty(Treatment[,2], n = 10)/10^6
plot(Treatment[,2]/1e6, Treatment[,4], type = "l",main = "NDC" ,axes=F, cex.main =2.5 , xlab = "", ylab="")
axis(2,cex.axis=2)
axis (1, cex.lab=2 ,cex.axis = 2,cex=2 , at =xl , labels = paste(xl, "Mbp", sep = ""), las=1)
abline(h= FSD , col="red")
abline(h= -FSD , col="red")
title( ylab = "Relative Coverage", xlab ="Genomic Position", cex.lab=2)
box()
dev.off()

## FINDING REGIONS OF SIGNIFICANT ENRICHMENT OR PEAKS##
library(R.utils)
T_SPD <- seqToIntervals(Treatment$POS[Treatment$NDC > FSD])
T_SPD <- as.data.frame(T_SPD)


## EXPORT LIST OF PEAKS AS A BED FILE##

fixSPD = function(x){                           
z <- as.data.frame(get(x))
z[,3] <-z[,2]
z[,2]<-z[,1]
z[,1] <- T[1,1]
names(z) <- c("chr.", "from","to")
z
}
x = "T_SPD"
assign(paste(x) ,fixSPD(x))
write.table(as.data.frame(get(x)), paste(x,".bed", sep=""), quote=F, sep="\t", row.names=F, col.names=F)

## INTERSECTING SPD WITH ANNOTATION FILE ## NEEDS BEDTOOLS ##
#Credit to Altuna Akalin # http://zvfak.blogspot.com/2011/02/calling-bedtools-from-r.html #

bedTools.2in = function(functionstring="bedIntersect",bed1,bed2,opt.string=""){
 #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F, sep="\t")
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

annotations <- read.table("annotation.gff", sep="\t") #import list of annotations
j <- bedTools.2in("bedtools intersect", annotations,get(x))
write.table(j, paste(x,"genes.txt", sep="_"), quote=F, sep="\t", row.names=F, col.names=F)


