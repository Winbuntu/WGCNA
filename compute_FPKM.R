library(data.table)


RC = as.data.frame(fread("RC.DRG.all.1year.tsv"))

MASK = as.data.frame(fread("MASK.table"))

RC.clean = RC[-match(MASK$Geneid,RC$Geneid),]

RC.clean.2 = RC.clean[,-19]

RC.to.FPKM <- function(x,fil.threshold = 2){
  
  x = as.data.frame(x)
  rownames(x) = x$Geneid
  gene.length = x$Length
  x.rc = x[,-c(1:6)]
  x.fpkm =  sweep(  (sweep(x.rc,1,gene.length,'/')*1000)   ,2,colSums(x.rc),"/")*1000000
  #boxplot( log2(x.fpkm+1) ,outline=FALSE )
  x.fpkm.fil = x.fpkm[ rowMeans(x.fpkm[,-13])>2  ,]
  
}

RC.clean.2.fpkm = RC.to.FPKM(RC.clean.2)
boxplot( log2(RC.clean.2.fpkm+1) ,outline=FALSE )








