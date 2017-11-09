library(data.table)

RC = as.data.frame(fread("DRG.all1year.gencode.vM15.RC.tsv"))

MASK = as.data.frame(fread("MASK.table"))

RC.clean = RC[-match(MASK$Geneid,RC$Geneid),]

RC.df = RC.clean[,7:43]

colSums(RC.df)
