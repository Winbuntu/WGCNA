library(data.table)


RC = as.data.frame(fread("RC.DRG.all.1year.tsv"))

MASK = as.data.frame(fread("MASK.table"))

RC.clean = RC[-match(MASK$Geneid,RC$Geneid),]

RC.df = RC.clean[,7:21]

row.names(RC.df) = RC.clean$Geneid

RC.df.clean = RC.df[,-13]

colnames(RC.df.clean)

colnames(RC.df.clean) = c("17E.5.5",
                          "17E.6.1",
                          "17E.6.11",
                          "17E.6.14",
                          "17E.6.16",
                          
                          "17E.6.18",
                          "17E.6.19",
                          "17E.6.27",
                          "17E.6.30",
                          "17E.6.34",
                          
                          "17E.6.36",
                          "17E.1.5",
                          "17E.3.5",
                          "17E.4.5")

write.csv(RC.df.clean,file = "RC.df.clean.csv")

#############################



