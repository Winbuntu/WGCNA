library(DESeq2)

RC.df.clean = read.csv("RC.df.clean.csv",row.names = 1)

cts <- as.matrix(RC.df.clean)


coldata <- data.frame(condition = c(
  "injury",
  "intact",
                                    "injury",
                                    "injury",
                                    "injury",
  
  
                                    "intact",
                                    "injury",
                                    "intact",
                                    "intact",
                                    "intact",
  
                                    "intact",
                                    "intact",
                                    "intact",
                                    "injury") ,
                      type = c(
                               "multiple cells",
                               "single cell",
                               "multiple cells",
                               "single cell",
                               "single cell",
                               
                               "single cell",
                               "single cell",
                               "single cell",
                               "single cell",
                               "single cell",
                               
                               "single cell",
                               "multiple cells",
                               "multiple cells",
                               "multiple cells"
                               ))
rownames(coldata) = colnames(cts)
#################################


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

dds <- dds[ rowSums(counts(dds)) > 1, ]

dds$condition <- relevel(dds$condition, ref="intact")
dds$condition


####################################


dds <- DESeq(dds)
res <- results(dds)

res

resOrdered <- res[order(res$padj),]

summary(res)

library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)

sum(resIHW$padj < 0.1, na.rm=TRUE)

resIHW.up = resIHW[intersect(which(resIHW$padj < 0.1),
                 which(resIHW$log2FoldChange > 0)) ,]

write.csv(as.data.frame(resIHW.up), 
          file="resIHW.up.csv")

####

resIHW.down = resIHW[intersect(which(resIHW$padj < 0.1),
                             which(resIHW$log2FoldChange < 0)) ,]

write.csv(as.data.frame(resIHW.down), 
          file="resIHW.down.csv")


###################


dds <- estimateSizeFactors(dds)
RC.df.clean.norm = counts(dds, normalized=TRUE)

RC.df.clean.norm.big = RC.df.clean.norm[rowMeans(RC.df.clean.norm) > 8,]


#################
library(FactoMineR)
pca.res = PCA(t(  log2(  RC.df.clean.norm.big +1)  ),graph = F)

PC1 <- as.numeric(pca.res$ind$coord[,1])
PC2 <- as.numeric(pca.res$ind$coord[,2])

sample.names = rownames(pca.res$ind$coord)

PCs <- data.frame(PC1,PC2,sample.names)

library(ggplot2)
library(ggrepel)

P<-ggplot(PCs, aes(PC1,PC2)) 
P +geom_point(shape= coldata$condition,size = 5,color = c("red","blue")[factor(coldata$type)] ) + 
  xlab(paste("PC1",as.character(round(pca.res$eig[,2][1],2)),"%")) + 
  ylab(paste("PC2",as.character(round(pca.res$eig[,2][2],2)),"%")) +
  ggtitle(   "PCA on RNAseq"    ) + 
  geom_text_repel(aes(label=sample.names) ,box.padding = unit(0.25, "lines"),
                  point.padding = unit(0.5, "lines")  )  + 
  theme_classic(base_size = 10) 


