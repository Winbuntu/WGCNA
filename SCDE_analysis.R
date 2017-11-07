library(scde)

RC.df.clean = read.csv("RC.df.clean.csv",row.names = 1)

#RC.df.clean = RC.df.clean[,-c(11,13,14)]

coldata # generated in DEseq.R 

#coldata = coldata[-c(11,13,14),]

sg = coldata$condition

names(sg) = rownames(coldata)

table(sg)

cd <- clean.counts(RC.df.clean, min.lib.size=1000,
                   min.reads = 1, min.detected = 5) ## changed here

o.ifm <- scde.error.models(counts = cd, 
                           groups = sg, 
                           n.cores = 1, 
                           threshold.segmentation = TRUE, 
                           save.crossfit.plots = FALSE, 
                           save.model.plots = FALSE, 
                           verbose = 1)

valid.cells <- o.ifm$corr.a > 0

table(valid.cells)

o.ifm <- o.ifm[valid.cells, ]

o.prior <- scde.expression.prior(models = o.ifm, 
                                 counts = cd, 
                                 length.out = 400, 
                                 show.plot = FALSE)



groups <- factor(sg)
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, 
                                    groups  =  groups, 
                                    n.randomizations  =  100, 
                                    n.cores  =  1, 
                                    verbose  =  1,
                                    return.posteriors = FALSE)


###############

head(ediff[order(ediff$Z, decreasing  =  TRUE), ],30)


write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "scde_results_filtered.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

scde.test.gene.expression.difference("Slc15a2", 
                                     models = o.ifm, 
                                     counts = cd, 
                                     prior = o.prior)
###################################


cd = RC.df.clean
cd <- cd[rowSums(cd>2)>3, ]

knn <- knn.error.models(as.matrix(cd), k = ncol(cd)/2, n.cores = 1, 
                        min.count.threshold = 2, min.nonfailed = 5, 
                        max.model.plots = 10)

# EVALUATION NOT NEEDED FOR SAKE OF TIME
varinfo <- pagoda.varnorm(knn, counts = cd, 
                          trim = 3/ncol(cd), 
                          max.adj.var = 5, 
                          n.cores = 1, plot = TRUE)

# normalize out sequencing depth as well
varinfo <- pagoda.subtract.aspect(varinfo, 
                                  colSums(cd[, rownames(knn)]>0))


library(FactoMineR)

pca.res = PCA(t(  varinfo$mat  ),graph = F)

PC1 <- as.numeric(pca.res$ind$coord[,1])
PC2 <- as.numeric(pca.res$ind$coord[,2])

sample.names = rownames(pca.res$ind$coord)

PCs <- data.frame(PC1,PC2,sample.names)

library(ggplot2)
library(ggrepel)

P<-ggplot(PCs, aes(PC1,PC2)) 
P +geom_point(shape = factor(coldata$condition)   ,size = 3) + 
  xlab(paste("PC1",as.character(round(pca.res$eig[,2][1],2)),"%")) + 
  ylab(paste("PC2",as.character(round(pca.res$eig[,2][2],2)),"%")) +
  ggtitle(   "PCA on RNAseq"    ) + 
  geom_text_repel(aes(label=sample.names) ,box.padding = unit(0.25, "lines"),
                  point.padding = unit(0.5, "lines")  )  + 
  theme_classic(base_size = 10) + ylim(c(-90,60)) + xlim(c(-140,90))

#############################################
# this is defined in compute_FPKM.R
# RC.clean.2.fpkm
good.genes = rownames(RC.clean.2.fpkm)




###########################
library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
#enableWGCNAThreads()
disableWGCNAThreads()

############

mydata=varinfo$mat;
dim(mydata)

mydata = mydata[match(good.genes,rownames(mydata)),]

#gene.names=names(sort(varinfo$arv,decreasing=T));

mydata.trans=t(mydata);

#######################################
library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
#enableWGCNAThreads()
disableWGCNAThreads()

############

dim(mydata)

#gene.names=names(sort(varinfo$arv,decreasing=T));

#mydata.trans=t(mydata);


n=5000;
datExpr=mydata.trans[,gene.names[1:n]];
SubGeneNames=gene.names[1:n];


powers = c(c(1:30));
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")

#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 22;

#calclute the adjacency matrix
adj= adjacency(datExpr,type = "signed", power = softPower);

TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) =SubGeneNames
dissTOM=1-TOM

geneTree = flashClust(as.dist(dissTOM),method="average");

plot(geneTree, xlab="", sub="",cex=0.3);


minModuleSize = 20;

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);

dynamicColors = labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

###########################################

datTraits = c(1,1,1,1,0,
              1,0,0,0,0,
              0,0,0,1)

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

############
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

pdf(file = "scde.pdf",width=11,height=30)
par(mfrow = c(1,1));
par(mar = c(6, 9, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()





# Define variable weight containing the weight column of datTrait
weight = data.frame(datTraits)

names(weight) = "weight"

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

#module.list = c("darkgrey","darkviolet","yellow","brown2","deeppink","midnightblue","pink4")
module.list = c("darkseagreen4","lavenderblush3","antiquewhite4","yellow4","ivory","saddlebrown","darkred","coral1")
b = sapply(module.list, function(x){which(dynamicColors == x)},simplify = T)

nm = length(module.list)

datExpr.selected = NULL

for(i in c(1:nm)){
  temp = datExpr[,b[[i]]]
  temp = t(temp)
  if(is.null(dim(datExpr.selected))){
    datExpr.selected = temp
  }else{
    datExpr.selected = rbind(datExpr.selected,temp)
    
  }
}
library(gplots)

palette.breaks <- seq(-3, 3, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)
heatmap.2(datExpr.selected, ColSideColors = c("red","blue")[datTraits+1] ,Rowv = NULL,trace = "none",density = "none",
          scale = "row",
          col=color.palette,
          breaks = palette.breaks,
          RowSideColors = unlist(mapply(rep,names(b),sapply(b,length)))    )


###################
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=  (SubGeneNames[which(dynamicColors==color)])
  write.table(module, paste("./scde_module/module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}


