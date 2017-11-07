FPKM = read.table("FPKM/genes.fpkm_table",header = T,row.names = 1)
head(FPKM)

dim(FPKM)
FPKM.clean = FPKM[,-13]

head(FPKM.clean)
par(mfrow=c(1,1))
boxplot(log2(FPKM.clean+1),outline=FALSE)


FPKM.fil = FPKM.clean[rowMeans(FPKM.clean) > 2,]


dim(FPKM.fil)


boxplot(log2(FPKM.fil+1),outline=FALSE)




head(FPKM.fil)

#######################################
library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
#enableWGCNAThreads()
disableWGCNAThreads()

############

mydata=log2(FPKM.fil+1)
dim(mydata)

gene.names = rownames(mydata)


mydata.trans=t(mydata);


datExpr=mydata.trans
SubGeneNames=gene.names

#######



powers = c(1:50)
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,
                      powerVector = powers,corFnc = cor,
                      corOptions = list(use = 'p'),
                      networkType = "signed")


par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##################


softPower = 18;

adj= adjacency(datExpr,type = "signed", power = softPower);

#Turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,
                          networkType = "signed", 
                          TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) = SubGeneNames
dissTOM = 1-TOM

######

#Hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#Plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="",cex=0.3);


####################################

minModuleSize = 20;

# Module identification using dynamic tree cut, you can also choose the hybrid method

dynamicMods = cutreeDynamic(dendro = geneTree,  
                            method="tree", 
                            minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, 
#                            distM = dissTOM, method="hybrid", deepSplit = 3, 
#                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#Get the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

## dynamicMods
##   0   1   2 
## 253 159  88

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

## dynamicColors
##      blue      grey turquoise 
##        88       253       159

plotDendroAndColors((geneTree),  dynamicColors,
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

datME=moduleEigengenes(datExpr,dynamicColors)$eigengenes
signif(cor(datME, use="p"), 2)

datTraits = c(1,1,1,1,0,
      1,0,0,0,0,
      0,0,0,1)


signif(cor(y,datME, use="p"),2)

p.values = corPvalueStudent(cor(y,datME, use="p"), nSamples = length(y))
par(mfrow = c(2,1));

barplot(signif(cor(y,datME, use="p"),2))
barplot(p.values,ylim = c(0,1.2))
#####################################

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

# Display the correlation values within a heatmap plot
pdf(file = "aa.pdf",width=11,height=30)
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
####################

# Define variable weight containing the weight column of datTrait
weight = data.frame(datTraits)

names(weight) = "weight"

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

module.list = c("darkred","pink","palevioletred2","plum2","brown2","firebrick4","orangered1","darkviolet","orangered3","salmon4")
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

heatmap.2(as.matrix( log2( diff.gene.exp+0.5 ) ),trace = "none",density = "none",
          Colv = as.dendrogram(complete.cluster),
          #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
          ColSideColors =  c("grey","red")[factor(type)],
          col=color.palette,
          breaks = palette.breaks,
          scale = c("row"),
          dendrogram = "both")
dev.off()

palette.breaks <- seq(-4, 4, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)
heatmap.2(datExpr.selected, ColSideColors = c("red","blue")[datTraits+1] ,Rowv = NULL,trace = "none",density = "none",
          scale = "row",
          col=color.palette,
          breaks = palette.breaks,
          RowSideColors = unlist(mapply(rep,names(b),sapply(b,length)))    )



module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=  convert_gene_name(SubGeneNames[which(dynamicColors==color)])
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}


