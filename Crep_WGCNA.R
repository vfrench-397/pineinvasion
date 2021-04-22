########################################################################################################
#Get all genes into WGCNA but remove the genes with low basemean values
#navigate into the correct directory
setwd("/project/bi594/Pine_invasion/")
library(DESeq2)
library(BiocManager)
library(MCMC.OTU)

load('Environment.4.19.21.RData')

#Merge count table with Genus IDs
countData<-read.csv("otu.csv", stringsAsFactors = FALSE)
colnames(countData)[1] = "Sample_ID"

count.trim <- purgeOutliers(countData,count.columns=2:282)
#Trimmed RV55 sample and a little under 230 taxa; adjust the arguments here; losing too many? too few? 
rownames(count.trim)=count.trim$cdat
count.trim$cdat <- NULL
class(count.trim)

t<- t(as.data.frame(lapply(count.trim,as.numeric)))
colnames(t) <- rownames(count.trim)
class(t) #t transposes data.frame into a matrix 
ncol(t)
nrow(t)

#Read in treatment data
samdf<- read.csv("variabletable_pi.csv")
rownames(samdf)[samdf$SampleID=='RV55'] #insert samples removed from purgeoutliers 
samdf.trim <- samdf[-c(21),]
treat=samdf.trim$site_code
g=data.frame(treat)
g
colData<- g #create coldata for DESeq 

class(colData)
str(colData)
nrow(colData)


dds<-DESeqDataSetFromMatrix(countData=t, colData=colData, design=~ treat) 

#diagdds = phyloseq_to_deseq2(ps.rarefied, ~ site + position)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(dds, fitType="local")

#dds<-DESeq(dds)


#res<- results(dds)
#filter for contigs with average(baseMean) >3 ****
#res3<-res[res$baseMean>3, ]
#dim(res) #281
#dim(res3) #56
#Didn't filter any counts 
#Not filtering by base mean because basemean filter formatted for gene counts, not OTU abundance data. Filtering previously with purgeOutliers function. 

# get rlog data (better transformation when size factors vary across samples)
rld <- rlogTransformation(dds, blind=FALSE, fitType="local") #This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
head(rld)
rld_wg=(assay(rld)) #Making matrix of rlogTransformation data 
head(rld_wg)
nrow(rld_wg)
#56
rldFiltered=(assay(rld))[(rownames((assay(rld))) %in% rownames(dds)),]
nrow(rldFiltered)
#56
write.csv(rldFiltered,file="Invasion_wgcna_allgenes.csv",quote=F,row.names=T)
#now we have our filtered data to take into WGCNA

#source("http://bioconductor.org/biocLite.R") #To download DESeq package (you can comment these lines out, they only need to be run once ever)
#biocLite("WGCNA")
#biocLite("flashClust")
#If using R version or greater you need to use BiocManager instead:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("WGCNA")
#^This should download both WGCNA and flashClust

####First part of tutorial:Data input and cleaning
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

dat=read.csv("Invasion_wgcna_allgenes.csv")
head(dat) 
rownames(dat)<-dat$X
head(dat)
dat$X=NULL
head(dat)
names(dat)
nrow(dat)
#56
datExpr0 = as.data.frame(t(dat))

#Don't run, no need for extra filtering? 
#gsg = goodSamplesGenes(datExpr0, verbose = 1); #verbose=1 default, change if we want more verbose data ? 
#gsg$allOK #if TRUE, no outlier taxa, if false run the script below

#if (!gsg$allOK)
#{if (sum(!gsg$goodGenes)>0)
 # printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
  #if (sum(!gsg$goodSamples)>0)
   # printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
  #datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
#}
#gsg=goodSamplesGenes(datExpr0, verbose = 1)
#gsg$allOK 
#dim(datExpr0) 
#264  #17 taxa excluded #probably don't filter 

### Outlier detection incorporated into trait measures. 
traitData= read.csv("Invasion_traits_WGCNA.csv", row.names=1)
traitData= read.csv("Invasion_traits_WGCNA.csv")
rownames(traitData)[traitData$Sample_ID=='RV55'] #insert samples removed from purgeoutliers 
trait.trim <- traitData[-c(21),]
rownames(trait.trim) <- trait.trim$Sample_ID
trait.trim$Sample_ID <- NULL
dim(trait.trim)
head(trait.trim)
names(trait.trim)

# Form a data frame analogous to expression data that will hold the clinical traits.
dim(datExpr0)
rownames(datExpr0)
# datTraits=allTraits
datTraits=trait.trim

table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
head(datExpr0)

#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0)) #type back to default, no direction in ITS counts
#Calculates (correlation or distance) network adjacency from given expression data or from a similarity.
# this calculates the whole network connectivity we choose signed because we care about direction of gene expression
k=as.numeric(apply(A,2,sum))-1 #Summing columns of adjacency matrix (-1 to account for self correlation) 
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average") #**** average 
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")
#no outliers 

# Remove outlying samples from expression and trait data
# remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
# datExpr=datExpr0[!remove.samples,]
# datTraits=datTraits[!remove.samples,]

save(datExpr0,datTraits, file="Invasion_Samples_Traits_ALL.RData")

################Moving on!  Network construction and module detection - this section can take a lot of time you might consider running it on a cluster for a larger dataset
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads use this in base R
allowWGCNAThreads() 
lnames = load(file="Invasion_Samples_Traits_ALL.RData")

#Figure out proper SFT
# Choose a set of (candidate) soft-thresholding powers
powers = c(seq(1, 90, by = 10), seq(100, 200, by = 10)); #default; Producing SFT.R.sq wayyyyy too small 
#may need to adjust these power values to hone in on proper sft value
#soft threshold powers are the power to which co-expression similarity is raised to calculate adjacency 
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="unsigned", verbose = 2) #want smallest value, closest to 0.9 (but still under)
#Printing table to help decide estimate for soft power threshold 
#Soft power threshold chosen based on SFT.R.sq being over .8 and mean k being below the hundreds 
#should always be less than 15 though for unsigned data 
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

#First, the user should ensure that variables (probesets, genes etc.) have not been filtered by 
#differential expression with respect to a sample trait. (with respect to traits we have not?)

#If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers
#(less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks) 
#and the mean connectivity remains relatively high (in the hundreds or above), 
#chances are that the data exhibit a strong driver that makes a subset of the 
#samples globally different from the rest. The difference causes high correlation 
#among large groups of genes which invalidates the assumption of the scale-free topology 
#approximation.

#If we want to keep the variation then according to FAQ unsigned networks with 20-30 samples should use power threshold of 8 

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower=8 #smallest value to plateau at ~0.85
adjacency=adjacency(datExpr0, power=softPower,type="unsigned") #must change method type here too!!
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "unsigned")
dissTOM= 1-TOM

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(10,6)
# pdf(file="dendrogram_thresh16.5_signed_1868.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
# dev.off()
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes

minModuleSize=2 #we only want large modules #set to 90 originally (How many taxa do we want? Arbitrary?)
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)
#5 modules 35 out of 56 taxa not in a module 

dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merg modules whose expression profiles are very similar
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 13)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_crep_nomerge.RData")

lnames = load(file = "Network_crep_nomerge.RData")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.6
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_0.6.RData")

###############Relating modules to traits and finding important genes
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Crep_Samples_Traits_ALL.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Network_signed_0.6.RData");
lnames = load(file = "Network_crep_nomerge.RData");
lnames

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)

#moduleColors
#black         blue       coral2     darkgrey  floralwhite       grey60        ivory 
#1338         2928         1167         1088          895          206          564 
#mediumorchid   orangered4       purple    steelblue 
#478         1557          615         1749 

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#represent module trait correlations as a heatmap
quartz()
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$pH7.5); #change Lipidrobust to your trait name
names(weight) = "pH7.5"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Gene-trait significance correlation plots
# par(mfrow=c(2,3))
module = "coral2"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for pH7.5",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Making VSD files by module for GO plot functions
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="coral2"]) #black  blue brown green  grey  pink   red 

c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size
table(moduleColors)
#moduleColors
#black         blue       coral2     darkgrey  floralwhite       grey60        ivory 
#1338         2928         1167         1088          895          206          564 
#mediumorchid   orangered4       purple    steelblue 
#478         1557          615         1749
head(c.vsd)
write.csv(c.vsd,"rlog_MMcoral2.csv",quote=F)

##############################heatmap of module expression with bar plot of eigengene, no resorting of samples...
#names(dis)
sizeGrWindow(8,7);
which.module="coral2" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset

#quartz()
# par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")
#this is a cool plot where you can see that genes in this module are upregulated in the pH7.5 treatment

##############################heatmap of module expression with bar plot of trait of interest by sample...
#here we just have binary traits, but if you have a continuous trait this code is cool
sizeGrWindow(8,7);
which.module="yellow" #pick module of interest
which.trait="fvfm" #change trait of interest here
datTraits=datTraits[order((datTraits$fvfm),decreasing=T),]#change trait of interest here

trait=datTraits[, paste(which.trait)]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset
genes=genes[rownames(datTraits),]

#quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(trait, col=which.module, main="", cex.main=2,
        ylab="fvfm",xlab="sample")#change trait of interest here

#Gene relationship to trait and important modules: Gene Significance and Module membership
allkME =as.data.frame(signedKME(t(dat), MEs))
head(allkME)
vsd=read.csv(file="rlog_MMcoral2.csv", row.names=1)
head(vsd)
gg=read.table("Crep454_iso2gene.tab", sep="\t")
head(gg)
library(pheatmap)

############################################
whichModule="coral2"
top=100

datME=MEs
vsd <- read.csv("Crep_wgcna_allgenes.csv", row.names=1)
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
summary(hubs)

gnames=c();counts=0
for(i in 1:length(hubs[,1])) {
  if (row.names(hubs)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(hubs)[i],2]
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
  } else { 
    gnames=append(gnames,i)
  }
} 
row.names(hubs)=gnames
length(hubs)

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#quartz()
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))

###fisher for GO
##########fisher of module vs whole dataset
library(WGCNA)
vsd <- read.csv("Crep_wgcna_allgenes.csv", row.names=1)
head(vsd)
options(stringsAsFactors=FALSE)
data=t(vsd)
allkME =as.data.frame(signedKME(data, MEs))

whichModule="coral2" # name your color and execute to the end

length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_fisher.csv",sep=""),quote=F)

modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

######--------------------end--------------------#######
