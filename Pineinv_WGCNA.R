########################################################################################################
###################Simmi

#Get all genes into WGCNA but remove the genes with low basemean values

#navigate into the correct directory
setwd("/project/bi594/Pine_invasion/")

library(WGCNA)
library(flashClust)
library(DESeq2)
library(BiocManager)
library(MCMC.OTU)

load('Environment.4.25.21.RData')

#Merge count table with Genus IDs
countData<-read.csv("otu.csv", stringsAsFactors = FALSE)
colnames(countData)[1] = "Sample_ID"

count.trim <- purgeOutliers(countData,count.columns=2:282, otu.cut= .000000000001)
#count.trim <- purgeOutliers(countData,count.columns=2:282)
#Trimmed RV55 sample and a little under 230 taxa; adjust the arguments here; losing too many? too few? 
rownames(count.trim)=count.trim$cdat
count.trim$cdat <- NULL
class(count.trim)

t<- t(as.data.frame(lapply(count.trim,as.numeric)))
colnames(t) <- rownames(count.trim)
class(t) #t transposes data.frame into a matrix 
sapply(t, is.numeric)
ncol(t)
nrow(t)
colnames(t) = c("PL1.1", "PL1.2", "PL1.3", "PL1.4", "PL1.5", "UN1.1", "UN1.2", "UN1.3","UN1.4","UN1.5", "INV1.1", "INV1.2", "INV1.3", "INV1.4", "PL2.1", "PL2.2", "PL2.3", "PL2.4", "PL2.5", "UN2.1", "UN2.3", "UN2.4", "UN2.5", "INV2.1", "INV2.2", "INV2.3", "INV2.4", "INV2.5")

#Read in treatment data
samdf<- read.csv("variabletable_pi.csv")
rownames(samdf)[samdf$SampleID=='RV55'] #insert samples removed from purgeoutliers 
samdf.trim <- samdf[-c(21),]
treat=samdf.trim$site_code
g=data.frame(treat)
g
colData<- g #create coldata for DESeq , **dont really need to do if running blind=TRUE

class(colData)
str(colData)
nrow(colData)


#**running deseq to filter super low counts**
dds<-DESeqDataSetFromMatrix(countData=t, colData=colData, design=~ treat) 

#diagdds = phyloseq_to_deseq2(ps.rarefied, ~ site + position)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(dds, fitType="local")

 
#Not filtering by base mean before rlog transforing the data because basemean filter formatted for gene counts, not OTU abundance data. Filtering previously with purgeOutliers function. 

# get rlog data (better transformation when size factors vary across samples)
rld <- rlogTransformation(dds, blind=TRUE, fitType="local")
#**should be blind=TRUE ? the benefit of wgcna is that it has no idea what our model was so normalize with this function**
#This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
head(rld)
rld_wg=(assay(rld)) #Making matrix of rlogTransformation data 
head(rld_wg)
nrow(rld_wg)
#281
rldFiltered=(assay(rld))[(rownames((assay(rld))) %in% rownames(dds)),]
nrow(rldFiltered)
#*filter by the ones removed from res3 function, to get rid of ones without the base mean, do we still do this if we didnt filter
#281

write.csv(rldFiltered,file="Invasion_wgcna_allgenes2.csv",quote=F,row.names=T)
#now we have our filtered data to take into WGCNA


####First part of tutorial:Data input and cleaning
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

dat=read.csv("Invasion_wgcna_allgenes2.csv")
head(dat) 
rownames(dat)<-dat$X
head(dat)
dat$X=NULL
head(dat)
names(dat)
nrow(dat)
#56
datExpr0 = as.data.frame(t(dat))


### Outlier detection incorporated into trait measures. 
#*make sure trait data is same as Expr0 and names are the same (linking things properly)
#traitData= read.csv("Invasion_traits_WGCNA.csv", row.names=1)
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
row.names(datTraits)= c("PL1.1", "PL1.2", "PL1.3", "PL1.4", "PL1.5", "UN1.1", "UN1.2", "UN1.3","UN1.4","UN1.5", "INV1.1", "INV1.2", "INV1.3", "INV1.4", "PL2.1", "PL2.2", "PL2.3", "PL2.4", "PL2.5", "UN2.1", "UN2.3", "UN2.4", "UN2.5", "INV2.1", "INV2.2", "INV2.3", "INV2.4", "INV2.5")

table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
head(datExpr0)

#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0)) #type back to default, no direction in ITS counts
#Calculates (correlation or distance) network adjacency from given expression data or from a similarity.
# this calculates the whole network connectivity we choose signed (this is the default) because we are dealing with OTU counts, not gene expression
k=as.numeric(apply(A,2,sum))-1 #Summing columns of adjacency matrix (-1 to account for self correlation) 
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average") #**** why did you change from default ("complete")?
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")
#no outliers 


save(datExpr0,datTraits, file="Invasion_Samples_Traits_ALL2.RData")

#################Victoria
######################### NETWORK CONSTRUCTION AND MODULE DETECTION ######################
################Moving on!  Network construction and module detection - this section can take a lot of time you might consider running it on a cluster for a larger dataset
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads use this in base R
allowWGCNAThreads() 
lnames = load(file="Invasion_Samples_Traits_ALL2.RData")

###################### SOFT POWER THRESHOLD #####################
#Figure out proper SFT
#**notes from tutorial: The function adjacency calculates the adjacency matrix from expression data.
#*Adjacency functions for both weighted and unweighted networks require the user to choose threshold parameters, for example by applying the approximate scale-free topology criterion 
#*general framework for soft thresholding that weighs each connection, suggested to look at scale-free topology index in WGCNA faq to decide soft threshold based on number of samples
#*dont want to include all possible data (not be too high) bc then its overfitting wgcna model
#*emphasizing network on stronger association/larger correlation coefficient and reduce the noise
# Choose a set of soft-thresholding powers *****what is this?? how do we choose these powers?
powers = c(seq(100, 190,by=10), seq(200,250, by=5))
powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))
; #may need to adjust these power values to hone in on proper sft value

# Choose a set of (candidate) soft-thresholding powers
#default; Producing SFT.R.sq wayyyyy too small 
#may need to adjust these power values to hone in on proper sft value
#soft threshold powers are the power to which co-expression similarity is raised to calculate adjacency 
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9 (but still under)
#Printing table to help decide estimate for soft power threshold 
#Soft power threshold chosen based on SFT.R.sq being over .8 and mean k being below the hundreds 
#should always be less than 15 though for signed data 
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

#First, the user should ensure that variables (probesets, genes etc.) have not been filtered by 
#differential expression with respect to a sample trait. (with respect to traits we have not?)

#If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers
#(less than 15 for signed or signed hybrid networks, and less than 30 for signed networks) 
#and the mean connectivity remains relatively high (in the hundreds or above), 
#chances are that the data exhibit a strong driver that makes a subset of the 
#samples globally different from the rest. The difference causes high correlation 
#among large groups of genes which invalidates the assumption of the scale-free topology 
#approximation.

#If we want to keep the variation then according to FAQ signed networks with 20-30 samples should use power threshold of 8 

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


#*looking at the scale dependence graph, make softPower (the soft threshold) the number that is under the 0.9 value suggested
#*if our scale dependence graph doesnt look great than may need to look online/ask for help

softPower=16 #smallest value to plateau at ~0.85
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#Calculates (correlation or distance) network adjacency from given expression data or from a similarity.

##################Corinne
############################### INITIAL DENDROGRAM ##############################
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average") #do we want to change method back to default?
sizeGrWindow(10,6)
# pdf(file="dendrogram_thresh16.5_signed_1868.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
# dev.off()
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes


#*each module gets a number and a size
#*assign color to each module, color represents clusters of coexpressed genes

minModuleSize=5 #we only want large modules #set to 90 originally (How many taxa do we want? Arbitrary?)
#dependent on size of transcriptome, size should be probably diff for ours. they used 90 bc doing GO enrichment but the base is 30
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
#DeepSplit: For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. The higher the value (or if TRUE), the more and smaller clusters will be produced.

table(dynamicMods)
#5 modules 35 out of 56 taxa not in a module 
#34 modules with 39 unaccounted 

dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

################ MERGING ##############################
#Merge modules whose expression profiles are very similar, 
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 16) #*will need to change this to the soft threshold decided earlier
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_invasion_0.7.RData")

lnames = load(file = "Network_invasion_0.7.RData")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.7 #*we can change the threshold of our height, merges them different based on value, dont want to overmerge things that arent that similar
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =1)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork0.7.pdf", width=20, height=20) 
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_0.7.RData")

#############################################################
#Relating modules to traits and finding important genes #I think skip because Taxa already assigned 
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
#lnames = load(file = "Network_invasion_nomerge.RData");
#The variable lnames contains the names of loaded variables.
#lnames
# Load network data saved in the second part.

#lnames = load(file = "Network_signed_0.6.RData"); #dont load for original non-merging look 
#lnames = load(file = "Network_invasion_nomerge.RData"); #don't load
#lnames

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)
#*can pick which colors you want to further explore

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

####################### MODULE HEATMAP #####################################
#represent module trait correlations as a heatmap

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

#*when see groups of modules that are all hot or all cold, they should be merged. this will have it more likely detect enrichment
#* can change MEDissThres= 0.7 based on this heatmap
#* can see tight correlations, what percentage of variation is explained by certain relationships and look for modules that are doing the same thing

######################### RELATING MODULES TO TRAITS #######################
#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$percC); #change Lipidrobust to your trait name
names(weight) = "percC" 
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

##################### MODULE CORRELATION PLOT ####################
#*how well things belong to module
#Gene-trait significance correlation plots
# par(mfrow=c(2,3))
module = "midnightblue" #*change this to the module we're going to look at
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for %C",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

########################## VSD FILES BY MODULE ###################### 
#Making VSD files by module - so we can tell what genera are in each module
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="midnightblue"]) #*change this also to the color of module we're looking at
#black  blue brown green  grey  pink   red

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
head(c.vsd)
write.csv(c.vsd,"rlog_MMmidnightblue.csv",quote=F)

#####################KNOW WHICH GENERA ARE IN WHICH MODULE###############
###fisher for GO
##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
vsd <- read.csv("Invasion_wgcna_allgenes.csv", row.names=1)
head(vsd)
options(stringsAsFactors=FALSE)
data=t(vsd)
allkME =as.data.frame(signedKME(data, MEs))

whichModule="midnightblue" #*name your color and execute to the end

length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which genera are in each module
yellow <- subset(inModule, module=="1")
yellow <- row.names(yellow)
yellow

grey <- subset(inModule, module=="1")
grey <- row.names(grey)
grey

midnightblue <- subset(inModule, module=="1")
midnightblue <- row.names(midnightblue)
midnightblue

###################### kMEs ###################
#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the yellow kMEs for only the genera in the yellow module
ykmeinput<- paste(yellow, sep=",")
yellowkme<- subset(modkME, rownames(modkME) %in% ykmeinput)
Genus <- rownames(yellowkme)
rownames(yellowkme) <- NULL
yellowkme <- cbind(Genus,yellowkme)

#plot the kMEs to show which genera fit best into the module
ggplot(data= yellowkme, aes(x=reorder(Genus, -kMEyellow), y=kMEyellow)) + 
  geom_bar(color= "black", fill="yellow", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Genus", y="kME")

#########
#make a subset of the midnightblue kMEs for only the genera in the midnightblue module
mbkmeinput<- paste(midnightblue, sep=",")
midnightbluekme<- subset(modkME, rownames(modkME) %in% mbkmeinput)
Genus <- rownames(midnightbluekme)
rownames(midnightbluekme) <- NULL
midnightbluekme <- cbind(Genus,midnightbluekme)

#plot the kMEs to show which genera fit best into the module
ggplot(data= midnightbluekme, aes(x=reorder(Genus, -kMEmidnightblue), y=kMEmidnightblue)) + 
  geom_bar(fill="midnightblue", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Genus", y="kME")


################
#make a subset of the grey kMEs for only the genera in the grey module
gkmeinput<- paste(grey, sep=",")
greykme<- subset(modkME, rownames(modkME) %in% gkmeinput)
Genus <- rownames(greykme)
rownames(greykme) <- NULL
greykme <- cbind(Genus,greykme)

#plot the kMEs to show which genera fit best into the module
ggplot(data= greykme, aes(x=reorder(Genus, -kMEgrey), y=kMEgrey)) + 
  geom_bar(fill="grey", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Genus", y="kME")

######################## HEATMAP OF GENERA ABUNDANCE BY SAMPLE ###################
##############################heatmap of module expression with bar plot of eigengene, no resorting of samples...
#names(dis)
sizeGrWindow(8,7);
which.module="midnightblue" #*change this also to the color of module we're looking at 
#pick module of interest
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

###################### HEATMAP OF GENERA ABUNDANCE ORDERED BY TRAIT ##################
##############################heatmap of module expression with bar plot of trait of interest by sample...
#*didnt do this in class but may want to create something that ranks phenotype with gene expression (as traits increase or decrease)
#here we just have binary traits, but if you have a continuous trait this code is cool
sizeGrWindow(8,7);
which.module="midnightblue" #pick module of interest
which.trait="percC" #change trait of interest here
datTraits=datTraits[order((datTraits$percN),decreasing=T),]#change trait of interest here

trait=datTraits[, paste(which.trait)]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset
genes=genes[rownames(datTraits),]

#quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(trait, col=which.module, main="", cex.main=2,
        ylab="%C",xlab="sample")#change trait of interest here

#*how well it belongs to module
#Gene relationship to trait and important modules: Gene Significance and Module membership
allkME =as.data.frame(signedKME(t(dat), MEs))
head(allkME)
vsd=read.csv(file="rlog_MMmidnightblue.csv", row.names=1)
head(vsd)
library(pheatmap)

################### TOP 100 HEAT MAP #########################
whichModule="midnightblue" #*color change
top=100
datME=MEs
vsd <- read.csv("Invasion_wgcna_allgenes.csv", row.names=1)
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
summary(hubs)

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#quartz()
setwd("/project/bi594/Pine_invasion/Figures/")
png(file="pheatmap_midnightblue.png", width=1000, height=1500)
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"kME",sep=""))
dev.off()
setwd("/project/bi594/Pine_invasion/")


######--------------------end--------------------#######

save.image(file='Environment.4.25.21.RData')
