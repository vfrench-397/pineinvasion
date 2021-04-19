#Pine invasion project
#Loading in env. 
setwd("/project/bi594/Pine_invasion/")
load('Environment4.13.21.RData')


setwd("/project/bi594/Pine_invasion/")

load(file, envir = parent.frame(), verbose = FALSE)fns <- list.files(path)
fns

#####################################
#The following tutorial is modified from:
#https://benjjneb.github.io/dada2/tutorial.html
#with edits by Carly D. Kenkel & Alizah Ali & Nicola Kriefall & Sarah Davies

# source("https://bioconductor.org/biocLite.R")
# biocLite("dada2")
# biocLite('vegan')
#####################################
library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

#Set path to unzipped, renamed fastq files

setwd("/project/bi594/Pine_invasion/rawdata/")
path <- "/project/bi594/Pine_invasion/rawdata/"

fns <- list.files(path)

#Let's make sure that all of our files are there
fns

################################
##### Trimming/Filtering #######
################################

fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R1", fastqs)]
fnRs <- fastqs[grepl("_R2", fastqs)]
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #getting sample names 

#########Visualize Raw data

#First, lets look at quality profile of R1 reads
##visualize quality score - this is where you can decide to visually create a quality cutoff of "X"
##Quality scores largely drop below 20 around 250 reads in the forward reads, and around 200 in the reverse reads
##x axis "cycle" = bp length (reads = no. basepairs)
plotQualityProfile(fnFs[c(1:9)])
plotQualityProfile(fnFs[c(10:18)])
plotQualityProfile(fnFs[c(19:27)])
plotQualityProfile(fnFs[c(28:35)])
plotQualityProfile(fnFs[c(36:44)])
plotQualityProfile(fnFs[c(45:53)])
plotQualityProfile(fnFs[c(54:59)])

plotQualityProfile(fnRs[c(1:9)])
plotQualityProfile(fnRs[c(10:18)])
plotQualityProfile(fnRs[c(19:27)])
plotQualityProfile(fnRs[c(28:35)])
plotQualityProfile(fnRs[c(36:44)])
plotQualityProfile(fnRs[c(45:53)])
plotQualityProfile(fnRs[c(54:59)])

#Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(path, "Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Find out if we need to remove primer
primerF<-grep("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "RV238_S210_L001_R1_001.fastq.gz")
primerF
#integer(0) - primer was not found, has already been filtered! Also confirmed in the SCC terminal 
primerR<-grep("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "RV238_S210_L001_R2_001.fastq.gz")
primerR

# Filter forward reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen= c(250,200), #end of single end reads = approx. 250 bp
                     maxN=0, #DADA does not allow Ns - filter out Ns (degenerate base)
                     maxEE= c(2,2), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     rm.phix=TRUE, #remove reads matching phiX genome
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

# reads.out = Number of reads remaining after filtering with quality scores 
#reads.out look good. Retaining more than 50% reads 
head(out)
tail(out)
#Must filter reads together to prevent mis-match sorting that can affect merging of forward and reverse reads 

#Filter reverse reads #
#outR <- filterAndTrim(fnRs, filtRs, truncLen= 200, #end of single end reads = approx. 300 bp
                      #maxN=0, #DADA does not allow Ns - filter out Ns (degenerate base)
                      #maxEE=1, #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                      #truncQ=2, 
                      #rm.phix=TRUE, #remove reads matching phiX genome
                      #compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE




#A word on Expected Errors vs a blanket quality threshold
#Take a simple example: a read of length two with quality scores Q3 and Q40, corresponding to error probabilities P=0.5 and P=0.0001. The base with Q3 is much more likely to have an error than the base with Q40 (0.5/0.0001 = 5,000 times more likely), so we can ignore the Q40 base to a good approximation. Consider a large sample of reads with (Q3, Q40), then approximately half of them will have an error (because of the P=0.5 from the Q2 base). We express this by saying that the expected number of errors in a read with quality scores (Q3, Q40) is 0.5.
#As this example shows, low Q scores (high error probabilities) dominate expected errors, but this information is lost by averaging if low Qs appear in a read with mostly high Q scores. This explains why expected errors is a much better indicator of read accuracy than average Q.

################################
##### Learn Error Rates #######
################################
#DADA2 learns its error model from the data itself by alternating estimation of the error rates and the composition of the sample until they converge on a jointly consistent solution (this is similar to the E-M algorithm)
#As in many optimization problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

setDadaOpt(MAX_CONSIST=30) #set higher to allow more cycles for convergence
errF <- learnErrors(filtFs, multithread=TRUE)
#101188750 bases in 404755 reads from 9 samples 
#output is list of 3, does this mean convergence after 3 rounds? 
errR <- learnErrors(filtRs, multithread=TRUE)
#105691000 bases in 528455 reads from 12 samples 
#Maximum cycles was set to 30, but Convergence was found after 4 rounds
#errF may take a long time 

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#error rates look relatively okay, not strictly linear but definitely a negative correlation 

################################
##### Dereplicate reads #######
################################
#Dereplication combines all identical sequencing reads into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)
# Name the derep-class objects by the sample names
names(derepFs) <-sample.names
names(derepRs) <-sample.names
################################
##### Infer Sequence Variants #######
################################

#Performing Joint sample inference 
#Must change some of the DADA options b/c original program optomized for ribosomal data, not ITS - from github, "We currently recommend BAND_SIZE=32 for ITS data." leave as default for 16S/18S
setDadaOpt(BAND_SIZE=32) #Keep set at 32 for ITS analysis 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]] #looking at number of sequence variants found; 
#158 sequence variants inferred from 6180 input sequences 
dadaRs[[1]]
#149 sequence variants inferred from 6460 input sequences 

#merging paired ends of forward and reverse reads for full sequences 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=FALSE)
head(mergers[[1]])
#originally not successfully merging, $maps were different lengths because sorted separately 
#maybe ask Sarah and James about this output? 


#construct sequence table
seqtab <- makeSequenceTable(mergers) #INPUT for WGCNA? 
dim(seqtab) #4887 sequence variants? 
head(seqtab) #gives sequence, then number of that sequence found in each sample

################################
##### Remove chimeras #######
################################
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# Identified 16 bimeras out of 1871 input sequences (~1.3%)
table(nchar(getSequences(seqtab.nochim))) #distribution of sequence lengths 

sum(seqtab.nochim)/sum(seqtab) #.9993371
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)
#For our sample, this ratio was 0.9998201, there was only 1 bimera

setwd('/project/bi594/Pine_invasion/')
write.csv(seqtab,file="pineinvasion_seqtab.csv")
write.csv(seqtab.nochim,file="pineinvasion_nochim.csv")

################################
##### Track Read Stats #######
################################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names 
head(track)
tail(track)

write.csv(track,file="ReadFilterStats_AllData_final.csv",row.names=TRUE,quote=FALSE)


################################
##### Assign Taxonomy #######
################################

#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. 
#DADA2 provides a native implementation of the RDP's naive Bayesian classifier. The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs the taxonomic assignments with at least minBoot bootstrap confidence.
#Here, I have supplied a modified version of the GeoSymbio ITS2 database listing more taxonomic info as phyloseq requires (Franklin et al. 2012)
#For example: GeoSymbio data (taken from "all clades" at https://sites.google.com/site/geosymbio/downloads):
#>A1.1
#modified version for phyloseq looks like this instead:
#>Symbiodinium; Clade A; A1.1

setwd('/project/bi594/Pine_invasion/')

taxa <- assignTaxonomy(seqtab.nochim, "/project/bi594/Pine_invasion/sh_general_release_dynamic_s_04.02.2020.fasta", multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL #this gets rid of all the sequences
#head(taxa) #prints full info, including sequence
head(taxa.print) #this just prints taxonomy; look like fungi? 
#minboot should be higher
#Obtain a csv file for the taxonomy so that it's easier to map the sequences for the heatmap.

unname(head(taxa, 30))
unname(taxa)

taxa<- sub('...', '', taxa)
write.csv(taxa, file="taxa.csv",row.name=TRUE,quote=FALSE)

#Now, save outputs so can come back to the analysis stage at a later point if desired
saveRDS(seqtab.nochim, file="final_seqtab_nochim.rds")
saveRDS(taxa, file="final_taxa_blastCorrected.rds")

#If you need to read in previously saved datafiles
seqtab.nochim <- readRDS("final_seqtab_nochim.rds")
taxa <- readRDS("final_taxa_blastCorrected.rds")



#everything above this was in dada2, now we're using phyloseq
################################
##### handoff 2 phyloseq #######
################################

setwd('/project/bi594/Pine_invasion/')

library('phyloseq')
library('Biostrings')
library('ggplot2')

#import dataframe holding sample information
#have your samples in the same order as the seqtab file in the rows, variables as columns
samdf<-read.csv("variabletable_pi.csv", header = TRUE, sep = ',')
rownames(samdf) <- samdf$SAMPLE #making rownames the same as sample names in seq.nochim to merge in phyloseq


#Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_names(samdf), 
               tax_table(taxa))
ps

#Collapse the data in the phyloseq object to genus-level IDs
glom <- tax_glom(ps, taxrank = 'Genus')
#Create a dataframe of the OTUs and taxonomic assignments at the genus level
otu<- data.frame(otu_table(glom))
tax_table <- data.frame(tax_table(glom))
#Replace the sequences in the OTU table with corresponding Genus-level IDs
Genus<- tax_table$Genus
colnames(otu) <- Genus
#the new OTU df will be the input for WGCNA
write.csv(otu,file="otu.csv")


#replace sequences with shorter names (correspondence table output below)
#ids<-taxa_names(ps)
#ids <- paste0("sq",seq(1, length(colnames(seqtab.nochim))))
#colnames(seqtab.nochim) <- ids

#Select the top 90 most abundant taxa
top90 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:90]
ps.top90 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top90 <- prune_taxa(top90, ps.top90)

#Obtain a csv file for the phyloseq data. Will give you the abundances for each sample and class. Useful for constructing the heatmap. Also, enables you to use ggplot, and construct more aesthetically pleasing plot.
psz <- psmelt(ps)
write.csv(psz, file="psz.csv")
psz90 <- psmelt(ps.top90)
write.csv(psz90, file="top90.csv")

#Functional assignments
##################### 
#Assign functional groups in Funguild and Fun^fun

#Format the input table so it can be run in Funguild
library(data.table)

funguild_input<- data.frame(taxa)
funguild_input$OTU <- row.names(funguild_input)  
funguild_input$Taxonomy <- paste(funguild_input$Kingdom, funguild_input$Phylum, funguild_input$Class, funguild_input$Order, funguild_input$Family, funguild_input$Genus, funguild_input$Species, sep =';')
funguild_input<- subset(funguild_input, select=-c(Kingdom, Phylum, Class, Order, Family, Genus, Species))

devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")
devtools::install_github("brendanf/FUNGuildR")

library(datastorr)
library(fungaltraits)
library(FUNGuildR)
library(jsonlite)

guilds<- funguild_assign(funguild_input)

############################


#Append the metadata to the phyloseq data
colnames(psz)[2] <- "SampleID" #rename the sample ID column so we can merge the two dataframes by this column
fulldf <- merge(psz, samdf, by="SampleID")
#do the same for the top 90 taxa
colnames(psz90)[2] <- "SampleID" #rename the sample ID column so we can merge the two dataframes by this column
fulldf90 <- merge(psz90, samdf, by="SampleID")

#Barplot of fungal abundance by forest type separated by Class
p <- ggplot(fulldf, aes(x = site_code, y=Abundance, fill=Class))+ 
  labs(x="Forest type", fill = "Class")
p + geom_bar(stat="identity", colour="black") +
  scale_x_discrete(labels=c('Invaded forest', 'Plantation', "Native forest"))

ggsave("Abundance.png", path = "/project/bi594/Pine_invasion/Figures/", width=10, height=6, dpi=300)

#Add guild annotations onto the abundance data
fullguild<- merge(fulldf, guilds, by="OTU")

p2 <- ggplot(fullguild, aes(x = site_code, y=Abundance, fill=trophicMode))+ 
  labs(x="Forest type", fill = "Class")
p2 + geom_bar(stat="identity", colour="black") +
  scale_x_discrete(labels=c('Invaded forest', 'Plantation', "Native forest"))



save.image(file='Environment.4.18.21.RData')
