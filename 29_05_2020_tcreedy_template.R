
# Introduction ------------------------------------------------------------

# By Thomas J. Creedy, 2020-04-03

# This script is a template for analysing metabarcoding data generated
# using the current Vogler Lab metabarcoding pipeline.

# I strongly recommend using these commands as a guide only and ensuring
# you understand what each command does and how to modify it to answer your
# specific questions.

# Not everything here will work for all data. Please try and figure out how
# to adapt these commands to your data - the purpose of this script is not
# to be a comprehensive universal analysis tool but to be a template to get
# you started.

# Setup -------------------------------------------------------------------
# Commands in this section clean up from any previous runs and implements
# any global settings.

# Remove all objects from the environment
rm(list = ls())

# Set global options
options(stringsAsFactors = F)


# Load libraries ----------------------------------------------------------
# Load the libraries (packages) you need for your data manipulation and 
# analysis. Remember, the order of loading is sometimes important.
library(stringr)
library(plyr)
library(dplyr)
library(reshape2)
library(vegan)
library(ggplot2)
library(tidyr)
library(picante)
library(ggpubr)
library(RAM)
library(BiodiversityR)
library(iNEXT)


# Load functions ----------------------------------------------------------
# Load any custom functions. These are assumed to be in the root of your
# working directory.

source("check_expected_richness.R")
source("rarexplore.R")

# Load data ---------------------------------------------------------------
# Load the three core data tables. If you have more than these data,
# great! Add loading commands for them here. I would generally suggest
# merging and/or reconfiguring them as necessary to be similar to the
# configuration of these tables.


# Load the reads table output by vsearch --usearch_global --tabbedout
reads_raw <- read.table("3pc_reads_mapV2.tsv", header = T, sep = "\t", comment.char = '')
#must remove all samples with FP in name as these are barcodes not MBC 
reads<-select(reads_raw, -matches(".FP."))
#make column 1 the row names and delete column 1 (could this be done before hand?)
MBC_reads <- data.frame(reads[,-1], row.names=reads[,1])
rm(reads)

# Load metadata table (change command if yours is a different format!)
# Note that here I'm assuming that the first column of the metadata
# table is the sample names.

metadata_raw <- read.csv("treedata_matt.csv", header=T)
metadata <- read.csv("treedata_matt.csv", header=T)

# If the sample names are in a different column, do this:

#metadata <- read.csv("metadata.csv")
#row.names(metadata) <- metadata$samplenames
#metadata <- subset(metadata, select = -samplenames)

# Load taxonomy table

# From SINTAX
taxonomy <- read.table("otus_midori_sintax.tsv", sep = "\t", row.names = 1)
colnames(taxonomy) <- c("taxonomy", "strand", "selected")

# Organise data -----------------------------------------------------------
# Here we make sure all the data corresponds properly to each other and 
# any reconfiguration or merging is done as needed.

# Generally, samples should be rows and species should be columns, so
# we transpose the usearch_global --tabbedout reads table

MBC_reads <- t(MBC_reads)

# If we have multiple metabarcoding samples that actually represent the 
# same ecological sample, here we might merge them together
#The rows from the same tray (i.e. same TX-YY) should be summed together
#All t1.1 should be summed, all t1.2 should be summed etc

#using gsub we can create a new character list that chops out all of _P0x using a Regex expression
#pattern - "_P[^_]* means match _ and P, ^match any character after P that is not in the set i.e 0 or higher number
#* means anything after match nothing after the last _
#the rownames MBC_reads takes the data from the row names in MBC_reads
samples_tree_trap<-gsub(pattern = "_P[^_]*", '', rownames(MBC_reads)) 
#the character list and puts it into a table that allows us to look at how many there are of each
#there should be 3's through out other than a the cocasional 2 for OTU delimitation calibration
samplenames<-table(samples_tree_trap)
#this summates the rows by sample_tree_trap, there all traps have been put together into their respective rows
MBC_reads<-rowsum(MBC_reads, group = samples_tree_trap)

#sapply takes the summated rows and divides each row by the number of rows that went in and then multiples by 3, returning a rounded number to 0
#it also transposes it at the end
MBC_reads<-sapply(samplenames, function(n){
 x<-(MBC_reads[n,]/samplenames[n])*3
 return(round(x,0))
}) %>% t

#The tree for each row should be identified, and this used to duplicate out the metadata table to have a row of the correct tree for each sample.
#extract rownames in mbc_reads, that meet the _T[^-]* requirements (_T and below, aka the tree number and trap number)
treenumbers<-str_extract(rownames(MBC_reads), "_T[^-]*") %>% gsub(pattern = "_", "", .)
#creates a dataframe of tree code and tree number next to it (column names)
temp <- data.frame(sample = rownames(MBC_reads), tree = str_extract(rownames(MBC_reads), "_T[^\\.]*") %>% gsub(pattern = "_T", "", .))
#changes tree value to an integer 
temp$tree<-as.integer(temp$tree)
#removes na values from tree column
temp <- temp[!is.na(temp$tree), ]
#merging extracted objects from temp into metadata
#by treeid and tree, all.x =T will add extra rows to the output, same with all.y but false in this case
metadata<-merge(metadata, temp, by.x = "treeid", by.y = "tree", all.x = T, all.y = F)
#make rownames the sample id
row.names(metadata) <- metadata$sample
#remove sample column as it is now the rownames
metadata <- subset(metadata, select = -sample)
#do something to colnames that do not equal "sample" confused on this one - TJC: this is exactly the same as the previous one, we
#just added it because the previous didn't seem to be working right. 
metadata <- metadata[, colnames(metadata)[! colnames(metadata) == "sample"]]
#should i remove the above line?

# Next, check that our metadata and reads table correspond. For both 
# tables, the sample names are in the row.names

# Are all MBC samples in metadata? If not, problem!!!
all(row.names(MBC_reads) %in% row.names(metadata))

# Make the metadata correspond to the reads. Takes rownames of metadata and puts it into the rownames of 
#MBC_reads
MBC_reads <- MBC_reads[row.names(metadata), ]
#to check
all(row.names(MBC_reads) %in% row.names(metadata))

# This keeps only metadata rows where the sample occurs in MBC, and 
# sorts the rows to match the MBC row order

###REORGANISE TAXONOMY###

# The data in the taxonomy table needs to be filtered and separated out
# As above
all(colnames(MBC_reads) %in% row.names(taxonomy))
taxonomy <- taxonomy[colnames(MBC_reads),]

# We create a detailed version of the taxonomy table with the scores
# This code is based on Yige Sun's rewriting of this step, thanks!
    # Add the OTUs as a column
taxonomy$otu <- row.names(taxonomy)
    # Separate the taxonomy column into multiple columns by the ',' character
taxdetailed <- separate(taxonomy, taxonomy, paste0('V', 1:(max(str_count(taxonomy$taxonomy, ','))+1)), 
                        sep = ",", fill = 'right', remove = T) %>% 
    # Turn the separate columns for taxa into rows
  pivot_longer(., cols = starts_with('V'), values_to = 'taxon', values_drop_na = T) %>%
    # Retain only the OTU and taxon columns
  select(otu, taxon) %>%
    # Separate the taxon column into three columns by colon and parentheses
  separate(., taxon, into = c('level', 'taxon', 'score', NA), sep = "[:()]")
  
# Record which taxa were selected
taxdetailed$selected <- taxdetailed$score == 1

# Rename the taxonomic levels
levelnames <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
names(levelnames) <- substr(levelnames, 1, 1)
#rewrites taxonomic level names
taxdetailed$level <- factor(levelnames[taxdetailed$level], levels = levelnames)

# Now we overwrite the original taxonomy table with the detailed one (creates a table of just 1s?????)
taxonomy <- dcast(taxdetailed, otu ~ level, value.var = "taxon")

# We may only wish to include taxonomy info if it's above a certain score
#taxonomy <- dcast(taxdetailed[taxdetailed$score > 0.65, ], otu ~ level, value.var = "taxon", fill = '')

# Set the row names
row.names(taxonomy) <- taxonomy$otu
taxonomy <- subset(taxonomy, select = -otu)    

# Finally, re-sort the taxonomy
taxonomy <- taxonomy[colnames(MBC_reads), ]
#taxonomy throughout for all OTUs... very nice
rm(taxdetailed)

# Filtering data ----------------------------------------------------------
# This is the most important step, where you remove irrelevant or likely
# incorrect data. Depending on your questions you could apply none, some
# or all of these steps.

#list orders of OTUs
orders <- unique(taxonomy$order)
# Taxonomic filtering - retain only OTUs of a certain taxon, e.g.:
taxonomy <- taxonomy[taxonomy$order == "Coleoptera", ]
#list arthropod orders
arthro_orders <- unique(taxonomy$order)
#
taxonomy <- taxonomy[taxonomy$order %in% arthro_orders, ]
# 1834 OTUs here in taxonomy
taxonomy<-na.omit(taxonomy)
# 1795 OTus here in taxonomy, 39 OTUs were na
MBC_reads <- MBC_reads[, rownames(taxonomy)]
#to 1795 OTUs here, 8309 (orginal before filtering to coleoptera) - 1795 = 6514 OTUs which are dropped
# i.e not coleoptera

# Read number filtering 
# Set counts less than a certain value to 0 - this isn't recommended
# because it doesn't take into account total reads per sample. It also
# generates difficulties with some of the accumulation tests later

#MBC_reads[MBC_reads < 2] <- 0

# Read proportion filtering
# Set counts that are less than a certain proportion of reads in a 
# sample to 0. The threshold here is 0.5%

threshold <- 0.0005
#MTH works through rows, if s is less than threshold, where s is equal to the s/sum of s, the value is given a 0
MBC_reads <- apply(MBC_reads, 1, function(s){
  s[s/sum(s) < threshold] <- 0
  return(s)
}) %>% t()
rm(threshold)
#1795 OTUs still

# After filtering you may have some OTUs or samples that now have 0 reads
# so these should be filtered out and the relevant lines in the metadata/
# taxonomy dropped as well

#MBC_reads <- MBC_reads[rowSums(MBC_reads) > 0, colSums(MBC_reads) > 0]
#Maybe this line needs editing. Go from 1795 to 63, that means there are 1733 OTUs with 0 reads
#metadata <- metadata[rownames(MBC_reads), ]
#taxonomy <- taxonomy[colnames(MBC_reads), ]

#BIG PROBLEM HERE, loose thousands of OTUs.

#MTH Which filtering mechanism shall I use? which would be the most appropiate for my question? 
#Should I retain max community composition or min? or perhaps somewhere in the middle?

#  Standardisation --------------------------------------------------------
# It is crucial with metabarcoding that you ensure your samples are 
# comparable, otherwise your analyses may be invalid. There are two ways
# we can consider whether samples are standard:
# They have received equivalent sampling depth
# They have accumulated a representative sample of a community
# These assumptions may vary depending on your question

#Which standarisation method should I use? which would be most appropriate for my question.

# Check the read numbers - have all samples been equally sequenced?
# Probably not!
#hist(rowSums(MBC_reads))

# Check the accumulation - have all samples recovered a good proportion
# of their expected OTU richness?
# Here we set the threshold proportion of recovery to 85%. This is 
# probably the lowest this should be.
jpeg("check_expectedrich.jpeg", width =480, height=480)
check <- check_expected_richness(MBC_reads, 0.85)
dev.off()

# For most community ecology, the best route is to remove samples that 
# did not reach the threshold. check_expected_richness returns a vector
# of samples that passed this as the first item

MBC_reads <- MBC_reads[check[[1]], ]

# An alternative to this standardisation is to rarefy. This process
# removes reads from samples until it reaches a threshold. This 
# standardises sampling depth but at the loss of real community data.
# It is more suitable for methodological analyses and is not recommended
# for community ecology.

# First lets review how many samples/OTUs might be lost at different 
# rarefaction target values:
#target_values = round(seq(min(rowSums(MBC_reads)), max(rowSums(MBC_reads)), length.out = 30))
#jpeg("rarexplore.jpeg", width =480, height=480)
#rarexplore(MBC_reads, target_values)
#dev.off()

# Pick a value that is a good trade-off of remaining samples and
# remaining OTUs, say 2500, and drop samples with below this value
#MBC_reads <- MBC_reads[rowSums(MBC_reads) >= 2500, ]
#MBC_reads <- rrarefy(MBC_reads, 2500)

# The most crucial part of standardisation is to convert your data into
# presence/absence. Read numbers are not true abundances and are rarely
# comparable between OTUs and are only comparable between samples if
# rarefaction has been done (which is rarely appropriate). The only time
# when you may retain read numbers is if you are interested in looking 
# at variation between one or a small number of similar OTUs in rarefied 
# samples.
MBC_reads[MBC_reads > 0] <- 1

# Finally, we need to clean up after standardisation
# Remove any OTUs which no longer have any reads
MBC_reads <- MBC_reads[, colSums(MBC_reads) > 0]
#474
# Remove any dropped samples/OTUs from metadata/taxonomy tables
metadata <- metadata[rownames(MBC_reads), ]
taxonomy <- taxonomy[colnames(MBC_reads), ]

#MTH for coleoptera using the check expected richness, MBC_reads has 63 OTUs
#______________________________________________________________
#EXPLORING DATA

#OPTIONAL create melted dataframe 
#MBC_reads_melt<-melt(MBC_reads)

#examine basic diversity indices
#OTU richness
specnumber(MBC_reads, groups = metadata$species)
specnumber(MBC_reads, groups = metadata$treeid)
specnumber(MBC_reads, groups = metadata$camp)
specnumber(MBC_reads, groups = metadata$GPS_alt)

#beta diveristy, needs fixing
#betadiver(MBC_reads, method = "j", groups = metadata$species)

#bar graphs of specnumber

#NMDs
reads_dist<-vegdist(MBC_reads, method = "jaccard")
NMDS.scree <- function(x){
  plot(rep(1,10), replicate(10, metaMDS(x, autotransform = F, k=1, trymax = 100)$stress), xlim = c(1,10),
       ylim = c(0, 0.30), xla ="# of Dimensions", ylab = "Stress", main ="NMDS stress plot")
  for(i in 1:10){
    points(rep(i +1,10), replicate(10, metaMDS(x, autotransform = F, k=i+1, trymax = 100)$stress))
  }
}

#run nmds 
jpeg("NMDsscreeplot.jpeg", width =480, height=480)
NMDS.scree(reads_dist)
dev.off()
#insuffient data, how have i lost thousands of OTUs?? I go from 8000 otus to 64, something wrong here
#pca on habitatvariables
#princomp(metadata)

#to build accumulation curves for each altitude 

#library(ggpubr)
#graph the altitudes of different trees, just exploration
#alt_bar<-ggplot(data=metadata, aes(x=treeid, y=GPS_alt, fill=GPS_alt)) +
  #geom_bar(stat="identity", position=position_dodge())+labs(x="tree number", y="tree altitude")+
  #scale_fill_continuous(name="trees and their altitude")
#alt_bar
#Chao estimates of species richness and bargraph for each tree
otuspecpooltree<-specpool(MBC_reads, metadata$treeid)
otuspecpooltree["treeid"]<-metadata_raw$treeid
otuspecpooltree["camp"]<-metadata_raw$camp
otuspecpooltree["scale_treealt"]<-scale(metadata_raw$GPS_alt)
otuspecpooltree["treealt"]<-metadata_raw$GPS_alt

#otuspecpooltreeplot<-(otuspecpooltree, aes(x=treeid, y=Species)) +
  #geom_bar(stat="identity", fill="steelblue")
#plot(otuspecpooltreeplot)

ggplot(otuspecpooltree, aes(x=treeid, y=chao)) + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=chao-chao.se, ymax=chao+chao.se))
#jpeg("chaotreeidbar", width = 480, height= 480, res = 300)

#Chao estimates of species richness and bargraph for each camp 
otuspecpoolcamp<-specpool(MBC_reads, metadata$camp)
otuspecpoolcamp["camp"]<-c("BC", "CA", "CP", "GU")

jpeg("chaobycampbar.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75, bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
ggplot(otuspecpoolcamp, aes(x=reorder(camp, -chao), y=chao)) + geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=chao-chao.se, ymax=chao+chao.se))
dev.off()

#alt_hist <-hist(metadata_raw$GPS_alt)
ggplot(otuspecpooltree, aes(x=scale_treealt, y=Species)) + geom_point() +geom_smooth(method=lm, se=FALSE)

jpeg("altitiude_OTU.jpeg", width =480, height=480)
ggplot(otuspecpooltree, aes(x=treealt, y=Species)) + geom_point() +geom_smooth(method=lm, se=F)
dev.off()

#altmod1<-glm(Species~treealt, family = poisson, data = otuspecpool)
#summary(altmod1)

#specaccume using vegan package
#Per tree spcies
#liquid amber
liquidamber_accum <- specaccum(MBC_reads, method = "random", permutations = 1000, subset=metadata$species=="Liquidambar styracaflua")
jpeg("liquidamber_accumVegan.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75,bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
plot(liquidamber_accum, ci.type = "poly", col = "blue", ci.col = "lightblue", 
     lwd = 2, ci.lty = 0, xlab = "number of traps", 
     ylab = "cumulative number of OTUs")
dev.off()
#pinus sp.
pinus_accum <- specaccum(MBC_reads, method = "random", permutations = 1000, subset=metadata$species=="Pinus sp.")
jpeg("pinus_accumVegan.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75,bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
plot(pinus_accum, ci.type = "poly", col = "blue", ci.col = "lightblue", 
     lwd = 2, ci.lty = 0, xlab = "number of traps", 
     ylab = "cumulative number of OTUs")
dev.off()

#___________________

#biodiveristyR package species accumulation curves using biodiveristyR package

metadata$trap.totals <- apply(MBC_reads,1,sum)
Accum.1 <- accumresult(MBC_reads, y=metadata, scale='trap.totals', method='exact', conditioned=TRUE, permutations = 1000)

#Per camp
jpeg("camp_accumBioDivR.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75,bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
accumcomp(MBC_reads, y=metadata, factor='camp', method='exact', legend=FALSE, conditioned=TRUE, xlab="trap")
dev.off()
#per tree species
jpeg("treesp_accumBioDivR.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75,bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
accumcomp(MBC_reads, y=metadata, factor='species', method='exact', legend=FALSE, conditioned=TRUE, xlab="trap")
dev.off()
#per altitudes
jpeg("peralt_accumBioDivR.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75,bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
accumcomp(MBC_reads, y=metadata, factor='GPS_alt', method='exact', legend=FALSE, conditioned=TRUE, xlab="trap")
dev.off()

#____________________________________
#iNEXT accumulation 
#Frequency raw

#need to split mbc_reads by camp into separate matrixs
#could use grep here instead

BC<-subset(MBC_reads,metadata$camp=="BC" )
BC_mat <- t(as.matrix(apply(BC,2,as.integer)))

CA<-subset(MBC_reads,metadata$camp=="CA" )
CA_mat <- t(as.matrix(apply(CA,2,as.integer)))

GU<-subset(MBC_reads,metadata$camp=="GU" )
Gu_mat <- t(as.matrix(apply(GU,2,as.integer)))

CP<-subset(MBC_reads,metadata$camp=="CP" )
CP_mat <- t(as.matrix(apply(CP,2,as.integer)))

camp_list = list(BC = BC_mat, CA = CA_mat, GU = Gu_mat, CP = CP_mat)

q0camp <- iNEXT(camp_list, q=0, datatype="incidence_raw", endpoint=500)
jpeg("q0campaccum_iNEXT.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75,bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
ggiNEXT(q0camp, facet.var = "none", type = 1)
dev.off()

#q0_1_2camp <- iNEXT(camp_list, q=c(0,1,2), datatype="incidence_raw", endpoint=500)
#ggiNEXT(q0_1_2camp, facet.var = "order", type = 1)






