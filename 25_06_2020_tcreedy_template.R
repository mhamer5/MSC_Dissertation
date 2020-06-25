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
library(tibble)
library(ape)
library(geosphere)
library(RAM)


# Load functions ----------------------------------------------------------
# Load any custom functions. These are assumed to be in the root of your
# working directory.

source("check_expected_richness.R")
source("rarexplore.R")
source("plot_adipart.R")
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
#TJC: yes, of course, but the idea is that this script reads the data from metabarcoding
# directly without having to do any manual or other adjustment to the data.
MBC_reads <- data.frame(reads[,-1], row.names=reads[,1])
rm(reads)

# Generally, samples should be rows and species should be columns, so
# we transpose the usearch_global --tabbedout reads table

MBC_reads <- t(MBC_reads)
# TJC: given you dropped samples above, you should also make sure to drop any OTUs
#  that have no reads now. 
MBC_reads <- MBC_reads[, colSums(MBC_reads)>0]


# Load metadata table (change command if yours is a different format!)
# Note that here I'm assuming that the first column of the metadata
# table is the sample names.

metadata_raw <- read.csv("treedata_matt.csv", header=T)
metadata_raw <- filter(metadata_raw, species!="Pinus sp.")
metadata <- read.csv("treedata_matt.csv", header=T)
climband<-c("2020", "2000", "2000", "2000", "2080", "2080", "2080","2020","2020","2020","2020","2020", "2000","2000","2000","2080","2080")
metadata<-cbind(metadata, climband)
# TJC: what's the point in having two objects with identical data?
# MTH: I use metadata_raw later on in exploratory stats

# If the sample names are in a different column, do this:

#metadata <- read.csv("metadata.csv")
#row.names(metadata) <- metadata$samplenames
#metadata <- subset(metadata, select = -samplenames)

# Load taxonomy table

# From SINTAX
taxonomy <- read.table("otus_midori_sintax.tsv", sep = "\t", row.names = 1)
colnames(taxonomy) <- c("taxonomy", "strand", "selected")


# Filtering low abundance OTU incidences ----------------------------------

# Read number filtering 
# Set counts less than a certain value to 0 - this isn't recommended
# because it doesn't take into account total reads per sample. It also
# generates difficulties with some of the accumulation tests later

#MBC_reads[MBC_reads < 2] <- 0

# Read proportion filtering
# Set counts that are less than a certain proportion of reads in a 
# sample to 0. 
# First generate a filtering matrix, which gives the proportion of 
# total reads in a sample for each OTU.
samplewise_p <- apply(MBC_reads, 1, function(s) s/sum(s) ) %>% t()

# Could then explore this matrix. At the simplest level, produce a
# histogram showing the distribution of proportion of per-OTU sample
# incidences that exceed a threshold value (here 0.5%, i.e. 0.005).
apply(samplewise_p, 2, function(o) sum(o >= 0.005) / sum(o > 0)) %>% hist()

# Could also find the maximum proportion for each OTU
maxprops <- apply(samplewise_p, 2, function(o) max(o)) %>% sort()

# Can use this to show how many OTUs would be lost completely for 
# different thresholds (because maxprops is sorted)
# E.g. what value would retain 95% of OTUs:
maxprops[floor((1 - 0.95) * length(maxprops))]
# E.g. what proportion of OTUs would be dropped if drop incidences less than 0.5%
sum(maxprops <= 0.000125) / length(maxprops)
#0.02392519

# An important consideration is that no sample should be left without any OTUs
# that occur only once, otherwise estimated richness cannot be calculated.
# Given that the minimum >0 proportion of any OTU within a sample is going to be
# for OTUs that only have 1 read, we can find the maximum threshold that fits
# this consideration
apply(samplewise_p, 1, function(s) min(s[s>0])) %>% min()
#8.805142e-06

# If this value is too low, you could pick a higher value and throw away samples 
# later. This isn't ideal though...

# So, pick a threshold. This is very data-dependent. Think carefully
threshold <- 0.000125

MBC_reads[samplewise_p < threshold] <- 0
MBC_reads <- MBC_reads[rowSums(MBC_reads) > 0, colSums(MBC_reads) > 0]
rm(samplewise_p)
#rm(threshold) ?

# Organise data -----------------------------------------------------------
# Here we make sure all the data corresponds properly to each other and 
# any reconfiguration or merging is done as needed.

# TJC: you need to drop the two samples that have "MT" in their name -
# don't worry, they would previously have been dropped later anyway but
# it now needs to be done here.
MBC_reads <- MBC_reads[! grepl("_MT", row.names(MBC_reads)), ]

# If we have multiple metabarcoding samples that actually represent the 
# same ecological sample, here we might merge them together
#The rows from the same tray (i.e. same TX-YY) should be summed together
#All t1.1 should be summed, all t1.2 should be summed etc

#using gsub we can create a new character list that chops out all of _P0x using a Regex expression
#pattern - "_P[^_]* means match _ and P, ^match any character after P that is not in the set i.e 0 or higher number
#* means anything after match nothing after the last _
#the rownames MBC_reads takes the data from the row names in MBC_reads
samples_tree_trap <- gsub(pattern = "_P[^_]*", '', rownames(MBC_reads)) 

# TJC: This section seems to be missing a part, it's not properly grouping 
# over trays since samples with different parts (i.e all, bulk, 
# coleoptera non-ms) are still in here. I can't remember now whether I gave
# you code for this or left it to you, in any case it can be done with the following:

samples_tree_trap <- gsub('^.*_T', 'T', samples_tree_trap)

#the character list and puts it into a table that allows us to look at how many there are of each
#there should be 3's through out other than a the occasional 2 for OTU delimitation calibration
samplenames<-table(samples_tree_trap)
#this summates the rows by sample_tree_trap, there all traps have been put together into their respective rows
MBC_reads<-rowsum(MBC_reads, group = samples_tree_trap)

# The below part can be removed now. The standardisation was only important for the filtering to be 
# appropriate, but now that the filtering has already happened we can ignore this.

#sapply takes the summated rows and divides each row by the number of rows that went in and then multiples by 3, returning a rounded number to 0
#it also transposes it at the end
#MBC_reads<-sapply(samplenames, function(n){
# x<-(MBC_reads[n,]/samplenames[n])*3
# return(round(x,0))
#}) %>% t

#The tree for each row should be identified, and this used to duplicate out the metadata table to have a row of the correct tree for each sample.
#extract rownames in mbc_reads, that meet the _T[^-]* requirements (_T and below, aka the tree number and trap number)

# TJC: this line is completely redundant, the object it creates isn't used anywhere.
#treenumbers<-str_extract(rownames(MBC_reads), "_T[^-]*") %>% gsub(pattern = "_", "", .)

#creates a dataframe of tree code and tree number next to it (column names)
# TJC: removed the '_' from the regex to fit with my added filter above line 121
temp <- data.frame(sample = rownames(MBC_reads), tree = str_extract(rownames(MBC_reads), "T[^\\.]*") %>% gsub(pattern = "T", "", .))
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

# Next, check that our metadata and reads table correspond. For both 
# tables, the sample names are in the row.names

# Are all MBC samples in metadata? If not, problem!!!
all(row.names(MBC_reads) %in% row.names(metadata))

# TJC: it doesn't take them and put them in. It takes the MBC_reads table and sorts the rows of that table
# according to the rows of the metadata table, then it overwrites the MBC_reads table.
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

#Remove Pinus trees
metadata<-metadata %>%
  rownames_to_column('sample') %>%
  filter(., species != "Pinus sp.") %>%
  column_to_rownames('sample')
MBC_reads <- MBC_reads[rownames(metadata), ]
# TJC: as always, when you remove samples, you should check to see if any
#  otus should now be removed too!
MBC_reads <- MBC_reads[, colnames(MBC_reads)]

taxonomy <- taxonomy[colnames(MBC_reads), ]

#list orders of OTUs
orders <- unique(taxonomy$order)
# Taxonomic filtering - retain only OTUs of a certain taxon, e.g.:

# TJC: since as you can see from your above command, there are NA values in 
# taxonomy$order, you need to ensure you also remove these
taxonomy <- taxonomy[taxonomy$order == "Coleoptera" & !is.na(taxonomy$order), ]
#list arthropod orders
arthro_orders <- unique(taxonomy$order)
#
taxonomy <- taxonomy[taxonomy$order %in% arthro_orders, ]
# 1834 OTUs here in taxonomy

# TJC: Running na.omit of a dataframe like this is really dangerous, because
#  you don't know whether you necessarily want to remove all rows with NAs,
#  they might be in columns that are fine and you're just blanket removing all of
#  them. Now we know there are some NAs in order and we rightfully want to remove
#  them, but we don't mind if they have NA in family or genus, do we? We just
#  care if they're Coleoptera. So much better to handle NAs explicitly, like I
#  did above on line 241
#taxonomy<-na.omit(taxonomy)
MBC_reads <- MBC_reads[, rownames(taxonomy)]
# TJC: now it's 1808 but that's a minor change, just retaining some OTUs which are
# Coleoptera but have unknown family/genus

# After filtering you may have some OTUs or samples that now have 0 reads
# so these should be filtered out and the relevant lines in the metadata/
# taxonomy dropped as well

MBC_reads <- MBC_reads[rowSums(MBC_reads) > 0, colSums(MBC_reads) > 0]

dim(MBC_reads)
# TJC: so now you have 1605 Coleoptera OTUS in 411 samples. Respectable. 

metadata <- metadata[rownames(MBC_reads), ]
taxonomy <- taxonomy[colnames(MBC_reads), ]


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
jpeg("check_expectedrich.jpg", width =600, height=800)
check <- check_expected_richness(MBC_reads, 0.85)
dev.off()

# For most community ecology, the best route is to remove samples that 
# did not reach the threshold. check_expected_richness returns a vector
# of samples that passed this as the first item

MBC_reads <- MBC_reads[check[[1]],]

# An alternative to this standardisation is to rarefy. This process
# removes reads from samples until it reaches a threshold. This 
# standardises sampling depth but at the loss of real community data.
# It is more suitable for methodological analyses and is not recommended
# for community ecology.

# First lets review how many samples/OTUs might be lost at different 
# rarefaction target values:
#target_values = round(seq(min(rowSums(MBC_reads)), max(rowSums(MBC_reads)), length.out = 30))
#jpeg("rarexplore.jpg", width =480, height=480)
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
#1605 OTUs in MBC_reads

# Finally, we need to clean up after standardisation
# Remove any OTUs which no longer have any reads

MBC_reads <- MBC_reads[rowSums(MBC_reads) > 0, colSums(MBC_reads) > 0]
#1587 OTUs in MBC_reads

#Remove any dropped samples/OTUs from metadata/taxonomy tables
metadata <- metadata[rownames(MBC_reads), ]
taxonomy <- taxonomy[colnames(MBC_reads), ]

#______________________________________________________________
#EXPLORING DATA

#OPTIONAL create melted dataframe 
#MBC_reads_melt<-melt(MBC_reads)

#examine basic diversity indices
#OTU richness
specnumber(MBC_reads, groups = metadata$species)
treeid<-specnumber(MBC_reads, groups = metadata$treeid)
specnumber(MBC_reads, groups = metadata$camp)
specnumber(MBC_reads, groups = metadata$GPS_alt)
#___________________________________
#Mantel Test 

dist.pres_abs = vegdist(MBC_reads, method = "jaccard")

#environmental vector (gpsalt or elevation)
elevation = metadata$GPS_alt
dist_elevation = dist(elevation, method = "euclidean")
#Mantel elvation test
reads_elevation = mantel(dist.pres_abs, dist_elevation, method = "spearman", permutations = 9999, na.rm = TRUE)

#distance matrix for ongitude and latitude
geo = data.frame(df$Longitude, df$Latitude)
d_geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)
#longitude and latitude 
long_lat = data.frame(metadata$longitude, metadata$latitude)
#haversine distance - The shortest distance between two points (i.e., the 'great-circle-distance' or 'as the crow flies'), 
#according to the 'haversine method'. 
#This method assumes a spherical earth, ignoring ellipsoidal effects
#Is there a better method that accounts for elevation changes?
d_geo = distm(long_lat, fun = distHaversine)
dist.geo = as.dist(d_geo)
abund_geo  = mantel(dist.pres_abs, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)


#_______________________________
#NMDs
#Thomas version so that it can take the raw community data
#NMDS.scree <- function(x, distance, autotransform, trymax){
#  plot(rep(1,5), replicate(5, metaMDS(x, distance = distance, autotransform = autotransform, 
    #                                  k=1, trymax = trymax)$stress), xlim = c(1,5),
    #   ylim = c(0, 0.4), xla ="# of Dimensions", ylab = "Stress", main ="NMDS stress plot euclidean @1000")
#  for(i in 1:7){
  #  points(rep(i +1,5), replicate(5, metaMDS(x, distance = distance, autotransform = autotransform, 
                                        #     k=i+1, trymax = trymax)$stress))
 # }
#}

#NMDS.scree(MBC_reads, 'jaccard', TRUE, 100)
#reads_dist<-vegdist(MBC_reads, method = "jaccard")
#CommunityNMDS_k5_jaccard <- metaMDS(reads_dist,
#                                    distance = "jaccard",
 #                                   k = 5,
  #                                  maxit = 99, 
   #                                 trymax = 100, 
          #                          wescores = T)
#CommunityNMDS_k3_euclidean$stress
#stressplot(CommunityNMDS_k3_euclidean)
#ordiplot(CommunityNMDS_k3_euclidean,type="p")
#ordihull(CommunityNMDS_k3_euclidean,groups=metadata$climband,draw="polygon", col="grey90",label=T, cex=0.5, main = "CommunityNMDS_k3_euclidean100trymax")
#_____________________________
#Chao exploration

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

jpeg("chaobytreebar.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75, bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
ggplot(otuspecpooltree, aes(x=reorder(treeid, -chao), y=chao)) + geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=chao-chao.se, ymax=chao+chao.se))
dev.off()

#specaccume using vegan package
#Per tree species
#liquid amber
liquidamber_accum <- specaccum(MBC_reads, method = "random", permutations = 1000, subset=metadata$species=="Liquidambar styracaflua")
jpeg("liquidamber_accumVegan.jpg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75,bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo"))
plot(liquidamber_accum, ci.type = "poly", col = "blue", ci.col = "lightblue", 
     lwd = 2, ci.lty = 0, xlab = "number of traps", 
     ylab = "cumulative number of OTUs")
dev.off()

#____________________________________________
### Additive diversity partitioning 
metadata <- metadata %>% rownames_to_column('sample')
sample_hierarchy <- metadata %>% subset(select = c("sample", "treeid", "camp"))
metadata <- metadata %>% column_to_rownames('sample')
adDivpart<-adipart(MBC_reads, sample_hierarchy, index=c("richness", "shannon", "simpson"), nsimul=99)
plot_adipart(adDivpart)

#_______________________________________________
#How does choa vary with altitude? 
#sp richness/chao?


#___________________
#beta dissimilarity 
#____________________
#GLM, elevation on species richness per tray and chao per tray
OTUrichness<-apply(MBC_reads,1,sum) 
GLMdataframe<-as.data.frame(cbind(OTUrichness, metadata$GPS_alt))
#colnames(GLMdataframe) <- c("OTUrich", "elevation", "camp")


OTUrich.alt.mod<-glm(formula=OTUrichness~metadata$GPS_alt, family = poisson, data=GLMdataframe)
summary(OTUrich.alt.mod)

propnull.dev<-(OTUrich.alt.mod$null.deviance - OTUrich.alt.mod$deviance)/OTUrich.alt.mod$null.deviance
propnull.dev
exp(3.377e-04)

# predict for a neat sequence of log elevation values
pred <- expand.grid(lgelevation = seq(7, 9, by=0.005249344))
head(pred)
tail(pred)
pred$fit <- predict(OTUrich.alt.mod, newdata=pred, type='response')
head(pred)
plot(GLMdataframe$OTUrichness ~ log(GLMdataframe$V2), data=GLMdataframe)
lines(fit ~ lgelevation, data=pred, col='red')

