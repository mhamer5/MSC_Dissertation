#!/usr/bin/env Rscript
options(stringsAsFactors = F)

suppressMessages(require(getopt))
suppressMessages(require(vegan))
suppressMessages(require(lme4))

# Set up options
spec <- matrix(c(
  'reads',        'r', 1, 'character',
  'metadata',     'm', 1, 'character',
  'paramversion', 'n', 1, 'numeric'
), byrow = T, ncol = 4)

# Read options
opt <- getopt(spec)

# Load data ---------------------------------------------------------------

# Load the reads table output by vsearch --usearch_global --tabbedout
reads <- t(read.table(opt$reads, header = T, sep = "\t", comment.char = '',
                      row.names = 1, check.names = F))

# Load metadata table (change command if yours is a different format!)
metadata <- read.csv(opt$metadata, row.names = 1)

# Filtering low abundance OTU incidences ----------------------------------

# Read proportion filtering
samplewise_p <- t(apply(reads, 1, function(s) s/sum(s) ))
threshold <- 0.000125

reads[samplewise_p < threshold] <- 0
reads <- reads[rowSums(reads) > 0, colSums(reads) > 0]
rm(samplewise_p, threshold)

# Organise data -----------------------------------------------------------

metadata <- metadata[row.names(reads), ]

# Filtering data ----------------------------------------------------------

#Remove Pinus trees
metadata <- metadata[metadata$species != "Pinus sp.",]

reads <- reads[row.names(metadata), ]

reads <- reads[rowSums(reads) > 0, colSums(reads) > 0]

metadata <- metadata[rownames(reads), ]

#  Standardisation --------------------------------------------------------

metadata <- cbind(metadata, t(estimateR(reads)))
metadata$p_obs_exp <- with(metadata, S.obs/S.chao1)
metadata <- metadata[metadata$p_obs_exp >= 0.85, ]
reads <- reads[row.names(metadata),]

reads[reads > 0] <- 1

reads <- reads[rowSums(reads) > 0, colSums(reads) > 0]
metadata <- metadata[rownames(reads), ]

# Mantel tests ------------------------------------------------------------

reads_jaccard <- vegdist(reads, 'jaccard')
perms <- 9999

manteltests <- list(
  mantel_elevation =  mantel(reads_jaccard, 
                             dist(metadata$DEM30_alt), 
                             method = "spearman", permutations = perms, na.rm = TRUE),
  mantel_spatial =  mantel(reads_jaccard,
                           dist(metadata[,c("utm_x", "utm_y")]),
                           method = "spearman", permutations = perms, na.rm = TRUE))

manteltests <- t(sapply(manteltests, function(x) c(stat = x$statistic, p = x$signif)))

# Additive diversity partitioning -----------------------------------------

metadata$sample <- row.names(metadata)
sample_hierarchy <-  subset(metadata, select = c("sample", "treeid", "camp"))
addpart <- adipart(reads, sample_hierarchy, index = "richness", nsimul = 9)
addpart <- cbind(stat = addpart$oecosimu$means - addpart$statistic, p = addpart$oecosimu$pval)
row.names(addpart) <- paste('adp', row.names(addpart), sep = '_')

# GLMM --------------------------------------------------------------------

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
metadata$treeid <- as.factor(metadata$treeid)
metadata$scalealt <- (metadata$DEM30_alt - mean(metadata$DEM30_alt))/sd(metadata$DEM30_alt)

m1 <- glmer(S.obs ~ scalealt + (1|treeid), data = metadata, family = poisson)
glmout <- coef(summary(m1))
phi <- overdisp_fun(m1)["ratio"]
glmout <- within(as.data.frame(glmout),
                 {`Std. Error` <- `Std. Error` * sqrt(phi)
                 `z value` <- Estimate / `Std. Error`
                 stat <- Estimate
                 `p` <- 2 * pnorm(abs(`z value`), lower.tail=FALSE)
                 })

glmout <- glmout[, c('stat',  'p')]
row.names(glmout) <- paste('glm', row.names(glmout), sep = '_')

# Bring together ----------------------------------------------------------

tests <- rbind(manteltests, addpart, glmout)
write.table(cbind(opt$paramversion, 
                  do.call('rbind', strsplit(row.names(tests), '_')), 
                  tests), stdout(), row.names = F, col.names = F, quote = F)
