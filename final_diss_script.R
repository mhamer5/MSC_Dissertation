rm(list=ls())
dev.off()
#Load packages----
library("tidyverse")
library("reshape")
library("splitstackshape")
library("stringr")
library("lessR")
library("rlang")
library("ggpubr")
#Set theme for plotting (necessary?)
theme_set(theme_pubr())

#Data Wrangling----------------------------------
stats_raw<-read.table("master_output.tsv")
names<-c("runid", "testid", "stat_name", "stat", "p_value")
colnames(stats_raw)<-names
stats_raw$runid = gsub("reads/", "", stats_raw$runid)

metadata_raw<-read.csv("runid_metadata.csv")
metadata_raw$repno<-rep(9, 73)
names<-c("runid", "minsize", "alpha", "readpercent", "repno")
colnames(metadata_raw)<-names
metadata_raw<-expandRows(metadata_raw, count = 5, count.is.col = TRUE, drop = TRUE)
all(metadata_raw[,1] %in% stats_raw[,1])
rm(names)

#control_stats <- as.data.frame(slice(stats_raw, 1:9))
#control_stats$minsize <- rep(8, 9)
#control_stats$alpha <- as.numeric(rep(2, 9))
#control_stats$readpercent <- as.factor(rep(100, 9))

#Merge dataframe, unreplicated dataframe for both denoising and readpercent runs 
data_merged<-Merge(stats_raw, metadata_raw, by="runid")
#extract alpha=X runs 
denoise_aX<-distinct(filter(data_merged, alpha == "X"))
#remove alpha==X from datamerged
data_merged<-filter(data_merged, !alpha == "X")
#change data to be correct types
data_merged$alpha<-as.numeric(data_merged$alpha) 
data_merged$minsize<-as.factor(data_merged$minsize)
data_merged<-distinct(data_merged)

#Delimiting Functions----------------------------------

#Function to filter for denoise runs and required statistic
denoiseruns_filtering<-function(., stattochoose){
  data_merged %>% 
  filter(!grepl("readpercent_", runid), 
         stat_name == {{ stattochoose }}) %>%
  distinct()
}
#Function to plot regression line (lm) of filtered denoise dataframe
regression_plot_denoise <- function(dataframe, xaxis, yaxis, factor, yaxislab){
  ggplot2::ggplot(dataframe, aes(x= {{ xaxis }},y= {{ yaxis }}, colour = {{ factor }})) + 
    geom_point(size=2,)+
    geom_smooth(method = lm, 
                se=T,
                linetype ="dashed",
                colour = "darkred") + 
    ylab({{ yaxislab }})
}
#Function to filter for readpercent runs and required stats
readpercent_filtering<-function(., stattochoose){
  data_merged %>% 
    filter(., !grepl("denoise_", runid)) %>%
    filter(stat_name == {{ stattochoose }} & !alpha == "X")
}
#Function to plot regression line (lm) of filtered readpercent dataframe
regression_plot_readpercent <- function(dataframe, xaxis, yaxis, yaxislab){
  ggplot2::ggplot(dataframe, aes(x= {{ xaxis }},y= {{ yaxis }})) + 
    geom_point(size=2,)+
    geom_smooth(method = lm, 
                se=T,
                linetype ="dashed",
                colour = "darkred") +
    ylab({{ yaxislab }})
}

#READPERCENT----------------------------------
#Mantel
#Elevation
filtered_Mantel_elevation<-readpercent_filtering(., "elevation")
READPERCENT_Mantel_elevation<-regression_plot_readpercent(filtered_Mantel_elevation, readpercent, stat, "Mantel_elevation_stat")
Mantel_elevation_lm_mod<-lm(formula=stat~readpercent, data = filtered_Mantel_elevation)
summary(Mantel_elevation_lm_mod)
rm(filtered_Mantel_elevation)
#spatial stat
filtered_Mantel_spatial<-readpercent_filtering(., "spatial")
READPERCENT_Mantel_spatial<-regression_plot_readpercent(filtered_Mantel_spatial, readpercent, stat, "Mantel_spatial_stat")
rm(filtered_Mantel_spatial)

#ADP
#alpha.1
filtered_ADP_alpha.1<-readpercent_filtering(., "alpha.1")
READPERCENT_ADP_alpha.1<-regression_plot_readpercent(filtered_ADP_alpha.1, readpercent, stat,  "ADP_alpha.1_stat")
ADP_alpha.1_lm_mod<-lm(formula=stat~readpercent, data = READPERCENT_ADP_alpha.1)
summary(ADP_alpha.1_lm_mod)
rm(filtered_ADP_alpha.1)
#alpha.2
filtered_ADP_alpha.2<-readpercent_filtering(., "alpha.2")
READPERCENT_ADP_alpha.2<-regression_plot_readpercent(filtered_ADP_alpha.2, readpercent, stat, "ADP_alpha.2_stat")
ADP_alpha.2_lm_mod<-lm(formula=stat~readpercent, data = filtered_ADP_alpha.2)
summary(ADP_alpha.2_lm_mod)
rm(filtered_ADP_alpha.2)
#beta.1
filtered_ADP_beta.1<-readpercent_filtering(., "beta.1")
READPERCENT_ADP_beta.1<-regression_plot_readpercent(filtered_ADP_beta.1, readpercent, stat,  "ADP_beta.1_stat")
rm(filtered_ADP_beta.1)
#beta.2
filtered_ADP_beta.2<-readpercent_filtering(., "beta.2")
READPERCENT_ADP_beta.2<-regression_plot_readpercent(filtered_ADP_beta.2, readpercent, stat, "ADP_beta.2_stat")
rm(filtered_ADP_beta.2)
#gamma 
filtered_ADP_gamma<-readpercent_filtering(., "gamma")
READPERCENT_ADP_gamma<-regression_plot_readpercent(filtered_ADP_gamma, readpercent, stat,  "ADP_gamma_stat")
rm(filtered_ADP_gamma)

#GLMM
#Intercept
filtered_GLMM_Intercept<-readpercent_filtering(., "(Intercept)")
READPERCENT_GLMM_Intercept<-regression_plot_readpercent(filtered_GLMM_Intercept, readpercent, stat,"GLMM_Intercept_stat")
rm(filtered_GLMM_Intercept)
#scalealt
filtered_GLMM_scalealt<-readpercent_filtering(., "scalealt")
READPERCENT_GLMM_scalealt<-regression_plot_readpercent(filtered_GLMM_scalealt, readpercent, stat, "GLMM_scalealt_stat")
rm(filtered_GLMM_scalealt)


figure_readpercent <- ggarrange(READPERCENT_Mantel_elevation, READPERCENT_Mantel_spatial, 
                                READPERCENT_ADP_alpha.1, READPERCENT_ADP_alpha.2, READPERCENT_ADP_beta.1,
                                READPERCENT_ADP_beta.2, READPERCENT_ADP_gamma, READPERCENT_GLMM_Intercept,
                                READPERCENT_GLMM_scalealt,
                   ncol = 3, nrow = 3) %>% 
                    ggexport(filename = "figure_readpercent.jpeg",
                    width = 1024, 
                    height = 768, 
                    verbose = T)

rm(READPERCENT_Mantel_elevation, READPERCENT_Mantel_spatial, 
   READPERCENT_ADP_alpha.1, READPERCENT_ADP_alpha.2, READPERCENT_ADP_beta.1,
   READPERCENT_ADP_beta.2, READPERCENT_ADP_gamma, READPERCENT_GLMM_Intercept,
   READPERCENT_GLMM_scalealt)


#DENOISING----------------------------------
#Mantel
#elevation stat
filtered_Mantel_elevation<-denoiseruns_filtering(., "elevation")
DENOISE_Mantel_elevation<-regression_plot_denoise(filtered_Mantel_elevation, alpha, stat, minsize, "Mantel_elevation_stat")
Mantel_elevation_mod<-lm(formula=stat~alpha+as.numeric(minsize), data=filtered_Mantel_elevation)
summary(Mantel_elevation_mod)
rm(filtered_Mantel_elevation)
#spatial stat
filtered_Mantel_spatial<-denoiseruns_filtering(., "spatial")
DENOISE_Mantel_spatial<-regression_plot_denoise(filtered_Mantel_spatial, alpha, stat, minsize, "Mantel_spatial_stat")
rm(filtered_Mantel_spatial)

#ADP
#alpha.1
filtered_ADP_alpha.1<-denoiseruns_filtering(., "alpha.1")
DENOISE_ADP_alpha.1<-regression_plot_denoise(filtered_ADP_alpha.1, alpha, stat, minsize, "ADP_alpha.1_stat")
rm(filtered_ADP_alpha.1)
#alpha.2
filtered_ADP_alpha.2<-denoiseruns_filtering(., "alpha.2")
DENOISE_ADP_alpha.2<-regression_plot_denoise(filtered_ADP_alpha.2, alpha, stat, minsize, "ADP_alpha.2_stat")
ADP_alpha.2_mod<-lm(formula=stat~alpha+as.numeric(minsize), data=filtered_ADP_alpha.2)
summary(ADP_alpha.2_mod)
rm(filtered_ADP_alpha.2)
#beta.1
filtered_ADP_beta.1<-denoiseruns_filtering(., "beta.1")
DENOISE_ADP_beta.1<-regression_plot_denoise(filtered_ADP_beta.1, alpha, stat, minsize, "ADP_beta.1_stat")
rm(filtered_ADP_beta.1)
#beta.2
filtered_ADP_beta.2<-denoiseruns_filtering(., "beta.2")
DENOISE_ADP_beta.2<-regression_plot_denoise(filtered_ADP_beta.2, alpha, stat, minsize, "ADP_beta.2_stat")
rm(filtered_ADP_beta.2)
#gamma 
filtered_ADP_gamma<-denoiseruns_filtering(., "gamma")
DENOISE_ADP_gamma<-regression_plot_denoise(filtered_ADP_gamma, alpha, stat, minsize, "ADP_gamma_stat")
rm(filtered_ADP_gamma)

#GLMM
#Intercept
filtered_GLMM_Intercept<-denoiseruns_filtering(., "(Intercept)")
DENOISE_GLMM_Intercept<-regression_plot_denoise(filtered_GLMM_Intercept, alpha, stat, minsize, "GLMM_Intercept_stat")
rm(filtered_GLMM_Intercept)
#scalealt
filtered_GLMM_scalealt<-denoiseruns_filtering(., "scalealt")
DENOISE_GLMM_scalealt<-regression_plot_denoise(filtered_GLMM_scalealt, alpha, stat, minsize, "GLMM_scalealt_stat")
rm(filtered_GLMM_scalealt)


figure_denoise <- ggarrange(DENOISE_Mantel_elevation, DENOISE_Mantel_spatial, 
                                DENOISE_ADP_alpha.1, DENOISE_ADP_alpha.2, DENOISE_ADP_beta.1,
                                DENOISE_ADP_beta.2, DENOISE_ADP_gamma, DENOISE_GLMM_Intercept,
                                DENOISE_GLMM_scalealt,
                                ncol = 3, nrow = 3,
                            common.legend = TRUE, legend = "bottom") %>% ggexport(filename = "figure_denoise.jpeg",
                                                                                  width = 1024, 
                                                                                  height = 768, 
                                                                                  verbose = T)
rm(DENOISE_Mantel_elevation, DENOISE_Mantel_spatial, 
   DENOISE_ADP_alpha.1, DENOISE_ADP_alpha.2, DENOISE_ADP_beta.1,
   DENOISE_ADP_beta.2, DENOISE_ADP_gamma, DENOISE_GLMM_Intercept,
   DENOISE_GLMM_scalealt)









