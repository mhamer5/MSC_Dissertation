rm(list=ls())
library("tidyverse")
library("reshape")
library("splitstackshape")
library("stringr")
library("lessR")
library("rlang")
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
#control_stats$readpercent <- rep(100, 9)

#Merge dataframe, unreplicated dataframe for both denoising and readpercent runs 
data_merged<-Merge(stats_raw, metadata_raw, by="runid")
#extract alpha=X runs 
denoise_aX<-distinct(filter(data_merged, alpha == "X"))
#remove alpha==X from datamerged
data_merged<-filter(data_merged, !alpha == "X")
#change data to be correct types
data_merged$alpha<-as.numeric(data_merged$alpha) 
data_merged$minsize<-as.factor(data_merged$minsize)

#Function to filter for denoise runs and required statistic
denoiseruns_filtering<-function(., stattochoose){
  data_merged %>% 
  slice(1:2674) %>% 
  filter(stat_name == {{ stattochoose }} & !alpha == "X") %>%
  distinct()
}
#Function to plot regression line (lm) of filtered denoise dataframe
regression_plot_denoise <- function(dataframe, xaxis, yaxis, factor, yaxislab){
  ggplot2::ggplot(dataframe, aes(x= {{ xaxis }},y= {{ yaxis }}, colour = {{ factor }})) + 
    geom_point(size=4,)+
    geom_smooth(method = lm, 
                se=T,
                linetype ="dashed",
                colour = "darkred") + 
    ylab({{ yaxislab }})
}

#Mantel
#elevation stat
filtered_Mantel_elevation<-denoiseruns_filtering(., "elevation")
regression_plot_denoise(filtered_alpha.1, alpha, stat, minsize, "Mantel_elevation_stat")
rm(filtered_Mantel_elevation)
#spatial stat
filtered_Mantel_spatial<-denoiseruns_filtering(., "spatial")
regression_plot_denoise(filtered_spatial, alpha, stat, minsize, "Mantel_spatial_stat")
rm(filtered_Mantel_spatial)

#ADP
#alpha.1
filtered_ADP_alpha.1<-denoiseruns_filtering(., "alpha.1")
regression_plot_denoise(filtered_ADP_alpha.1, alpha, stat, minsize, "ADP_alpha.1_stat")
rm(filtered_ADP_alpha.1)
#alpha.2
filtered_ADP_alpha.2<-denoiseruns_filtering(., "alpha.2")
regression_plot_denoise(filtered_ADP_alpha.2, alpha, stat, minsize, "ADP_alpha.2_stat")
rm(filtered_ADP_alpha.2)
#beta.1
filtered_ADP_beta.1<-denoiseruns_filtering(., "beta.1")
regression_plot_denoise(filtered_ADP_alpha.1, alpha, stat, minsize, "ADP_beta.1_stat")
rm(filtered_ADP_beta.1)
#beta.2
filtered_ADP_beta.2<-denoiseruns_filtering(., "beta.2")
regression_plot_denoise(filtered_ADP_beta.2, alpha, stat, minsize, "ADP_beta.2_stat")
rm(filtered_ADP_beta.2)
#gamma 
filtered_ADP_gamma<-denoiseruns_filtering(., "gamma")
regression_plot_denoise(filtered_ADP_gamma, alpha, stat, minsize, "ADP_gamma_stat")
rm(filtered_ADP_gamma)

#GLMM
#Intercept
filtered_GLMM_Intercept<-denoiseruns_filtering(., "(Intercept)")
regression_plot_denoise(filtered_GLMM_Intercept, alpha, stat, minsize, "GLMM_Intercept_stat")
rm(filtered_GLMM_Intercept)
#scalealt
filtered_GLMM_scalealt<-denoiseruns_filtering(., "scalealt")
regression_plot_denoise(filtered_GLMM_scalealt, alpha, stat, minsize, "GLMM_scalealt_stat")
rm(filtered_GLMM_scalealt)




















#need a function to filter for readpercent runs, do not dereplicate rows


#Exploratory stats---------
#ReadPercent 
readpercent_regression<-function(filtered_dataframe, stat_name, ylabel){
  dplyr::filtered_dataframe<-stats_raw %>% slice(298:657) %>% filter(stat_name == stat_name)
  dplyr::filtered_dataframe$readpercent <- as.numeric(str_extract(mantel_elevation$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
  dplyr::filtered_dataframe<-filter(control_stats, stat_name == stat_name)  %>%  union(filtered_dataframe, control_stats)
  
  ggplot2::ggplot(filtered_dataframe, aes(x=readpercent, y= stat)) + geom_point()+
    geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
    labs(y="ylabel")
}


#Mantel
readpercent_regression(mantel_elevation, elevation, mantel_elevation_stat)

mantel_elevation<-stats_raw %>% slice(298:657) %>% filter(stat_name == c("elevation"))
mantel_elevation$readpercent <- as.numeric(str_extract(mantel_elevation$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
mantel_elevation<-filter(control_stats, stat_name == "elevation") %>% union(mantel_elevation, control_stats)

ggplot(mantel_elevation, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="mantel_elevation_stat")

ggplot(mantel_elevation, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="mantel_elevatation_pvalue")

rm(mantel_elevation)

mantel_spatial<-stats_raw %>% slice(298:657) %>% filter(stat_name == c("spatial"))
mantel_spatial$readpercent <- as.numeric(str_extract(mantel_spatial$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
mantel_spatial<-filter(control_stats, stat_name == "spatial") %>% union(mantel_spatial, control_stats)

ggplot(mantel_spatial, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="mantel_spatial_stat")

ggplot(mantel_spatial, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="mantel_elevatation_pvalue")

rm(mantel_spatial)

#GLMM

GLMM_intercept<-stats_raw %>% slice(298:657) %>% filter(stat_name == "(Intercept)")
GLMM_intercept$readpercent <- as.numeric(str_extract(GLMM_intercept$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
GLMM_intercept<-filter(control_stats, stat_name == "(Intercept)") %>% union(GLMM_intercept, control_stats)

ggplot(GLMM_intercept, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="GLMM_intercept_stat")

ggplot(GLMM_intercept, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="GLMM_intercept_pvalue")

rm(GLMM_intercept)

GLMM_scalealt<-stats_raw %>% slice(298:657) %>% filter(stat_name == "scalealt")
GLMM_scalealt$readpercent <- as.numeric(str_extract(GLMM_scalealt$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
GLMM_scalealt<-filter(control_stats, stat_name == "scalealt") %>% union(GLMM_scalealt, control_stats)

ggplot(GLMM_scalealt, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="GLMM_scalealt_stat")

ggplot(GLMM_scalealt, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="GLMM_scalealt_pvalue")

rm(GLMM_scalealt)

#Additive Diversity Partitioning
#alpha.1
ADP_alpha.1<-stats_raw %>% slice(298:657) %>% filter(stat_name == "alpha.1")
ADP_alpha.1$readpercent <- as.numeric(str_extract(ADP_alpha.1$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
ADP_alpha.1<-filter(control_stats, stat_name == "alpha.1") %>% union(ADP_alpha.1, control_stats)

ggplot(ADP_alpha.1, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_alpha.1_stat")

ggplot(ADP_alpha.1, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_alpha.1_pvalue")
rm(ADP_alpha.1)
#alpha.2
ADP_alpha.2<-stats_raw %>% slice(298:657) %>% filter(stat_name == "alpha.2")
ADP_alpha.2$readpercent <- as.numeric(str_extract(ADP_alpha.2$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
ADP_alpha.2<-filter(control_stats, stat_name == "alpha.2") %>% union(ADP_alpha.2, control_stats)

ggplot(ADP_alpha.2, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_alpha.2_stat")

ggplot(ADP_alpha.2, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_alpha.2_pvalue")
rm(ADP_alpha.2)
#gamma
ADP_gamma<-stats_raw %>% slice(298:657) %>% filter(stat_name == "gamma")
ADP_gamma$readpercent <- as.numeric(str_extract(ADP_gamma$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
ADP_gamma<-filter(control_stats, stat_name == "gamma") %>% union(ADP_gamma, control_stats)

ggplot(ADP_gamma, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="gamma_stat")

ggplot(ADP_gamma, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="gamma_pvalue")
rm(ADP_gamma)
#beta.1
ADP_beta.1<-stats_raw %>% slice(298:657) %>% filter(stat_name == "beta.1")
ADP_beta.1$readpercent <- as.numeric(str_extract(ADP_beta.1$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
ADP_beta.1<-filter(control_stats, stat_name == "beta.1") %>% union(ADP_beta.1, control_stats)

ggplot(ADP_beta.1, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="beta.1_stat")

ggplot(ADP_beta.1, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="beta.1_pvalue")
rm(ADP_beta.1)
#beta.2
ADP_beta.2<-stats_raw %>% slice(298:657) %>% filter(stat_name == "beta.2")
ADP_beta.2$readpercent <- as.numeric(str_extract(ADP_beta.2$runid, pattern = "_(.*)_") %>% gsub(pattern ="_","",.))
ADP_beta.2<-filter(control_stats, stat_name == "beta.2") %>% union(ADP_beta.2, control_stats)

ggplot(ADP_beta.2, aes(x=readpercent, y= stat)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="beta.2_stat")

ggplot(ADP_beta.2, aes(x=readpercent, y= p_value)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="beta.2_pvalue")
rm(ADP_beta.2)

#Could I turn the blocks into a single function that would plot the results for me when I parse it 
#stats_raw$stat_name?








