rm(list=ls())
library("ggplot2")
library("tidyverse")
library("reshape")
library("splitstackshape")
library("stringr")
library("lessR")
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

control_stats <- as.data.frame(slice(stats_raw, 1:9))
control_stats$readpercent <- rep(100, 9)
#Exploratory stats---------
#ReadPercent 
readpercent_regression<-function(
  
)
#Mantel
)

  
  
  
  
  
  
  
  
  
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

#Denoising----------------------

merged<-Merge(stats_raw, metadata_raw, by="runid")
merged$alpha<-as.numeric(alpha)
#mantel
#elevation
mantel_elevation<-merged %>% filter(stat_name == "elevation" & !alpha == "X") %>% slice(1:297)
mantel_elevation<-mantel_elevation[!duplicated(mantel_elevation),]
mantel_elevation$alpha<-as.numeric(mantel_elevation$alpha)
mantel_elevation$minsize<-as.factor(mantel_elevation$minsize)

ggplot(mantel_elevation, aes(x=alpha, y= stat, colour=minsize)) + geom_point()+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="mantel_elevation_stat")

#build a graph based on minsize only by alpha.2
ADP_alpha.2_X<-merged %>% filter(stat_name == "alpha.2" & alpha == "X") %>% slice(1:297)
ADP_alpha.2_X<-ADP_alpha.2_X[!duplicated(ADP_alpha.2_X),]
ADP_alpha.2_X$minsize<-as.factor(ADP_alpha.2_X$minsize)

graph_ADP_alpha.2_X<-ggplot(ADP_alpha.2_X, aes(x=minsize, y= stat, colour=runid)) + geom_point(size=4,)+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_alpha.2_X")

#ADP_alpha.1
ADP_alpha.1<-merged %>% filter(stat_name == "alpha.1" & !alpha == "X") %>% slice(1:297)
ADP_alpha.1<-ADP_alpha.1[!duplicated(ADP_alpha.1),]
ADP_alpha.1$alpha<-as.numeric(ADP_alpha.1$alpha)
ADP_alpha.1$minsize<-as.factor(ADP_alpha.1$minsize)

ggplot(ADP_alpha.1, aes(x=alpha, y= stat, colour=minsize)) + geom_point(size=4,)+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_alpha.1")

#ADP_alpha.2
ADP_alpha.2<-merged %>% filter(stat_name == "alpha.2" & !alpha == "X") %>% slice(1:297)
ADP_alpha.2<-ADP_alpha.2[!duplicated(ADP_alpha.2),]
ADP_alpha.2$alpha<-as.numeric(ADP_alpha.2$alpha)
ADP_alpha.2$minsize<-as.factor(ADP_alpha.2$minsize)

graph_ADP_alpha.2<-ggplot(ADP_alpha.2, aes(x=alpha, y= stat, colour=minsize)) + geom_point(size=4,)+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_alpha.2")

#How can I combine ADP_alpha.2 and ADP_alpha.2_X?

#ADP_gamma
ADP_gamma<-merged %>% filter(stat_name == "alpha.2" & !alpha == "X") %>% slice(1:297)
ADP_gamma<-ADP_gamma[!duplicated(ADP_gamma),]
ADP_gamma$alpha<-as.numeric(ADP_gamma$alpha)
ADP_gamma$minsize<-as.factor(ADP_gamma$minsize)

ggplot(ADP_gamma, aes(x=alpha, y= stat, colour=minsize)) + geom_point(size=4,)+
  geom_smooth(method=lm, se=T, linetype="dashed", color="darkred") + 
  labs(y="ADP_gamma")







