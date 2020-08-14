#This could do with a lot of cleaning up 
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
library("grid")
library("patchwork")
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

control_stats <- as.data.frame(slice(stats_raw, 1:9))
control_stats$minsize <- rep(8, 9)
control_stats$alpha<-rep(2,9)
control_stats$readpercent<-rep(100,9)
#Merge dataframe, unreplicated dataframe for both denoising and readpercent runs 
data_merged<-Merge(stats_raw, metadata_raw, by="runid")
#extract alpha=X runs 
denoise_aX<-filter(data_merged, alpha == "X")
#change X to 0
denoise_aX$alpha<-0
#remove alpha==X from datamerged
data_merged<-filter(data_merged, !alpha == "X")
#change data to be correct types
data_merged$alpha <- as.numeric(as.character(data_merged$alpha))
data_merged<-union(data_merged,denoise_aX)
data_merged<-distinct(data_merged)
denoise_runs<-data_merged[grepl("denoise_", data_merged$runid),]
denoise_runs<-union(denoise_runs, control_stats)
readpercent_runs<-data_merged[grepl("readpercent_", data_merged$runid),]
readpercent_runs<-union(readpercent_runs,control_stats)


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
#filtered_Mantel_elevation<-readpercent_filtering(., "elevation")
#READPERCENT_Mantel_elevation<-regression_plot_readpercent(filtered_Mantel_elevation, readpercent, stat, "Mantel_elevation_stat")
#Mantel_elevation_lm_mod<-lm(formula=stat~readpercent, data = filtered_Mantel_elevation)
#summary(Mantel_elevation_lm_mod)
#rm(filtered_Mantel_elevation)
#spatial stat
#filtered_Mantel_spatial<-readpercent_filtering(., "spatial")
#READPERCENT_Mantel_spatial<-regression_plot_readpercent(filtered_Mantel_spatial, readpercent, stat, "Mantel_spatial_stat")
#rm(filtered_Mantel_spatial)

#ADP
#alpha.1
#filtered_ADP_alpha.1<-readpercent_filtering(., "alpha.1")
#READPERCENT_ADP_alpha.1<-regression_plot_readpercent(filtered_ADP_alpha.1, readpercent, stat,  "ADP_alpha.1_stat")
#ADP_alpha.1_lm_mod<-lm(formula=stat~readpercent, data = READPERCENT_ADP_alpha.1)
#summary(ADP_alpha.1_lm_mod)
#rm(filtered_ADP_alpha.1)
#alpha.2
#filtered_ADP_alpha.2<-readpercent_filtering(., "alpha.2")
#READPERCENT_ADP_alpha.2<-regression_plot_readpercent(filtered_ADP_alpha.2, readpercent, stat, "ADP_alpha.2_stat")
#ADP_alpha.2_lm_mod<-lm(formula=stat~readpercent, data = filtered_ADP_alpha.2)
#summary(ADP_alpha.2_lm_mod)
#rm(filtered_ADP_alpha.2)
#beta.1
#filtered_ADP_beta.1<-readpercent_filtering(., "beta.1")
#READPERCENT_ADP_beta.1<-regression_plot_readpercent(filtered_ADP_beta.1, readpercent, stat,  "ADP_beta.1_stat")
#rm(filtered_ADP_beta.1)
#beta.2
#filtered_ADP_beta.2<-readpercent_filtering(., "beta.2")
#READPERCENT_ADP_beta.2<-regression_plot_readpercent(filtered_ADP_beta.2, readpercent, stat, "ADP_beta.2_stat")
#rm(filtered_ADP_beta.2)
#gamma 
#filtered_ADP_gamma<-readpercent_filtering(., "gamma")
#READPERCENT_ADP_gamma<-regression_plot_readpercent(filtered_ADP_gamma, readpercent, stat,  "ADP_gamma_stat")
#rm(filtered_ADP_gamma)

#GLMM
#Intercept
#filtered_GLMM_Intercept<-readpercent_filtering(., "(Intercept)")
#READPERCENT_GLMM_Intercept<-regression_plot_readpercent(filtered_GLMM_Intercept, readpercent, stat,"GLMM_Intercept_stat")
#rm(filtered_GLMM_Intercept)
#scalealt
#filtered_GLMM_scalealt<-readpercent_filtering(., "scalealt")
#READPERCENT_GLMM_scalealt<-regression_plot_readpercent(filtered_GLMM_scalealt, readpercent, stat, "GLMM_scalealt_stat")
#rm(filtered_GLMM_scalealt)


#figure_readpercent <- ggarrange(READPERCENT_Mantel_elevation, READPERCENT_Mantel_spatial, 
#                                READPERCENT_ADP_alpha.1, READPERCENT_ADP_alpha.2, READPERCENT_ADP_beta.1,
 #                               READPERCENT_ADP_beta.2, READPERCENT_ADP_gamma, READPERCENT_GLMM_Intercept,
  #                              READPERCENT_GLMM_scalealt,
   #                ncol = 3, nrow = 3) %>% 
    #                ggexport(filename = "figure_readpercent.jpeg",
     ##               width = 1024, 
       #             height = 768, 
        #            verbose = T)

#rm(READPERCENT_Mantel_elevation, READPERCENT_Mantel_spatial, 
 #  READPERCENT_ADP_alpha.1, READPERCENT_ADP_alpha.2, READPERCENT_ADP_beta.1,
  # READPERCENT_ADP_beta.2, READPERCENT_ADP_gamma, READPERCENT_GLMM_Intercept,
   #READPERCENT_GLMM_scalealt)


#DENOISING----------------------------------
#Mantel
#elevation stat
#filtered_Mantel_elevation<-denoiseruns_filtering(., "elevation")
#DENOISE_Mantel_elevation<-regression_plot_denoise(filtered_Mantel_elevation, alpha, stat, minsize, "Mantel_elevation_stat")
#Mantel_elevation_mod<-lm(formula=stat~alpha+as.numeric(minsize), data=filtered_Mantel_elevation)
#summary(Mantel_elevation_mod)
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
#ADP_alpha.2_mod<-lm(formula=stat~alpha+as.numeric(minsize), data=filtered_ADP_alpha.2)
#summary(ADP_alpha.2_mod)
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

#Investigating more closely----
#Denoising---- 
#Elevation
#8 plots, each plot with alpha on x axis, four row 2 column lay out 
#y axis = stat left hand column is elevation and right is spatial stat
#each row is a different value of minsize 
#do this with the p values as well, make the plot character close if not signifcant, open if signifcant for each point
#change the limits of the y axis to make there be a 0 i the y axis

#GLMM - this has been very hard to plot, stat values are massively different to each other 
denoises_controls_GLM_Intercept<-denoise_runs %>%
  filter(stat_name == "(Intercept)")

denoises_controls_GLM_scalealt<-denoise_runs %>%
  filter(stat_name == "scalealt")

g1<-ggplot(subset(denoises_controls_GLM_Intercept,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_point(data=denoises_controls_GLM_Intercept[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_GLM_Intercept[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  geom_line()+
  theme(axis.text=element_text(size=5)) +
  ylab("GLM Intercept")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(0,4.5)
g2<-g1+facet_grid(~minsize)
g2
g3<-ggplot(subset(denoises_controls_GLM_scalealt,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_point(data=denoises_controls_GLM_scalealt[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_GLM_scalealt[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  geom_line()+
  theme(axis.text=element_text(size=5)) +
  ylab("GLM Scaled Altitude")+
  xlab("Denoise Alpha")+
  theme_gray(base_size = 10) +
  ylim(0,0.1) 
g4<-g3+facet_grid(~minsize)
g5<-g2/g4
g5 %>% ggexport(filename = "figure_denoise_GLM.jpeg",
               width = 1024, 
                height = 768, 
                verbose = T,
                res = 250)
rm(g1,g2, g3, g4, g5,denoises_controls_GLM_Intercept,denoises_controls_GLM_scalealt)


#Mantel - this one is easier, values are very similar in size, unlike GLM stat values
denoises_controls_Mental_elevation<-denoise_runs %>%
  filter(stat_name == "elevation")

denoises_controls_Mental_spatial<-denoise_runs %>%
  filter(stat_name == "spatial")

g1<-ggplot(subset(denoises_controls_Mental_elevation,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_line()+
  geom_point(data=denoises_controls_Mental_elevation[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_Mental_elevation[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  theme(axis.text=element_text(size=5)) +
  ylab("Mantel elevation")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(0,0.75)
g2<-g1+facet_grid(~minsize)

g3<-ggplot(subset(denoises_controls_Mental_spatial,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_line()+
  geom_point(data=denoises_controls_Mental_spatial[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_Mental_spatial[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  theme(axis.text=element_text(size=5)) +
  ylab("Mantel Spatial")+
  xlab("Denoise Alpha")+
  theme_gray(base_size = 10) +
  ylim(0,0.75) 
g4<-g3+facet_grid(~minsize)
g5<-g2/g4
g5 %>% ggexport(filename = "figure_denoise_Mantel.jpeg",
                width = 1024, 
                height = 768, 
                verbose = T,
                res = 250)
rm(g1,g2, g3, g4, g5, denoises_controls_Mental_elevation,denoises_controls_Mental_spatial)

#ADP 

denoises_controls_adp_alpha.1<-denoise_runs %>%
  filter(stat_name == "alpha.1")

denoises_controls_adp_alpha.2<-denoise_runs %>%
  filter(stat_name == "alpha.2")

g1<-ggplot(subset(denoises_controls_adp_alpha.1,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_point(data=denoises_controls_adp_alpha.1[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_adp_alpha.1[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  geom_line()+
  theme(axis.text=element_text(size=5)) +
  ylab("alpha.1")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(-8,0)
g2<-g1+facet_grid(~minsize)

g3<-ggplot(subset(denoises_controls_adp_alpha.2,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_point(data=denoises_controls_adp_alpha.2[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_adp_alpha.2[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  geom_line()+
  theme(axis.text=element_text(size=5)) +
  ylab("alpha.2")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(0,400) 
g4<-g3+facet_grid(~minsize)
g5<-ggarrange(g2,g4, ncol =1, nrow=2)
g5
denoises_controls_adp_beta.1<-denoise_runs %>%
  filter(stat_name == "beta.1")

denoises_controls_adp_beta.2<-denoise_runs %>%
  filter(stat_name == "beta.2")

g6<-ggplot(subset(denoises_controls_adp_beta.1,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_point(data=denoises_controls_adp_beta.1[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_adp_beta.1[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  geom_line()+
  theme(axis.text=element_text(size=5)) +
  ylab("Beta.1")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(0,450)
g7<-g6+facet_grid(~minsize)

g8<-ggplot(subset(denoises_controls_adp_beta.2,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_point(data=denoises_controls_adp_beta.2[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_adp_beta.2[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  geom_line()+
  theme(axis.text=element_text(size=5)) +
  ylab("Beta.2")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(0,150) 
g9<-g8+facet_grid(~minsize)
g10<-ggarrange(g7,g9, ncol =1, nrow=2)
g10
denoises_controls_adp_gamma<-denoise_runs %>%
  filter(stat_name == "gamma")

g11<-ggplot(subset(denoises_controls_adp_gamma,alpha!=0), aes(x=alpha,y=stat))+
  geom_point(shape=1, size=1)+
  geom_point(data=denoises_controls_adp_gamma[29:32, ], aes(x=alpha, y=stat), colour="red", size=1.5)+
  geom_point(data=denoises_controls_adp_gamma[33, ], aes(x=alpha, y=stat), colour="blue", size=1.5)+
  geom_line()+
  theme(axis.text=element_text(size=5)) +
  ylab("Gamma")+
  xlab("Denoise Alpha")+
  theme_gray(base_size = 10) +
  ylim(0,450)
g12<-g11+facet_grid(~minsize)

figure_denoise_ADP<-ggarrange(g5, g10, g12, ncol =1, nrow=3)
figure_denoise_ADP

figure_denoise_ADP%>%ggexport(filename = "figure_denoise_ADP.jpeg",
                verbose = T,
                width = 1500,
                height = 2000,
                res = 250)

rm(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12, figure_denoise_ADP)

#ReadPercentage----
#control_stats <- as.data.frame(slice(stats_raw, 1:9))
#control_stats$minsize <- rep(8, 9)
#control_stats$alpha <- as.numeric(rep(2, 9))
#control_stats$readpercent <- rep(100, 9)
#control_stats<-control_stats %>% select(runid,testid,stat_name, stat, p_value, readpercent)

#GLM
read_percent_GLM_intercept<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "(Intercept)") 

read_percent_GLM_scalealt<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "scalealt") 

g1<-ggplot(read_percent_GLM_intercept, aes(x=readpercent,y=stat))+
  geom_point(size=2)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  ylab("Intercept")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(0,4)

g2<-ggplot(read_percent_GLM_scalealt, aes(x=readpercent,y=stat))+
  geom_point(size=2)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  xlab("Read Percentage (%)")+
  ylab("Scaled Altitude")+
  theme_gray(base_size = 10) +
  ylim(0,0.06)

g3<-ggarrange(g1, g2,
         ncol = 1, nrow= 2,
         font.label = list(size = 8))
g3
g3%>%ggexport(filename = "figure_readpercent_GLM.jpeg",
              width = 1024, 
              height = 768, 
              verbose = T,
              res = 250)
rm(g1,g2,g3,read_percent_GLM_scalealt,read_percent_GLM_scalealt)
#Mantel
read_percent_MANTEL_elevation<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "elevation") 

read_percent_MANTEL_spatial<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "spatial") 

g1<-ggplot(read_percent_MANTEL_elevation, aes(x=readpercent,y=stat))+
  geom_point(size=1)+
  geom_point(data=read_percent_MANTEL_elevation[40, ], aes(x=readpercent, y=stat), colour="bue", size=1.5)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  ylab("Elevation")+
  xlab("")+
  theme_gray(base_size = 10) +
  ylim(0,0.75)

g2<-ggplot(read_percent_MANTEL_spatial,  aes(x=readpercent,y=stat))+
  geom_point(size=1)+
  geom_point(data=read_percent_MANTEL_spatial[40, ], aes(x=readpercent, y=stat), colour="blue", size=1.5)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  xlab("Read Percentage (%)")+
  ylab("Spatial")+
  theme_gray(base_size = 10) +
  ylim(0,0.75)

g3<-ggarrange(g1, g2,
              ncol = 1, nrow= 2,
              font.label = list(size = 8))
g3%>%ggexport(filename = "figure_readpercentage_Mantel.jpeg",
           width = 1024, 
           height = 768, 
           verbose = T,
           res = 250)

rm(g1,g2,g3,read_percent_MANTEL_elevation,read_percent_MANTEL_spatial)

#ADP
#Define graphs
read_percent_ADP_alpha.2<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "alpha.2")
read_percent_ADP_alpha.1<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "alpha.1")
read_percent_ADP_beta.1<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "beta.1")
read_percent_ADP_beta.2<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "beta.2")
read_percent_ADP_gamma<-readpercent_runs %>%
  select(runid,testid,stat_name, p_value, stat, readpercent) %>%
  filter(stat_name == "gamma")

alpha.2<-ggplot(read_percent_ADP_alpha.2, aes(x=readpercent,y=stat))+
  geom_point(size=1)+
  geom_point(data=read_percent_ADP_alpha.2[41, ], aes(x=readpercent, y=stat), colour="blue", size=1.5)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  xlab("Read Percentage (%)")+
  ylab("Alpha.2")+
  theme_gray(base_size = 10) +
  ylim(0,200)
alpha.1<-ggplot(read_percent_ADP_alpha.1, aes(x=readpercent,y=stat))+
  geom_point(size=1)+
  geom_point(data=read_percent_ADP_alpha.1[41, ], aes(x=readpercent, y=stat), colour="blue", size=1.5)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  xlab("Read Percentage (%)")+
  ylab("Alpha.1")+
  theme_gray(base_size = 10)+
  ylim(-5,0)
beta.1<-ggplot(read_percent_ADP_beta.1, aes(x=readpercent,y=stat))+
  geom_point(size=1)+
  geom_point(data=read_percent_ADP_beta.1[41, ], aes(x=readpercent, y=stat), colour="blue", size=1.5)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  xlab("Read Percentage (%)")+
  ylab("Beta.1")+
  theme_gray(base_size = 10) +
  ylim(0,200)
beta.2<-ggplot(read_percent_ADP_beta.2, aes(x=readpercent,y=stat))+
  geom_point(size=1)+
  geom_point(data=read_percent_ADP_beta.2[41, ], aes(x=readpercent, y=stat), colour="blue", size=1.5)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  xlab("Read Percentage (%)")+
  ylab("Beta.2")+
  theme_gray(base_size = 10) + 
  ylim(0,50)
gamma<-ggplot(read_percent_ADP_gamma, aes(x=readpercent,y=stat))+
  geom_point(size=1)+
  geom_point(data=read_percent_ADP_gamma[41, ], aes(x=readpercent, y=stat), colour="blue", size=1.5)+
  geom_smooth(method = lm, 
              se=T,
              linetype ="dashed",
              colour = "darkred")+
  theme(axis.text=element_text(size=5.8)) +
  xlab("Read Percentage (%)")+
  ylab("Gamma")+
  theme_gray(base_size = 10) +
  ylim(0,250)
gamma
adp_readpercent<-(alpha.1|alpha.2)/(beta.1|beta.2)/gamma +
  plot_layout(widths=1)
adp_readpercent
adp_readpercent%>%ggexport(filename = "figure_readpercentadp.jpeg",
              width = 1024, 
              height = 768, 
              verbose = T,
              res = 250)

rm(list=ls())
dev.off()
















