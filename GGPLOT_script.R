# GGPLOT figures

setwd("File pathway")

################
library(lattice)

library("plotly")
library(ggplot2)
library("extrafont")
library("wesanderson")
library("scales")
library(nlme)
library(plyr)
library(rlang)


########### Clean ggplot graphs

cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color = 'black'))

# Figure egg size averaged by species vs. mean temperature
numeggs <- read.csv("MeanLifehistoryCompleteDatabase.csv", header =T) #123 spp / 241 obs. Averaged
numeggs <- read.csv("FecundityHHT.csv", header =T) #23 spp. / 36 obs

# plotting Size vs SST
ggplot(numeggs, aes(sstmean, (volume), alpha=0.4, color=Symbionts))+
  ylab('Egg size (mm3)')+
  geom_point(size =2.5)+
  scale_y_continuous(trans = log_trans())+
  cleanup

# plotting Size vs Fecundity
ggplot(numeggs, aes(Eggs_per_area, (eggsize), alpha=0.4, color=Symbionts))+
  ylab('Egg size (mm3)')+
  geom_point(size =2.5)+
  geom_hline(yintercept = 52169051.82, color="pink")+
  scale_x_continuous(trans = log_trans())+
  scale_y_continuous(trans = log_trans())+
  cleanup

