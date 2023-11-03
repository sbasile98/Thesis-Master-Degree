library(dplyr)
library(ggplot2)
library(reshape2)

setwd("C:/Users/simon/OneDrive/Desktop")

res <- read.csv("Grid_1155_585.csv", sep = ";", header = FALSE)
resML <- read.csv("Grid_ML_1155_583.csv", sep = ";", header = FALSE)

resFIN <- merge(x=res, y=resML, by="V1")

colnames(resFIN)[colnames(resFIN) == 'V1'] = "Generation"
colnames(resFIN)[colnames(resFIN) == 'V2.x'] = "Fitness"
colnames(resFIN)[colnames(resFIN) == 'V2.y'] = "FitnessML"
colnames(resFIN)[colnames(resFIN) == 'V3.y'] = "Time"

#ggplot(data=subset(resFIN, resFIN$Generation >= 0 & resFIN$Generation <= 50 & resFIN$Fitness <= 1800)) + 
  #geom_line(aes(x=Generation, y=Fitness, color="EA")) + 
  #geom_line(aes(x=Generation, y=FitnessML, color="EAML"))

ggplot(data=subset(resFIN, resFIN$Generation >= 150 & resFIN$Generation <= 350 & resFIN$Fitness <= 10000)) + 
  geom_line(aes(x=Time, y=Fitness, color="EA")) + 
  geom_line(aes(x=Time, y=FitnessML, color="EAML"))