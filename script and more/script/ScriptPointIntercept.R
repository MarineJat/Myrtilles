setwd("D:/Utilisateurs/alois/Documents/Etudiant/2024Stage_Master/Jattiot/dataANDanalyses") 

library(stringr)
library(dplyr)
library(ggplot2)

rm(list=ls(all=TRUE))

options(na.action = "na.fail")

e<-read.delim("PointIntercept.csv",sep=";")
head(e)
names(e)

#renaming column to have shorter names
e$adminNames<-as.factor(word(e$administrativeNames,-1))
e$shortloc<-as.factor(word(e$locality, 2,sep = fixed("|")))
e$NEWloc<-as.factor(ifelse(e$experiment=="Browsed",
                           str_sub(e$siteNumber,start=1,end=-2),str_sub(e$siteNumber,start=1,end=-3))) # use of the short code (as siteNumber) for LOCALITY
e$siteNumber<-as.factor(e$siteNumber)
e$experiment<-as.factor(e$experiment)
e$year<-as.factor(e$year)
e$SiteStation<-as.factor(paste(e$siteNumber,e$stationNumber,sep="_"))
e$SiteStationYear<-as.factor(paste(e$siteNumber,e$stationNumber,e$year,sep="_"))

e2<-e[,c("NaTron_datasetName","adminNames","shortloc","NEWloc","siteNumber","SiteStation","SiteStationYear","experiment","year",
      "individualCount","scientificName")]
head(e2)

e3<-droplevels(e2[e2$scientificName!="Pin height",])
dim(e3)
head(e3)

evac<-droplevels(e3[e3$scientificName=="Vaccinium myrtillus",])
head(evac)

## trial for trondelag with only vaccinium
evacT<-evac[evac$NaTron_datasetName=="SUSTHERB Tr\xf8ndelag",]
dev.off() ## annule les parametrages graphiques
  hist(evacT$individualCount)

hist(log(evacT$individualCount))
par(mfrow=c(2,1))
hist(log(evacT$individualCount)[evacT$experiment=="Unbrowsed"],main="Unbrowsed",col=3)
hist(log(evacT$individualCount)[evacT$experiment=="Browsed"],main="Browsed",col="darkgreen")

