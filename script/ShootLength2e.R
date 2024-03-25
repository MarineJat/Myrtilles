setwd("~/Desktop/script and more/script")

library(stringr)
library(dplyr)
library(ggplot2)
library(plyr)
library(visreg)
#
library(nlme)
library(MASS)
library(lme4)
library(MuMIn)
library(corrplot)

rm(list=ls(all=TRUE))

options(na.action = "na.fail")

d<-read.delim("summerBrowsing.csv",sep=";")
head(d)
names(d)
#write.table(d[d$shoot.length..centimeter.>60&!is.na(d$shoot.length..centimeter.),],"possibleoutliers.csv",sep=";")


#renaming column to have shorter names
d$adminNames<-as.factor(word(d$administrativeNames,-1))
d$shortloc<-as.factor(word(d$locality, 2,sep = fixed("|")))
d$NEWloc<-as.factor(ifelse(d$experiment=="Browsed",
                           str_sub(d$siteNumber,start=1,end=-2),str_sub(d$siteNumber,start=1,end=-3)))
# use of the short code (as siteNumber) for LOCALITY
d$siteNumber<-as.factor(d$siteNumber)
d$Region<-as.factor(mapvalues(d$NaTron_datasetName,from=unique(d$NaTron_datasetName),c("Hedmark","Telemark","Tingvoll","Trondelag")))
d$longi<-as.numeric(gsub(",",".",d$decimalLongitude))
d$lati<-as.numeric(gsub(",",".",d$decimalLatitude))
d$exp<-as.factor(d$experiment)
d$year<-as.factor(d$year)
d$yearc<-as.numeric(as.character(d$year)) ## year as a continuous variable
d$recNum<-d$recordNumber
d$recNum_0<-ifelse(d$recNum==0,1,0) # indicator of "0" on recordnumber
d$recNum_NA<-ifelse(is.na(d$recNum)==T,1,0) # indicator of "NA" on recordnumber

d$nbframe_occCOUNT<-d$number.of.frames.per.occurrence..count.
d$nbframe<-d$number.of.frames.per.plot..count
d$nbframe_NA<-ifelse(is.na(d$nbframe)==T,1,0) # indicator of NA on nbframe per plotcount
d$nbframe_0<-ifelse(d$nbframe==0,1,0) # indicator of "0" on nbframe per plotcount
d$nbflowers<-d$number.of.flowers..count.
d$BUB<-as.factor(d$browsed..boolean.) ## a factor with yes and no
d$BUB[!(d$BUB%in%c("Yes","No"))]<-NA
d$BUB<-droplevels(d$BUB)
d$BUB01<-as.numeric(as.character(mapvalues(d$BUB,c("No","Yes"),c(0,1))))
# a binary variable with 0=UB and 1 =B
d$BUB01_NA<-ifelse(is.na(d$BUB01)==T,1,0) # indicator of NA on BUB01
d$nbleaves<-as.numeric(d$number.of.leaves..count.)
d$shootlg<-d$shoot.length..centimeter.
d$shootlg_NA<-ifelse(is.na(d$shootlg)==T,1,0) # indicator of NA on shootlg
d$shootlg_0<-ifelse(d$shootlg==0,1,0) # indicator of 0 on shootlg
#
## creation of new composite variables
d$SiteStation<-as.factor(paste(d$siteNumber,d$stationNumber,sep="_"))
d$SiteStationYear<-as.factor(paste(d$siteNumber,d$stationNumber,d$year,sep="_"))
d$SiteYear<-as.factor(paste(d$siteNumber,d$year,sep="_"))


## choose the region
#d<-droplevels(d[d$Region=="Tingvoll",])
#d<-droplevels(d[d$Region=="Telemark",])
#d<-droplevels(d[d$Region=="Hedmark",])

d<-droplevels(d[d$Region=="Trondelag",])


# looking at the file to remove the columns with NA
summary(d)

## datafile with all modalities per station per site per year
fic0<-expand.grid(year=unique(d$year),stationNumber=unique(d$stationNumber),NEWloc=unique(d$NEWloc),exp=unique(d$experiment))
fic0$exp<-mapvalues(fic0$exp,from=c("Browsed","Unbrowsed"),c("B","UB"))
fic0$siteNumber<-paste(fic0$NEWloc,fic0$exp,sep="")
fic0$SiteStationYear<-as.factor(paste(fic0$siteNumber,fic0$stationNumber,fic0$year,sep="_"))
dim(fic0)
## for Trondelag, fic0 has 2400 lines= 8 years, 30 sites (15 localities x 2 experimental modalities B /UB) and 10 stations per site

## datafile with all modalities per site per year
fic1<-expand.grid(year=unique(d$year),NEWloc=unique(d$NEWloc),exp=unique(d$experiment))
fic1$exp<-mapvalues(fic1$exp,from=c("Browsed","Unbrowsed"),c("B","UB"))
fic1$siteNumber<-paste(fic1$NEWloc,fic1$exp,sep="")
fic1$SiteYear<-as.factor(paste(fic1$siteNumber,fic1$year,sep="_"))
dim(fic1)
## for Trondelag, fic1 has 240 lines= 8 years, 30 sites (15 localities x 2 experimental modalities B /UB)

# creating a smaller file with columns of interest and short column names

d2<-d[,c("Region","NEWloc","siteNumber","SiteYear","SiteStationYear",
         "exp","year","yearc","scientificName","recNum","recNum_0","recNum_NA",
         "BUB01","BUB01_NA","nbframe","nbframe_NA","nbframe_0","shootlg","shootlg_0","shootlg_NA")]
str(d2)

# creating a file with only vaccinium
dvac<-droplevels(d2[d2$scientificName=="Vaccinium myrtillus",])
head(dvac,25)
summary(dvac)

### NB : distribution of the nb of record (should be 10...)
table(dvac$recNum) # goes up to 36
plot(table(dvac$recNum))
table(dvac$recNum[dvac$nbframe>124],dvac$nbframe[dvac$nbframe>124])

##  verifying that the nb of unique modalities
length(unique(dvac$SiteStationYear)) ## there are less than 2400
length(unique(dvac$SiteStation))
length(unique(dvac$SiteYear))
length(unique(dvac$siteNumber))
length(unique(dvac$NEWloc))

str(dvac)
summary(dvac)


### aggregating the values of RESPONSE variable at the station level, so with "SiteStationYear"

## 1. for shoot length

temp_mean<-with(dvac,aggregate(shootlg,list(SiteStationYear=dvac$SiteStationYear),mean,na.rm=TRUE))
colnames(temp_mean)<-c("SiteStationYear","meanshootlg")
mergefile1<-merge(fic0,temp_mean,by="SiteStationYear",all.x=T)
head(mergefile1,50)
dim(mergefile1)
dim(mergefile1[is.na(mergefile1$meanshootlg)==T,])
table(is.na(temp_mean$meanshootlgth))

# 2. for nbframe_plotCOUNT

temp_frame<-with(dvac,aggregate(nbframe,list(SiteStationYear=dvac$SiteStationYear),max,na.omit=T))
colnames(temp_frame)<-c("SiteStationYear","nbframe")
mergefile2<-merge(mergefile1,temp_frame,by="SiteStationYear",all.x=T)
head(mergefile2)
dim(mergefile2)
dim(mergefile2[is.na(mergefile2$nbframe)==T,])

# 3. for BUB

temp_BUB<-with(dvac,aggregate(BUB01,list(SiteStationYear=dvac$SiteStationYear),sum,na.rm=TRUE))
colnames(temp_BUB)<-c("SiteStationYear","BUBsum")
mergefile3<-merge(mergefile2,temp_BUB,by="SiteStationYear",all.x=T)
head(mergefile3)
dim(mergefile3)
dim(mergefile3[is.na(mergefile3$BUBsum)==T,])

#4. number of myrtillus individuals measured in each station (= quadrat) a given year
maxpied<-aggregate(dvac$recNum,list(SiteStationYear=dvac$SiteStationYear),max)
colnames(maxpied)<-c("SiteStationYear","maxpied")
mergefile4<-merge(mergefile3,maxpied,by="SiteStationYear",all.x=T)
dim(mergefile4)
head(mergefile4)

mf<-mergefile4
str(mf)
mf$propBUB<-mf$BUBsum/mf$maxpied
head(mf)
hist(mf$propBUB)
hist(log(mf$propBUB/(1-mf$propBUB)))

hist(mf$nbframe)
hist(log(mf$nbframe))

names(mf)
summary(mf)
mf<-droplevels(mf)

## beware: when nbframe=125 and maxpied = 0 or NA, this means that there was no myrtillus
## we create an indicator for these cases

droplevels(dvac[dvac$SiteStationYear=="SLUB_7_2008",])
na.omit(mf[mf$nbframe>125,])
mf[mf$nbframe==125&is.na(mf$nbframe)==F,]
dim(mf[mf$nbframe==125&is.na(mf$nbframe)==F,])
mf125<-plot(mf$nbframe,mf$maxpied)



## Paired boxplot
### Creation of 8 boxplots representing the "shoot length" according to the treatment "B" "UB" according to the year
library(dplyr)
library(ggplot2)

### 2008 

# Filter data for 2008
df_2008 <- subset(mergefile4, year == "2008")

# Create the ggplot with color="NEWloc"
ggplot(df_2008, aes(x=exp, y=meanshootlg, colour=NEWloc))+
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2008", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 

# Create the ggplot with color="exp"
ggplot(df_2008, aes(x=exp, y=meanshootlg, colour=exp))+
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2008", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 
  
### 2010 
df_2010 <- subset(mergefile4, year == "2010")
ggplot(df_2010, aes(x=exp, y=meanshootlg, colour=exp))+ # possible also with colour="NEWloc"
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2010", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 

### 2012 
df_2012 <- subset(mergefile4, year == "2012")
ggplot(df_2012, aes(x=exp, y=meanshootlg, colour=exp))+ # possible also with colour="NEWloc"
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2012", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 

### 2014
df_2014 <- subset(mergefile4, year == "2014")
ggplot(df_2014, aes(x=exp, y=meanshootlg, colour=exp))+ # possible also with colour="NEWloc"
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2014", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 

### 2016 
df_2016 <- subset(mergefile4, year == "2016")
ggplot(df_2016, aes(x=exp, y=meanshootlg, colour=exp))+ # possible also with colour="NEWloc"
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2016", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 

### 2019
df_2019 <- subset(mergefile4, year == "2019")
ggplot(df_2010, aes(x=exp, y=meanshootlg, colour=exp))+ # possible also with colour="NEWloc"
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2019", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 

### 2021
df_2021 <- subset(mergefile4, year == "2021")
ggplot(df_2021, aes(x=exp, y=meanshootlg, colour=exp))+ # possible also with colour="NEWloc"
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2021", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 

### 2023
df_2023 <- subset(mergefile4, year == "2023")
ggplot(df_2023, aes(x=exp, y=meanshootlg, colour=exp))+ # possible also with colour="NEWloc"
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length in 2021", x = "Treatment", y = "Mean shoot length") +
  theme_classic() 
         
### Creation of 2 boxplots representing the "shoot length" as a function of years, respectively for "B" and "UB"
# Browsed
df_B <- subset(mergefile4, exp == "B")
ggplot_B <- ggplot(mergefile4, aes(x=year, y=meanshootlg, colour=year))+  # colour = NEWloc also
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length Browsed", x = "Years", y = "Mean shoot length") +
  theme_classic() 


# Unbrowsed
df_UB <- subset(mergefile4, exp == "UB")
ggplot_UB <- ggplot(mergefile4, aes(x=year, y=meanshootlg, colour=year))+  # colour = NEWloc also
  geom_point()+
  theme(legend.position="none")+
  geom_boxplot(alpha=0)+
  geom_line(aes(group=NEWloc), colour="grey70")+
  labs(title = "Mean shoot length Unbrowsed", x = "Years", y = "Mean shoot length") +
  theme_classic() 

plot(ggplot_B)
plot(ggplot_UB)



