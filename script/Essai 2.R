#### SCRIPT STAGE M2 ONLY VACCINIUM ####


############################# Chargement des packages ##########################

library(stringr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(vegan)
library(cluster)
library(ca)
library(MASS)
library(caMx)
library(FactoMineR)
library(ggmosaic)

############################# Chargement base de données #######################

rm(list=ls(all=TRUE))

setwd("~/Desktop/script and more/script")

d<-read.csv("summerBrowsing.csv",sep=";") # d = sheet summerBrowsing

# Renaming the qualitative variable in the "region" (= "NaTron_datasetName") column --> does not work
d <- d %>% mutate(NaTron_datasetName = ifelse(NaTron_datasetName == "SUSTHERB Tr<bf>ndelag", "SUSTHERB Trondelag", NaTron_datasetName))
#?????????#

# Show first rows of table to check
head(d)

# creating a file with only vaccinium
dvac<-droplevels(d[d$scientificName=="Vaccinium myrtillus",])

# Show first rows of table to check
head(dvac)

############################# EXPLORATORY STATISTICS ###########################


################################################################################
################################################################################
### SHOOT LENGTH
################################################################################
# 1. Differences between the treatment (Browsed, Unbrowsed) ?

# Creation of a dataframe in order to have the 3 columns with the 2 treatments et the shoot lenght
bubshootlength <- dvac %>% select(`experiment`, `shoot.length..centimeter.`)

# Renommer les colonnes 
bubshootlength <- bubshootlength %>% 
  rename(Treatment = experiment, ShootLength = `shoot.length..centimeter.`)

# Création du boxplot avec ggplot2
ggplot(bubshootlength, aes(x = Treatment, y = ShootLength)) +
  geom_boxplot(fill = "lightblue", color = "blue") +
  labs(title = "Comparison of Shoot Length",
       x = "Parameter",
       y = "Shoot Length (cm)")

### Analysis : The shoot length are higher when it's unbrowsed so it's similar to what Hallvard has found

################################################################################
# 2. Study of correlation between some studied variables and response variables 
# Select numeric variables for correlation analysis
tablongilati <- dvac %>%
  select(decimalLatitude, decimalLongitude, `shoot.length..centimeter.`)

# Check that the subset is correct
head(tablongilati)

# Calculate the correlation matrix
cor_matrix <- cor(tablongilati, use = "complete.obs")

# Check the structure of cor_matrix
str(cor_matrix)

# Create the heatmap
ggplot(data = as.data.frame(as.table(cor_matrix)), 
       aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Heatmap of Correlation between Variables",
       x = "Variables",
       y = "Variables")

################????????????????????????? MARCHE PAS 

### Analysis : According to the correlation heatmap, we can see that the longitude and the latitude both present a correlation with the shoot lenght
###           However, the correlation seems to be low around (0.1)


# 2.2 Linear regression between shoot length over the years for each region
# Creation of a new table with the 3 columns
tabregion = dvac%>%select("NaTron_datasetName", "year", "shoot.length..centimeter.")

# Remove missing values (NA) and LM
model <- lm(shoot.length..centimeter. ~ year * NaTron_datasetName, data = na.omit(tabregion))

# Representation of the LM with ggplot2
ggplot(na.omit(tabregion), aes(x = year, y = shoot.length..centimeter.)) +
  geom_point(aes(color = NaTron_datasetName)) +  # Colored dots by region
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Regression line
  facet_wrap(~NaTron_datasetName, scales = "free_y") +  # Facets for each group
  labs(title = "Linear regression between shoot length over the years for each region",
       x = "Years",
       y = "Shoot Lenghts",
       color = "Regions")  # Color legend

### Question : comment faire apparaître "SUSTHERB Trondelag en légende" ?

### Analysis : The 3 regions of interest, Hedmark, Telemark and Trondelag, show an increase in shoot lengths over the years (browsed or unbrowsed)
###           While Tingvoll shows a decrease in shoot length over time


# 2.3 Dispersion of Shoot Length Over Years by Experiment
# Creation of a new table with the 3 columns
tabBUB <- dvac %>% select("experiment", "year", "shoot.length..centimeter.")

# Remove missing values (NA) and LM
model1 <- lm(shoot.length..centimeter. ~ year * experiment, data = na.omit(tabBUB))

# Creating a scatter plot with "browsed" and "unbrowsed" colored differently
ggplot(na.omit(tabBUB), aes(x = year, y = shoot.length..centimeter., color = experiment)) +
  geom_point() +
  labs(title = "Dispersion of Shoot Length Over Years by Experiment",
       x = "Years",
       y = "Shoot Length",
       color = "Experiment")


################################################################################
# 3. Statistics to see if the variables "site number" (qualitative), "locality" (qualitative), "regions" (qualitative), "experiment" (qualitative), "station number" (from 1 to 10) and years "year" are correlated or not with the variable to be explained "shoot length"?

# 3.1 Pearson correlation (for quantitative variables)
# 3.1.1 Shoot length / Years
cor(dvac$shoot.length..centimeter., dvac$year, method="pearson")
### Analysis : A Pearson correlation coefficient of 0.1344572 indicates a weak positive correlation between the 2 variables. 
###           This means that when the value of one variable increases, the value of the other variable tends to increase as well, but the relationship is weak.

# 3.1.1 Shoot length / station number
complete_shootStation <- dvac[complete.cases(dvac$shoot.length..centimeter., dvac$stationNumber), ] # Exclude observations with missing values in both variables
correlation_shootStation <- cor(complete_shootStation$shoot.length..centimeter., complete_shootStation$stationNumber, method="pearson") # Calculate correlation
correlation_shootStation
### Analysis : A Pearson correlation coefficient of 0.00653058 indicates a very weak positive correlation between the 2 variables. 
###           This means that when the value of one variable increases, the value of the other variable tends to increase as well, but the relationship is very weak.


# 3.2 ANOVA or Kruskal-Wallis test (for qualitative variables)
# We use ANOVA for qualitative variables like "site number", "locality", "regions", and "experiment". 
# The null hypothesis of ANOVA is that the group means are equal

# ANOVA1 : shoot length and site number
# We replace the missing values with the average of the variable shoot.length..centimeter. (or any other imputation method you prefer) before creating the model
# Replace missing values with mean
dvac$shoot.length..centimeter.[is.na(dvac$shoot.length..centimeter.)] <- mean(dvac$shoot.length..centimeter., na.rm = TRUE)

# Check if there are any missing values left
sum(is.na(dvac$shoot.length..centimeter.))

# Convert siteNumber to factor
dvac$siteNumber <- as.factor(dvac$siteNumber)

# Check number of levels in siteNumber
print(table(dvac$siteNumber))

# Rerun the linear regression model
model_site <- lm(shoot.length..centimeter. ~ siteNumber, data = dvac)

# Rerun the ANOVA
anova_ShootSite <- anova(model_site)
print(anova_ShootSite)

### Analysis : The very low p-value (< 2.2e-16 < p=0.05) suggests that the effect of siteNumber on shoot.length..centimeter. is statistically significant
###            We can reject H0 according to which the group means are equal.
###            The siteNumber variable has a significant impact on the dependent variable shoot.length..centimeter.

# See which sites have an impact on the shoot length centimeter
summary(model_site)
### Analysis : some of them but not all


# ANOVA2 : shoot length and locality
modele_locality <- lm(shoot.length..centimeter. ~ locality, data = dvac)
anova_locality <- anova(modele_locality)
print(anova_locality)
### Analysis : significant effect *** of locality on shoot.length..centimeter

# See which localities have an impact on the shoot length centimeter
summary(modele_locality)
### Analysis : almost all the localities have an impact on the shoot length


# ANOVA3 : shoot length and regions
modele_NaTron_datasetName <- lm(shoot.length..centimeter. ~ NaTron_datasetName, data = dvac)
anova_NaTron_datasetName <- anova(modele_NaTron_datasetName)
print(anova_NaTron_datasetName)
### Analysis : significant effect *** of regions on shoot.length..centimeter

# See which regions have an impact on the shoot length centimeter
summary(modele_NaTron_datasetName)
### Analysis : we can see that Tingvoll doesn't really have an impact on the shoot length centimeter
### ?????????? il manque Hedmark


# ANOVA4 : shoot length and experiment
modele_experiment <- lm(shoot.length..centimeter. ~ experiment, data = dvac)
anova_experiment <- anova(modele_experiment)
print(anova_experiment)
### Analysis : significant effect *** of experiment on shoot.length..centimeter


# ANOVA5 : shoot length and station number
modele_stationNumber <- lm(shoot.length..centimeter. ~ stationNumber, data = dvac)
anova_stationNumber <- anova(modele_stationNumber)
print(anova_stationNumber)
### Analysis : low effect * of station number on shoot.length..centimeter


# ANOVA6 : shoot length and year
modele_year <- lm(shoot.length..centimeter. ~ year, data = dvac)
anova_year <- anova(modele_year)
print(anova_year)
### Analysis : significant effect *** of year on shoot.length..centimeter


# 3.3 Multiple linear regression (for variables mixtes) 
# Purpose: evaluate the combined effect of several independent variables on the dependent variable
modele_multiple <- lm(shoot.length..centimeter. ~ siteNumber + locality + NaTron_datasetName + experiment + stationNumber + year, data = dvac)
summary(modele_multiple)
### Analysis : significant effect *** of some site number on shoot.length..centimeter but a lot of NA
###           Very low p-value (p-value: < 2.2e-16) suggests that the overall model is statistically significant.
#?????????#




#### à changer shoot length par frames et vérifier les tableaux dataframe refaire avec % de Lila 



################################################################################
################################################################################
### NUMBER OF FRAMES PER PLOT (number.of.frames.per.plot..count.)
# Nb de carrés de 5cm x 5 cm necessaires pour arriver à mesurer 10 plants de myrtille. 
# Cela donne une indication de la densité des plants de myrtille.
# DENSITY
################################################################################

# 1. Differences between the treatment (Browsed, Unbrowsed) ?

# Creation of a dataframe in order to have the 3 columns with the 2 treatments and the nb of frames per plot
bubframes <- data.frame(
  Treatment = rep(c("Browsed", "Unbrowsed"), each = 50),
  Frames = c(rnorm(50, mean = 20, sd = 5),
             rnorm(50, mean = 25, sd = 5))
)

# Creation dof the boxplot with ggplot2
ggplot(bubframes, aes(x = experiment, y = number.of.frames.per.plot..count.)) +
  geom_boxplot(fill = "lightblue", color = "blue") +
  labs(title = "Comparison of Number of frames per plot",
       x = "Parameter",
       y = "Frames per plot")

### Analysis : 

################################################################################
# 2. Study of correlation between some studied variables and response variables 

# 2.1 Heatmap correlation for numeric variables (longitude and latitude) 

# Data
Longilati <- data.frame(decimalLatitude,decimalLongitude,number.of.frames.per.plot..count.)

# Calcul of the matrix correlation 
cor_matrix <- cor(Longilati)

# Creation of the correlation heatmap with ggplot2
ggplot(data = reshape2::melt(cor_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Correlation Heatmap",
       x = "Variables",
       y = "Variables")

### Analysis : According to the correlation heatmap, we can see that the longitude and the latitude both present a correlation with the shoot lenght
###           However, the correlation seems to be low around (0.1)

################################################################################
# 3. Statistics to see if the variables "site number" (qualitative), "locality" (qualitative), "regions" (qualitative), "experiment" (qualitative), "station number" (from 1 to 10) and years "year" are correlated or not with the variable to be explained "Nb of frames per plot"

# 3.1. Nb of frames per plot / station number
complete_framesStation <- d[complete.cases(d$number.of.frames.per.plot..count., d$stationNumber), ] # Exclude observations with missing values in both variables
correlation_framesStation <- cor(complete_framesStation$number.of.frames.per.plot..count., complete_framesStation$stationNumber, method="pearson") # Calculate correlation
correlation_framesStation
### Analysis : A Pearson correlation coefficient of 0.02 indicates a weak positive correlation between the 2 variables. 
###           This means that when the value of one variable increases, the value of the other variable tends to increase as well, but the relationship is very weak.


# 3.2. ANOVA or Kruskal-Wallis test (for qualitative variables)

# ANOVA1 : nb of frames per plot and site number
# We replace the missing values with the average of the variable number.of.frames.per.plot..count.(or any other imputation method you prefer) before creating the model
# Replace missing values with mean
d$number.of.frames.per.plot..count.[is.na(d$number.of.frames.per.plot..count.)] <- mean(d$number.of.frames.per.plot..count., na.rm = TRUE)

# Check if there are any missing values left
sum(is.na(d$number.of.frames.per.plot..count.))

# Convert siteNumber to factor
d$siteNumber <- as.factor(d$siteNumber)

# Check number of levels in siteNumber
print(table(d$siteNumber))

# Rerun the linear regression model
model_site <- lm(number.of.frames.per.plot..count. ~ siteNumber, data = d)

# Rerun the ANOVA
anova_FramesSite <- anova(model_site)
print(anova_FramesSite)

### Analysis : The very low p-value (< 2.2e-16 < p=0.05) suggests that the effect of siteNumber on nb of frames per plot is statistically significant
###            We can reject H0 according to which the group means are equal.
###            The siteNumber variable has a significant impact on the dependent variable shoot.length..centimeter.

# See which sites have an impact on the nb of frames per plot
summary(model_site)
### Analysis : some of them but not all


# ANOVA2 : nb of frames per plot and locality
modele_locality <- lm(number.of.frames.per.plot..count. ~ locality, data = d)
anova_locality <- anova(modele_locality)
print(anova_locality)
### Analysis : significant effect *** of locality on shoot.length..centimeter

# See which localities have an impact on the nb of frames per plot
summary(modele_locality)
### Analysis : almost all the localities have an impact on the nb of frames per plot


# ANOVA3 : nb of frames per plot and regions
modele_NaTron_datasetName <- lm(number.of.frames.per.plot..count. ~ NaTron_datasetName, data = d)
anova_NaTron_datasetName <- anova(modele_NaTron_datasetName)
print(anova_NaTron_datasetName)
### Analysis : significant effect *** of regions on the nb of frames per plot

# See which regions have an impact on the nb of frames per plot
summary(modele_NaTron_datasetName)
### Analysis : all the regions have an impact on the nb of frames per plot
### ?????????? il manque Hedmark


# ANOVA4 : nb of frames per plot and experiment
modele_experiment <- lm(number.of.frames.per.plot..count. ~ experiment, data = d)
anova_experiment <- anova(modele_experiment)
print(anova_experiment)
### Analysis : significant effect *** of experiment on the nb of frames per plot


# ANOVA5 : nb of frames per plot and station number
modele_stationNumber <- lm(number.of.frames.per.plot..count. ~ stationNumber, data = d)
anova_stationNumber <- anova(modele_stationNumber)
print(anova_stationNumber)
### Analysis : Strong effect of station number on the nb of frames per plot


# ANOVA6 : nb of frames per plot and year
modele_year <- lm(number.of.frames.per.plot..count. ~ year, data = d)
anova_year <- anova(modele_year)
print(anova_year)
### Analysis :  the year is not significant on the nb of frames per plot


# 3.3 Multiple linear regression (for variables mixtes) 
# Purpose: evaluate the combined effect of several independent variables on the dependent variable
modele_multiple <- lm(number.of.frames.per.plot..count. ~ siteNumber + locality + NaTron_datasetName + experiment + stationNumber + year, data = d)
summary(modele_multiple)
### Analysis : significant effect *** of some site number on the nb of frames per plot but a lot of NA
###           Very low p-value (p-value: < 2.2e-16) suggests that the overall model is statistically significant.
#?????????#

################################################################################
### Exploratory graphs

# 4. Preliminary charts
# 4.1. Box plots to compare "number.of.frames.per.plot..count." between the different modalities of “experiment”
boxplot(number.of.frames.per.plot..count. ~ experiment, data = d, main = "Boxplot by Experiment")
### Analysis : the number of squares necessary to have a minimum of 10 individuals increases when the study plot is grazed
# this shows the impact of herbivory because more squares are needed to have a greater density of blueberries

# 4.2. Scatter plot to visualize the relationship between "stationNumber" and "number.of.frames.per.plot..count."
plot(d$stationNumber, d$number.of.frames.per.plot..count., main = "Scatter Plot", xlab = "Station Number", ylab = "Number of Frames")
### Analysis: Approximately the same dispersion (same number of plots needed to obtain at least 10 individuals) for stations 1 to 10 per site (exclosure or control). 
# However, it is observed that the dispersion is generally higher on station 1.

# 5. Graphs for quantitative variables
# 5.1. Boxplot for "number.of.frames.per.plot..count." based on "experiment."
ggplot(d, aes(x = experiment, y = number.of.frames.per.plot..count.)) +
  geom_boxplot() +
  labs(title = "Distribution of the number of plots per experiment")
### Analysis: This graph confirms 1.1 (the number of plots increases with herbivory). To have the same density, more plots are needed.

# 5.2. Boxplot for "number.of.frames.per.plot..count." based on "locality."
ggplot(d, aes(x = locality, y = number.of.frames.per.plot..count.)) +
  geom_boxplot() +
  labs(title = "Distribution of the number of plots per locality") +
  
  # Add option for text size on the x-axis
  theme(axis.text.x = element_text(size = 8))  # Choose the size that suits you
### Analysis: This graph is good, but the x-axis text size needs to be reduced
### ???????

# 6. Graphs for cross relationships
# 6.1. Contingency table to visualize relationships between two qualitative variables.
ggplot(data = d, aes(x = product(experiment, siteNumber))) +
  geom_mosaic(aes(weight = value)) +
  labs(title = "Contingency table between entre experiment et siteNumber")

# 6.2. Scatterplot for "number.of.frames.per.plot..count." based on a qualitative variable
ggplot(d, aes(x = locality, y = number.of.frames.per.plot..count., color = experiment)) +
  geom_point(na.rm = FALSE) +
  labs(title = "Relationship between locality and number of plots")

################################################################################
# 7. Exploratory Stats Tests

# 7.1. Linear regression
lm_result <- lm(number.of.frames.per.plot..count. ~ siteNumber + locality + NaTron_datasetName + experiment + stationNumber + year, data = d)
summary(lm_result)
### Analysis: The model seems to be significant (p-value < 2.2e-16 < 0.05).
# Examining the residuals would provide additional information about the model quality (not done).
# The "siteNumber" variables seem to be specific indicators. Some sites have positive coefficients, while others have negative coefficients. This suggests significant variations between sites.
# Issue: The "locality" variables are mentioned as NA (not available) - appropriate for the model?
## ???????????
# Adjusted R-squared is approximately 0.3455, suggesting that the model explains about 34.55% of the total variation in the dependent variable.

# 7.2. Chi-square test for independence between "experiment" and "siteNumber."
chisq_result <- chisq.test(d$experiment, d$siteNumber)
summary(chisq_result)
### Analysis: p = 1 > 0.05, no significant association between the "experiment" and "siteNumber" variables.







