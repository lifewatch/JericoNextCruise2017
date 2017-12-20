#script adapted from:
# Name: Amanda Winegardner, Xavier Giroux-Bougard, Bérenger Bourgeois, Emmanuelle Chrétien and Monica Granados 
# Date: March 2016
# Description: Multivariate analyses II
# Dataset: File names are "DoubsSpe.csv", "DoubsEnv.csv"
# Notes: Material in R script obtained from 
#        Borcard D., Gillet F. et Legendre P., 2011. Numerical Ecology with R. Springer.
#

#And adapted from Reinhoud De blok
#***************************************************************************************#
#http://qcbs.ca/wiki/r_workshop10
# 0. Loading the data ####

rm(list=ls())
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("plyr")
install.packages("rJava")
install.packages("venneuler")

# To do so, you need to install the package "devtools" and to perform these few steps
install.packages("devtools")
library(devtools)

# Install mvpart, MVPARTwrap and rdaTest from github
devtools::install_github("cran/mvpart")
devtools::install_github("cran/MVPARTwrap")
devtools::install_github("philippec/fonctions_R_git/rdaTest")

# Load packages
library(vegan)
library(mvpart)
library(MVPARTwrap)
library(rdaTest)
library(labdsv)
library(MASS)
library(plyr)
library(stats)
library(magrittr)
library(tibble)
# Species community data frame (FCM clusters per station): 

df<-read.csv("jerico_LW_FCM_ExampleforRDA", header = TRUE, dec = ".", sep = ",")
names(df)
df$StationCode<-as.factor(df$StationCode)
levels(df$StationCode)
#define that OriginalStationName must be taken as the rownames
df<- df %>% remove_rownames %>% column_to_rownames(var="StationCode")
#make two subsets: species per station & environmental parameter per station
#start with checking the column names
names(df)
spe<-subset(df, select=c(10, 13:143))
pig<-subset(df, select=c(151:170))
env<-subset(df, select=c(145:149, 171:173, 179, 181,185:187, 189))
#check for NA values
pig<-pig[-20,]#row 20 is full of NA values, delete this row in all dataframes
env<-env[-20,]
spe<-spe[-20,]

#check the new dataframe; df is too big for R so write it to a .csv file
write.csv(env, "LWenv.csv")
write.csv(pig, "LWpigm.csv")
write.csv(spe, "LWspe.csv")
####eerst eens kijken of je clusters gecorreleerd zijn

##in case you wanne remove a row, can happen with this code:
#df <- df[-8,] # site number 8 contains no species and must be removed from your dataframe

#***************************************************************************************#

# 1.Explore the data ####

## 1.1 Species data ####
names(spe)
dim(spe)
str(spe)
head(spe)
summary(spe) 

### Species distibution, all species confounded

(ab <- table(unlist(spe)))
barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=grey(5:0/5))

### Number of absences
sum(spe==0)

### Proportion of zeros in the community data set
sum(spe==0)/(nrow(spe)*ncol(spe))

## Apply an appropriate transformation for community composition data

?decostand # this function provide some standardiztion methods for community domposition data (see method options)

# Hellinger transformation = standardisation
spe.hel<-decostand(spe, method="hellinger", na.rm=TRUE)#ignore missing values in row or column standardization
pig.hel<-decostand(pig, method="hellinger", na.rm=TRUE)
###
mosthighlycorrelated <- function(mydataframe,numtoreport)
{
  # find the correlations
  cormatrix <- cor(mydataframe)
  # set the correlations on the diagonal or lower triangle to zero,
  # so they will not be reported as the highest ones:
  diag(cormatrix) <- 0
  cormatrix[lower.tri(cormatrix)] <- 0
  # flatten the matrix into a dataframe for easy sorting
  fm <- as.data.frame(as.table(cormatrix))
  # assign human-friendly names
  names(fm) <- c("First.Variable", "Second.Variable","Correlation")
  # sort and print the top n correlations
  head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
}
correlationspe<-as.data.frame(mosthighlycorrelated(spe.hel, 10000))
#a lot of clusters correlate 100% with each other, due to the large number of clusters. Probably the clusters vary in  a similar way over the stations 
correlationpig<-as.data.frame(mosthighlycorrelated(pig.hel, 10))

## 1.2 Environmental data ####
names(env)
dim(env)
str(env)
head(env)
summary(env)
pairs(env, main="Bivariate Plots of the Environmental Data" )#collinearity of explanatory variables is present!
str(env)
# Standardization of the 11 environmental data
#decostand acts similar to scale()
env.z <- decostand(env, method="standardize", na.rm=TRUE)#ignore missing values in row or column standardization
apply(env.z, 2, mean) # the data are now centered (means~0)
apply(env.z, 2, sd)   # the data are now scaled (standard deviations=1)

# Note: explanatory variables expressed in different units
#       should be standardized before computing distances measures and ordination analysis. 
#       It is however not necessary for multivariate tree analysis.


#***************************************************************************************#
# 2. Canonical Analysis
## 2.1. Redundancy Analysis (RDA) ####

### Prepare the data

#delete column with only NA 
correlationenv<-as.data.frame(mosthighlycorrelated(env.z, 10))
names(env.z)
#env.z<-subset(env.z, select=c(1,3,5:10, 12:14))

### Run the RDA of all exlanatory variables of env on species abundances
?rda

spe.rda <- rda(spe.hel~., data=env.z, na.action = na.exclude) # the dot "." means that all variables of env.z will be included.
pig.rda<- rda(pig.hel~., data=env.z, na.action = na.exclude) 
### Extract the results
summary(spe.rda, display=NULL)
#These results contain the proportion of variance of Y explained by the X variables (constrained proportion, 80.7% here), the unexplained variance of Y (unconstrained proportion, 19.29% here) and then summarize the eigenvalues, the proportions explained and the cumulative proportion of each canonical axis (each canonical axis = each constraining variable, in this case, the environmental variables from env). 
summary(pig.rda, display=NULL)#75.8% explained variance
### Select the significant explanatory variables by forward selection
?ordiR2step
ordiR2step(rda(spe.hel~1, data=env.z, na.action = na.exclude), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=10000)
##no env variable can explain the variance in data!!



ordiR2step(rda(pig.hel~1, data=env.z, na.action = na.exclude), scope= formula(pig.rda), direction= "stepwise", R2scope=TRUE, pstep=10000)
#Nutr_NO2 is an explanatory variable for pigments

env.signif <- subset(env.z, select = c("Nutr_NO2"))
pig.rda.signif<-rda(pig.hel~., data=env.signif)
### Extract the results
summary(pig.rda.signif)#in this case, nutr_NO2 only explains 17% of the variance.

### Calculate the adjusted R^2 of the RDA
(R2adj <- RsquareAdj(pig.rda.signif)$adj.r.squared)
#only 12% is explained
### Test the significance of the model and the significance of axis
?anova.cca
anova(pig.rda, step=1000)#model is not significant
anova(pig.rda.signif, step=1000, by="axis") #model is significant

### Plot the results
#### quick plots scaling 1 and 2
windows()
plot(pig.rda.signif, scaling=1, main="Triplot RDA (scaling 1)")#euclidean distance
windows()
plot(pig.rda.signif, scaling=2, main="Triplot RDA (scaling 2)")#distance based on their correlation

#### advanced plots scaling 1
windows()
plot(pig.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(pig.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(pig.rda.signif, display="species", choices=c(1), scaling=1),
       scores(pig.rda.signif, display="species", choices=c(2), scaling=1),
       col="black",length=0)
text(scores(pig.rda.signif, display="species", choices=c(1), scaling=1),
     scores(pig.rda.signif, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(pig.rda.signif, display="species", scaling=1)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(pig.rda.signif, display="bp", choices=c(1), scaling=1),
       scores(pig.rda.signif, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(pig.rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(pig.rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(pig.rda.signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#### advanced plots scaling 2
windows()
plot(spe.rda.signif, scaling=2, main="Triplot RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spe.rda.signif, display="species", choices=c(1), scaling=2)*2,
       scores(spe.rda.signif, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(spe.rda.signif, display="species", choices=c(1), scaling=2)*2.1,
     scores(spe.rda.signif, display="species", choices=c(2), scaling=2)*2.1,
     labels=rownames(scores(spe.rda.signif, display="species", scaling=2)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(spe.rda.signif, display="bp", choices=c(1), scaling=2),
       scores(spe.rda.signif, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(spe.rda.signif, display="bp", choices=c(1), scaling=2)+0.05,
     scores(spe.rda.signif, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(spe.rda.signif, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  

####Multivariate regression tree is also another option which is more robust




# 5. Other useful ordination methods####

?cca # Constrained Correspondence Analysis (CCA)

?metaMDS # Nonmetric Multidimesional Scaling

?CCorA # Canonical Correlation Analysis

help(coinertia, package=ade4) # Coinertia Analysis 

help(mfa, package=ade4) # Multiple Factorial Analysis

https://r-forge.r-project.org/R/?group_id=195  # packages AEM and PCNM contain functions to perform spatial analysis 

