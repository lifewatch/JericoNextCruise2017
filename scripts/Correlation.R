#Comparing data obtained from different sensors

#Script adapted from: http://little-book-of-r-for-multivariate-analysis.readthedocs.io/en/latest/src/multivariateanalysis.html

#combine data from the various sensors into one matrix (df)
#column 1 needs to contain all stations

#If you want to compare different variables that have different units, are very different variances, it is a good idea to first standardise the variables.
df <- as.data.frame(scale(df[2:14])) #similar standardisation to decostand function of vegan

library("car")
scatterplotMatrix(df[2:6])#compare between sensors in column 2 to 6
plot(df$S4, df$S5)#compare between two sensors (e.g. sensor 4 and 5)
text(df$S4, df$S5, df$S1, cex=0.7, pos=4, col="red")#Stations are also labelled

sapply(df[2:14],mean)
sapply(df[2:14],sd)

### calculating correlations for multivariate data ####
#keep in mind that correlations ONLY describe a linear relationship!
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

mosthighlycorrelated(df[2:14], 10)#top10 most correlated sensorscouples
library(Hmisc)
rcorr(as.matrix(df[2:14]), type="pearson")#linear correlation
rcorr(as.matrix(df[2:14]), type="spearman")#rand
library(PerformanceAnalytics)
chart.Correlation(df[2:14])
library(corrplot)
xx<-cor(df[2:14])
corrplot(xx, type="upper", order="hclust")
