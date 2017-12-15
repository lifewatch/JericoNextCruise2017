rm(list=ls())
#set directory to this folder

library(dplyr)
#### start with reading all the datafile into R ####
#check column names and delete columns which are not relevant to the final complete combined dataset

#station locations, Underway and CTD data at 3m, pigment and zooplankton data
df1<-read.csv("Cruise_vliz_data.csv", header = TRUE, dec = ".", sep = ",")
names(df1)
df1$stationcode<-as.factor(df1$stationcode)
#delete last row since it only contains NA values
df1<-df1[-44,]

#FCM data from RWS
df2<-read.csv("Copy of SimonStevin_discrete FCM phytoplankon size fractions_RWS (00000002).csv", header = TRUE, dec = ".", sep = ",")
names(df2)
df2$stationcode<-as.factor(df2$stationcode)
#delete "filenames" and "AnalysisDate" columns from this dataset as this is not relevant for the complete dataset
df2<-select(df2, -filenames)
df2<-select(df2, -AnalysisDate)
names(df2)
#change column names to the standard format
colnames(df2)[1] <- "pico_conc_RWS"
colnames(df2)[2] <- "nano_conc_RWS"
colnames(df2)[3] <- "micro_conc_RWS"
colnames(df2)[4] <- "Synechoccocus_conc_RWS"
colnames(df2)[5] <- "pico_totFLR_RWS"
colnames(df2)[6] <- "nano_totFLR_RWS"
colnames(df2)[7] <- "micro_totFLR_RWS"
colnames(df2)[8] <- "synechoccocus_totFLR_RWS"


#FCM data from CNRS-LOG
df3<-read.csv("Copy of FCM and fluoroProbe cnrs.csv", header = TRUE, dec = ".", sep = ",")
df3$stationcode<-as.factor(df3$stationcode)
names(df3)
#delete column "Institut" since the institute is mentioned in the column names of the phytoplankton size classes
df3<-select(df3, -Institut)
names(df3)
#change column names to the standard names
colnames(df3)[3] <- "pico_conc_CNRS"
colnames(df3)[4] <- "nano_conc_CNRS"
colnames(df3)[5] <- "micro_conc_CNRS"
colnames(df3)[2] <- "Synechoccocus_conc_CNRS"
colnames(df3)[7] <- "pico_totFLR_CNRS"
colnames(df3)[8] <- "nano_totFLR_CNRS"
colnames(df3)[9] <- "micro_totFLR_CNRS"
colnames(df3)[6] <- "synechoccocus_totFLR_CNRS"


#FCM data from VLIZ-Ugent
df4<-read.csv("fractionate_vliz_FCM.csv", header = TRUE, dec = ".", sep = ";")
names(df4)
df4$stationcode<-as.factor(df4$stationcode)
#Delete columns "StartDate_station" and "EndDate_station
df4<-select(df4, -StartDate_station)
df4<-select(df4, -EndDate_station)
#change column names to the standard names
colnames(df4)[6] <- "pico_conc_VLIZ"
colnames(df4)[7] <- "nano_conc_VLIZ"
colnames(df4)[8] <- "micro_conc_VLIZ"
colnames(df4)[9] <- "macro_conc_VLIZ"
colnames(df4)[5] <- "pico_totFLR_VLIZ"
colnames(df4)[2] <- "nano_totFLR_VLIZ"
colnames(df4)[3] <- "micro_totFLR_VLIZ"
colnames(df4)[4] <- "macro_totFLR_VLIZ"
names(df4)

#Chemtaxdata according to Muylaert et al for the North Sea
df5<-read.csv("CHEMTAX_Muylaert_corrected.csv", header = TRUE, dec = ".", sep = ";")
names(df5)
colnames(df5)[1] <- "stationcode"
df5$stationcode<-as.factor(df5$stationcode)

#Fluoroprobe data from CNRS-LOG
df6<-read.csv("cnrs_fluoroprobe.csv", header = TRUE, dec = ".", sep = ",")
names(df6)
df6$stationcode<-as.factor(df6$stationcode)
#add institute name to the column Fluoroprobe
colnames(df6)[3] <- "FluoroProbeCNRS"
df6<-select(df6, -Institut)

#FRRF data from CNRS-LOG
df7<-read.csv("Synth-FRRF-TabStatC-VLIZ17.csv", header = TRUE, dec = ".", sep = ",")
names(df7)
#change codestation to stationcode
colnames(df7)[3] <- "stationcode"
df7$stationcode<-as.factor(df7$stationcode)
#need to select euphoticlayer only in the StatParam column for the complete station dataset
#check levels of StatParam
levels(df7$StatParam)#we need "meanoneuphoticlayer"
df7<-filter(df7, StatParam == "meanoneuphoticlayer")
#delete column codestatistic
df7<-select(df7, -codestatistic)

#WISP data from RWS
df8<-read.csv("Copy of wisp_RWS.csv", header = TRUE, dec = ".", sep = ",")
names(df8)
df8$stationcode<-as.factor(df8$stationcode)

#merge all datasets based on stationcode

df12 = full_join(df1, df2, by = "stationcode")
names(df12)
df123 = full_join(df12, df3, by = "stationcode")
names(df123)
df1234= full_join(df123, df4, by = "stationcode")
names(df1234)
df12345= full_join(df1234, df5, by = "stationcode")
names(df12345)
df123456= full_join(df12345, df6, by = "stationcode")
names(df123456)
df1234567= full_join(df123456, df7, by = "stationcode")
names(df1234567)
df12345678= full_join(df1234567, df8, by = "stationcode")
names(df12345678)

#add extra parameters to this dataset

#add distance to shore
library(sdmpredictors)
layercodes <- c("BO2_curvelmean_ss", "MS_biogeo05_dist_shore_5m" )
#"BO2_curvelmean_ss" = Mean surface current velocity (m/s)
#"MS_biogeo05_dist_shore_5m" = distance to shore unit Kilometers, resolution of used datalayer is 5arcminutes
env <- load_layers(layercodes, equalarea = FALSE)
df9 <- extract(env, cbind(df12345678$Longitude, df12345678$Latitude))

df9<-as.data.frame(df9)
#merge together
df123456789 <- cbind(df12345678,df9)
names(df123456789)

#add bottom depth to the dataset
#Bathymetry data can be retrieved from [Emodnet](http://portal.emodnet-bathymetry.eu/)
library(raster)
r <- raster("emodnet_mean.asc")
bathy <- data.frame(coordinates(r))
bathy$z <- as.vector(r)
bathy$x <- round(bathy$x, digits = 2)
bathy$y <- round(bathy$y, digits = 2)
require(dplyr)
bathy <- summarise(group_by(bathy, x, y),
                   z = mean(z))
#The bathymetry can also be selected for the specific station points
df123456789$x <- round(df123456789$Longitude, digits = 2)
df123456789$y <- round(df123456789$Latitude, digits = 2)
df123456789 <- left_join(df123456789, bathy)
df123456789 <- dplyr::select(df123456789, -x, -y)
#rename z to depth #always check if it is the correct column number!
names(df123456789)
colnames(df123456789)[130] <- "Depth"
#the unit of depth is meter
#write this complete cruise data set to a csv file
write.csv(df123456789, "Cruise_data_2017_V9_complete.csv") 

#after the export, column 1 needs to be removed as this contains the R rownumbers of the df123456789 which do not mean anything. 
