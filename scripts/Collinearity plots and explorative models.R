-+#### Install packages ####
  require(mgcv) 
  require(lattice)
  
  #### Dataset####
  File <- c("Cruise_data_2017_V8_complete.csv")
  dataframe <-  read.csv(File, header = TRUE)
  dataframe <-  as.data.frame (dataframe)
  head(dataframe)
  dim(dataframe)
  
  ####Generate Cleveland dotplots####
  dotchart(dataframe$Nutr_NO3,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$Nutr_SiO2,labels=dataframe$stationcode,cex=.7)
  
  dotchart(dataframe$Pigm_Fucoxanthin,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$Pigm_Chlorophyll_c3,labels=dataframe$stationcode,cex=.7)
  
  dotchart(dataframe$CTD_Salinity,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$UW_Temperature,labels=dataframe$stationcode,cex=.7)
  
  dotchart(dataframe$nano_conc_vliz,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$nano_concRWS,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$Nano_concentration_CNRS,labels=dataframe$stationcode,cex=.7)
  
  dotchart(dataframe$micro_conc_vliz,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$micro_concRWS,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$Micro_concentration_CNRS,labels=dataframe$stationcode,cex=.7)
  dotchart(dataframe$Zoopl_Calanoidea,labels=dataframe$stationcode,cex=.7)
  
  dataframe.selection <- dataframe[,-c(2:5)]
  head(dataframe.selection)
  
  
  #### k means clustering####
  d <- dist(dataframe, method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward") 
  plot(fit) # display dendogram
  groups <- cutree(fit, k=5) # cut tree into 5 clusters
  # draw dendogram with red borders around the 5 clusters 
  rect.hclust(fit, k=4, border="red")
  
  d <- dist(dataframe.selection, method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward") 
  plot(fit) # display dendogram
  groups <- cutree(fit, k=5) # cut tree into 5 clusters
  # draw dendogram with red borders around the 5 clusters 
  rect.hclust(fit, k=4, border="red")
  
  
  #### Delete 'outlier' location####
  test <- NULL
  test <- which(dataframe$stationcode=="37")
  dataframe <- dataframe[-test,]
  
  test <- NULL
  test <- which(dataframe$stationcode=="38")
  dataframe <- dataframe[-test,]
  
  #test <- NULL
  #test <- which(dataframe$stationcode=="2")
  #dataframe <- dataframe[-test,]
 
  ####Correlations####
   panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  
  
  
  pairs(~Nutr_NH4 + Nutr_NO2 + Nutr_NO3	+ Nutr_NOX + Nutr_PO4	+ Nutr_SiO2, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Nutrient Matrix", na.action = na.omit, cex.labels = 2)
  
  # Pigm_Pheophorbide_a,Pigm_19hexanoylooxyfucoxanthine	 , Pigm_Echinenone, Pigm_Pheophytin_a hadden geen data en dus niet meegenomen in deze pairplotr uit.
  pairs(~Pigm_Chlorophyll_c3	+ Pigm_Chlorophyll_c2	+	Pigm_Chlorophyllide_a	+	Pigm_Peridine	+	Pigm_19butanoyloxyfucoxanthin	+	Pigm_Fucoxanthin	+	Pigm_Neoxanthine	+	Pigm_Prasinoxanthin	+	Pigm_Diadinoxanthin	+	Pigm_Antheraxanthin	+	Pigm_Alloxanthin	+	Pigm_Diatoxanthin	+	Pigm_Zeaxanthin	+	Pigm_Chlorophyll_b	+	Pigm_Chlorophyll_a	+		Pigm_BetaCarotene, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Pigment Matrix", na.action = na.omit, cex.labels = 1.5)
  
  pairs(~CTD_Pressure + CTD_Temperature +	CTD_Conductivity +	CTD_Obs +	CTD_Altimeter +	CTD_Par +	CTD_FlECO.AFL +	CTD_Sbeox +	CTD_Sbeox0PS +	CTD_Salinity +	CTD_SvCM +	CTD_Density, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="CTD Matrix", na.action = na.omit, cex.labels = 1.5)
  
  #Zoopl_Branchiopoda,Zoopl_Mollusca, Zoopl_Mysidae, Zoopl_Porcellanidae.zoe weinig data, dus niet meegenomen.
  pairs(~Zoopl_Amphipoda	 + Zoopl_Annelida +	Zoopl_Anomura +	Zoopl_Appendicularia +	Zoopl_Brachyura.zoe +	Zoopl_Calanoidea +	Zoopl_Caridae.zoe +	Zoopl_Chaetognatha +	Zoopl_Cirripeda.cypris +	Zoopl_Cirripeda.nauplius +	Zoopl_Cnidaria +	Zoopl_Ctenophora +	Zoopl_Cumacea +	Zoopl_Echinodermata +	Zoopl_Harpacticoida +	Zoopl_Noctiluca +	Zoopl_Pisces.egg +	Zoopl_Pisces.larvae, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Zooplankton Matrix", na.action = na.omit, cex.labels = 1.5)
  
  pairs(~Turb_Secchi + UW_id_nav	+ UW_OctansHeading	+ UW_OdomDepth200khz	+ UW_CourseOverGround	+ UW_SpeedOverGround	+ UW_OdomDepth33khz	+ UW_NavDepth50khz	+ UW_Speedlog	+ UW_FLRTchla	+ UW_Salinity	+ UW_Temperature	+ UW_RelativeHumidity	+ UW_WindDirection	+ UW_WindSpeed	+ UW_AirPressure	+ UW_SoundVelocity	+ UW_WaterFlow, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Underway Matrix", na.action = na.omit, cex.labels = 1.5)
  
  
  pairs(~pico_conc_vliz +	nano_conc_vliz +	micro_conc_vliz +	macro_conc_vliz +	nano_tflr_vliz +	micro_tflr_vliz +	macro_tflr_vliz +	pico_tflr_vliz, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Phytoplankton Matrix 1", na.action = na.omit, cex.labels = 1.5)
  
  
  pairs(~pico_conc_vliz +	nano_conc_vliz +	micro_conc_vliz +	macro_conc_vliz +	nano_tflr_vliz +	micro_tflr_vliz +	macro_tflr_vliz +	pico_tflr_vliz + Synechococcus_concentration_CNRS +	Pico_concentration_CNRS	+ Nano_concentration_CNRS +	Micro_totFLRed_CNRS	+ Synechococcus_totTFLRed_CNRS	+ Pico_totFLR_CNRS	+ Nano_totFLRed_CNRS	+ Micro_totFLRed_CNRS	, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Phytoplankton Matrix VLIZ-CNRS", na.action = na.omit, cex.labels = 1.5)
  
  pairs(~pico_conc_vliz +	nano_conc_vliz +	micro_conc_vliz +	macro_conc_vliz +	nano_tflr_vliz +	micro_tflr_vliz +	macro_tflr_vliz +	pico_tflr_vliz 	+ pico_concRWS	+ nano_concRWS	+ micro_concRWS	+ Synechoccocus_concRWS	+pico_totFLRRWS	+ nano_totFLRRWS	+ micro_totFLRRWS	+ synechoccocus_totFLRRWS	, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Phytoplankton Matrix VLIZ-RWS", na.action = na.omit, cex.labels = 1.5)
  
  pairs(~ Synechococcus_concentration_CNRS +	Pico_concentration_CNRS	+ Nano_concentration_CNRS +	Micro_totFLRed_CNRS	+ Synechococcus_totTFLRed_CNRS	+ Pico_totFLR_CNRS	+ Nano_totFLRed_CNRS	+ Micro_totFLRed_CNRS	+ pico_concRWS	+ nano_concRWS	+ micro_concRWS	+ Synechoccocus_concRWS	+pico_totFLRRWS	+ nano_totFLRRWS	+ micro_totFLRRWS	+ synechoccocus_totFLRRWS	, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Phytoplankton Matrix RWS-CNRS", na.action = na.omit, cex.labels = 1.5)
  
    pairs(~Pigm_Chlorophyll_c3	+ Pigm_Chlorophyll_c2	+	Pigm_Chlorophyllide_a	+	Pigm_19butanoyloxyfucoxanthin + Pigm_Fucoxanthin	+	Pigm_Alloxanthin	+	Pigm_Zeaxanthin	+	Pigm_Chlorophyll_b	+	Pigm_Chlorophyll_a	+ Pigm_Diadinoxanthin +	Fv.Fm	+ Sigma	+ NPQ	+ NSV	+ rETR	+ RCII	+ JVPII.Sigma	+ aLHII	+ JVPII.aLHII	, data=dataframe,
        lower.panel=panel.smooth, upper.panel=panel.cor, 
        pch=20, main="Selectie Pigment  & FRRF Matrix", na.action = na.omit, cex.labels = 1.5)
    
    pairs(~Nutr_NH4 + Nutr_NO2 + Nutr_NO3	+ Nutr_NOX + Nutr_PO4	+ Nutr_SiO2 + Fv.Fm	+ Sigma	+ NPQ	+ NSV	+ rETR	+ RCII	+ JVPII.Sigma	+ aLHII	+ JVPII.aLHII, data=dataframe,
          lower.panel=panel.smooth, upper.panel=panel.cor, 
          pch=20, main="Nutrient Matrix", na.action = na.omit, cex.labels = 1.5)
 
     
  ####Simple plots####
  plot(dataframe$MS_biogeo05_dist_shore_5m, dataframe$Zoopl_Calanoidea)
  plot(dataframe$MS_biogeo05_dist_shore_5m, log10(dataframe$Zoopl_Calanoidea))
  
  plot(dataframe$MS_biogeo05_dist_shore_5m, log10(dataframe$Diatoms_Reinhoud))
  
  plot(dataframe$Diatoms_Reinhoud, dataframe$Zoopl_Calanoidea)
  plot(log10(dataframe$Diatoms_Reinhoud), log10(dataframe$Zoopl_Calanoidea))

  
  ####Exploratory models####
  ##predict zooplankton based on nutrients
  M <- gam(log10(Zoopl_Calanoidea) ~ s(MS_biogeo05_dist_shore_5m, k = 4) + s(Nutr_NO3, k = 4) + s(Nutr_SiO2, k = 4) + s(CTD_Salinity, k = 4),data=dataframe)
  summary(M)
  anova(M)
  plot(M,page=1)
  plot(resid(M))
  AIC(M)
  gam.check(M) 
  
  
  #predict zooplankton based on nutrients and latitude and longitude
  M <- gam(log10(Zoopl_Calanoidea) ~ s(MS_biogeo05_dist_shore_5m, k = 4) + s(Nutr_NO3, k = 4) + s(Nutr_SiO2, k = 4) + s(CTD_Salinity, k = 4) + s(Longitude, Latitude, k = 10),data=dataframe)
  summary(M)
  anova(M)
  plot(M,page=1)
  plot(resid(M))
  AIC(M)
  gam.check(M)
  
  plot.gam(M,select=5,main= "",scheme = 2, xlim =c(-1,5), ylim = c(49.5,53)) # don't know how to change the colors.
  
  map('worldHires','UK', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','France', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','Belgium', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','Netherlands', add = TRUE, fill =TRUE, col = "gray90")
  box()
  points(dataframe$Longitude, dataframe$Latitude, pch = 20)
  
  
  #predict nanoplankton based on nutrients
  M <- gam(log10(nano_conc_vliz) ~ s(MS_biogeo05_dist_shore_5m, k = 4) + s(Nutr_NO3, k = 4) + s(Nutr_SiO2, k = 4) + s(CTD_Salinity, k = 4),data=dataframe)
  summary(M)
  anova(M)
  plot(M,page=1)
  plot(resid(M))
  AIC(M)
  gam.check(M)
  
  
  #predict nanoplankton based on nutrients and latitude and longitude
  M <- gam(log10(nano_conc_vliz) ~ s(MS_biogeo05_dist_shore_5m, k = 4) + s(Nutr_NO3, k = 4) + s(Nutr_SiO2, k = 4) + s(CTD_Salinity, k = 4) + s(Longitude, Latitude, k = 10),data=dataframe)
  summary(M)
  anova(M)
  plot(M,page=1)
  plot(resid(M))
  AIC(M)
  gam.check(M) 
  
  plot.gam(M,select=5,main= "",scheme = 2, xlim =c(-1,5), ylim = c(49.5,53)) # don't know how to change the colors.
  
  map('worldHires','UK', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','France', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','Belgium', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','Netherlands', add = TRUE, fill =TRUE, col = "gray90")
  box()
  points(dataframe$Longitude, dataframe$Latitude, pch = 20)
  
  
  #Multitrofic model (use phytoplankton concentration to predict calanoids)
  M <- gam(log10(Zoopl_Calanoidea) ~ s(nano_conc_vliz, k = 4) + s(Nutr_NO3, k = 4) + s(Nutr_SiO2, k = 4) + s(CTD_Salinity, k = 4) + s(Longitude, Latitude, k = 10),data=dataframe)
  summary(M)
  anova(M)
  plot(M,page=1)
  plot(resid(M))
  AIC(M)
  gam.check(M)
  
  plot.gam(M,select=5,main= "",scheme = 2, xlim =c(-1,5), ylim = c(49.5,53)) # don't know how to change the colors.
  
  map('worldHires','UK', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','France', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','Belgium', add = TRUE, fill =TRUE, col = "gray90")
  map('worldHires','Netherlands', add = TRUE, fill =TRUE, col = "gray90")
  box()
  points(dataframe$Longitude, dataframe$Latitude, pch = 20)
  