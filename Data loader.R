### This program brings in the forcings from the 2024 SWAN data set and the 2024
### ERA5 data set (comprising AirT and SeaT). The forcings variables are as follows:
###   XWindv indicates wind velocity W to E (ie strength of west wind) -- 1 hour
###   YWindv indicates wind velocity S to N (ie strength of south wind) -- 1 hour
###   aice is the sea ice fraction (relative square covering w/ ice) (%) -- 1 hour
###   Hsig is significant wave height (m) -- 1 hour
###   RTpeak is peak wave period (sec.) -- 1 hour
###   Tm01 is mean absolute wave period (sec.) -- 1 hour
###   Dspr is width of directional distribution (deg.) -- 1 hour
###   WavePower is ___? -- 1 hour
###   AirT (merged in below) is air temperature (Kelvin) -- 6 hours
###   SeaT (merged in below) is sea temperature (Kelvin) -- 6 hours
### It also brings in the erosion data.

library(purrr)
library(RColorBrewer)
library(lubridate)

### Bring in 1950-2004 SWAN forcings data files and combine into a single frame.
setyrs_frcs <- 1950:2004
frame <- data.frame()
for(i in seq_along(setyrs_frcs)){
  #fra <- read.table(paste0("E:\\My Drive\\Project support\\Coastal erosion\\Kaktovik forcings\\2024 SWAN data\\statistics_data_1950-2004\\SWAN_output\\SWAN_data_WP_",setyrs_frcs[i],".PNTF05K")) #Windows
  fra <- read.table(paste0("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Coastal erosion/Kaktovik forcings/2024 SWAN data/statistics_data_1950-2004/SWAN_output/SWAN_data_WP_",setyrs_frcs[i],".PNTF05K")) #MacOS
  frame <- rbind(frame,fra)
}
colnames(frame) <- c("Time","Tsec","Xp","Yp","XWindv","YWindv","aice","Hsig",
                     "RTpeak","Tm01","Dspr","WavePower")

### Format dates and place into new column: `stamp`
year <- as.numeric(substr(frame$Time,start=1,stop=4))
month <- as.numeric(substr(frame$Time,start=5,stop=6))
day <- as.numeric(substr(frame$Time,start=7,stop=8))
hour <- floor(frame$Tsec/60/60) %% 24 #Tsec in Frame reports seconds since midnight, April 30/May 1
minute <- ((frame$Tsec/60/60 - hour)*60) %% 1440
frame$stamp <- as.POSIXct(paste0(year,"-",month,"-",day," ",hour,":",minute),format="%Y-%m-%d %H:%M",tz="GMT")

### Bring in 2004-2018 SWAN forcings and tack onto frame
#frame2 <- readRDS("E:\\My Drive\\Project support\\Coastal erosion\\Summer 2023\\frame") #Windows
frame2 <- readRDS("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Coastal erosion/Summer 2023/frame") #MacOS
frame2 <- frame2[year(frame2$stamp)>2004,]
frame2 <- frame2[minute(frame2$stamp)==0,]
frame[((nrow(frame)+1):(nrow(frame)+nrow(frame2))),c("Time","Tsec","Xp","Yp","XWindv",
                                                     "YWindv","aice","Hsig","RTpeak",
                                                     "Tm01","Dspr","stamp","AirT",
                                                     "SeaT")] <-
  frame2[,c("Time","Tsec","Xp","Yp","XWindv","YWindv","aice","Hsig","RTpeak","Tm01",
                                                    "Dspr","stamp","AirT","SeaT")]
frame <- frame[!is.na(frame$stamp),]

### Bring in ERA5 forcings and format dates
fr <- data.frame()
for(i in seq_along(setyrs_frcs)){
  #f <- read.table(paste0("E:\\My Drive\\Project support\\Coastal erosion\\Kaktovik forcings\\2024 SWAN data\\statistics_data_1950-2004\\ERA5\\temp2m_sst_ERA5_",setyrs_frcs[i],".txt")) #Windows
  f <- read.table(paste0("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Coastal erosion/Kaktovik forcings/2024 SWAN data/statistics_data_1950-2004/ERA5/temp2m_sst_ERA5_",setyrs_frcs[i],".txt")) #MacOS
  fr <- rbind(fr,f)
}
colnames(fr) <- c("Time","AirT","SeaT")
fr$stamp <- as.POSIXct(paste0(as.character(floor(fr$Time)),round((fr$Time-floor(fr$Time))*100)),
                       tz="GMT",format="%Y%m%d%H")

### Merge ERA5 forcings into frame
ind <- match(fr$stamp,frame$stamp)
frame[ind,c("AirT","SeaT")] <- fr[,c("AirT","SeaT")]

### Bring in erosion
setyrs_ero <- c(1950,1955,1969,1975,1979,2000,2004,2006,2007,2008,2009,2011,2014,2015,
            2016,2017,2018,2019)
slidyrs_ero <- c(1950,1955,1969,1975,1979,2000,2004,2006,2007,2008,2009,2010,2011,2012,
             2014,2015,2016,2017,2018,2019)
table <- data.frame()
for(i in 1:17){
  if(i != 7){
    #tab <- read.table(paste0("E:\\My Drive\\Project support\\Coastal erosion\\Kaktovik erosion\\Kaktovik_bluff_",setyrs_ero[i],"to",setyrs_ero[i+1],".txt"),
    #                    header=TRUE,sep=",") #Windows
    tab <- read.table(paste0("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Coastal erosion/Kaktovik erosion/Kaktovik_bluff_",setyrs_ero[i],"to",setyrs_ero[i+1],".txt"),
                        header=TRUE,sep=",") #MacOS
  } else {
    #tab <- read.table(paste0("E:\\My Drive\\Project support\\Coastal erosion\\Kaktovik erosion\\Kaktovik_bluff_",setyrs_ero[i],"to",setyrs_ero[i+1],".txt"),
    #                  header=TRUE,sep="\t") #Windows
    tab <- read.table(paste0("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Coastal erosion/Kaktovik erosion/Kaktovik_bluff_",setyrs_ero[i],"to",setyrs_ero[i+1],".txt"),
                        header=TRUE,sep="\t") #MacOS
  }
  if(i == 15){
    tab <- rbind(rbind(tab[1:which(tab$TransectID == 36),],rep(NA,ncol(tab))),
                 tab[(which(tab$TransectID == 36) + 1):nrow(tab),]) 
    ### Insert a new row into tab for transect 36's missing right endpoint (7/13/2017)
    tab[which(tab$TransectID == 36)+1,] <- c(rep(NA,5),36,NA,NA,"07/13/2017",rep(NA,28))
    rownames(tab) <- NULL
  }
  table <- rbind(table,tab[c("IntersectX","IntersectY","TransectID","Uncertainty",
                             "NSM","EPR","EPRunc","Feature","ShorelineID","Date_",
                             "Year_")])
  assign(paste0("table_",setyrs_ero[i+1]),tab)
}
table$Feature <- as.factor(table$Feature)
table$TransectID <- as.factor(table$TransectID)

### Format dates
table$ShorelineID <- as.POSIXct(table$ShorelineID,format="%m/%d/%Y",tz="GMT") #Assume UTC time
table$ShorelineID_doy <- as.numeric(strftime(table$ShorelineID,format="%j"))
table$ShorelineID <- year(table$ShorelineID)
table$stamp <- as.POSIXct(paste0(table$ShorelineID,table$ShorelineID_doy,"120000"),
                          format="%Y%j%H%M%S",tz="GMT") #Readings are presumed to happen at noon
table$Date_ <- substr(table$Date_,7,10)
table$Date_[table$Date_=='004'] <- 2004
table$Date_[table$Date_=='006'] <- 2006

#Remove duplicate stamp-transect observations
table <- table[!duplicated(table[, c("stamp", "TransectID")]),]

### Merge erosion into frame and save converted-to-numeric frame
rownames(table) <- NULL
frame[,c("ShorelineID_doy",paste0("EPR_",1:161),paste0("EPRunc_",1:161),paste0("Feature_"),1:161)] <- #,"setyr")] <- 
  cbind(NA,matrix(rep(NA,161*3*nrow(frame)),nrow(frame),161*3),NA)
rownames(frame) <- NULL
pb <- txtProgressBar()
for(i in 1:161){
  tab <- table[table$TransectID == i,]
  frame_ind <- match(tab$stamp,frame$stamp)
  frame_ind <- frame_ind[!is.na(frame_ind)]
  tab_ind <- match(frame$stamp,tab$stamp)
  tab_ind <- tab_ind[!is.na(tab_ind)]
  frame[unique(frame_ind[which(!is.na(frame_ind))]),"ShorelineID_doy"] <- tab[tab_ind,"ShorelineID_doy"]
  frame[frame_ind,paste0("EPR_",i)] <- tab[tab_ind,"EPR"]
  frame[frame_ind,paste0("EPRunc_",i)] <- tab[tab_ind,"EPRunc"]
  frame[frame_ind,paste0("Feature_",i)] <- tab[tab_ind,"Feature"]
  setTxtProgressBar(pb,i/161)
}
close(pb)
frame_names <- colnames(frame)[colnames(frame) != "stamp"]
frame_temp <- as.data.frame(sapply(frame[,frame_names],as.numeric))
frame <- cbind(frame_temp,frame$stamp)
colnames(frame) <- c(colnames(frame)[1:500],"stamp")

#saveRDS(frame,file="E:\\My Drive\\Project support\\Coastal erosion\\Spring 2025\\frame") #Windows
saveRDS(frame,file="/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Coastal erosion/Spring 2025/frame") #MacOS
