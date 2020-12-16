### This code opens a .gdb file and converts the information to a .csv file with 
### the same format of the Marine Cadastre 2015 .csv files
### AIS data can be downloaded from https://marinecadastre.gov/ais/
###
### Note 1: Pre-2014 data does not include VesselName, IMO, or CallSign for vessel protection. 
###         The variables are in the .csv file will be NULL.
### Note 2: Pre-2014 data MMSI numbers are encoded, the MMSIs here cannot be traced to the actual vessel. 
### Note 3: The Vessel and Voyage layers do not contain spatial information and therefore will return a warning when being read, 
###         warnings have been suppressed for those two lines.
### Note 4: The broadcast layer takes a long time to read, depending on the file, can take up to ~15 minute
###
### 16 December 2020 Vanessa ZoBell



# Variables to change: 
# idir: the directory of your .gdb file
# outDir: the directory where you want the .csv file to be saved
# dsn: the file name of your .gdb file

idir = 'G:/Ch.4_AmbientNoise/AmbientNoiseCINMS/AIS_ALL/MarineCadastre'  #director where your .gdb file lives
setwd(idir)
outDir = 'G:/Ch.4_AmbientNoise/AmbientNoiseCINMS/AIS_ALL/MarineCadastre' #where you want the CSV file to write to
dsn = 'Zone9_2014_01.gdb' #the name of your file


library(gdalUtils)
library(rgdal)
library(dplyr)
library(tidyverse)




year = strsplit(dsn, '_')
year = year[[1]]
year = year[2]

if (year > 2012){
  Vessel = paste(gsub('.gdb', '', dsn),'_Vessel', sep = "")
  Voyage = paste(gsub('.gdb', '', dsn), '_Voyage', sep = "")
  Broadcast = paste(gsub('.gdb', '', dsn), '_Broadcast', sep = "")
  
  options(warn=-1)
  vsl = readOGR(dsn = dsn, layer = Vessel, dropNULLGeometries = FALSE)
  vyg = readOGR(dsn = dsn, layer = Voyage, dropNULLGeometries = FALSE)
  options(warn=0)
  broad = readOGR(dsn = dsn, layer = Broadcast, dropNULLGeometries = FALSE)
} else if (year < 2013) {
  options(warn=-1)
  vsl = readOGR(dsn = dsn, layer = 'Vessel', dropNULLGeometries = FALSE)
  vyg = readOGR(dsn = dsn, layer = 'Voyage', dropNULLGeometries = FALSE)
  options(warn=0)
  broad = readOGR(dsn = dsn, layer = 'Broadcast', dropNULLGeometries = FALSE)
}


broadALL = cbind(data.frame(broad@data), data.frame(broad@coords))

broad2vsl = left_join(broadALL, vsl, by = "MMSI")
broad2vyg = left_join(broad2vsl, vyg, by = "VoyageID")
dataALL = select(broad2vyg, MMSI.x, BaseDateTime, coords.x2, 
                 coords.x1, SOG, COG, Heading, Name, IMO, CallSign, 
                 VesselType, Status, Length, Width, Draught, Cargo)

dataALL = rename(dataALL, 'MMSI' = 'MMSI.x', 'LAT' = 'coords.x2', 'LON' = 'coords.x1', 'VesselName' = 'Name', 'Draft' = 'Draught')


write.csv(dataALL, paste(outDir, gsub('.gdb', '.csv', dsn), sep = "/"), row.names = FALSE)


# depending on how big your .gdb file, it can take ~1 - 20 minute to run...