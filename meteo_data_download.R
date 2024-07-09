# R script

# download meteorological data (telemetric data) from repository of IMWM-NRI (pl: IMGW)
# website: https://danepubliczne.imgw.pl/datastore
# author: Janusz Godziek


library(climate)
library(sf)
library(rgdal)
library(downloader)
library(stringr)


#install.packages('downloader')

#repair R to read polish characters
Sys.setlocale("LC_CTYPE", "english")


#climate package docs
# https://www.rdocumentation.org/packages/climate/versions/1.0.4

?meteo_imgw()

pth = 'C:/Users/PC COMPUTER/Desktop/JG/SCIENCE/BIOFORCLIM_DR/Meteo_wind_regime/data/imgw/precip/'

#download data from imgw repository (climate station)
imgw_clim <- meteo_imgw(interval = 'hourly', rank = 'climate', year = 1950:2022,
                        station =  'ZAWOJA', col_names = 'short')

write.csv(imgw_clim, paste0(pth, 'Zawoja_clim.csv')) #save to csv

#download precipitation data
imgw_precip <- meteo_imgw(interval = 'daily', rank = 'precip', year = 1990:2022, #error podczas pobierania danych z 1985 r. (2022-07-25)
                          station =  c('ŚMIETANOWA', 'STAŃCOWA', 
                                       'LIPNICA WIELKA', 'ZAWOJA'), col_names = 'short')

write.csv(imgw_precip, paste0(pth, 'BG_precip.csv')) #save to csv

?meteo_imgw

#################
#built-in data about meteo stations
imgw_stations <- imgw_meteo_stations

#read csv with station names
imgw_st_names <- read.csv('D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/metadane_stacje/wykaz_stacji.csv')

#convert to spatial object
imgw_stations <- st_as_sf(imgw_stations, coords = c('X', 'Y'))
#imgw_stations <- merge(imgw_stations, imgw_st_names, by.x = 'id', by.y = 'id') #merge data
imgw_stations <- st_set_crs(imgw_stations, 'EPSG:4326')

st_write(imgw_stations,
         'D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/metadane_stacje/metadane_stacje.shp',
         driver = 'ESRI Shapefile', append = F)


#download data from IMGW telemetric stations

#create months sequence
months <- seq(1:12) #create sequence 1-12
months <- str_pad(months, 2, pad = "0") #add 0s for numbers 1-9

#create years sequence
years <- seq(from = 2008, to = 2009)
length(years)

#generate list of months (yyyy-mm)
ym <- list() #create empty list

for (y in years){
  for (mth in months){
    m_temp <- paste0(y, '-', mth) #paste year and month
    ym <- append(ym, list(m_temp)) #append to list (for each 'y' turn of for loop)
  }
}

ym <- unlist(ym) #convert list to vector
ym

#download data from IMGW repository by year
for (i in ym){
  m <- i #month (yyyy-mm)
  m2 <- str_replace(m, '-', '_') #change - to _
  y <- substr(i, 1, 4)  #extract year from ym(i)
  zip_url <- paste0('https://danepubliczne.imgw.pl/datastore/getfiledown/Arch/Telemetria/Meteo/', 
                    y,'/Meteo_', m, '.ZIP')
  
  #download file to temporary memory
  temp <- tempfile() #create temp file
  downloader::download(url = zip_url, destfile = temp, mode = 'wb') #download zip from url to temp directory
  
  setwd('D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/dane_telem') #set working directory
  #dir.create(m2) #create folder
  unzip(temp, exdir = paste0('./', m2)) #extract data to folder
  unlink(temp) #remove the temp file
}

#functions to download data -> different versions (for .ZIP and .zip)
#ZIP
imgw_telem <- function(yr_begin, yr_end){
  
  #create months sequence
  months <- seq(1:12) #create sequence 1-12
  months <- str_pad(months, 2, pad = "0") #add 0s for numbers 1-9
  
  #create years sequence
  years <- seq(from = yr_begin, to = yr_end)
  
  #generate list of months (yyyy-mm)
  ym <- list() #create empty list
  
  for (y in years){
    for (mth in months){
      m_temp <- paste0(y, '-', mth) #paste year and month
      ym <- append(ym, list(m_temp)) #append to list (for each 'y' turn of for loop)
    }
  }
  
  ym <- unlist(ym) #convert list to vector
  
  for (i in ym){
    m <- i #month (yyyy-mm)
    m2 <- str_replace(m, '-', '_') #change - to _
    y <- substr(i, 1, 4)  #extract year from ym(i)
    zip_url <- paste0('https://danepubliczne.imgw.pl/datastore/getfiledown/Arch/Telemetria/Meteo/', 
                      y,'/Meteo_', m, '.ZIP')
    
    #download file to temporary memory
    temp <- tempfile() #create temp file
    downloader::download(url = zip_url, destfile = temp, mode = 'wb') #download zip from url to temp directory
    
    setwd('D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/dane_telem') #set working directory
    #dir.create(m2) #create folder
    unzip(temp, exdir = paste0('./', m2)) #extract data to folder
    unlink(temp) #remove the temp file
  }
  
}
#zip
imgw_telem2 <- function(yr_begin, yr_end){
  
  #create months sequence
  months <- seq(1:12) #create sequence 1-12
  months <- str_pad(months, 2, pad = "0") #add 0s for numbers 1-9
  
  #create years sequence
  years <- seq(from = yr_begin, to = yr_end)
  
  #generate list of months (yyyy-mm)
  ym <- list() #create empty list
  
  for (y in years){
    for (mth in months){
      m_temp <- paste0(y, '-', mth) #paste year and month
      ym <- append(ym, list(m_temp)) #append to list (for each 'y' turn of for loop)
    }
  }
  
  ym <- unlist(ym) #convert list to vector
  
  for (i in ym){
    m <- i #month (yyyy-mm)
    m2 <- str_replace(m, '-', '_') #change - to _
    y <- substr(i, 1, 4)  #extract year from ym(i)
    zip_url <- paste0('https://danepubliczne.imgw.pl/datastore/getfiledown/Arch/Telemetria/Meteo/', 
                      y,'/Meteo_', m, '.zip')
    
    #download file to temporary memory
    temp <- tempfile() #create temp file
    downloader::download(url = zip_url, destfile = temp, mode = 'wb') #download zip from url to temp directory
    
    setwd('D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/dane_telem') #set working directory
    #dir.create(m2) #create folder
    unzip(temp, exdir = paste0('./', m2)) #extract data to folder
    unlink(temp) #remove the temp file
  }
  
}

imgw_telem(2010, 2021) #.ZIP
imgw_telem2(2010, 2021) #.zip



#download data for single month
m <- '2021-03' #month
m2 <- str_replace(m, '-', '_') #change - to _
y <- substr(m, 1, 4)  #extract year from ym(i)
zip_url <- paste0('https://danepubliczne.imgw.pl/datastore/getfiledown/Arch/Telemetria/Meteo/', 
                  y, '/Meteo_', m, '.zip') #this link works
zip_url


setwd('D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/dane_telem') #set working directory
downloader::download(url = zip_url, destfile = paste0('./zip/', m2, '.zip'), mode = 'wb') #download zip from url to temp directory
unzip(zipfile = paste0('./zip/', m2, '.zip'), 
      exdir = './zip') #extract data to folder

?unzip

#download file to temporary memory
temp <- tempfile() #create temp file
downloader::download(url = zip_url, destfile = temp) #download zip from url to temp directory

setwd('D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/dane_telem') #set working directory
dir.create(m2) #create folder
unzip(temp, exdir = paste0('./', m2)) #extract data to folder
unlink(temp) #remove the temp file


#read file with telemetric data
telem <- read.csv('D:/GIS/DANE_PRZESTRZENNE/Dane_IMGW/dane_telem/2022_01/B00202A_2022_01.csv', sep = ';')


#unzip examples
unzip('Meteo_2022-04.zip', exdir = './2022_04') #extract data to folder
unzip('Meteo_2022-04.zip', list = T) #list files located in zip






