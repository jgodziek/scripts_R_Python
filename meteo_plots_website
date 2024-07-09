# R SCRIPT

# script used for plotting meteorological data from Automatic Weather Stations (monthly plots)
# website: https://klimatbgpn.us.edu.pl/index.php/aktualne-pomiary-asm/
# author: Janusz Godziek



library(gstat)
library(ggplot2)
library(lubridate)
library(xts)
library(writexl)
library(openair)
library(plotly)
library(scales)
library(clifro)
library(ggpubr)
library(ggplotify)
library(cowplot)
library(modeest)
library(gt)
library(gridExtra)
library(magick)
library(patchwork)
library(padr)
library(zoo)
library(magrittr)
library(dplyr)
library(stringr)

#install.packages("padr")

#set language to polish
Sys.setlocale("LC_CTYPE", "polish")


#data for one month - ANALYSIS

wykres_mies <- function(plik_dane, #ścieżka dostępu do pliku z danymi
                        stacja_pl, #nazwa stacji po polsku
                        stacja, #nazwa stacji bez polskich znaków
                        wys, #wysokość stacji (razem ze skrótem m n.p.m.)
                        okres, #rok miesiąc w formacie YYYY-mm
                        pth) #ścieżka dostępu do zapisu pliku png
{
  
  #read data, set column names -> Dropbox
  dane <- read.table(plik_dane, sep = ",", skip = 4, 
                     col.names = c("TIMESTAMP","RECORD","BAT","TCP","TC","RH","WD","WS"))
  
  #create date column from timestamp
  dane$DATE <- as.Date(dane$TIMESTAMP, "%Y-%m-%d")
  #move 'DATE' column to the first place in the table
  dane <- dane[,c(ncol(dane),1:(ncol(dane)-1))]
  
  #extract month
  #dane$MONTH <- as.yearmon(dane$TIMESTAMP) #format 'maj 2022'
  dane$MONTH <- format(as.Date(dane$DATE, format="%d/%m/%Y"),"%Y-%m") #format '2022-05'
  #set 'MONTH' column to the first place
  dane <- dane[,c(ncol(dane),1:(ncol(dane)-1))]
  
  #convert TIMESTAMP column from character to POSIXct class
  dane$TIMESTAMP <- as.POSIXct(dane$TIMESTAMP,tz="UTC")
  
  #converts data to numeric (multiple columns)
  cols.num <- c("BAT","TCP","TC","RH","WD","WS")
  dane[cols.num] <- sapply(dane[cols.num],as.numeric)
  sapply(dane, class)
  
  dane_1month <- dane[dane$MONTH == okres, ] #select data for one month
  
  stacja_info <- data.frame(a = c("stacja", "wysokość", "miesiąc"),
                            b = c(stacja_pl, wys, okres))
  
  #remove rows with NA
  #! removes also columns with complete data! - need to be worked out on other way
  dane_1month_noNA <- na.omit(dane_1month)
  
  #podstawowe statystyki
  #temperatura
  Tmax = round(max(dane_1month_noNA$TC), 2)
  Tmin = round(min(dane_1month_noNA$TC), 2)
  Tmean = round(mean(dane_1month_noNA$TC), 2)
  ampT = Tmax - Tmin
  
  #wilgotność
  RHmin = round(min(dane_1month_noNA$RH), 2)
  RHmean = round(mean(dane_1month_noNA$RH), 2)
  
  #prędkość wiatru
  WSmax = round(max(dane_1month_noNA$WS), 2)
  WSmean = round(mean(dane_1month_noNA$WS), 2)
  
  #extract datetime of max temperature 
  x <- dane_1month_noNA %>%
    slice(which.max(TC))
  Tmax_time <- format(x[1, "TIMESTAMP"], "%Y-%m-%d %H:%M")
  
  #extract datetime of min temperature 
  x <- dane_1month_noNA %>%
    slice(which.min(TC))
  Tmin_time <- format(x[1, "TIMESTAMP"], "%Y-%m-%d %H:%M")
  
  #extract datetime of min relative humidity 
  x <- dane_1month_noNA %>%
    slice(which.min(RH))
  RHmin_time <- format(x[1, "TIMESTAMP"], "%Y-%m-%d %H:%M")
  
  #extract datetime of max wind speed
  x <- dane_1month_noNA %>%
    slice(which.max(WS))
  WSmax_time <- format(x[1, "TIMESTAMP"], "%Y-%m-%d %H:%M")
  
  
  #bulid dataframe with characteristic values
  dane_char <- data.frame(a=c("maksymalna prędkość wiatru", "średnia prędkość wiatru", 
                              "temperatura maksymalna", "temperatura minimalna", 
                              "średnia temperatura", "amplituda temperatur", 
                              "minimalna wilgotność względna", "średnia wilgotność względna"),
                          wartosc = c(paste(WSmax, "m/s"), paste(WSmean, "m/s"), 
                                      paste(Tmax, "°C"), paste(Tmin, "°C"), 
                                      paste(Tmean, "°C"), paste(ampT, "°C"),
                                      paste(RHmin, "%"), paste(RHmean, "%")),
                          czas = c(WSmax_time, "-", Tmax_time, Tmin_time, "-", "-", 
                                   RHmin_time, "-"))
  
  
  # construct table with ggtexttable from ggpubr package         
  stacja_info2 <- ggtexttable(stacja_info, rows = NULL, cols = c("", ""),
                              theme = ttheme("mRedWhite", base_size = 20))
  
  stacja_info2 <- stacja_info2 %>%
    tab_add_border(linetype = 1, linewidth = 1, linecolor = "gray35")
  stacja_info2 <- stacja_info2 + theme(
    plot.margin = unit(c(0.2,5,0.5,0.2),"cm")) +
    annotate("text", x = 0.5, y = 0.5, label= "UŚ & BgPN", 
             angle = 20, size = 12, fontface = "bold", alpha = 0.2)
  
  
  # construct table with ggtexttable from ggpubr package         
  dane_char2 <- ggtexttable(dane_char, rows = NULL, cols = c("", "wartość", "data i godzina (UTC)"),
                            theme = ttheme("mBlue"))
  
  #add border
  dane_char2 <- dane_char2 %>%
    tab_add_border(linetype = 1, linewidth = 1, linecolor = "gray35")
  
  dane_char2 <- dane_char2 + theme(plot.margin = unit(c(0.2,4,0.5,0.2),"cm")) +
    annotate("text", x = 0.5, y = 0.5, label= "UŚ & BgPN", 
             angle = 20, size = 20, fontface = "bold", alpha = 0.2)
  
  #plot data v1b - temperature [geom_line] -> OK
  Ptemp_1month <- ggplot(dane_1month, aes(x = TIMESTAMP, y = TC)) +
    geom_line(colour = "mediumblue", size = 0.3) +
    theme_gray() + 
    scale_x_datetime(date_breaks = "1 day", date_labels = "%d", minor_breaks = "1 day") + 
    scale_y_continuous(breaks = seq(-40, 50, by=2)) +
    theme(
      axis.text.x=element_text(size=10, face = "bold", vjust = 1),
      axis.text.y=element_text(size=10, face = "bold"),
      axis.title.x = element_blank(), #remove xaxis title
      axis.title.y = element_text(size=10, face = "bold", vjust = 3),
      panel.border = element_rect(colour = "gray35", fill=NA, size=1),
      plot.margin = unit(c(0.2,0.2,0.5,1),"cm")
    ) + ylab("temperatura [°C]") + geom_hline(yintercept = 0, size = 0.5, col = "gray35") +
    annotate("text", x= as.POSIXct(dane_1month$TIMESTAMP[nrow(dane_1month)/2]), 
             y=(Tmax-Tmin)/2, label= "UŚ & BgPN", angle = 20, size = 25, 
             fontface = "bold", alpha = 0.2)
  
  
  #plot data v2 - relative humidity [geom_col]
  Prhum_1month <- ggplot(dane_1month, aes(x = TIMESTAMP, y = RH)) +
    geom_col(fill = 'mediumblue', width = 600, alpha = 0.3) +
    theme_gray() +
    scale_x_datetime(date_breaks = "1 day", date_labels = "%d", minor_breaks = "1 day") + 
    scale_y_continuous(breaks = seq(0, 100, by=10)) +
    labs(caption = 'wyk. J. Godziek') +
    theme(
      axis.text.x=element_text(size=10, face = "bold", vjust = 1),
      axis.text.y=element_text(size=10, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=10, face = "bold", vjust = 3),
      panel.border = element_rect(colour = "gray35", fill=NA, size=1),
      plot.margin = unit(c(0.2,0.2,0.5,1),"cm")
    ) + ylab("wilgotność względna [%]") +
    annotate("text", x= as.POSIXct(dane_1month$TIMESTAMP[nrow(dane_1month)/2]), 
             y=50, label= "UŚ & BgPN", angle = 10, size = 25, 
             fontface = "bold", alpha = 0.2) +
    labs(caption = 'wyk. J. Godziek') +
    theme(plot.caption = element_text(size=12, face="italic", vjust = 0.5))
  
  
  
  #plot data v3 - wind speed [geom_line]
  Pws_1month <- ggplot(dane_1month, aes(x = TIMESTAMP, y = WS)) +
    geom_line(colour = "darkred", size = 0.1) +
    theme_gray() +
    scale_x_datetime(date_breaks = "1 day", date_labels = "%d", minor_breaks = "1 day") + 
    scale_y_continuous(breaks = seq(0, 100, by=0.5)) +
    theme(
      axis.text.x=element_text(size=10, face = "bold", vjust = 1),
      axis.text.y=element_text(size=10, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=10, face = "bold", vjust = 3),
      panel.border = element_rect(colour = "gray35", fill=NA, size=1),
      plot.margin = unit(c(0.2,0.2,0.5,1),"cm")
    ) + ylab("prędkość wiatru [m/s]") +
    annotate("text", x= as.POSIXct(dane_1month$TIMESTAMP[nrow(dane_1month)/2]), 
             y=WSmax/2, label= "UŚ & BgPN", angle = 20, size = 25, 
             fontface = "bold", alpha = 0.2)
  
  
  #windrose (library clifro) -> OK
  Pwd_1month_Clifro <- windrose(speed = dane_1month$WS, direction = dane_1month$WD,
                                n_directions = 36,
                                ggtheme = "gray", col_pal = "Set1", speed_cuts = c(100, 101, 102)) #enables producing rose for one wind speed interval
  Pwd_1month_Clifro <-Pwd_1month_Clifro + theme(legend.position = "none", #hide legend of wind speed
                                                panel.border = element_rect(colour = "gray35", fill=NA, size=1), #add frame
                                                axis.text.x=element_text(size=12, face = "bold", vjust = 1),
                                                axis.text.y=element_text(size=10, face = "bold"),
                                                plot.margin = unit(c(0.2,0.2,0.5,0.2),"cm"))
  #add title
  Pwd_1month_Clifro <- Pwd_1month_Clifro + ggtitle('kierunki wiatru') + theme(
    plot.title = element_text(size=8, face = "bold", vjust = 1, hjust = 0.5)) +
    annotate("text", x = 0, 
             y = 0, label= "UŚ & BgPN", angle = 45, size = 12, 
             fontface = "bold", alpha = 0.2)
  
  
  #library patchwork -> merge all plots and tables
  plot_meteoData <- (stacja_info2 + dane_char2 + Pwd_1month_Clifro) / 
    Pws_1month / Ptemp_1month / Prhum_1month
  
  #format strings (for file-saving purposes)
  okres2 <- gsub('-', '_', okres)
  nazwaPliku <- c(paste0(stacja, '_', okres2))
  
  #export plot
  #170 dpi
  ggsave(plot_meteoData, width = 2000, height = 2000, dpi = 170, units = "px",
         device = "png", limitsize = F, filename = paste(nazwaPliku, ".png", sep = ""),
         path = pth)
  
}


#function generating vector of months in the 'yyyy-mm' format (for full years)
ym_seq <- function(yr_begin, #year to start (format: yyyy)
                   yr_end, ##year to finish (format: yyyy)
                   mth_begin, #month to start (format: month number, eg. 2 [February])
                   mth_end) #month to finish (format: month number, eg. 9 [September])
{
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
  ym <- ym[!(ym %in% head(ym, (mth_begin-1)))] #remove first elements
  ym <- head(ym, -(12-mth_end)) #remove last elements
  return(ym) #return value of ym as an output of function
}

#month data for iteration
mth <- ym_seq(2023, 2024, 8, 2) #sometimes works, sometimes not -> check the result!!!
mth

#mth <- c('2023-03', '2023-04', '2023-05', '2023-06', '2023-07')

#path to data file
data_pth <- 'D:/JG/SCIENCE/BIOFORCLIM_DR/stacje_meteo/BabiaG/_DANE/2024_03_04'

#path to save output plots
save_pth <- 'D:/JG/SCIENCE/BIOFORCLIM_DR/stacje_meteo/BabiaG/materialy_stronawww/wykresy/2024_03_20'

#data frame to iterate through stations
stations <- data.frame(name = c('Stonow', 'Plaj', 'Perc'),
                       name_pl = c('Stonów', 'Płaj', 'Perć'),
                       alt = c('807 m n.p.m.', '1222 m n.p.m.', '1415 m n.p.m.'))


#automatic generation of monthly plots
for (s in 1:nrow(stations)){ #iterate through stations
  for (m in mth){ #for each station iterate through months
    wykres_mies(plik_dane = paste0(data_pth, '/', stations[s,1], '_10min.dat'),
                stacja_pl = stations[s, 2],
                stacja = stations[s, 1],
                wys = stations[s, 3],
                okres = m,
                pth = save_pth)
  }
}



#repair the data

setwd(data_pth)
getwd()


data <- read.table('Plaj_10min.dat', sep = ",", skip = 4, 
                   col.names = c("TIMESTAMP","RECORD","BAT","TCP","TC","RH","WD","WS"))

library(data.table)
sel_data <- data[data$TIMESTAMP %like% "2024-02", ] #select records from February
summary(sel_data)


sel_data$TC <- replace(sel_data$TC, sel_data$TC == -100, NA)
sel_data$RH <- replace(sel_data$RH, sel_data$RH == -100, NA)

#save
write.csv(sel_data, 'Plaj_2024_02_repaired.dat')



#Stonow: 807, Płaj: 1222, Perc: 1415
#function to plot monthly data
wykres_mies(plik_dane = paste0(data_pth, "/Plaj_10min_repaired.dat"),
            stacja_pl = 'Płaj',
            stacja = 'Plaj',
            wys = '1222 m n.p.m.',
            okres = '2024-02',
            pth = save_pth)






