# HOT_RAIN performs Rhythmicity Analysis Incorporating Nonparameteric methods to determine minima for 
# picophytoplankton abundance and maxima for phytoplankton Qc at Station ALOHA (Thaben and Westermark, 2014).
#
# Started:  22/Dec/2021 Annette Hynes, UW
# Modified: 03/Jan/2022 Repeat rain analysis for negative data to get trough as well as peak times
# Modified:	20/Jun/2023 Streamline data access and re-run using latest data available

library(tidyverse)
library(viridis)
library(readxl)
library(rain)
#renv::activate("~/Desktop/renvtest/popcycle/")
library(popcycle)

###############
### SeaFlow ###
###############

# Download SeaFlow data from Zenodo (doi.org/10.5281/zenodo.7154076)
url <- "https://zenodo.org/record/7154076/files/SeaFlow_dataset_v1.5.xlsx" 
file_name <- tempfile()
try(download.file(url,file_name,method="curl"))
if (is.na(file.size(file_name))) download.file(url, file_name,method="auto")

# Read SeaFlow data
all_SF <- readxl::read_xlsx(file_name) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")) %>%
  mutate(diam_croco = as.numeric(diam_croco),
         Qc_croco = as.numeric(Qc_croco))

# Select data around Station aloha location
lat <- 22.75; lon <- -158 
aloha <- all_SF %>% dplyr::filter(lat >= 22.25 & lat <= 23.25 & lon <= -157.5 & lon >= -158.5) %>%
  dplyr::filter(time < as.POSIXct("2021-02-01", tz = "UTC"))
rm(all_SF)  # clear some memory

# Convert from short to long format for easier plotting
aloha_long <- aloha %>% 
	tidyr::pivot_longer(cols = -c(time, lat, lon, depth, cruise), 
		names_to = c('.value','pop'), names_sep = "_") %>%
	arrange(time)

# Specific language:  use euk instead of picoeuk
group.colors <- c(prochloro = viridis::viridis(4)[1],
                  synecho = viridis::viridis(4)[2], 
                  euk = viridis::viridis(4)[3],
                  croco = viridis::viridis(4)[4])

aloha_long <- aloha_long %>% 
	mutate(pop = case_when(pop == 'picoeuk' ~ 'euk', TRUE ~ pop),
		pop = factor(pop, levels = names(group.colors)), Program = "SeaFlow")

# Date conversion
aloha_long <- aloha_long %>% 
	mutate(local_time = lubridate::with_tz(time, tzone = 'HST'),
		year = lubridate::year(local_time),
        month = lubridate::month(local_time),
        yday = lubridate::yday(local_time),
        hour = lubridate::hour(local_time))

# Split double cruises by bumping late month cruises to the next month
aloha_long <- aloha_long %>% 
	mutate(month = case_when(cruise == "KM1709" ~ month + 1,
		cruise == "KM2011" ~ month + 1, TRUE ~ month))

#################################
### NONPARAMETRIC RHYTHMICITY ###
#################################

# Use only cruises that have at least 48 hours worth of data in order to resolve a 24-hr period.

ALOHA_48 <- aloha_long %>%
    dplyr::group_by(cruise) %>%
    dplyr::mutate(duration = as.numeric(difftime(dplyr::last(local_time), dplyr::first(local_time), units = 'hours'))) %>%
    dplyr::filter(duration >= 48 & !cruise %in% c('KOK1606', 'MGL1704', 'KM1906')) # Exclude there-and-back-again cruises that seem long

hour_list <- format(lubridate::floor_date(ALOHA_48$local_time, unit="hour"), "%Y-%m-%dT%H")
time <- as.POSIXct(hour_list, format = "%Y-%m-%dT%H", tz = "HST")
ALOHA_48$local_hour <- lubridate::hour(ALOHA_48$local_time)
ALOHA_48$local_datehour <- time

ALOHA_48_hr <- ALOHA_48 %>%
    group_by(local_datehour, pop, cruise) %>%
    summarize_all(funs(mean), na.rm = TRUE) 
   
ALOHA_48_hr <- ALOHA_48_hr %>%
  dplyr::group_by(cruise) %>%
  dplyr::mutate(consec_day = as.integer(date - first(date)))
  
cruise_list <- unique(ALOHA_48_hr$cruise)
param_list <- c('abundance', 'Qc')
phyto_list <- unique(ALOHA_48_hr$pop)

# For peaks, run normally.  For troughs, multiply de-trended data by -1.
rain_file <- paste0(save_path, '../HOT_rain_results_', date(), '.csv')
big_rain <- NULL

#big_rain <- read_csv(rain_file)
for (phyto in phyto_list){
    print(phyto)
    for(voyage in cruise_list2){
        print(voyage)
        df <- ALOHA_48_hr[which(ALOHA_48_hr$pop == phyto & ALOHA_48_hr$cruise == voyage), c('local_datehour', param_list, "local_hour", "consec_day")]      
        df_full <- df[, c("local_datehour", param)] %>%
            ungroup() %>%
            complete(local_datehour = seq(min(local_datehour), max(local_datehour), by = "hour"))   # pad the empty hours
        x <- as.numeric((df_full$local_datehour - df_full$local_datehour[1])/3600)  #   Times by hour from first time point
        for (param in param_list){
            print(param)
            y <- dplyr::pull(df_full, param)
            this_lm <- coef(lm(y ~ x))
            data_white <- (y - (this_lm[2]*x + this_lm[1]))#*(-1)   # Remove long-term trends by subtracting the linear regression
            if (param == "abundance"){
              data_white <- (1-)*data_white # Target minima instead of maxima for abundance
            }
            dfw <- data.frame(x, data_white)
           
            dat_rain <- rain(data_white, deltat=1, period=24, nr.series=1, peak.border=c(0.3,0.7), method = 'independent', na.rm = TRUE, verbose=TRUE)
            dat_rain$cruise <- voyage
            dat_rain$pop <- phyto
            dat_rain$param <- param
            dat_rain$first_peak <- df_full$local_datehour[1] + 3600*dat_rain$phase
            dat_rain$peak_hour <- lubridate::hour(dat_rain$first_peak)
            
            big_rain <- rbind(big_rain, dat_rain)
            write.csv(big_rain, file = rain_file, row.names = FALSE) # Can take a long time.  Save results in case R is disrupted.
        }   # end param loop
    }   # end cruise loop
}   # end pop loop
