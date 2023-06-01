# HOT_figures puts together final figures of the Station ALOHA SeaFlow data for publication. 
#
# Started: 25/May/2023  Annette Hynes, UW

library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)
library(suncalc)
library(googlesheets4)
library(viridis)
library(latex2exp)
library(readxl)
library(ggpubr)
library(lubridate)
library(marmap)
library(oce)
library(knitr)
library(broom)
library(ggh4x)
library(grDevices)
renv::activate("~/Desktop/renvtest/popcycle/")
library(popcycle)

#################
### LOAD DATA ###
#################

save_path <- '~/Documents/SeaFlow/HOT/'
file_name <- paste0(save_path, 'Data/SeaFlow_dataset_v1.5.xlsx') # From Zenodo (10.5281/zenodo.7154076)
all_SF <- readxl::read_xlsx(file_name)
lat <- 22.75; lon <- -158 # Station ALOHA location
aloha <- dplyr::filter(all_SF, lat >= 22.25 & lat <= 23.25 & lon <= -157.5 & lon >= -158.5)
rm(all_SF)  # clear some memory

aloha$time <- as.POSIXct(aloha$time, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
aloha <- dplyr::filter(aloha, time < as.POSIXct("2021-02-01", tz = "UTC"))

# Convert from short to long format for easier plotting
abundance <- aloha %>%
  tidyr::pivot_longer(
    cols = starts_with("abundance_"),
    names_to = "pop",
    names_prefix = "abundance_",
    values_to = "abundance",
    values_drop_na = FALSE
  ) %>%
  dplyr::select(time, lat, lon, cruise, pop, abundance)

aloha$diam_croco <- as.numeric(aloha$diam_croco)
diam <- aloha %>%
  tidyr::pivot_longer(
    cols = starts_with("diam_"),
    names_to = "pop",
    names_prefix = "diam_",
    values_to = "diam",
    values_drop_na = FALSE
  ) %>%
  dplyr::select(time, lat, lon, cruise, pop, diam)

aloha$Qc_croco <- as.numeric(aloha$Qc_croco)
Qc <- aloha %>%
  tidyr::pivot_longer(
    cols = starts_with("Qc_"),
    names_to = "pop",
    names_prefix = "Qc_",
    values_to = "Qc",
    values_drop_na = FALSE
  ) %>%
  dplyr::select(time, lat, lon, cruise, pop, Qc)

biomass <- aloha %>%
  tidyr::pivot_longer(
    cols = starts_with("biomass_"),
    names_to = "pop",
    names_prefix = "biomass_",
    values_to = "biomass",
    values_drop_na = FALSE
  ) %>%
  dplyr::select(time, lat, lon, cruise, pop, biomass)

df_list <- list(abundance, diam, Qc, biomass)
aloha_long <- df_list %>% purrr::reduce(full_join, by = c("time", "lat", "lon", "cruise", "pop"))

# Specific language:  use euk instead of picoeuk
aloha_long$pop <- as.factor(gsub('picoeuk', 'euk', aloha_long$pop))
group.colors <- c(prochloro = viridis::viridis(4)[1],
                  synecho = viridis::viridis(4)[2], euk = viridis::viridis(4)[3],
                  croco = viridis::viridis(4)[4])
aloha_long$pop <- factor(aloha_long$pop, levels = names(group.colors))

date_HST <- format(aloha_long$time, format = "%Y-%m-%d %H:%M:%S", tz = 'HST')
aloha_long$local_time <- as.POSIXct(date_HST, tz = "HST")

aloha_long$year <- as.factor(lubridate::year(aloha_long$local_time))
aloha_long$month <- as.factor(lubridate::month(aloha_long$local_time))
aloha_long$hour <- as.factor(lubridate::hour(aloha_long$local_time))

mo_list <- c('KM1709', 'KM2011')    # Split double cruises by bumping late month cruises to the next month
ind_mo <- which(aloha_long$cruise %in% mo_list)
aloha_long$month[ind_mo] <- as.factor(as.numeric(as.character(aloha_long$month[ind_mo])) + 1)
aloha_long$day <- as.factor(lubridate::yday(aloha_long$local_time))

#################
### ABUNDANCE ###
#################

### INFLUX ###
influx <- read.csv(paste0(save_path, "Data/HOT_influx.csv"), skip = 5, header = F)  # Data from HOT-DOGS
in_head <- read.csv(paste0(save_path, "Data/HOT_influx.csv"), skip = 2, nrows = 1)
colnames(influx) <- colnames(in_head)
  
influx[influx == -9] <- NA # Flagged values
in_date <- stringr::str_pad(influx$date, 6, side = "left", pad = 0)
in_time <- stringr::str_pad(influx$time, 6, side = "left", pad = 0)
t_influx <- as.POSIXct(paste0(in_date, ' ', in_time), format = '%m%d%y %H%M%S', tz = 'UTC')
influx$DateTime <- t_influx

# Convert abundances from 10E5 cells/mL to cells/uL
influx$peuk <- as.numeric(as.character(influx$ebact))*(100)
influx$peuk[which(influx$peuk < -1)] <- NA
influx$syn <- as.numeric(as.character(influx$sbact))*(100)
influx$syn[which(influx$syn < -1)] <- NA
influx$pro <- as.numeric(as.character(influx$pbact))*(100)
influx$pro[which(influx$pro < -1)] <- NA

post_influx <- influx[, c('DateTime', 'pro', 'syn', 'peuk')]

influx_db <- NULL
for (i in seq(1, nrow(post_influx))){
    pro <- data.frame(time = post_influx$DateTime[i], abundance = post_influx$pro[i], pop = 'prochloro')
    syn <- data.frame(time = post_influx$DateTime[i], abundance = post_influx$syn[i], pop = 'synecho')
    pico <- data.frame(time = post_influx$DateTime[i], abundance = post_influx$peuk[i], pop = 'euk')
    influx_db <- rbind(influx_db, pro, syn, pico)
}

influx_db$month <- lubridate::month(influx_db$time)
influx_db$year <- lubridate::year(influx_db$time)

# Split double cruises by bumping late month cruises to the next month for KM2011/HOT323
mo_list_date <- c(as.POSIXct("2020-09-26 13:19:36", tz = "UTC"))    
ind_mo <- which(influx_db$time %in% mo_list_date)
influx_db$month[ind_mo] <- influx_db$month[ind_mo] + 1

influx_db$Amin <- NA # Blank variable for error bars
influx_db$Amax <- NA
influx_db$pop <- factor(influx_db$pop, levels = c('prochloro', 'synecho', 'euk', 'croco'))
influx_db$Program <- "HOT"

### SeaFlow daily means ###
ALOHA_in <- aloha_long[, c("time", "lat", "lon", "cruise", "pop", 'abundance', 'diam', 'Qc', 'biomass', 'local_time', 'year', 'month', 'day')]
all_ALOHA_day <- ALOHA_in %>%
        group_by(day, month, year, pop) %>%
        select(!cruise) %>%
        summarise(across(everything(), .f = list(mean = mean, sd = sd), na.rm = TRUE))

all_ALOHA_day$Amin <- all_ALOHA_day$abundance_mean - all_ALOHA_day$abundance_sd # values for error bars
all_ALOHA_day$Amax <- all_ALOHA_day$abundance_mean + all_ALOHA_day$abundance_sd # values for error bars
all_ALOHA_day$year <- as.numeric(as.character(all_ALOHA_day$year))
all_ALOHA_day$month <- as.numeric(as.character(all_ALOHA_day$month))

SF_day <- as.data.frame(all_ALOHA_day[, c('abundance_mean', 'Amin', 'Amax', 'pop', 'month', 'year')])
SF_day$Program <- 'SeaFlow'
SF_day$abundance <- SF_day$abundance_mean

### Plot abundances together together ###

col_seaflow <- viridis::viridis(9)[5]
col_HOT <- '#FFFFFFFF'

prog.colors <- c(SeaFlow =  viridis::viridis(9)[5], HOT = '#FFFFFFFF')

SF_day <- SF_day[, c('abundance', 'Amin', 'Amax', 'pop', 'month', 'year', 'Program')]
all_HOT <- influx_db[, c('abundance', 'Amin', 'Amax', 'pop', 'month', 'year', 'Program')]
all_day <- rbind(SF_day, all_HOT)
all_day$Program <- as.factor(all_day$Program)

# abundance: averaged daily, faceted by month
# Use only the HOT points that have SeaFlow data

all_day$mo_yr <- paste0(month.abb[all_day$month], ' ', all_day$year)
SF_moyr <- unique(all_day$mo_yr[which(as.character(all_day$Program) == 'SeaFlow')])   # List of month-year pairs with SeaFlow data
inclusive <- all_day[which(all_day$mo_yr %in% SF_moyr), ]         # Allow only HOT data with SeaFlow counterpart

fig_name <- paste0(save_path, "HOT_abundance_daily_point_on_station_sd_SF.pdf")
pdf(fig_name, width = 15, height = 8)
g <- ggplot2::ggplot(inclusive, aes(x = year, y = abundance, fill = Program)) +
    ggplot2::geom_linerange(ggplot2::aes(x = year, ymin = Amin, ymax = Amax), color = col_seaflow) +
    ggplot2::geom_point(pch = 21, size = 3) +
    ggplot2::facet_grid(rows = vars(pop), cols = vars(month), scales = 'free_y') +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::scale_fill_manual(values = prog.colors) +
    ggplot2::labs(y = unname(latex2exp::TeX('Abundance (10$^6$ cells L$^{-1}$)')), x = 'Year')
print(g)
dev.off()

### Abundance property-property plot ###
# Hourly means for abundance

ALOHA_in <- aloha_long[, c("time", "year", "month", "day", "hour", "cruise", "pop", 'abundance')]

abund_hr <- ALOHA_in %>%
  dplyr::group_by(year, day, hour, pop, cruise) %>%
  dplyr::summarize(abund_SF = mean(abundance, na.rm = TRUE))

influx_db$day <- lubridate::yday(influx_db$time)
influx_db$hour <- lubridate::hour(influx_db$time)
#abund_hot$hour <- 12  # Calibration is based on 12:00 GMT (02:00 HST)

AvA <- merge(influx_db, abund_hr)
AvA$pop <- factor(AvA$pop, levels = c('prochloro', 'synecho', 'euk')) #get the order I want
droplevels(AvA) # drop croco

influx.colors <- c(prochloro = viridis::viridis(4)[1],
                  synecho = viridis::viridis(4)[2], euk = viridis::viridis(4)[3])   # No croco

# plot
p <- list()
i <- 1
for (phyto in c("prochloro", "synecho", "euk")){
  this_df <- dplyr::filter(AvA, pop == phyto)
  lim <- max(max(max(this_df$abundance, na.rm = TRUE), max(this_df$abund_SF)), 6)   # FR wants to spread the Syn axes a little
  print(lim)
  p[[i]] <- ggplot2::ggplot(this_df, aes(x = abundance, y = abund_SF, fill = pop)) +
    ggplot2::geom_abline(slope = 1, intercept = 0) + 
    ggplot2::geom_point(pch = 21, size = 3, alpha = 0.5) +
    scale_fill_manual(values = influx.colors) +
    ggplot2::coord_equal(xlim = c(0, lim), ylim = c(0, lim)) + 
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::labs(x = '', y = '') 
  i <- i + 1
}

# Get single legend 
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs,function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

my_leg <- get_legend(p[[1]])
x <- grid::textGrob(unname(latex2exp::TeX('Influx Abundance (10$^6$ cells L$^{-1}$)')), gp = gpar(fontsize=18))
y <- grid::textGrob(unname(latex2exp::TeX('SeaFlow Abundance (10$^6$ cells L$^{-1}$)')), rot = 90, gp = gpar(fontsize=18))

fig_name <- paste0(save_path, "HOT_abundance_SF_vs_Influx.pdf")
pdf(fig_name, width = 15, height = 5)
gridExtra::grid.arrange(p[[1]] + theme(legend.position="none"), p[[2]] + theme(legend.position="none"), 
  p[[3]] + theme(legend.position="none"), my_leg, nrow = 1, widths=c(2, 2, 2, 1), 
  left = y, bottom = x)

dev.off()

### Mean + 2 sd to define bloom ###

bloom <- inclusive %>%
    group_by(Program, pop) %>%
    mutate(abund_mean = mean(abundance, na.rm = TRUE), abund_sd = sd(abundance, na.rm = TRUE))

bloom_sf <- subset(bloom, Program == "SeaFlow")
bloom_sf$plus_2sd <- bloom_sf$abund_mean + 2*bloom_sf$abund_sd

g <- ggplot2::ggplot(bloom, aes(x = year, y = abundance, fill = Program)) +
    ggplot2::geom_linerange(ggplot2::aes(x = year, ymin = Amin, ymax = Amax), color = col_seaflow) +
    ggplot2::geom_point(pch = 21, size = 3) +
    ggplot2::geom_hline(data = bloom_sf, aes(yintercept = abund_mean, linetype = "Mean"), color = 'black') +
    ggplot2::geom_hline(data = bloom_sf, aes(yintercept = plus_2sd, linetype = '2.0 sd'), color = 'firebrick3') +
    scale_linetype_manual(name = "Mean + 2 sd", values = c(1, 2),
        guide = guide_legend(override.aes = list(color = c("firebrick3", "black")))) +
    ggplot2::facet_grid(rows = vars(pop), cols = vars(month), scales = 'free_y') +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::scale_fill_manual(values = prog.colors) +
    ggplot2::labs(y = unname(latex2exp::TeX('Abundance (10$^6$ cells L$^{-1}$)')), x = 'Year', title = 'Define blooms by abundance standard deviation')

fig_name <- paste0(save_path, "HOT_abundance_bloom_sd.pdf")
pdf(fig_name, width = 15, height = 8)
    print(g)
dev.off()

####################
### HOURLY MEANS ###
####################

ALOHA_in <- aloha_long[, c("time", "lat", "lon", "cruise", "pop", 'abundance', 'diam', 'Qc', 'biomass', 'local_time')]

# Mid-month sunrise and sunset
midmonth <- as.Date(paste0(lubridate::year(ALOHA_in$local_time), '-', lubridate::month(ALOHA_in$local_time), '-15'), tz = 'HST')
sun_mid <- suncalc::getSunlightTimes(date = midmonth, lat = lat, lon = lon, keep = c('sunrise', 'sunset'), tz = "HST")
sunrise_mid <- lubridate::hour(sun_mid$sunrise) + lubridate::minute(sun_mid$sunrise)/60

# START AT DAWN -- Do this before getting hourly means.  Subtracting sunrise after led to weird hour rounding resulting in some hours missing and some duplicated.

sunrise_seconds <- sunrise_mid*3600
suntime <- ALOHA_in$local_time - sunrise_seconds
ALOHA_in$sundate <- lubridate::as_date(suntime)
ALOHA_in$sunhour <- lubridate::hour(suntime)
ALOHA_in$hour <- lubridate::hour(ALOHA_in$time)

all_ALOHA_hr <- ALOHA_in %>%
        dplyr::group_by(sunhour, sundate, pop, cruise) %>%
        dplyr::summarize_all(mean, na.rm = TRUE)

all_ALOHA_hr$date <- lubridate::as_date(all_ALOHA_hr$local_time)

all_ALOHA_hr$year <- as.numeric(as.character(lubridate::year(all_ALOHA_hr$local_time)))
all_ALOHA_hr$month <- as.numeric(as.character(lubridate::month(all_ALOHA_hr$local_time)))
ind_mo <- which(all_ALOHA_hr$cruise %in% mo_list) # Split multiple cruises within a month
all_ALOHA_hr$month[ind_mo] <- as.numeric(as.character(all_ALOHA_hr$month[ind_mo])) + 1
all_ALOHA_hr$day <- as.numeric(as.character(lubridate::day(all_ALOHA_hr$local_time)))

# Get mid-month sunrise, sunset again.  Don't trust the hourly means
midmonth <- as.Date(paste0(all_ALOHA_hr$year, '-', all_ALOHA_hr$month, '-15'), tz = 'HST')
sun_mid <- suncalc::getSunlightTimes(date = midmonth, lat = lat, lon = lon, keep = c('sunrise', 'sunset'), tz = "HST")
sunrise_mid <- lubridate::hour(sun_mid$sunrise) + lubridate::minute(sun_mid$sunrise)/60
sunset_mid <- lubridate::hour(sun_mid$sunset) + lubridate::minute(sun_mid$sunset)/60
all_ALOHA_hr$sunrise_mid <- sunrise_mid
all_ALOHA_hr$sundusk_mid <- sunset_mid - sunrise_mid
all_ALOHA_hr$sundate <- as.factor(all_ALOHA_hr$sundate)

all_ALOHA_hr$pop <- factor(all_ALOHA_hr$pop, levels = c('prochloro', 'synecho', 'euk', 'croco' )) #get the order I want

# Different subsets
cruise_cut <- c('KOK1606', 'MGL1704', 'KM1906')   # These cruises pass through Station ALOHA so quickly, they are not helpful for diel patterns
all_ALOHA_hr_no_outlier <- subset(all_ALOHA_hr, !(cruise %in% cruise_cut))
ind_285 <- which(all_ALOHA_hr_no_outlier$cruise == "KOK1608" & all_ALOHA_hr_no_outlier$pop == "croco")  # Croco diameter is very strange and high during this cruise
all_ALOHA_hr_no_outlier <- all_ALOHA_hr_no_outlier[-ind_285, ]
all_ALOHA_hr_no_outlier$year <- as.factor(all_ALOHA_hr_no_outlier$year)

#######################
### DIAMETER AND QC ###
#######################

# Remove low abundance points first.
all_ALOHA_hr_robust <- subset(all_ALOHA_hr_no_outlier, abundance > 0.02) # 0.048 = About 30 cells per 3-min file

# Get number of days in each population:
all_ALOHA_hr_robust$n_days <- NA
phyto_list <- c("prochloro", 'synecho', 'euk', 'croco')

for(phyto in phyto_list){
    ind <- which(all_ALOHA_hr_robust$pop == phyto)
    this_phyto <- subset(all_ALOHA_hr_robust, pop == phyto)
    nday <- length(unique(droplevels(this_phyto$sundate)))
    all_ALOHA_hr_robust$n_days[ind] <- nday
}

# Draw with dual axes

coef <- (0.261*((4*pi/3)^0.860))/(2^(3*0.860))  # relationship between Qc and diameter in Ribalet et al, 2019
expo <- 3*0.860

med_diam <- all_ALOHA_hr_robust %>%
    dplyr::group_by(pop) %>%
    summarize(median = median(diam, na.rm = TRUE))

g1 <- ggplot(all_ALOHA_hr_robust, aes_string(x = 'sunhour', y = 'diam', group = 'sundate')) +
    geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    geom_line(alpha = 0.25, color = prog.colors['SeaFlow'], size = 2) +
    theme_bw(base_size = 20) +
    scale_y_continuous(sec.axis = sec_axis(~ coef*.^expo, labels = NULL)) +
    facet_wrap(vars(pop), ncol = 1, scales = 'free_y') +
    geom_text(aes(label = paste0('n = ', all_ALOHA_hr_robust$n_days), x = -Inf, y = Inf), hjust = 0, vjust = 1.5, check_overlap = T, size = 6) +
    ggplot2::labs(y = unname(latex2exp::TeX('Diameter ($\\mu$m)')), x = 'Hours since Dawn')

g2 <- ggplot(all_ALOHA_hr_robust, aes_string(y = 'diam')) +
    geom_histogram(alpha = 0.5, fill = "grey30") +
    geom_hline(data = med_diam, aes(yintercept = median)) +
    theme_bw(base_size = 20) +
    scale_y_continuous(labels = NULL, sec.axis = sec_axis(~ 1000*coef*.^expo, name = unname(latex2exp::TeX('Carbon Quota (fg C cell$^{-1}$)')))) +
    ggplot2::theme(legend.position = 'none') +
    facet_wrap(vars(pop), ncol = 1, scales = 'free_y') +
    ggplot2::labs(x = 'No. hours', y = '')

fig_name <- paste0(save_path, 'HOT_diameter_Qc_summary.pdf')
pdf(fig_name, width = 8, height = 12)
    grid.arrange(g1, g2, ncol = 2)
dev.off()

###################
### RHYTHMICITY ###
###################

# Rain results

rain_file1 <- paste0(save_path, 'Data/HOT_rain_results_2021_12_25.csv')  # peaks (relevant for Qc)
big_rain_peak <- read.csv(rain_file1)
rain_file2 <- paste0(save_path, 'Data/HOT_rain_results_2022_01_03_troughs.csv')   # troughs (relevant for abundance)
big_rain_trough <- read.csv(rain_file2)

big_rain <- rbind(subset(big_rain_peak, param == 'Qc'), subset(big_rain_trough, param == 'abundance'))

t <- as.POSIXct(big_rain$first_peak, format = '%Y-%m-%d %H:%M:%S', tz = "HST")
pk_hr <- lubridate::hour(t)
big_rain$peak_hour <- pk_hr
sig <- big_rain$pVal <= 0.05    # Significant p-values--are they periodic or not?
big_rain$periodic <- sig
big_rain$pop <- gsub("picoeuk", 'euk', big_rain$pop)
big_rain$pop <- factor(big_rain$pop, levels = names(group.colors))

cruise_cut <- c('HOT317', 'KOK1606', 'KOK1608', 'MGL1704', 'KM1906', 'KOK1515')   # These cruises either have some outlier values or pass through Station ALOHA so quickly, they are not helpful for diel patterns
big_rain <- subset(big_rain, !(cruise %in% cruise_cut))

# Group into 2-hour blocks
TOD_2 <- cut(big_rain$peak_hour, breaks = c(-0.5, 1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5),
    labels = c( '00:00 - 01:59', '02:00 - 03:59', '04:00 - 05:59', '06:00 - 07:59', '08:00 - 09:59', '10:00 - 11:59', '12:00 - 13:59', '14:00 - 15:59', '16:00 - 17:59', '18:00 - 19:59', '20:00 - 21:59', '22:00 - 23:59'))
TOD_2 <- factor(TOD_2, levels = c('Aperiodic', '00:00 - 01:59', '02:00 - 03:59', '04:00 - 05:59', '06:00 - 07:59', '08:00 - 09:59', '10:00 - 11:59', '12:00 - 13:59', '14:00 - 15:59', '16:00 - 17:59', '18:00 - 19:59', '20:00 - 21:59', '22:00 - 23:59'))
TOD_2[which(big_rain$periodic == FALSE)] <- 'Aperiodic'
big_rain$TOD_2 <- TOD_2

# abundance
rain_trough <- subset(big_rain, param == "abundance")
n.cruise <- length(unique(rain_trough$cruise))

count_trough <- rain_trough %>%
    group_by(pop, TOD_2) %>%
    summarise(n = n())

fig_name <- paste0(save_path, "HOT_rain_abundance_bar_2hr_nocolor.pdf")
pdf(fig_name, width = 8, height = 10)
g <- ggplot2::ggplot(rain_trough, aes(TOD_2, fill = TOD_2)) +
    geom_rect(aes(xmin = -Inf, xmax = 4.5, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    geom_rect(aes(xmin = 11, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    ggplot2::geom_bar(aes(y = (..count..)/n.cruise), alpha = 1, colour = 'black', fill = 'white', show.legend = FALSE) +
    ggplot2::theme_bw(base_size = 22) +
    #ggplot2::scale_fill_manual(values = day_colors) +
    ggplot2::scale_x_discrete(labels = c('NA', '1', '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', '23')) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
    ggplot2::facet_wrap(vars(pop), ncol = 1) +
    ggpubr::rotate_x_text() +
    ggplot2::labs(x = 'Abundance Minima Hour (HST)', y = 'Percent of Cruises')
print(g)
dev.off()

# QC
rain_peak <- subset(big_rain, param == "Qc")
n.cruise <- length(unique(rain_peak$cruise))

rain_peak$TOD_2 <- factor(rain_peak$TOD_2, levels = c('Aperiodic', '00:00 - 01:59', '02:00 - 03:59', '04:00 - 05:59', '06:00 - 07:59', '08:00 - 09:59', '10:00 - 11:59', '12:00 - 13:59', '14:00 - 15:59', '16:00 - 17:59', '18:00 - 19:59', '20:00 - 21:59', '22:00 - 23:59'))

fig_name <- paste0(save_path, "HOT_rain_Qc_bar_2hr_nocolor.pdf")
pdf(fig_name, width = 8, height = 10)
g <- ggplot2::ggplot(rain_peak, aes(TOD_2, fill = TOD_2)) +
    geom_rect(aes(xmin = -Inf, xmax = 4.5, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    geom_rect(aes(xmin = 11, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    ggplot2::geom_bar(aes(y = (..count..)/n.cruise), alpha = 1, colour = 'black', fill = 'white', show.legend = FALSE) +
    ggplot2::theme_bw(base_size = 22) +
    #ggplot2::scale_fill_manual(values = day_colors, name = "Time of Day (HST)", drop = FALSE) +
    ggplot2::scale_x_discrete(labels = c('NA', '1', '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', '23'), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
    ggplot2::facet_wrap(vars(pop), ncol = 1) +
    ggpubr::rotate_x_text() +
    ggplot2::labs(x = 'Carbon Quota Maxima Hour (HST)', y = 'Percent of Cruises')
print(g)
dev.off()

###############
### BIOMASS ###
###############

# Get particulate carbon from HOT
HOT_POC <- read.csv(paste0(save_path, "Data/HOT_PC_ATP.csv"), skip = 6, header = F)
POC_head <- read.csv(paste0(save_path, "Data/HOT_PC_ATP.csv"), skip = 3, nrows = 1)
colnames(HOT_POC) <- colnames(POC_head)

HOT_POC$cruise <- paste0('HOT', substr(as.character(HOT_POC$botid), 1, 3))
HOT_POC[HOT_POC == -9] <- NA
HOT_POC$date <- str_pad(HOT_POC$date, 6, side = "left", pad = 0)
HOT_POC$time <- str_pad(HOT_POC$time, 6, side = "left", pad = 0)
HOT_POC$date <- as.POSIXct(paste0(HOT_POC$date, " ", HOT_POC$time), format = "%m%d%y %H%M%S", tz = "UTC")

HOT_POC$day <- lubridate::yday(HOT_POC$date)
HOT_POC$year <- lubridate::year(HOT_POC$date)
HOT_POC$month <- lubridate::month(HOT_POC$date) 

# Split double cruises by bumping late month cruises to the next month for KM2011/HOT323, HOT319 borders Jan and Feb
mo_list_date <- c(as.POSIXct("2020-01-31 02:14:22", tz = "UTC"), as.POSIXct("2020-09-28 00:59:58", tz = "UTC"), as.POSIXct("2020-09-29 01:07:07", tz = "UTC"))    
ind_mo <- which(HOT_POC$date %in% mo_list_date)
HOT_POC$month[ind_mo] <- HOT_POC$month[ind_mo] + 1

HOT_Clive <- HOT_POC %>%  
    dplyr::group_by(cruise) %>%
    dplyr::summarize(PC = mean(pc, na.rm = TRUE), ATP = mean(atp, na.rm = TRUE), month = mean(month), year = mean(year))

PC_only <- HOT_POC[!is.na(HOT_POC$pc), ]
HOT_Clive <- merge(HOT_Clive, PC_only[, c("date", "cruise")])   # Carry the time stamp for PC

# Estimate percentage live particulate carbon from ATP based on Henderikx-Freitas et al, 2021 and Christian & Karl, 1994.
# Note that PC and ATP are taken on different days
HOT_Clive$pc_live_250 <- (HOT_Clive$ATP*250/HOT_Clive$PC)/(12.01*10^3)  # percentage of live carbon and convert to g ATP/g C
HOT_Clive$pc_live_400 <- HOT_Clive$ATP*400/HOT_Clive$PC/(12.01*10^3)
HOT_Clive$pc_live_150 <- (HOT_Clive$ATP*150/HOT_Clive$PC)/(12.01*10^3) 

umolpkg2ugpL <- 12.01*1036/(10^3)   # Convert umol C/kg to ug C/L.
HOT_Clive$C_live_150 <- HOT_Clive$PC*umolpkg2ugpL*HOT_Clive$pc_live_150
HOT_Clive$C_live_250 <- HOT_Clive$PC*umolpkg2ugpL*HOT_Clive$pc_live_250
HOT_Clive$C_live_400 <- HOT_Clive$PC*umolpkg2ugpL*HOT_Clive$pc_live_400
HOT_Clive <- na.omit(HOT_Clive)

# Compare per-cruise biomass estimates between methods
all_ALOHA_cruise <- aloha_long %>%
    group_by(pop, month, year) %>%
    summarize(mean_biomass = mean(biomass, na.rm = TRUE))

all_ALOHA_cruise$month <- as.numeric(as.character(all_ALOHA_cruise$month))
all_ALOHA_cruise$year <- as.numeric(as.character(all_ALOHA_cruise$year))

total_ALOHA_cruise <- all_ALOHA_cruise %>%
    group_by(month, year) %>%
    summarize(total_biomass = sum(mean_biomass, na.rm = TRUE))

PC_all <- merge(HOT_Clive, total_ALOHA_cruise, by = c('month', 'year'))
PC_all$SF2PC <- PC_all$total_biomass/PC_all$C_live_250
PC_all$bloom <- FALSE
ind_b1 <- which(PC_all$month == 7 & PC_all$year == 2016)	# Eukaryote blooms: July 2016, Aug 2016, Aug 2017, Aug 2019
ind_b2 <- which(PC_all$month == 8 & PC_all$year == 2016)
ind_b3 <- which(PC_all$month == 8 & PC_all$year == 2017)
ind_b4 <- which(PC_all$month == 8 & PC_all$year == 2019)
ind_bloom <- c(ind_b1, ind_b2, ind_b3, ind_b4)
PC_all$bloom[ind_bloom] <- TRUE

SF2PC <- PC_all %>%
	dplyr::group_by(bloom) %>%
	dplyr::summarize(SF2PC_mean = mean(SF2PC), SF2PC_sd = sd(SF2PC), SF2PC_min = min(SF2PC), SF2PC_max = max(SF2PC))

#########################
### C FIXATION: DAILY ###
#########################
# Take hourly means only from daylight hours.  Use data table "robust" with Qc outlier cruises removed and abundance >= 0.02.
# Remove days with fewer than 6 hours between first and last daily point.
sun <- suncalc::getSunlightTimes(date = as.Date(all_ALOHA_hr_robust$local_time, tz = "HST"), lat = lat, lon = lon, keep = c('sunrise', 'sunset'), tz = "HST")
all_ALOHA_hr_robust$sunrise <- sun$sunrise
all_ALOHA_hr_robust$sunset <- sun$sunset

all_ALOHA_hr_robust$light <- "night"
ind_day <- which(all_ALOHA_hr_robust$local_time < all_ALOHA_hr_robust$sunset & all_ALOHA_hr_robust$local_time > all_ALOHA_hr_robust$sunrise)
all_ALOHA_hr_robust$light[ind_day] <- "day"
all_ALOHA_hr_robust$daylength <- as.numeric(difftime(all_ALOHA_hr_robust$sunset, all_ALOHA_hr_robust$sunrise, units = "hours"))
all_ALOHA_light <- subset(as.data.frame(all_ALOHA_hr_robust), light == "day")
all_ALOHA_light$hour <- as.numeric(lubridate::hour(all_ALOHA_light$local_time))    # Otherwise the linear model treats it as categorical

Qc_mean <- all_ALOHA_light %>%
    dplyr::group_by(date, pop) %>%
    dplyr::summarize(Qc_mean = mean(Qc, na.rm = TRUE), data_day = last(sunhour) - first(sunhour), n = length(Qc))#, gap = diff(sunhour))

# Estimate specific growth via exponential growth in cellular carbon
ALOHA_lm_day <- all_ALOHA_light %>%
    dplyr::group_by(date, pop) %>%
    tidyr::nest() %>%
    dplyr::mutate(
        #model = purrr::map(data, ~lm(Qc ~ hour, data = .)),
        model = purrr::map(data, ~lm(log(Qc) ~ sunhour, data = .)),
        tidied = purrr::map(model, broom::tidy)
    ) %>%
    tidyr::unnest(tidied)

ALOHA_lm_day_df <- as.data.frame(ALOHA_lm_day[, c("pop", 'date', 'term', 'estimate', 'std.error', 'p.value')])
exp_slope <- subset(ALOHA_lm_day_df, term == 'sunhour')
    colnames(exp_slope) <-  c("pop", "date", "term", "r", "r_se", "r_p")
exp_int <- subset(ALOHA_lm_day_df, term == '(Intercept)')
    colnames(exp_int) <-  c("pop", "date", "term", "log_Qc0", "log_Qc0_se", "log_Qc0_p")

tform_exp <- merge(subset(exp_slope, select = -c(term)), subset(exp_int, select = -c(term)), by = c('pop', 'date'))
tform_exp <- merge(Qc_mean, tform_exp)

tform_exp$r_p[which(tform_exp$r_p >= 0.05)] <- NA   # Remove insignificant linear regressions
tform_exp <- tform_exp[!is.na(tform_exp$r_p), ]   # Remove flagged p.values
tform_exp$r[which(tform_exp$r <= 0)] <- NA   # Remove negative C fixation estimates
tform_exp <- tform_exp[!is.na(tform_exp$r), ]  # Remove lines with NA C fixation estimate
tform_exp$data_day[which(tform_exp$data_day < 6)] <- NA    # Tag days when first and last data points are less than 6 hours apart
tform_exp <- tform_exp[!is.na(tform_exp$data_day), ]  # Remove lines with short data days

# Detect outliers

r_stats <- tform_exp %>%
    group_by(pop) %>%
    summarize(r_mean = mean(r), r_sd = sd(r), r_med = median(r), r_IQR = IQR(r), 
              r_25 = quantile(r, probs = c(0.25)), r_75 = quantile(r, probs = c(0.75)),
              r_min = min(r), r_max = max(r))

r_stats$mean_plus_6sd <- r_stats$r_mean + 6*r_stats$r_sd

tform_exp$r[which(tform_exp$r > 0.3)] <- NA   # Remove exceptionally high rates
tform_exp <- tform_exp[!is.na(tform_exp$r), ]

tform_exp$rmax <- tform_exp$r + tform_exp$r_se
tform_exp$rmin <- tform_exp$r - tform_exp$r_se

tform_exp$year <- lubridate::year(tform_exp$date)
tform_exp$month <- lubridate::month(tform_exp$date)
# Split double cruises by bumping late month cruises to the next month for KM2011/HOT323
mo_list_date <- c("2020-09-27", "2020-09-28", "2020-09-29")    
ind_mo <- which(as.character(tform_exp$date) %in% mo_list_date)
tform_exp$month[ind_mo] <- tform_exp$month[ind_mo] + 1

fig_name <- paste0(save_path, "HOT_C_specific_growth_daily_sd.pdf")
pdf(fig_name, width = 15, height = 8)
g <- ggplot2::ggplot(tform_exp, aes(x = year, y = r)) +
    ggplot2::geom_linerange(ggplot2::aes(x = year, ymin = rmin, ymax = rmax), color = col_seaflow) +
    ggplot2::geom_point(pch = 21, size = 3, fill = col_seaflow) +
    ggplot2::facet_grid(rows = vars(pop), cols = vars(month)) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::labs(y = unname(latex2exp::TeX('Net scatter-based specific growth rate (hr$^{-1}$)')), x = 'Year')
print(g)
dev.off()

##########################
### Primary production ###
##########################

# Station ALOHA Primary Productivity array
HOT_pp_read <- read.csv(paste0(save_path, 'Data/HOT_primary_production.txt'), skip = 3, na.strings = "-9")
units <- HOT_pp_read[1, ]
HOT_pp <- HOT_pp_read[-1, 2:15]
HOT_pp <- as.data.frame(sapply(HOT_pp, as.numeric))
HOT_pp$date <- stringr::str_pad(HOT_pp$date, 6, pad = "0")
HOT_pp$stime <- stringr::str_pad(HOT_pp$stime, 4, pad = "0")
HOT_pp$etime <- stringr::str_pad(HOT_pp$etime, 4, pad = "0")
d_HOT <- as.POSIXct(paste0(HOT_pp$date), format = '%y%m%d', tz = 'HST')
dt1_HOT <- as.POSIXct(paste0(HOT_pp$date, ' ', HOT_pp$stime), format = '%y%m%d %H%M', tz = 'HST')
dt2_HOT <-  as.POSIXct(paste0(HOT_pp$date, ' ', HOT_pp$etime), format = '%y%m%d %H%M', tz = 'HST')
HOT_pp$Date <- d_HOT
HOT_pp$date_start <- dt1_HOT
HOT_pp$date_end <- dt2_HOT

HOT_pp$l12[which(HOT_pp$l12 == -9)] <- NA
HOT_pp$year <- lubridate::year(HOT_pp$Date)
HOT_pp$day <- lubridate::yday(HOT_pp$Date)
HOT_pp$month <- lubridate::month(HOT_pp$Date)
# Split double cruises by bumping late month cruises to the next month for KM2011/HOT323
mo_list_date <- c("2020-09-26")   
ind_mo <- which(as.character(HOT_pp$Date) %in% mo_list_date)
HOT_pp$month[ind_mo] <- HOT_pp$month[ind_mo] + 1

# C fixation from SeaFlow

abund_dawn <- all_ALOHA_light %>%
    group_by(pop, cruise, date) %>%
    summarise(abundance_dawn = first(abundance), Qc_dawn = first(Qc), Qc_min = min(Qc))

SF_PP <- merge(abund_dawn, tform_exp, by = c('date', 'pop'))  # Note:  date is local date.  It keeps sequential hours of daylight together

sun <- suncalc::getSunlightTimes(date = as.Date(SF_PP$date, tz = "HST"), lat = lat, lon = lon, keep = c('sunrise', 'sunset'), tz = "HST")
daylength <- as.double(sun$sunset - sun$sunrise, units = "hours")
SF_PP$daylength <- daylength

SF_PP$year <- lubridate::year(SF_PP$date) 
SF_PP$day <- lubridate::yday(SF_PP$date) 
SF_PP$PP_exp <- SF_PP$abundance_dawn*(exp(SF_PP$log_Qc0))*(exp(SF_PP$r*SF_PP$daylength) - 1) # Estimate net primary production
SF_PP$pop <- factor(SF_PP$pop, levels = c('croco', 'euk', 'synecho', 'prochloro'))    # Reverse order to put most numerous on bottom

# mg/m^3 = ug/L = pg/uL

# Cruise mean

PP_cruise <- HOT_pp %>%
    group_by(month, year) %>%
    summarize(mean_PP = mean(l12, na.rm = TRUE))

SF_cruise <- SF_PP %>%
    group_by(pop, month, year) %>%
    summarize(mean_PP = mean(PP_exp, na.rm = TRUE))

# Plot biomass and PP together
# HOT
PCH <- HOT_Clive[, c('month', 'year')]
PCH$param <- 'Biomass'
PCH$value <- HOT_Clive$C_live_250
PCH$min <- HOT_Clive$C_live_150
PCH$max <- HOT_Clive$C_live_400

PPH <- PP_cruise[, c('month', 'year')]
PPH$param <- 'Productivity'
PPH$value <- PP_cruise$mean_PP
PPH$min <- NA
PPH$max <- NA

C_HOT <- rbind(PCH, PPH)
C_HOT$param <- factor(C_HOT$param, levels = c('Biomass', 'Productivity'))
C_HOT$mo_yr <- paste0(month.abb[C_HOT$month], ' ', C_HOT$year)
SF_moyr <- unique(paste0(month.abb[aloha_long$month], ' ', aloha_long$year))   # List of month-year pairs with SeaFlow data
C_HOT_incl <- C_HOT[which(C_HOT$mo_yr %in% SF_moyr), ]         # Allow only HOT data with SeaFlow counterpart

### SeaFlow cruise means ###
BMSF <- all_ALOHA_cruise[, c('pop', 'month', 'year')]
BMSF$param <- 'Biomass'
BMSF$value <- all_ALOHA_cruise$mean_biomass

PPSF <- SF_cruise[, c('pop', 'month', 'year')]
PPSF$param <- 'Productivity'
PPSF$value <- SF_cruise$mean_PP

C_SF <- rbind(BMSF, PPSF)
C_SF$param <- factor(C_SF$param, levels = c('Biomass', 'Productivity'))
C_SF$pop <- factor(C_SF$pop, levels = c('croco', 'euk', 'synecho', 'prochloro'))

fig_name <- paste0(save_path, "HOT_all_carbon_stacked_per_cruise.pdf")
pdf(fig_name, width = 15, height = 8)
p <- ggplot2::ggplot(C_SF, aes(x = year)) +
    ggplot2::geom_bar(aes(y = value, fill = pop), alpha = 0.5, color = NA, position = 'stack', stat = "identity") +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::scale_fill_manual(values = group.colors) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "population")) +
    ggplot2::geom_point(data = C_HOT_incl, aes(y = value, color = 'HOT'), fill = 'white', pch = 21, size = 3, alpha = 1) +
    ggplot2::geom_linerange(data = C_HOT_incl, aes(x = year, ymin = min, ymax = max, color = 'HOT')) +
    ggplot2::scale_color_manual(name = "", values = c("HOT" = "black")) +
    ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
    ggplot2::facet_grid(cols = vars(month), rows = vars(param), scales = 'free_y') +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::labs(x = 'Year', y = unname(latex2exp::TeX('Productivity ($\\mu$g C L$^{-1}$ d$^{-1}$) or Biomass ($\\mu$g C L$^{-1}$)')))
print(p)
dev.off()


################################
### Crocosphaera N2 fixation ###
################################

croco_A <- dplyr::filter(all_ALOHA_day[, c("day", "month", "year", "pop", "abundance_mean", "biomass_mean")], pop == "croco")
croco_PP <- dplyr::filter(SF_PP[, c("date", "pop", "r", "year", "month", "PP_exp")], pop == "croco")
croco_PP$day <- lubridate::yday(croco_PP$date)
croco_Nfix <- merge(croco_A, croco_PP)

croco_Nfix$bloom <- "No"
ind1 <- which(croco_Nfix$month == 8 & croco_Nfix$year == 2016)
ind2 <- which(croco_Nfix$month == 8 & croco_Nfix$year == 2017)
ind3 <- which(croco_Nfix$month == 8 & croco_Nfix$year == 2019)
ind4 <- which(croco_Nfix$month == 9 & croco_Nfix$year == 2019)
croco_Nfix$bloom[c(ind1, ind2, ind3, ind4)] <- "Yes"

Npcell <- 7.3/(0.4*10^6)	# mol N cell^-1 d^-1 (Wilson et al, 2017)
Npcell_max <- Npcell + 1.5/(0.4*10^6)	# sd of croco N2 fixation is 1.5 nmol N L^-1 d^-1
Npcell_min <-  Npcell - 1.5/(0.4*10^6)	

croco_Nfix$Nfix <- croco_Nfix$abundance_mean*(10^6)*Npcell
croco_Nfix$Nmax <- croco_Nfix$abundance_mean*(10^6)*Npcell_max
croco_Nfix$Nmin <- croco_Nfix$abundance_mean*(10^6)*Npcell_min

C2N <- 6.6 	# Wilson et al, 2017
g2mol <- 12.01	# Atomic weight of carbon
ng2ug <- 10^3	# convert ng to ug
croco_Nfix$Csupp <- croco_Nfix$Nfix*C2N*g2mol/(ng2ug)	# Amount of carbon biomass that can be supported by N2 fixed (in ug C L^-1)
croco_Nfix$Csupp_max <- croco_Nfix$Nmax*C2N*g2mol/(ng2ug)
croco_Nfix$Csupp_min <- croco_Nfix$Nmin*C2N*g2mol/(ng2ug)

g1 <- ggplot(croco_Nfix) +
  geom_abline(slope = 1, intercept = 0) +
  ggplot2::geom_linerange(ggplot2::aes(x = PP_exp, ymin = Csupp_min, ymax = Csupp_max, color = bloom)) +
  ggplot2::geom_point(aes(x = PP_exp, y = Csupp, fill = bloom), pch = 21, size = 3) +
  labs(x = unname(latex2exp::TeX("Crocosphaera scatter-based production ($\\mu$g C L$^{-1}$ d$^{-1}$)")), y = unname(latex2exp::TeX("Carbon biomass supported by Crocosphaera N$_2$ fixation ($\\mu$g C L$^{-1}$ d$^{-1}$)"))) +
  theme_bw(base_size = 18)

g2 <- ggplot(croco_Nfix) +
  geom_abline(slope = 1, intercept = 0) +
  ggplot2::geom_linerange(ggplot2::aes(x = biomass_mean, ymin = Csupp_min, ymax = Csupp_max, color = bloom)) +
  ggplot2::geom_point(aes(x = biomass_mean, y = Csupp, fill = bloom), pch = 21, size = 3) +
  labs(x = unname(latex2exp::TeX("Crocosphaera biomass ($\\mu$g C L$^{-1}$)")), y = unname(latex2exp::TeX("Carbon biomass supported by Crocosphaera N$_2$ fixation ($\\mu$g C L$^{-1}$ d$^{-1}$)"))) +
  theme_bw(base_size = 18)

png(file = paste0(save_path, "Croco_N2_fix_support_C.png"), width = 6000, height = 3200, res = 300)
t <- textGrob("Crocosphaera carbon supported by N2 fixation", gp = gpar(fontsize = 18))
grid.arrange(g1, g2, ncol = 2, top = t)
dev.off()
