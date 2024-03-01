# HOT_figures puts together final figures of the Station aloha SeaFlow data for publication. 
#
# Started: 25/May/2023  Annette Hynes, UW
 
library(tidyverse)
library(popcycle)
library(grid)
library(gridExtra)
library(ggpubr)
library(suncalc)
library(viridis)
library(latex2exp)
library(readxl)
library(marmap)
library(oce)
library(broom)
library(grDevices)
library(dunn.test)

###############
### SeaFlow ###
###############

setwd("/R")

# Download SeaFlow data from Zenodo (doi.org/10.5281/zenodo.7154076)
url <- "https://zenodo.org/record/7154076/files/SeaFlow_dataset_v1.5.xlsx" 
file_name <- tempfile()
try(download.file(url, file_name, method = "auto"))

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

pop.labels <- c("Prochlorococcus", "Synechococcus", "Small eukaryotes", "Crocosphaera")
names(pop.labels) <- c("prochloro", "synecho", "euk", "croco")

aloha_long <- aloha_long %>% 
  mutate(pop = case_when(pop == 'picoeuk' ~ 'euk', TRUE ~ pop),
         pop = factor(pop, levels = names(group.colors)), Program = "SeaFlow")

# Local time
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

### Calculate daily means ###
aloha_long_mean <- aloha_long %>%
  select(!c(cruise, depth)) %>%
  group_by(Program, pop, year, yday, month) %>%
  summarise_all(list(mean = mean, sd = sd), na.rm = TRUE) %>%
  rename(abundance = abundance_mean,
         DateTime = time_mean) %>%
  mutate(date = paste(year, month))

##############
### Influx ###
##############
influx_file <- '../Data/HOT_influx.csv' # Data from HOT-DOGS 
influx <- read_csv(influx_file, skip = 2, col_types = cols("c","c","c","d","d","d","d","c"))[-1,-c(8)]     # abundances are in 10E5 cells/mL

# Date conversion
influx <- influx %>% 
	mutate(DateTime = strptime(paste(date,time), format = "%m%d%y %H%M%S", tz = 'HST'),
		year = lubridate::year(DateTime),
		month = lubridate::month(DateTime),
		yday = lubridate::yday(DateTime),
		hour = lubridate::hour(DateTime),
		date = paste(year, month))

# Convert abundances from 10E5 cells/mL to cells/uL
post_influx <- influx %>% 
	mutate(pbact = case_when(pbact < 0 ~ NA, TRUE ~ pbact * 100),
		sbact = case_when(sbact < 0 ~ NA, TRUE ~ sbact * 100),
    ebact = case_when(ebact < 0 ~ NA, TRUE ~ ebact * 100)) %>%
  rename(prochloro = pbact, synecho = sbact, euk = ebact)

# Only keep influx data after SeaFlow started (December, 2014)
influx_keep <- post_influx %>% filter(DateTime > lubridate::as_date('2014-12-01'))

# Convert from short to long format for easier plotting
influx_long <- influx_keep %>% 
	tidyr::pivot_longer(cols = c(prochloro, synecho, euk), names_to = "pop", 
	  values_to = "abundance") %>%
	mutate(pop = factor(pop, levels = names(group.colors)), abundance_sd = NA,
    Program = "HOT")

### Use only the HOT points that have SeaFlow data
influx_long <- influx_long[influx_long$date %in% unique(aloha_long_mean$date),]

################################
### MERGE SeaFlow and Influx ###
################################
all_day <- bind_rows(influx_long %>% select(Program, DateTime, year, month, yday, pop, abundance, abundance_sd), 
	aloha_long_mean %>% select(Program, DateTime, year, month, yday, pop, abundance, abundance_sd)) %>%
  mutate(Program = factor(Program, levels = c("SeaFlow", "HOT")))

###################
### DAILY MEANS ###
###################

### abundance: averaged daily, faceted by month ###
prog.colors <- c(SeaFlow =  viridis::viridis(9)[5], HOT = '#FFFFFFFF')

g <- all_day %>% 
  arrange(Program) %>%
  ggplot2::ggplot(aes(x = year, y = abundance, fill = Program)) +
  ggplot2::geom_linerange(aes(ymin = abundance - abundance_sd, ymax =  abundance + abundance_sd)) +
  ggplot2::geom_point(size = 3, pch = 21) +
  ggplot2::facet_grid(pop ~ month, scales = 'free_y', labeller = labeller(pop = pop.labels)) +
  ggplot2::theme_bw(base_size = 18) +
  ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggplot2::scale_fill_manual(values = prog.colors) +
  ggplot2::labs(y = unname(latex2exp::TeX('Abundance (10$^6$ cells L$^{-1}$)')), x = 'Year')

fig_name <- paste0("../Figures/HOT_abundance_daily_point_on_station_sd_SF.pdf")
pdf(fig_name, width = 15, height = 9)
print(g)
dev.off()

### Abundance property-property plot ###
# Hourly means for abundance
all_day_p <- all_day %>%
  pivot_wider(names_from = Program, values_from = c(abundance, abundance_sd), 
    id_cols = c(year, yday, pop))

# plot
p <- list()
lim <- c(350, 6, 8)
i <- 1
for (phyto in c("prochloro", "synecho", "euk")){
  p[[i]] <- all_day_p %>%
    filter(pop == phyto) %>%
    ggplot(aes(abundance_HOT, abundance_SeaFlow)) + 
    geom_linerange(aes(ymin = abundance_SeaFlow - abundance_sd_SeaFlow, ymax =  abundance_SeaFlow + abundance_sd_SeaFlow)) +
    geom_point(size = 3, pch = 21, fill = "white")+
    geom_abline(slope = 1, intercept = 0, col = "red3", lty= 2) +
    xlim(0, lim[i]) +
    ylim(0, lim[i]) + 
    theme_bw(base_size = 18) +
    labs(x = "", y = "", title = pop.labels[[i]])
  i <- i + 1
}

fig_name <- "../Figures/HOT_abundance_SF_vs_Influx.pdf"
pdf(fig_name, width = 12, height = 4)
fig <- gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]],  
	nrow = 1,widths=c(2, 2, 2, 1),
    left =  grid::textGrob(unname(latex2exp::TeX('SeaFlow abundance (10$^6$ cells L$^{-1}$)')), rot = 90), 
    bottom = grid::textGrob(unname(latex2exp::TeX('Influx abundance (10$^6$ cells L$^{-1}$)'))))
dev.off()

### Mean + 2 sd to define bloom ###

bloom <- all_day %>%
    group_by(Program, pop) %>%
    mutate(abund_mean = mean(abundance, na.rm = TRUE), abund_sd = sd(abundance, na.rm = TRUE))

bloom_sf <- subset(bloom, Program == "SeaFlow")
bloom_sf$plus_2sd <- bloom_sf$abund_mean + 2*bloom_sf$abund_sd

g <- ggplot2::ggplot(bloom, aes(x = year, y = abundance, fill = Program)) +
    ggplot2::geom_linerange(ggplot2::aes(x = year, ymin = abundance - abundance_sd, ymax = abundance + abundance_sd)) +
    ggplot2::geom_point(pch = 21, size = 3) +
    ggplot2::geom_hline(data = bloom_sf, aes(yintercept = abund_mean, linetype = "Mean"), color = 'black') +
    ggplot2::geom_hline(data = bloom_sf, aes(yintercept = plus_2sd, linetype = '2.0 sd'), color = 'firebrick3') +
    scale_linetype_manual(name = "Mean + 2 sd", values = c(1, 2),
        guide = guide_legend(override.aes = list(color = c("firebrick3", "black")))) +
    ggplot2::facet_grid(rows = vars(pop), cols = vars(month), scales = 'free_y', labeller = labeller(pop = pop.labels)) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::scale_fill_manual(values = prog.colors) +
    ggplot2::labs(y = unname(latex2exp::TeX('Abundance (10$^6$ cells L$^{-1}$)')), x = 'Year', title = 'Define blooms by abundance standard deviation')

fig_name <- "../Figures/HOT_abundance_bloom_sd.pdf"
pdf(fig_name, width = 15, height = 9)
    print(g)
dev.off()

####################
### HOURLY MEANS ###
####################

# find sunrise and sunset
# START AT DAWN -- Do this before getting hourly means.  Subtracting sunrise after led to weird hour rounding resulting in some hours missing and some duplicated.
# Mid-month sunrise and sunset
midmonth <- as.Date(paste0(lubridate::year(aloha_long$local_time), '-', lubridate::month(aloha_long$local_time), '-15'), tz = 'HST')
sun_mid <- suncalc::getSunlightTimes(date = midmonth, lat = lat, lon = lon, keep = c('sunrise', 'sunset'), tz = "HST")
sunrise_mid <- lubridate::hour(sun_mid$sunrise) + lubridate::minute(sun_mid$sunrise)/60

all_aloha_hr <-  aloha_long %>%
  mutate(suntime = local_time - (sunrise_mid * 3600),
         sundate = lubridate::as_date(suntime),
         sunhour = lubridate::hour(suntime)) %>%
  dplyr::group_by(Program, cruise, pop, sunhour, sundate) %>%
  dplyr::summarize_all(mean, na.rm = TRUE) %>%
  mutate(date = lubridate::as_date(local_time)) %>%
  arrange(time)

# Split double cruises by bumping late month cruises to the next month
all_aloha_hr <- all_aloha_hr %>% 
	mutate(month = case_when(cruise == "KM1709" ~ month + 1, cruise == "KM2011" ~ month + 1, TRUE ~ month), 
		month = stringr::str_pad(month, 2, side = "left", pad = 0))
                    
### Curation
cruise_cut <- c('KOK1606', 'MGL1704', 'KM1906')   # These cruises pass through Station aloha so quickly, they are not helpful for diel patterns
all_aloha_hr_no_outlier <- subset(all_aloha_hr, !(cruise %in% cruise_cut))

ind_285 <- which(all_aloha_hr_no_outlier$cruise == "KOK1608" & all_aloha_hr_no_outlier$pop == "croco")  # Croco diameter is very strange and high during this cruise
all_aloha_hr_no_outlier <- all_aloha_hr_no_outlier[-ind_285, ]
all_aloha_hr_no_outlier$year <- as.factor(all_aloha_hr_no_outlier$year)

# Remove low abundance points first.
all_aloha_hr_robust <- subset(all_aloha_hr_no_outlier, abundance > 0.02) # 0.048 = About 30 cells per 3-min file

### Plot  DIAMETER and QC during time of day

# Get number of days in each population:
all_aloha_hr_robust <- all_aloha_hr_robust %>% 
  group_by(pop) %>%
  mutate(n_days = length(unique(sundate)))

# Draw with dual axes

coef <- (0.261*((4*pi/3)^0.860))/(2^(3*0.860))  # relationship between Qc and diameter in Ribalet et al, 2019
expo <- 3*0.860

med_diam <- all_aloha_hr_robust %>%
    dplyr::group_by(pop) %>%
    summarize(median = median(diam, na.rm = TRUE))

g1 <- all_aloha_hr_robust %>%
    ggplot(aes(x = sunhour, y = diam, group = sundate)) +
    geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    geom_line(alpha = 0.25, color = prog.colors['SeaFlow'], linewidth = 2) +
    theme_bw(base_size = 20) +
    scale_x_continuous(breaks = seq(0, 24, by = 6), expand = c(0, 0)) +
    scale_y_continuous(sec.axis = sec_axis(~ coef*.^expo, labels = NULL)) +
    facet_wrap(vars(pop), ncol = 1, scales = 'free_y', labeller = labeller(pop = pop.labels)) +
    geom_text(aes(label = paste0(' n = ', n_days), x = -Inf, y = Inf), hjust = 0, vjust = 1.5, check_overlap = T, size = 6) +
    ggplot2::labs(y = unname(latex2exp::TeX('Diameter ($\\mu$m)')), x = 'Hours since sunrise')

g2 <- all_aloha_hr_robust %>%
    ggplot(aes(y = diam, group = sundate)) +
    geom_histogram(alpha = 0.5, fill = "grey30") +
    geom_hline(data = med_diam, aes(yintercept = median)) +
    theme_bw(base_size = 20) +
    scale_y_continuous(labels = NULL, sec.axis = sec_axis(~ 1000*coef*.^expo, name = unname(latex2exp::TeX('Carbon quota (fg C cell$^{-1}$)')))) +
    ggplot2::theme(legend.position = 'none') +
    facet_wrap(vars(pop), ncol = 1, scales = 'free_y', labeller = labeller(pop = pop.labels)) +
    ggplot2::labs(x = 'No. hours', y = '')

fig_name <- '../Figures/HOT_diameter_Qc_summary.pdf'
pdf(fig_name, width = 8, height = 12)
    gridExtra::grid.arrange(g1, g2, ncol = 2)
dev.off()

###################
### RHYTHMICITY ###
###################

# Rain results

rain_file <- '../Data/HOT_rain_results_2023_08_14.csv'  # peaks (relevant for Qc)
big_rain <- read.csv(rain_file)

sig <- big_rain$pVal <= 0.05    # Significant p-values--are they periodic or not?
big_rain$periodic <- sig
big_rain$pop <- factor(big_rain$pop, levels = names(group.colors))

# allow for aperiodic
TOD <- factor(big_rain$peak_hour, levels = c('Aperiodic', as.character(0:23)))
TOD[which(big_rain$periodic == FALSE)] <- 'Aperiodic'
big_rain$TOD <- TOD

# abundance
rain_trough <- subset(big_rain, param == "abundance")
n.cruise <- length(unique(rain_trough$cruise))

count_trough <- rain_trough %>%
    group_by(pop, TOD) %>%
    summarise(n = n())

fig_name <- "../Figures/HOT_rain_abundance_bar_nocolor.pdf"
pdf(fig_name, width = 8, height = 10)
g <- ggplot2::ggplot(rain_trough, aes(TOD, fill = TOD)) +
    geom_rect(aes(xmin = -Inf, xmax = 8, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    geom_rect(aes(xmin = 20, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    ggplot2::geom_bar(aes(y = after_stat(count)/n.cruise), alpha = 1, colour = 'black', fill = 'white', show.legend = FALSE) +
    ggplot2::theme_bw(base_size = 22) +
    ggplot2::scale_x_discrete(labels = c('NA', as.character(0:23)), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
    ggplot2::facet_wrap(vars(pop), ncol = 1, labeller = labeller(pop = pop.labels)) +
    ggpubr::rotate_x_text() +
    ggplot2::labs(x = 'Abundance minima hour (HST)', y = 'Percent of cruises')
print(g)
dev.off()


# QC
rain_peak <- subset(big_rain, param == "Qc")
n.cruise <- length(unique(rain_peak$cruise))

count_peak <- rain_peak %>%
  group_by(pop, TOD) %>%
  summarise(n = n())

fig_name <- "../Figures/HOT_rain_Qc_bar_nocolor.pdf"
pdf(fig_name, width = 8, height = 10)
g <- ggplot2::ggplot(rain_peak, aes(TOD, fill = TOD_2)) +
    geom_rect(aes(xmin = -Inf, xmax = 8, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    geom_rect(aes(xmin = 20, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
    ggplot2::geom_bar(aes(y = after_stat(count)/n.cruise), alpha = 1, colour = 'black', fill = 'white', show.legend = FALSE) +
    ggplot2::theme_bw(base_size = 22) +
    ggplot2::scale_x_discrete(labels = c('NA', as.character(0:23)), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
    ggplot2::facet_wrap(vars(pop), ncol = 1, labeller = labeller(pop = pop.labels)) +
    ggpubr::rotate_x_text() +
    ggplot2::labs(x = 'Carbon quota maxima hour (HST)', y = 'Percent of cruises')
print(g)
dev.off()


###############
### BIOMASS ###
###############

# Get particulate carbon from HOT
HOT_POC <- read.csv(paste0("../Data/HOT_PC_ATP.csv"), skip = 6, header = F)
POC_head <- read.csv(paste0("../Data/HOT_PC_ATP.csv"), skip = 3, nrows = 1)
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

umolpkg2ugpL <- 12.01*1024/(10^3)   # Convert umol C/kg to ug C/L using seawater density at Station ALOHA.
HOT_Clive$PC_ugpL <- HOT_Clive$PC*umolpkg2ugpL
HOT_Clive$C_live_150 <- HOT_Clive$PC_ugpL*HOT_Clive$pc_live_150
HOT_Clive$C_live_250 <- HOT_Clive$PC_ugpL*HOT_Clive$pc_live_250
HOT_Clive$C_live_400 <- HOT_Clive$PC_ugpL*HOT_Clive$pc_live_400

HOT_Clive <- na.omit(HOT_Clive)

# Compare per-cruise biomass estimates between methods
all_aloha_cruise <- aloha_long %>%
    group_by(pop, month, year) %>%
    summarize(mean_biomass = mean(biomass, na.rm = TRUE))

all_aloha_cruise$month <- as.numeric(as.character(all_aloha_cruise$month))
all_aloha_cruise$year <- as.numeric(as.character(all_aloha_cruise$year))

total_aloha_cruise <- all_aloha_cruise %>%
    group_by(month, year) %>%
    summarize(total_biomass = sum(mean_biomass, na.rm = TRUE))

PC_all <- merge(HOT_Clive, total_aloha_cruise, by = c('month', 'year'))
PC_all$SF2PC_tot <- PC_all$total_biomass/PC_all$PC_ugpL
PC_all$SF2PC_live <- PC_all$total_biomass/PC_all$C_live_250
PC_all$bloom <- FALSE
ind_b1 <- which(PC_all$month == 7 & PC_all$year == 2016)	# Eukaryote blooms: July 2016, Aug 2016, Aug 2017, Aug 2019
ind_b2 <- which(PC_all$month == 8 & PC_all$year == 2016)
ind_b3 <- which(PC_all$month == 8 & PC_all$year == 2017)
ind_b4 <- which(PC_all$month == 8 & PC_all$year == 2019)
ind_bloom <- c(ind_b1, ind_b2, ind_b3, ind_b4)
PC_all$bloom[ind_bloom] <- TRUE

SF2PC <- PC_all %>%
	dplyr::group_by(bloom) %>%
	dplyr::summarize(live_mean = mean(SF2PC_live), live_sd = sd(SF2PC_live), live_min = min(SF2PC_live), live_max = max(SF2PC_live),
	                 tot_mean = mean(SF2PC_tot), tot_sd = sd(SF2PC_tot), tot_min = min(SF2PC_tot), tot_max = max(SF2PC_tot))


#########################
### C FIXATION: DAILY ###
#########################
# Take hourly means only from daylight hours.  Use data table "robust" with Qc outlier cruises removed and abundance >= 0.02.
# Remove days with fewer than 6 hours between first and last daily point.
sun <- suncalc::getSunlightTimes(date = as.Date(all_aloha_hr_robust$local_time, tz = "HST"), lat = lat, lon = lon, keep = c('sunrise', 'sunset'), tz = "HST")
all_aloha_hr_robust$sunrise <- sun$sunrise
all_aloha_hr_robust$sunset <- sun$sunset

all_aloha_hr_robust$light <- "night"
ind_day <- which(all_aloha_hr_robust$local_time < all_aloha_hr_robust$sunset & all_aloha_hr_robust$local_time > all_aloha_hr_robust$sunrise)
all_aloha_hr_robust$light[ind_day] <- "day"
all_aloha_hr_robust$daylength <- as.numeric(difftime(all_aloha_hr_robust$sunset, all_aloha_hr_robust$sunrise, units = "hours"))
all_aloha_light <- subset(as.data.frame(all_aloha_hr_robust), light == "day")
all_aloha_light$hour <- as.numeric(lubridate::hour(all_aloha_light$local_time))    # Otherwise the linear model treats it as categorical

Qc_mean <- all_aloha_light %>%
    dplyr::group_by(date, pop) %>%
    dplyr::summarize(Qc_mean = mean(Qc, na.rm = TRUE), data_day = last(sunhour) - first(sunhour), n = length(Qc))#, gap = diff(sunhour))

# Estimate specific growth via exponential growth in cellular carbon
aloha_lm_day <- all_aloha_light %>%
    dplyr::group_by(date, pop) %>%
    tidyr::nest() %>%
    dplyr::mutate(
        model = purrr::map(data, ~lm(log(Qc) ~ sunhour, data = .)),
        tidied = purrr::map(model, broom::tidy)
    ) %>%
    tidyr::unnest(tidied)

aloha_lm_day_df <- as.data.frame(aloha_lm_day[, c("pop", 'date', 'term', 'estimate', 'std.error', 'p.value')])
exp_slope <- subset(aloha_lm_day_df, term == 'sunhour')
    colnames(exp_slope) <-  c("pop", "date", "term", "r", "r_se", "r_p")
exp_int <- subset(aloha_lm_day_df, term == '(Intercept)')
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

tform_exp$rmax <- tform_exp$r + tform_exp$r_se
tform_exp$rmin <- tform_exp$r - tform_exp$r_se

tform_exp$year <- lubridate::year(tform_exp$date)
tform_exp$month <- lubridate::month(tform_exp$date)

# Split double cruises by bumping late month cruises to the next month for KM2011/HOT323
mo_list_date <- c("2020-09-27", "2020-09-28", "2020-09-29")    
ind_mo <- which(as.character(tform_exp$date) %in% mo_list_date)
tform_exp$month[ind_mo] <- tform_exp$month[ind_mo] + 1

fig_name <- "../figures/HOT_C_specific_growth_daily_sd.pdf"
pdf(fig_name, width = 15, height = 9)
g <- ggplot2::ggplot(tform_exp, aes(x = year, y = r)) +
    ggplot2::geom_linerange(ggplot2::aes(x = year, ymin = rmin, ymax = rmax), color = prog.colors[1]) +
    ggplot2::geom_point(pch = 21, size = 3, fill = prog.colors[1]) +
    ggplot2::facet_grid(rows = vars(pop), cols = vars(month), labeller = labeller(pop = pop.labels)) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::labs(y = unname(latex2exp::TeX('Net scatter-based cellular growth rate (h$^{-1}$)')), x = 'Year')
print(g)
dev.off()

##########################
### Primary production ###
##########################

# Station aloha Primary Productivity array
HOT_pp_read <- read.csv(paste0('../Data/HOT_primary_production.txt'), skip = 3, na.strings = "-9")
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

abund_dawn <- all_aloha_light %>%
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
PCH$value <- HOT_Clive$PC_ugpL

PPH <- PP_cruise[, c('month', 'year')]
PPH$param <- 'Productivity'
PPH$value <- PP_cruise$mean_PP

C_HOT <- rbind(PCH, PPH)
C_HOT$param <- factor(C_HOT$param, levels = c('Biomass', 'Productivity'))
C_HOT$mo_yr <- paste0(month.abb[C_HOT$month], ' ', C_HOT$year)
SF_moyr <- unique(paste0(month.abb[aloha_long$month], ' ', aloha_long$year))   # List of month-year pairs with SeaFlow data
C_HOT_incl <- C_HOT[which(C_HOT$mo_yr %in% SF_moyr), ]         # Allow only HOT data with SeaFlow counterpart

### SeaFlow cruise means ###
BMSF <- all_aloha_cruise[, c('pop', 'month', 'year')]
BMSF$param <- 'Biomass'
BMSF$value <-  all_aloha_cruise$mean_biomass

PPSF <- SF_cruise[, c('pop', 'month', 'year')]
PPSF$param <- 'Productivity'
PPSF$value <- SF_cruise$mean_PP

C_SF <- rbind(BMSF, PPSF)
C_SF$param <- factor(C_SF$param, levels = c('Biomass', 'Productivity'))
C_SF$pop <- factor(C_SF$pop, levels = c('croco', 'euk', 'synecho', 'prochloro'))

fig_name <- "../Figures/HOT_all_carbon_stacked_per_cruise.pdf"
pdf(fig_name, width = 15, height = 9)
p <- ggplot2::ggplot(C_SF, aes(x = year)) +
    ggplot2::geom_bar(aes(y = value, fill = pop), alpha = 0.5, color = NA, position = 'stack', stat = "identity") +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::scale_fill_manual(values = group.colors, labels = pop.labels) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Population")) +
    ggplot2::geom_point(data = C_HOT_incl, aes(y = value, color = 'HOT'), fill = 'white', pch = 21, size = 3, alpha = 1) +
    ggplot2::scale_color_manual(name = "", values = c("HOT" = "black")) +
    ggplot2::scale_x_continuous(breaks=seq(2014, 2021, 1), labels=c("2014", '',  "2016", '', '2018', '', '2020', ''), minor_breaks = seq(2015, 2021, 2)) +
    ggplot2::facet_grid(cols = vars(month), rows = vars(param), scales = 'free_y') +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::labs(x = 'Year', y = unname(latex2exp::TeX('Productivity ($\\mu$g C L$^{-1}$ d$^{-1}$) or biomass ($\\mu$g C L$^{-1}$)')))
print(p)
dev.off()


################################
### Crocosphaera N2 fixation ###
################################

croco_A <- dplyr::filter(aloha_long_mean[, c("yday", "month", "year", "pop", "abundance", "biomass_mean")], pop == "croco")
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

croco_Nfix$Nfix <- croco_Nfix$abundance*(10^6)*Npcell
croco_Nfix$Nmax <- croco_Nfix$abundance*(10^6)*Npcell_max
croco_Nfix$Nmin <- croco_Nfix$abundance*(10^6)*Npcell_min

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

png(file = "../Figures/Croco_N2_fix_support_C.png", width = 6000, height = 3200, res = 300)
t <- grid::textGrob("Crocosphaera carbon supported by N2 fixation", gp = gpar(fontsize = 18))
grid.arrange(g1, g2, ncol = 2, top = t)
dev.off()

###################
### VARIABILITY ###
###################

### Abundance using lag ###

data_h <- aloha_long %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(time, unit = "hours")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(abundance)/mean(abundance, na.rm = TRUE)), gap = as.numeric(diff(time))) %>%
  dplyr::mutate(resolution = "hourly")
data_h <- dplyr::filter(data_h, gap < 5)   # exclude gaps >= 5 hours

data_d <- aloha_long %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(time, unit = "days")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(abundance)/mean(abundance, na.rm = TRUE)), gap = as.numeric(diff(time))) %>%
  dplyr::mutate(resolution = "daily")
data_d <- dplyr::filter(data_d, gap < 3)  # exclude gaps >= 3 days

data_m <- aloha_long %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(time, "months")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(abundance)/mean(abundance, na.rm = TRUE)), gap = as.numeric(diff(time))) %>%
  dplyr::mutate(resolution  = "monthly")
data_m <- dplyr::filter(data_m, gap < 45)  # exclude gaps >= 45 days

data_s <- aloha_long %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(time, "3 months")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(abundance)/mean(abundance, na.rm = TRUE)), gap = as.numeric(diff(time))) %>%
  dplyr::mutate(resolution = "seasonal")
data_s <- dplyr::filter(data_s, gap < 100)  # exclude gaps >= 100 days

data_y <- aloha_long %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(time, "years")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(abundance)/mean(abundance, na.rm = TRUE)), gap = as.numeric(diff(time))) %>%
  dplyr::mutate(resolution = "annual")

data_variance <- bind_rows(data_h, data_d, data_m, data_s, data_y) %>%
  mutate(resolution = factor(resolution, levels = c("hourly","daily","monthly","seasonal","annual")),
         pop = factor(pop, levels = c("prochloro","synecho","euk","croco")))

data_variance <- dplyr::filter(data_variance, !is.na(var))

box_stats <- data_variance %>%
  group_by(pop, resolution) %>%
  summarize(n = dplyr::n(), min_var = min(var, na.rm = TRUE), max_var = max(var, na.rm = TRUE), 
            q25 = quantile(var, probs = 0.25, na.rm = TRUE), q50 = quantile(var, probs = 0.5, na.rm = TRUE), 
            q75 = quantile(var, probs = 0.75, na.rm = TRUE))

box_stats <- box_stats %>% 
  group_by(pop, resolution) %>%
  mutate(IQR = q75 - q25, outlier1 = q25 - IQR*1.5, ymin = max(min_var, outlier1),
         outlier2 = q75 + IQR*1.5, ymax = min(max_var, outlier2)) # calculate whiskers

box_stats <- box_stats %>%
  group_by(pop) %>%
  mutate(offset = 0.1*max(ymax))

p1 <- ggplot(box_stats) +
  geom_boxplot(aes(x = resolution, ymin = ymin, lower = q25, middle = q50, upper = q75, ymax = ymax), stat = "identity", fill = "darkgrey", alpha = 0.5) +
  geom_text(aes(x = resolution, y = ymax + offset, label = n)) +
  ylab("Fold change in cell abundance") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid( pop ~ ., scales = "free_y", labeller = labeller(pop = pop.labels))

png(paste0("../Figures/Variability_gap.png"), width = 1000, height = 2000, res = 300)
print(p1)
dev.off()

### statistical test 

# Test normality using Shapiro-Wilk
var_norm <- data_variance %>%
  group_by(pop, resolution) %>%
  do(tidy(shapiro.test(.$var)))
# Not normal:  Pro-seasonal, euk-annual, croco-annual

pop_list <- unique(var_norm$pop)

# Multiple comparisons:  Kruskall-Wallis test, then pairwise Dunn rank sum test (and Wilcoxon for comparison though Dunn is considered more appropriate)

for(phyto in pop_list){
  print(phyto)
  this_pop <- subset(data_variance, pop == phyto)
  KW <- kruskal.test(x = this_pop$var, g = this_pop$resolution)
  print(KW)
  print(KW$p.value)
  if(KW$p.value < 0.01){
    PW <- pairwise.wilcox.test(x = this_pop$var, g = this_pop$resolution, p.adjust.method = "BH")
    print(PW)
    D <- dunn.test::dunn.test(x = this_pop$var, g = this_pop$resolution, method = "bh", alpha = .01)
    print(D)
  }
}

### Fold change abundance versus specific growth


data_d_dated <- aloha_long %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(time, unit = "days")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(abundance)/mean(abundance, na.rm = TRUE)), gap = as.numeric(diff(time)), date = time[1:length(time)-1]) %>%
  dplyr::mutate(resolution = "daily")
data_d_dated <- dplyr::filter(data_d_dated, gap < 3)  # exclude gaps >= 3 days

growth_fold_abund <- merge(data_d_dated, tform_exp, by = c("date", "pop"))

linreg <- growth_fold_abund %>%
  dplyr::nest_by(pop) %>%
  mutate(model = list(lm(var ~ r, data = data))) 

Rsq <- linreg %>% 
  summarise(rsq = summary(model)$r.squared)

g <- ggplot(growth_fold_abund, aes(x = r, y = var)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth") +
  ylab("Fold change in cell abundance") +
  xlab(latex2exp::TeX('Net scatter-based cellular growth rate (h$^{-1}$)')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid( pop ~ ., scales = "free_y", labeller = labeller(pop = pop.labels))

png(paste0("../Figures/growth_fold_abundance.png"), width = 1000, height = 2000, res = 300)
print(g)
dev.off()
    
### Specific Growth ###

data_d <- tform_exp %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(date, unit = "days")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(n)/mean(n, na.rm = TRUE)), gap = as.numeric(diff(date))) %>%
  dplyr::mutate(resolution = "daily")
data_d <- dplyr::filter(data_d, gap < 3)  # exclude gaps >= 3 days

data_m <- tform_exp %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(date, "months")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(n)/mean(n, na.rm = TRUE)), gap = as.numeric(diff(date))) %>%
  dplyr::mutate(resolution  = "monthly")
data_m <- dplyr::filter(data_m, gap < 45)  # exclude gaps >= 45 days

data_s <- tform_exp %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(date, "3 months")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(n)/mean(n, na.rm = TRUE)), gap = as.numeric(diff(date))) %>%
  dplyr::mutate(resolution = "seasonal")
data_s <- dplyr::filter(data_s, gap < 100)  # exclude gaps >= 100 days

data_y <- tform_exp %>% 
  dplyr::group_by(pop, time = lubridate::floor_date(date, "years")) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  dplyr::summarize(var = abs(diff(n)/mean(n, na.rm = TRUE)), gap = as.numeric(diff(date))) %>%
  dplyr::mutate(resolution = "annual")

data_variance <- bind_rows(data_d, data_m, data_s, data_y) %>%
  mutate(resolution = factor(resolution, levels = c("daily","monthly","seasonal","annual")),
         pop = factor(pop, levels = c("prochloro","synecho","euk","croco")))

data_variance <- dplyr::filter(data_variance, !is.na(var))

box_stats <- data_variance %>%
  group_by(pop, resolution) %>%
  summarize(n = dplyr::n(), min_var = min(var, na.rm = TRUE), max_var = max(var, na.rm = TRUE), 
            q25 = quantile(var, probs = 0.25, na.rm = TRUE), q50 = quantile(var, probs = 0.5, na.rm = TRUE), 
            q75 = quantile(var, probs = 0.75, na.rm = TRUE))

box_stats <- box_stats %>% 
  group_by(pop, resolution) %>%
  mutate(IQR = q75 - q25, outlier1 = q25 - IQR*1.5, ymin = max(min_var, outlier1),
         outlier2 = q75 + IQR*1.5, ymax = min(max_var, outlier2)) # calculate whiskers

box_stats <- box_stats %>%
  group_by(pop) %>%
  mutate(offset = 0.1*max(ymax))

p2 <- ggplot(box_stats) +
  geom_boxplot(aes(x = resolution, ymin = ymin, lower = q25, middle = q50, upper = q75, ymax = ymax), stat = "identity", fill = "darkgrey", alpha = 0.5) +
  geom_text(aes(x = resolution, y = ymax + offset, label = n)) +
  ylab(" Fold change in cellular growth rate") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid( pop ~ ., labeller = labeller(pop = pop.labels)) 

png(paste0("../Figures/Variability_growth_gap.png"), width = 1000, height = 2000, res = 300)
print(p2)
dev.off()

### statistical test 

# Test normality using Shapiro-Wilk
var_norm <- data_variance %>%
  group_by(pop, resolution) %>%
  do(tidy(shapiro.test(.$var)))
# Not normal:  Error

pop_list <- unique(data_variance$pop)

# Multiple comparisons:  Kruskall-Wallis test

for(phyto in pop_list){
  print(phyto)
  this_pop <- subset(data_variance, pop == phyto)
  KW <- kruskal.test(x = this_pop$var, g = this_pop$resolution)
  print(KW)
  print(KW$p.value)
 # if(KW$p.value < 0.01){
    PW <- pairwise.wilcox.test(x = this_pop$var, g = this_pop$resolution, p.adjust.method = "BH")
    print(PW)
    D <- dunn.test::dunn.test(x = this_pop$var, g = this_pop$resolution, method = "bh", alpha = .01)
    print(D)
  #}
}

# No significant differences among distributions

###################
### CHLOROPHYLL ###
###################

# Get Chla from HOT
HOT_Chl <- read.csv(paste0("../Data/HOT_Chl.csv"), skip = 6, header = F)
Chl_head <- read.csv(paste0("../Data/HOT_Chl.csv"), skip = 3, nrows = 1)
colnames(HOT_Chl) <- colnames(Chl_head)

HOT_Chl$cruise <- paste0('HOT', substr(as.character(HOT_Chl$botid), 1, 3))
HOT_Chl[HOT_Chl == -9] <- NA
HOT_Chl$date <- str_pad(HOT_Chl$date, 6, side = "left", pad = 0)
HOT_Chl$time <- str_pad(HOT_Chl$time, 6, side = "left", pad = 0)
HOT_Chl$date <- as.POSIXct(paste0(HOT_Chl$date, " ", HOT_Chl$time), format = "%m%d%y %H%M%S", tz = "UTC")

HOT_Chl$day <- lubridate::yday(HOT_Chl$date)
HOT_Chl$year <- lubridate::year(HOT_Chl$date)
HOT_Chl$month <- lubridate::month(HOT_Chl$date) 

# Split double cruises by bumping late month cruises to the next month for KM2011/HOT323, HOT319 borders Jan and Feb
mo_list_date <- c(as.POSIXct("2020-01-31 02:14:22", tz = "UTC"), as.POSIXct("2020-09-28 00:59:58", tz = "UTC"), as.POSIXct("2020-09-29 01:07:07", tz = "UTC"))    
ind_mo <- which(HOT_Chl$date %in% mo_list_date)
HOT_Chl$month[ind_mo] <- HOT_Chl$month[ind_mo] + 1

HOT_Chl_mo <- HOT_Chl %>%  
  dplyr::group_by(cruise) %>%
  dplyr::summarize(chla = mean(chl, na.rm = TRUE), month = mean(month), year = mean(year))

# Merge with SeaFlow biomass

Chl_all <- merge(HOT_Chl_mo, total_aloha_cruise, by = c('month', 'year'))

Chl_all$bloom <- FALSE
ind_b1 <- which(Chl_all$month == 7 & Chl_all$year == 2016)	# Eukaryote blooms: July 2016, Aug 2016, Aug 2017, Aug 2019
ind_b2 <- which(Chl_all$month == 8 & Chl_all$year == 2016)
ind_b3 <- which(Chl_all$month == 8 & Chl_all$year == 2017)
ind_b4 <- which(Chl_all$month == 8 & Chl_all$year == 2019)
ind_bloom <- c(ind_b1, ind_b2, ind_b3, ind_b4)
Chl_all$bloom[ind_bloom] <- TRUE

g <- ggplot(Chl_all) +
  geom_point(aes(x = total_biomass, y = chla), pch = 21, size = 3, fill = "white") +
  geom_abline(slope = 1/128, intercept = 0) +
  geom_abline(slope = 1/80, intercept = 0, linetype = 2) +
  geom_abline(slope = 1/176, intercept = 0, linetype = 2) +
  labs(x = latex2exp::TeX('Biomass measured by SeaFlow ($\\mu$g C L$^{-1}$)'), y = latex2exp::TeX('Fluorometric chlorophyll a ($\\mu$g C L$^{-1}$)')) +
  theme_bw(base_size = 18)

png(paste0("../Figures/Chla_vs_biomass.png"), width = 2000, height = 2000, res = 300)
  print(g)
dev.off()
