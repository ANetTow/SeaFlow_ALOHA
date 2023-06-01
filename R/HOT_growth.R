# Statistical tests on specific growth rates

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

##############################
### C-specific growth rate ###
##############################

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
all_ALOHA_hr_robust <- subset(all_ALOHA_hr_no_outlier, abundance > 0.02) # 0.048 = About 30 cells per 3-min file

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

#########################
### Statistical tests ###
#########################

### Normality (Shapiro-Wilk's test)
growth_norm <- tform_exp %>%
    group_by(pop) %>%
    do(tidy(shapiro.test(.$r)))

# Pro and picoeuk deviate from normality (p < 0.05) so use a nonparameteric test.

pop_list <- unique(tform_exp$pop)

### Compare among populations

# Multiple comparisons:  Kruskall-Wallis test, then pairwise Dunn rank sum test (and Wilcoxon for comparison though Dunn is considered more appropriate)

KW <- kruskal.test(x = tform_exp$r, g = tform_exp$pop)
    print(KW)
    print(KW$p.value)
    if(KW$p.value < 0.01){
        PW <- pairwise.wilcox.test(x = tform_exp$r, g = tform_exp$pop, p.adjust.method = "BH")
        print(PW)
        D <- dunn.test::dunn.test(x = tform_exp$r, g = tform_exp$pop, method = "bh", alpha = .01)
        print(D)
    }

### Compare among seasons (within population)
tform_exp$season <- NA
tform_exp$season[tform_exp$month %in% c(12, 1, 2)] <- 'winter'
tform_exp$season[tform_exp$month %in% c(3, 4, 5)] <- 'spring'
tform_exp$season[tform_exp$month %in% c(6, 7, 8)] <- 'summer'
tform_exp$season[tform_exp$month %in% c(9, 10, 11)] <- 'fall'
tform_exp$season <- as.factor(tform_exp$season)

# Multiple comparisons:  Kruskall-Wallis test, then pairwise Dunn rank sum test (and Wilcoxon for comparison though Dunn is considered more appropriate)

for(phyto in pop_list){
    print(phyto)
    this_pop <- subset(tform_exp, pop == phyto)
    KW <- kruskal.test(x = this_pop$r, g = this_pop$season)
    print(KW)
    print(KW$p.value)
    if(KW$p.value < 0.01){
        PW <- pairwise.wilcox.test(x = this_pop$r, g = this_pop$season, p.adjust.method = "BH")
        print(PW)
        D <- dunn.test::dunn.test(x = this_pop$r, g = this_pop$season, method = "bh", alpha = .01)
        print(D)
    }
}

############
### PLOT ###
############

col_8 <- viridis::inferno(8)
season_cols <- c(winter = col_8[3], spring = col_8[4], summer = col_8[5], fall = col_8[7])
tform_exp$season <- factor(tform_exp$season, levels = c('winter', 'spring', 'summer', 'fall'))

# box & whisker
fig_name <- paste0(save_path, 'HOT_seasonal_C_growth_boxplot.png')
png(fig_name, width = 600, height = 800)
g <- ggplot(tform_exp, aes(x = season, y = r, fill = season)) +
    geom_boxplot() +
    ggplot2::scale_fill_manual(values = season_cols) +
    facet_wrap(vars(pop), ncol = 1) +
    labs(y = latex2exp::TeX('Specific growth rate (h$^{-1}$)')) +
    theme_bw(base_size = 18)
print(g)
dev.off()
