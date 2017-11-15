######################################################################
## Do the data i/o for the meningococcal data. Code for this was taken
## from Benjamin and then modified
##
## Author: Michael Höhle <http://www.math.su.se/~hoehle>
## Date:   2017-09-15 (fimo)
######################################################################

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(magrittr)
library(scales)

Sys.setlocale("LC_TIME", "en_US.UTF-8")

## Menigococcal data from surveillance
library(surveillance)
data(imdepi)

## Aggregate across event types, age and gender
meningo_cases <- imdepi$events@data %>% as_tibble %>% mutate(count = 1L)

## Grab coordinates
meningo_coords <- coordinates(imdepi$events)
meningo_cases %<>%
  mutate(x_coord = meningo_coords[, 1], y_coord = meningo_coords[, 2])

## Add dates and remove unneeded columns
meningo_cases %<>%
  mutate(day = ceiling(time),
         date = as.Date("2001-12-31") + lubridate::days(day),
         year = as.integer(lubridate::year(date)),
         month = as.integer(lubridate::month(date))) %>%
  select(-eps.t, -eps.s, -.obsInfLength, -.sources, -.bdist, -.influenceRegion)

## Munge district-level data
district_data <- imdepi$stgrid %>% as_tibble %>%
  mutate(population = as.integer(popdensity * area))

tile_location <- tibble(tile = unique(district_data$tile),
                        location = 1:length(unique(district_data$tile)))
district_data <- left_join(district_data, tile_location, by = "tile")

## Get the total population (population is constant across time)
total_pop <- sum((district_data %>% filter(BLOCK == 1))$population)

## Monthly counts with covariates
meningo_monthly <- meningo_cases %>% select(tile, BLOCK) %>%
  mutate(count = 1L) %>%
  group_by(tile, BLOCK) %>%
  summarize(count = n()) %>%
  ungroup %>%
  right_join(district_data %>% select(tile, location, BLOCK, area, popdensity, population),
             by = c("tile", "BLOCK")) %>%
  mutate(count = ifelse(is.na(count), 0L, count)) %>%
  mutate(total_pop = total_pop) %>%
  rename(time = BLOCK) %>%
  arrange(location, time) %>%
  mutate(year = 2002 + floor((time - 1) / 12),
         month = ifelse(time %% 12 == 0, 12, time %% 12),
         date = as.Date(paste(year, month, "01", sep = "-")))

## Extract dates
dates <- (meningo_monthly %>% filter(tile == "01001") %>% select(date))$date


## Store result to be used further on.
saveRDS(meningo_monthly, file = "meningo_monthly.rds")
