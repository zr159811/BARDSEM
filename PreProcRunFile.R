# Zachary Roman
# Run file for Romand and Brandt (2020)
# < ZacharyJoseph.Roman@uzh.ch >

library(rgdal)
library(spdep)
library(sf)
library(foreign)
library(haven)
library(ggplot2)
library(rstan)

# Allows rstan to parallelize MCMC chains across physical cores
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# SOUTHERN HOMICIDE DATA
dat1 <- read.dta("data/homicide1990.dta12")


# Shape FILE OGR
southMap <- readOGR("data/shapeData/south.shp",verbose = F)

# SHAPE FILE SF
map_sf<- read_sf("data/shapeData/south.shp")

# The /data/ folder contains codebooks associated with the data
# See Roman and Brandt (2020) for official citations on the data collection.



# UCR Arrests extension data
# http://doi.org/10.3886/E108164V3.
# https://www.openicpsr.org/openicpsr/project/108164/version/V3/view?path=/openicpsr/108164/fcr:versions/V3/county_ucr_offenses_known_1960_2017_dta.zip&type=file

# UCR DATA
dat2 <- read_dta("UCR_arrests/county_ucr_offenses_known_yearly_1960_2017.dta")
colnames(dat1)[1] <- "id"
dat2 <- dat2[dat2$year == "1990",]
dat2$fips_state_county<- as.numeric(dat2$fips_state_county)


map_sf2 <- merge(x = map_sf,y = dat2,by.x = "FIPSNO",by.y = "fips_state_county")

# Full data: Missing 2 counties
map_sf3 <- merge(x = map_sf2,y = dat1,by.x = "FIPSNO",by.y = "fips")

# Identify missing cases to manufacture W from OGR without them
southMap2 <- southMap[as.numeric(as.character(southMap$FIPS)) %in% map_sf3$FIPSNO,]

# Make W matrix
nbList = poly2nb(southMap2)
W = nb2mat(nbList, style="C")

# Normalize  W matrix
W <- ifelse(W > 0,1,0)
W_c <- W / rowSums(W)

neighbors_sf <- as(nb2lines(nbList, coords = coordinates(southMap2)), 'sf')
neighbors_sf <- st_set_crs(neighbors_sf, st_crs(x = southMap2))

map_sf4 <- st_as_sf(x = map_sf3)




# # # # # # # # # # # # # # # # #
# PLOTS 

# Subsetting map data to NC for example state level plot
map_sfX <-  map_sf4[map_sf4$state %in% c("north carolina","south carolina"),]
smapXX <- southMap2[southMap2$STATE_NAME  %in% c("North Carolina","South Carolina"),]
nbListX = poly2nb(smapXX)
neighbors_sfX <- as(nb2lines(nbListX, coords = coordinates(smapXX)), 'sf')
neighbors_sfX <- st_set_crs(neighbors_sfX, st_crs(x = smapXX))


ggplot() +
  geom_sf(data = map_sfX) +
  geom_sf(data = neighbors_sfX, lty = 1, size = 0.5) +
  geom_point(aes(x = coordinates(smapXX)[,1],
                 y = coordinates(smapXX)[,2])) +
  coord_sf() +
  xlab(" ") +
  ylab(" ") +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank()) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



# # # # # # # # # # # # # #
# Neighborhood connection plot

ggplot() +
  geom_sf(data = map_sf4) +
  geom_sf(data = neighbors_sf, lty = 1, size = 0.4) +
  geom_point(aes(x = coordinates(southMap2)[,1],
                 y = coordinates(southMap2)[,2]), size = 0.4) +
  coord_sf() +
  xlab("Longitude") +
  ylab("Latitude")   +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank()) +   theme_minimal() 


# # # # # # # # # # # # # #
# Homicide rate chloro-plot

ggplot() +
  geom_sf(data = map_sf4, aes(fill = scale(as.numeric(hrate)))) +
  scale_fill_continuous(low = "black",high = "red") +
  theme_minimal() +
  coord_sf() +
  xlab("Longitude") +
  ylab("Latitude")   


# # # # # # # # # # # # # # # # #
# Pre-Proc and Model

varZ <- c("sname",
          "cname",
          "cfips",
          "gini",
          "population",
          "actual_robbery_total",       #F1
          "actual_burg_total",          #F1
          "actual_theft_total",         #F1
          "actual_mtr_veh_theft_total", #F1
          "unemployment",  #F2
          "ln_income",     #F2
          "poverty",       #F2
          "edeprivation",  #F2
          "hrate",                #Eta
          "officers_assaulted",   #Eta
          "actual_assault_total", #Eta
          "actual_rape_total")    #Eta

dat <- map_sf4[,varZ]

#Process UCR data

dat$pop1000 <- dat$population/1000


# Rates reflect occurrences per 1,000 people
dat$ln_actual_robbery_total <- dat$actual_robbery_total/dat$pop1000
dat$ln_actual_burg_total <- dat$actual_burg_total/dat$pop1000
dat$ln_actual_theft_total <- dat$actual_theft_total/dat$pop1000
dat$ln_actual_mtr_veh_theft_total <- dat$actual_mtr_veh_theft_total/dat$pop1000
dat$ln_officers_assaulted <- dat$officers_assaulted/dat$pop1000
dat$ln_actual_assault_total <- dat$actual_assault_total/dat$pop1000
dat$ln_actual_rape_total <- dat$actual_rape_total/dat$pop1000


# # # # # # # #
# BARDSEM 

rt <- stanc("SEM_ENDLAG_impacts.stan")
sm <- stan_model( stanc_ret = rt , verbose = FALSE)

x <- data.frame(scale(dat$unemployment), # F1
                scale(dat$gini),         # F1
                scale(dat$poverty),      # F1
                scale(dat$edeprivation), # F1
                scale(dat$ln_actual_burg_total),         # F2
                scale(dat$ln_actual_theft_total),        # F2
                scale(dat$ln_actual_mtr_veh_theft_total))# F2


y <- data.frame(scale(dat$hrate), # ETA
                scale(dat$ln_actual_robbery_total),
                scale(dat$ln_officers_assaulted),
                scale(dat$ln_actual_assault_total))

datenstan <- list(N = nrow(dat),
                  Kx = ncol(x),
                  Ky =  ncol(y),
                  y =  y,
                  x = x,
                  W = W_c)

fit <- sampling(sm,
                data = datenstan,
                chains = 6,
                iter = 4000,
                warmup = 2000)

parz <-c("b1",
         "ly",
         "lx",
         "tx",
         "ty",
         "phi[1,2]",
         "phi[2,1]",
         "phi[2,2]",
         "phi[1,1]",
         "rho")

summary(fit,pars = parz)$summary







