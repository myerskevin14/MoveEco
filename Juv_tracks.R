################## Prep. workspace ###############################

# load packages relevant to this script:
library( tidyverse ) 
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( sf ) 
library(purrr)
library(lubridate)
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )


load("TracksWorkspace.RData")

setwd("C:/Users/KEVINMYERS14/Documents/EEB 627/MoveEco/GOEA_move")

#load movement data
dataraw <- read_csv ("GOEA_move.csv")
table(dataraw$unit)
#load eagle data
records <- read.csv("JuvRec.csv")
table(records$unit)
#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read( "~/EEB 627/MoveEco//Data/BOPNCA_Boundary.shp" )


##############

#######################################################################
######## cleaning data ###############################################
# Data cleaning is crucial for accurate analysis # 
# Trapping records provide info on when individuals were fitted #
# with transmitters.#
colnames( records )

head(records)
#convert date.time column to correct POSIXct format using lubridate
records$fledge <- dmy(records$fledge_date)

records$disp <- dmy(records$disp_date)

# Add a day so that we can ignore records from the trapping day #
# and start only with  those from the following day:
records$StartDate <- records$fledge + lubridate::days(20)
records$EndDate <- records$disp - lubridate::days(2) ####Error: Incompatible classes: <character> + <Period> ###
 ##########################################################################################
#convert to day of year (1-366)
records$StartDay <- lubridate::yday( records$StartDate )
records$EndDay <- lubridate::yday( records$EndDate )
#check 
head( records); dim( records)

unique(records$unit)
head(records)
###################################################################
# Clean GPS data

# Data provided by blm did not include hdop, vdop, or time_to_fix so I cannot filter data by quality of the fixes. 

#start by creating a new dataframe to store cleaned location records:
datadf <- dataraw 
#which columns do we have?
colnames( datadf )
head(datadf)
table(datadf$unit)
# Filter to remove inaccurate locations
datadf <- datadf %>% 
  dplyr::filter( latitude > 0 ) %>% 
  #remove superfluous columns
 dplyr::select(-...1 )

 #view
head( datadf ); dim( datadf )
#How many rows did we remove?
# Answer: 1 (...1)
#
dim( dataraw ) - dim( datadf )
# What % of data did we loose?
# Answer: 45860-1964
1964/45860
# 4% of data

# We also need to set a time column containing date and time information #
# in POSIX format (as required by amt)#
# We rely on lubridate for this. If you haven't used lubridate before #
# go here: https://cran.r-project.org/web/packages/lubridate/vignettes/lubridate.html
# to learn more about how to easily manipulate time and dates in R #

# Data are stored in year, month, day, hour, minute and some have seconds while others do not. You can either remove or add seconds depending on if it matters to you. 

# To ADD seconds
datadf$date.time <- as.character(datadf$date.time)
datadf$date.time <- ifelse(nchar(datadf$date.time) == 16, 
                           paste0(datadf$date.time, ":00"), 
                           datadf$date.time)

head(datadf)
table(datadf$unit)
# To REMOVE seconds
 #datadf$date.time <- sub("(:\\d{2})$", "", datadf$date.time)  # Remove seconds if present
 #datadf$date <- ymd_hm(datadf$date.time, tz = "UTC")

# Convert 'date.time' to POSIXct in UTC
datadf <- datadf %>%
  mutate(date = ymd_hms(date.time, tz = "UTC")) %>%
 filter(!is.na(date))


datadf$date <- ymd_hms(datadf$date.time, tz = "UTC")

# Convert from UTC to MST (Mountain Standard Time)
datadf$date <- with_tz(datadf$date, tz = "MST")

# and create new column where we convert it to posixct
datadf$ts <- as.POSIXct( datadf$date )
#view
head( datadf ); dim( datadf )

# # check if any data are missing
all( complete.cases( datadf ) )
# # none so we can move on

# we also add year, month, and day of year information using lubridate
datadf <- datadf %>% 
  mutate( year = lubridate::year(date),
          mth = lubridate::month(date),
          wk = lubridate::week(date),
          jday = lubridate::yday(date) )
table(datadf$unit)


# We need to remove records for fixes that were recorded before the #
# units were fitted to the animals so we append relevant information #
# from the records dataframe. We do that by combining datadf to records df#
colnames(records)
colnames(datadf)


datadf <- records %>%  
  dplyr::select( unit, territory, sex, band, mass, hatch_date, 
                 fledge_date, StartDay, EndDay ) %>% 
  right_join( datadf, by = "unit" )
#view
head( datadf ); dim( datadf )

# filtered by 20 days after fledging and 2 days before dispersal
datadf <- datadf %>% 
  group_by( unit ) %>%
  #group_by(year) %>%
  dplyr::filter( jday >= StartDay ) %>% 
  dplyr::filter( jday <= EndDay ) %>% 
  #dplyr::filter( mth >= 5 & mth <= 8 ) %>% 
  ungroup() # You need to ungroup for code to work later!!! # and possibly adjust months considered the "breeding season"

table(datadf$unit)
datadf %>% dplyr::filter( unit == 132255 ) %>% tail()
#view
tail( datadf ); dim( datadf )

#3##### UP TO HERE@@@@@@@@

# serial IDs are cumbersome so we create a new individual ID column:
datadf <- datadf %>%
  group_by(unit, year) %>%
  mutate(id = cur_group_id()) %>% ungroup()

# Now we've cleaned it temporally and removed errors, but now we need to clean it spatially.

##################################################################
### Define coordinate system and projection for the data ######
# location data were recorded using WGS84 in lat long #
# We use the epsg code to define coordinate system for our sites #
# How? Google epsg WGS84 # First result should  take you here: #
# https://spatialreference.org/ref/epsg/wgs-84/ 
# where you can find that epgs = 4326 for this coordinate system #
# If you are not familiar with geographic data, vector, rasters & #
# coordinate systems go here: 
# https://bookdown.org/robinlovelace/geocompr/spatial-class.html #
# to learn more. #

# For amt, crs need to be provided using sp package so:
crsdata <- 4326# sp::CRS( "+init=epsg:4326" )
# We also want to transform the lat longs to easting and northings #
# using UTM. For this we need to know what zone we are in. Go: #
# http://www.dmap.co.uk/utmworld.htm
# We choose zone 11:
crstracks <- sf::st_crs( NCA_Shape )#sp::CRS( "+proj=utm +zone=11" )
#We convert the NCA shapefile to the same projection as our tracks
#NCA_Shape <- sf::st_transform( NCA_Shape, crstracks )
# We are now ready to make tracks using atm package
#We first check sample size #
table(datadf$unit)
# How many individuals have we dropped so far?
# 0 but some lost a significant amount of points 

# We can also get an idea of the data collection for each individual
# by plotting histograms
#sampling duration
ggplot( datadf, aes( x = jday, group = id ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( ) +
  facet_wrap( ~ id )

# What do the histograms tell you about the nature of the data #
# Sample size, intensity for different individuals? #
# Answer:  Some individuals have much more data than others and it varies by jday

#######################################################################
###### Creating tracks, calculating step lengths and turning angles ###
####              for all individuals at the same time:           #####
########################################################################
# A track is an amt object that creates a path based on the groupings that you care about.
head(datadf)
str(datadf)
#amt requires us to turn data into tracks for further analyses.
trks <- datadf %>% 
  #make track. Note you can add additional columns to it (lat, lon, and ts are the names from your df) You have to select the other columns that you want to include.
  amt::make_track( .y = latitude, .x = longitude, .t = ts, 
                   #define columns that you want to keep, relabel if you need:
                   id = id, territory = territory, unit = unit,
                   sex = sex, year = year, mth = mth, wk =wk, jday = jday,
                   #assign correct crs
                   crs = crsdata )

# Reproject to UTM to convert lat lon to easting northing:
#trks <- amt::transform_coords( trks, crstracks )
trks <- amt::transform_coords( trks, crs_to = crstracks )

#Turn into a tibble list by grouping and nest by individual IDs:
trks <- trks %>%  amt::nest( data = c(-"id" ))
#view
trks

# Creates a tibble of tibbles. Each row represents 1 bird and all their points. Should all have same number of columns, but different number of rows based on how many points were collected.

# We plot overall paths for each individual:
for( i in 1:dim(trks)[1]){
  a <- as_sf_points( trks$data[[i]] ) %>% 
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$id[i]), subtitle = paste0('Year:', trks$year[i]) ) +
    geom_sf(data = NCA_Shape, inherit.aes = FALSE ) +
    geom_sf() 
  print(a)
} 

# Which ones have migration paths?
# Answer:  almost all 
#
# Any ideas on how to remove migration data?
# Answer:  Remove by time? 
# Here we rely on NCA polygon, removing records that exist East of the #
# NCA. We can extra the extent of a polygon: many territories are outside of NCA boundaries to start with 

# assign crs
NCA_Shape <- st_transform(NCA_Shape, crs = crstracks)

# Create a 10 km buffer. Might want to make bigger?
NCA_buff <- st_buffer(NCA_Shape, dist = 50000)

sf::st_bbox( NCA_buff )

# Extract the coordinates of the bounding box
xmin <- as.numeric(st_bbox(NCA_buff)$xmin)
xmax <- as.numeric(st_bbox(NCA_buff)$xmax)
ymin <- as.numeric(st_bbox(NCA_buff)$ymin)
ymax <- as.numeric(st_bbox(NCA_buff)$ymax)

# filter points to bounding box
trks <- trks %>% 
  mutate(
    breeding = map(data, ~ filter(., x_ >= xmin & x_ <= xmax & y_ >= ymin & y_ <= ymax))
  )

#we check that it worked for our breeding data
for( i in 1:dim(trks)[1]){
  a <- as_sf_points( trks$breeding[[i]] ) %>% 
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$id[i]), subtitle = paste0('Year:', trks$year[i])  ) +
    geom_sf(data = NCA_Shape, inherit.aes = FALSE ) +
    geom_sf() 
  print(a)
} 
#still need to look at a couple of individuals in more detail

# # Note we created two other groups of tibbles for the breeding season
# sl = step length: How much were the animals moving. (you can use this code for turning angle too)
# # Plot step lengths
for( i in 1:dim(trks)[1]){
  a <-  steps( trks$breeding[[i]] ) %>% # if you change formats you need to carry columns because amt will get rid of it.
    #a <-  steps( trks$migrating[[i]] ) %>%
    mutate( jday = lubridate::yday( t1_ ) ) %>%
    group_by( jday ) %>%
    summarise( sl_ = log( sum(sl_) ) ) %>% # May have to log if step lengths are really long. log ( sum(sl_)). if you don't want to log it sl_ = sl-
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$id[i]), subtitle = paste0('Year:', trks$year[i])  ) +
    geom_line( aes( y = sl_, x = jday) )
  print(a)
}

# Some analyses require that your data not be correlated (kde) and others that they're equally spaced. 
# Use MEDIAN not MEAN
# You can resample data to get a more even fix rate.
# We focus on breeding season data:
# Estimate sampling rate for each individual by looping through 
# data using purr function map

# you might need a different sampling rate for each bird if the birds. Compare the one that you have the most data for at finer and coarser scales. If the coarse and fine scales are similar then you can lump them. If not you have to separate them.
sumtrks <- trks %>%  summarize( 
  map( breeding, amt::summarize_sampling_rate ) )
#view
sumtrks[[1]]

trks
# resample at 3 hours. Might have to section out data.
# Add tibbles with added step lengths calculated by bursts from #
# breeding season data:
trks.all <- trks %>% 
  mutate( red =  map( breeding, function(x) x %>%  
                        track_resample( rate = hours(3), 
                                        tolerance = minutes(30) ) ),
          steps = map( breeding, function(x) x %>%  
                         track_resample( rate = hours(3), 
                                         tolerance = minutes(30)) %>% 
                         steps_by_burst() ) )
#view
trks.all


#note that this creates a new set of tibbles called steps - that uses
# the breeding season data

# We can now unnest the dataframes of interest
#Starting with all breeding season data
trks.pfdp <- trks.all %>% dplyr::select( id, red ) %>%
  unnest( cols = red )
################################################################################
# Checking points per year, month, week

# Group by Alias and year, then summarize to get the count
points_by_unit_year <- trks.pfdp %>%
  group_by(unit, year) %>%
  summarize(count = n())

# Check
print(points_by_unit_year)

# Group by Alias and month
pts_mth <- trks.pfdp %>%
  group_by(unit, year, mth) %>%
  summarize(rows_per_month = n()) 

# Check
print(pts_mth)

# Group by week
pts_wk <- trks.pfdp %>%
  group_by(unit, year, wk) %>%
  summarize(rows_per_week = n()) 

# Check
table(pts_wk$unit)


# Get ranges
range(pts_mth$rows_per_month)
range(pts_wk$rows_per_week)
# 
#############################################################
########## step lengths and turning angles  ##################
#########################
# Step lengths don't match up well with my data. The two different types of transmitters and taking data at different intervals would probably mean I need a 5+ hour step length. 

#the step dataframes resampled at 5sec intervals
trks.steps <- trks.all %>% dplyr::select( id, steps ) %>%
  unnest( cols = steps )
head( trks.steps )

# We can plot step lengths by:
trks.steps %>%   
  ggplot(.) +
  #geom_density( aes( x = sl_, fill = as.factor(burst_)), alpha = 0.4 ) +
  geom_histogram( aes( x = sl_ ) ) +
  xlab("Step length" ) + 
  #ylim( 0, 0.01 ) + xlim(0, 2000 ) +
  theme_bw( base_size = 19 )  +
  theme( legend.position = "none" ) +
  facet_wrap( ~id, scales = 'free_y' )
#What does the plot tell us about the step lengths traveled by the individual?
# Answer: fairly consistent across individuals. 
#############################################################################
# Saving relevant objects and data ---------------------------------
#save breeding season data (not thinned)
write_rds( trks.pfdp, "trks.pfdp")
#save breeding season data (turned into steps)
write_rds( trks.steps, "trks.steps" )


#save workspace in case we need to make changes
save.image( "TracksWorkspace.RData" )

########## end of save #########################
############### END OF SCRIPT ########################################