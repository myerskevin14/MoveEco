#################################################################
# Script developed by Jen Cruz to clean and format location data #
# We also convert cleaned location data to tracks using atm package   # 
# We rely heavily on amt getting started vignette here:       #
# https://cran.r-project.org/web/packages/amt/vignettes/p1_getting_started.html#
#                                                               #
# Data are Prairie Falcon locations collected during Spring/Summer #
# of 2021 at Morley Nelson Birds of Prey NCA.                      #
# Data were collected for multiple individuals and at #
# different frequencies including 2 sec intervals when the individuals#
# were moving (every 2-3 days), and 30min fixes otherwise to define #
# breeding season range. # Frequency shifted to hourly once #
# individuals left their breeding grounds. #
#################################################################

################## Prep. workspace ###############################

# Install new packages from "CRAN" repository if you don't have them. # 
install.packages( "tidyverse" ) #actually a collection of packages 
install.packages( "amt" )
#trying to install amt directly from github
# install.packages( "devtools" )
# devtools::install_github("jmsigner/amt")
install.packages( "sf" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation and plotting
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt ) #creating tracks from location data
library( sf ) #handling spatial data
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

# Set working directory. This is the path to your Rstudio folder for this 
# project. If you are in your correct Rstudio project then it should be:
getwd()
# if so then:
workdir <- getwd()

# set path to where you can access your data #
# Note that the path will be different for your.#
#datapath <- "Z:/Common/PrairieFalcons/"
datapath <- "Data/"
#import GPS data# 
# Fixes are stored as separate CSV files for each individual
## We therefore create a function that imports multiple files at once:
load_data <- function( path ){
  # extract all file names in the folder
  myfiles <- dir( path, pattern = '\\.csv', full.names = TRUE )
  for( i in 1:length(myfiles) ){
    mydata <- read.csv( file = myfiles[i], 
              #remove white spaces  
              strip.white =TRUE, 
              #include column headers
              header = TRUE, 
              # read the serial column as a character instead of number:
              colClasses = c("serial" = "character") ) 
  # create df for first file and append rows for other files
   ifelse( i == 1,
           df <- mydata, 
           df <- bind_rows(df, mydata) )
  } 
  #return df for all individuals
  return( df )
}

#apply function to import all files as list of databases:
dataraw <- load_data( paste0(datapath, '2021/') )
#Note that the files are all in a subdirectory
head(dataraw)
# Import trapping records with details of when radiotrackers were 
# fitted to the individuals
records <- read.csv( file = paste0( datapath,"survey_0.csv" ),
                     #replaces those values with NA
                     na.strings = c(""," ","NA"), 
                     # include column headings
                     header = TRUE )
#check
head( records ); dim( records )

#import polygon of the NCA as sf spatial file:
#NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")
NCA_Shape <- sf::st_read( paste0( datapath, "BOPNCA_Boundary.shp") )
##############

#######################################################################
######## cleaning data ###############################################
# Data cleaning is crucial for accurate analysis # 
# Trapping records provide info on when individuals were fitted #
# with transmitters.#
colnames( records )
# we keep transmitter id, date and sex
records <- records %>% dplyr::select( Telemetry.Unit.ID, Sex, 
                                      Date.and.Time )
#view
records
#convert date to correct format using lubridate
records$StartDate <- lubridate::mdy_hms( records$Date.and.Time, 
                                    tz = "MST" )
# Add a day so that we can ignore records from the trapping day #
# and start only with  those from the following day:
records$StartDate <- records$StartDate + lubridate::days(1)
#convert to day of year
records$StartDay <- lubridate::yday( records$StartDate )
#unit IDs were missing starting number of their serial number #
# we append those so we can match it to the GPS serial IDs:
records$serial <- paste0( '894608001201',records$Telemetry.Unit.ID )

#check 
head( records); dim( records)

# Using serial ID unique to each individual found in df, add territory column
# In Territory column each serial ID is linked to its corresponding territory
records <- records %>%
  mutate(territory = case_when(
    endsWith(serial, "47221") ~ "SG",
    endsWith(serial, "47775") ~ "CRW",
    endsWith(serial, "47874") ~ "SDTP",
    endsWith(serial, "48120") ~ "PR_II",
    endsWith(serial, "46751") ~ "HHGS_DS",
    endsWith(serial, "46983") ~ "HHGS_US",
    endsWith(serial, "47197") ~ "Mac",
    endsWith(serial, "48229") ~ "CRW_new",
    endsWith(serial, "48377") ~ "CFR",
  ))

unique(records$territory)

###################################################################
# Clean GPS data
# GPS units often provide information on the quality of the fixes they #
# obtained.#
# The units from Cellular track technologies provide HDOP, VDOP and #
# time to fix information # 
# Start by viewing what those look like in the dataset #

hist( dataraw$vdop, breaks = 50 )
hist( dataraw$hdop, breaks = 50 )
hist( dataraw$time_to_fix )

# Remove 2D fixes and fixes where HDOP or VDOP ≥10 following #
# D’eon and Delparte (2005).#
# Also those where time to fix > 20min or with 0 satellites:

#start by creating a new dataframe to store cleaned location records:
datadf <- dataraw 
#which columns do we have?
colnames( datadf )
# Filter to remove inaccurate locations
datadf <- datadf %>% dplyr::filter( hdop < 10 ) %>%
  dplyr::filter( vdop < 10 ) %>%
  dplyr::filter( time_to_fix <= 20 ) %>% 
  dplyr::filter( nsats > 0 ) %>%
  dplyr::filter( lat > 0 ) %>% 
  #remove superfluous columns
  dplyr::select( -inactivity, -geo, -data_voltage, -solar_current, 
                 -solar_charge )

#view
head( datadf ); dim( datadf )
#How many rows did we remove?
## Answer: 28836, 274664 left 
# 303500-28836
dim( dataraw ) - dim( datadf )
dim(dataraw)
# What % of data did we loose?
# Answer: almost 10%
#28836/303500
# We also need to set a time column containing date and time information #
# in POSIX format (as required by amt)#
# We rely on lubridate for this. If you haven't used lubridate before #
# go here: https://cran.r-project.org/web/packages/lubridate/vignettes/lubridate.html
# to learn more about how to easily manipulate time and dates in R #
# Data are stored in year, month, day, hour, minute, second format in our data. 
# We define correct format with lubridate 
datadf$date <- lubridate::ymd_hms( datadf$GPS_YYYY.MM.DD_HH.MM.SS,
              tz = "UTC" )
datadf$date <- lubridate::with_tz( datadf$date, tz = "MST" )
# and create new column where we convert it to posixct
datadf$ts <- as.POSIXct( datadf$date )
#view
head( datadf ); dim( datadf )

# # check if any data are missing
all( complete.cases( datadf ) )
#TRUE # none so we can move on

# we also add month, week, and day of year information using lubridate
datadf <- datadf %>% 
  mutate( mth = lubridate::month(date),
          wk = lubridate::week(date),
          jday = lubridate::yday(date) )

# We need to remove records for fixes that were recorded before the #
# units were fitted to the animals so we append relevant information #
# from the records dataframe. We do that by combining datadf to records df#
datadf <- records %>%  dplyr::select( serial, territory, Sex, StartDay ) %>% 
  right_join( datadf, by = "serial" )
#view
head( datadf ); dim( datadf )
# Then using StartDay to filter records, removing those that occurred#
#  earlier when unit was turned on, but not fitted to animal #
datadf <- datadf %>% 
  group_by( serial ) %>% 
  dplyr::filter( jday > StartDay ) %>% ungroup()
#view
head( datadf ); dim( datadf )
# serial IDs are cumbersome so we create a new individual ID column:
datadf$id <- group_indices( datadf, serial )

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
crsdata <- 4326#
# We also want to transform the lat longs to easting and northings #
# using UTM. For this we need to know what zone we are in. Go: #
# http://www.dmap.co.uk/utmworld.htm
# We choose zone 11:
#note we can also use the crs of the NCA polygon: 
crstracks <- sf::st_crs( NCA_Shape )#sp::CRS( "+proj=utm +zone=11" )
# We are now ready to make tracks using atm package
#We first check sample size #
table( datadf$id, datadf$wk )
# How many individuals have we dropped so far?
# 1
# We can also get an idea of the data collection for each individual
# by plotting histograms
#sampling duration
ggplot( datadf, aes( x = jday, group = id ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( ) +
  facet_wrap( ~ id )
#speeds travelled
ggplot( datadf, aes( x = speed, group = id ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( ) +
  facet_wrap( ~ id )

# What do the histograms tell you about the nature of the data #
# Sample size, intensity for different individuals? #
## Answer: When in the year that the animals were data was collected 
#
#Why is the first bar on the speed histograms so tall?
# Answer: Data from transporting the unit to the field before deployment
#
# do we need to remove data based on these?
## Answer: No, outliers are not present 
#
#######################################################################
###### Creating tracks, calculating step lengths and turning angles ###
####              for all individuals at the same time:           #####
########################################################################
#amt requires us to turn data into tracks for further analyses.
trks <- datadf %>% 
  #make track. Note you can add additional columns to it
  amt::make_track( .y = lat, .x = lon, .t = ts, 
    #define columns that you want to keep, relabel if you need:
    id = id, territory = territory,
    sex = Sex, mth = mth, wk = wk,
    jday = jday, speed = speed, alt = alt, 
    #assign correct crs
    crs = crsdata )

# Reproject to UTM to convert lat lon to easting northing:
#trks <- amt::transform_coords( trks, crstracks )
trks <- amt::transform_coords( trks, crs_to = crstracks )
#Turn into a tibble list by grouping and nest by individual IDs:
trks <- trks %>%  amt::nest( data = -"id" )
#view
trks

# Remember we have multiple types of data including detailed data for flights #
# 3 times a week, 30min fixes during the day, then hourly fixes during #
# migration. We start by focusing on data during breeding season. #
# That means we need to remove migration locations.
# How do we know when individuals started migrating North?
# We plot overall paths for each individual:
for( i in 1:dim(trks)[1]){
  a <- as_sf_points( trks$data[[i]] ) %>% 
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$id[i]) ) +
    geom_sf(data = NCA_Shape, inherit.aes = FALSE ) +
    geom_sf() 
  print(a)
} 
# Which ones have migration paths?
# Answer: all of them
#
# Any ideas on how to remove migration data?
# Answer: decide a date cut off for sampling, limit by NCA shape as buffer.
# 
# Here we rely on NCA polygon, removing records that exist East of the #
# NCA. We can extra the extent of a polygon:
sf::st_bbox( NCA_Shape )
#Then use the Eastern-most coordinate to filter out data 
xmax <- as.numeric(st_bbox(NCA_Shape)$xmax) #627081.5
#Then use the Northern-most coordinate to filter out data 
ymax <- as.numeric(st_bbox(NCA_Shape)$ymax) + 10000 #627081.5

#subset those tracks less than as breeding and those > as migrating:
trks <- trks %>% mutate(
  breeding = map( data, ~ filter(., x_ < xmax ) ) )

trks <- trks %>% mutate(
  breeding = map( breeding, ~ filter(., y_ < ymax ) ) )

#some individuals come back to overwinter at the NCA and so #
# we need to remove those records as well #
# we do so using month column to remove anything after June
trks <- trks %>% mutate(
  breeding = map( breeding, ~filter(., mth < 7 ) ),  
  migrating = map( data, ~filter(., mth > 6 ) ),
  locals = map( migrating, ~ filter(., x_ < xmax ) ),
  locals = map(locals, ~ filter(., y_ < ymax ) ) 
)

# We focus on breeding season data for visualization as that is the 
# one of interest in latest weeks. 
for( i in 1:dim(trks)[1]){
  a <- as_sf_points( trks$locals[[i]] ) %>%  #change breeding to locals for wintering data
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$id[i]) ) +
    geom_sf(data = NCA_Shape, inherit.aes = FALSE ) +
    geom_sf() 
  print(a)
} 

# For homework plot the locals instead.#what do you note? are they all overwintering
# at the NCA? which ones are? List individuals here:
#Answer: 9,8,3,1,6,5,
#
# Despite us setting a sampling (fix) rate for our transmitters, bad weather, 
# thick canopy etc can cause fix attempts to fail. Our fix rate may therefore 
# not be exactly what we set it for. If we want measures of distance (step lengths),
# or we want to use this data for  AKDE, SSFs, iSSFs or HMMs (movement models) in discrete 
# time, we need fix rates to be equally spaced. The first step to do this is to
# estimate sampling rate for each individual by looping through 
# data using purr function map
sumtrks <- trks %>%  summarize( 
  map( breeding, amt::summarize_sampling_rate ) )
#view
sumtrks[[1]]
# This plots the actual sampling rate for each individual separately. #
# Look at the median. What is it?  (units reported on last column)
# Answer: 5 second tolerance on either side for high resolution rates.
#
# What about min and max? Those values give you an idea of the breaks in 
# your data. The median allows you to work out the most common fix rate. 
# Probably the one of interest for a lot of questions. #
# we will choose two sampling rates 30min (median value) for Range analysis #
# and RRFs and 5sec for movement questions, SSFs, iSSFs. 


trks.all <- trks %>% mutate(
  # Here we take breeding season data and resample at 5 seconds, allowing +- 4sec:
  highres = map( breeding, function(x) x %>%  
          track_resample( rate = seconds(5), 
                           tolerance = seconds(4) ) ),
  #Now repeat the process but with 30 min sampling rate 
  red = map(breeding, function( x ) x %>%  
               track_resample( rate = minutes(30),
               tolerance = minutes(5) ) ) ) 
#view
trks.all
#note that this creates a new set of tibbles called steps - that uses
# the breeding season data

#Now unnest the dataframes of interest
#Pull out the 5sec breeding season data
trks.breed <- trks.all %>% dplyr::select( id, highres ) %>% 
  unnest( cols = highres ) 
tail( trks.breed )

#you can also plot them once you unnested the resampled points for fine tune cleaning 
trks.breed %>% 
  #dplyr::filter( jday < 178 ) %>% 
as_sf_points( . ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf( aes( colour = as.factor(jday) ) ) +
  #plot separate for each individual
  facet_wrap( ~id )

#we remove remaining migrating tracks
trks.breed <- trks.breed %>% dplyr::filter( jday < 178 )

# Now for 30min intervals
trks.thin <- trks.all %>% dplyr::select( id, red ) %>% 
  unnest( cols = red ) %>% dplyr::filter( jday < 178 )
tail( trks.thin )

#Also migration data:
trks.mig <- trks.all %>% dplyr::select( id, migrating ) %>% 
  unnest( cols = migrating ) 
head( trks.mig )

### repeat here creating the unnesting for locals! ################################################################################


trks.locals <- trks.all %>% dplyr::select( id, highres ) %>% 
  unnest( cols = highres ) 
tail( trks.locals )

#you can also plot them once you unnested the resampled points for fine tune cleaning 
trks.locals %>% 
  #dplyr::filter( jday < 178 ) %>% 
  as_sf_points( . ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf( aes( colour = as.factor(jday) ) ) +
  #plot separate for each individual
  facet_wrap( ~id )

#we remove remaining migrating tracks
trks.locals <- trks.locals %>% dplyr::filter( jday < 178 )

# Now for 30min intervals
trks.thin <- trks.all %>% dplyr::select( id, red ) %>% 
  unnest( cols = red ) %>% dplyr::filter( jday < 178 )
tail( trks.thin )

## remove all individuals that did not overwinter at NCA using #
# the plot to determine which ones to keep.
## Answer: Based on plots, keep 9,8,3,1,6,5
#
#############################################################################
# Saving relevant objects and data ---------------------------------
#save breeding season data (not thinned)
write_rds( trks.breed, "Data/trks.breed")

#save breeding season data (thinned)
write_rds( trks.thin, "Data/trks.thin" )

#save migration data (unthinned)
write_rds( trks.mig, "Data/trks.mig" )

#save workspace in case we need to make changes
save.image( "TracksWorkspace.RData" )

########## end of save #########################
############### END OF SCRIPT ########################################
