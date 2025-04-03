##################################################################
# Script developed by Jen Cruz and modified by Ashley Santiago to estimate 
# ranges using AKDE     
# For this script we rely on Fleming et al.(2015) Ecology 96(5):1182-1188#
# We use ctmm first, and then use amt                           #
# For instructions on how to use ctmm directly check out:       #
# https://cran.r-project.org/web/packages/ctmm/vignettes/variogram.html #
# https://cran.r-project.org/web/packages/ctmm/vignettes/akde.html #
###################################################################

################## prep workspace ###############################

# Clean your workspace to reset your R environment. #
rm( list = ls() )
#install.packages( "ctmm" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library(dplyr)
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library(knitr)
library( amt )
library(ggplot2)
library( sf )
library(ctmm )#for more detailed functionality 
#####################################################################
## end of package load ###############
load("AKDEworkspace.RData")
###################################################################
#### Load or create data -----------------------------------------

#load cleaned data:
#download full data for breeding season monitoring of golden eagles at the Birds of Prey NCA
trks.pfdp <- read_rds( "trks.pfdp" )

#Look at data 
head(trks.pfdp)
unique(trks.pfdp$year)

#check class
class( trks.pfdp )
#check that the crs was correctly imported 
get_crs( trks.pfdp )

#order the ids numerically
trks.pfdp <- arrange(trks.pfdp, id)
###############################################################
##### Estimate ranges using AKDE continuous-time movement model:#
################################################################

# We start by plotting points for each individual:
# This is simply to view our data and make sure it looks okay
view(trks.pfdp %>% filter(id == 7))

trks.pfdp %>%
  group_by( id ) %>% 
  as_sf_points() %>% 
  summarise( do_union = FALSE ) %>%
  st_cast("LINESTRING" ) %>% 
  #uncommment the line below if you want to view one individual
  filter( id == 13 ) %>% 
  ggplot(., aes( color = as.factor(id)
  ) ) +
  theme_bw( base_size = 15 ) + 
  geom_sf() 

####################################################

####################################################
##### ALL individuals using ctmm     ###############
###################################################

# 1) Plot variograms for all individuals

# extract names for individuals first into an object
ids <- ( trks.pfdp$id )
unique(ids)

### RUN variograms  

#create objects to store results
svf.t <- list() #storing the variogram data
ctmm.t <- list()
xlimz <- c(0,36 %#% "hour" )
#set plot parameters
par( mfrow = c(3,3))
#loop through all individuals 
for( i in unique(ids) ){
  #print progress
  print( i )
  # extract data for individual i
  t <- trks.pfdp %>% filter( id == i )
  #convert to ctmm object and add to list
  ctmm.t[[i]] <- as_telemetry( t )
  #Calculate empirical variograms:
  svf.t[[i]] <- variogram( ctmm.t[[i]] )
  #plot variograms for each individual
  plot( svf.t[[i]], xlim =  xlimz )
}

# Note how they are unique for each individual. This supports 
# the authors suggestions to estimate unique movement models
# for each individual

#create and object to store results 
m.best <- list()

# Had issues with the for loop if I don't run this code after creating the list. Got an error about null values, but this fixed it.
for( i in unique(ids)) {
  m.best[[i]] <- 1
}

# use AIC to pick a best model from the model choices #
# we also plot the empirical variograms vs the model results #
# loop through each individual 
for( i in unique(ids)){  # can't use 1:length(ids) for b/c we no longer have all individuals
  print( i )
  #use empirical variogram estimated in the previous step 
  guess <- ctmm.guess(data = ctmm.t[[i]], variogram = svf.t[[i]],
                      interactive = FALSE )
  #guess function is what we want to guess the bandwidth
  # and compare fit using AIC to select the top model
  m.best[[i]] <- ctmm.select( ctmm.t[[i]], guess, verbose = TRUE,
                              trace = 2 )
  #select here in the ctmm.select is comparing the models
  #view summary output for model comparison for each individual
  print(summary( m.best[[i]] ))
}

head(m.best)
# OUF is multiple central foraging areas
# use individual names to replace those in the list: 

akde.pfdp <- list()

for(i in unique(ids)) {
  akde.pfdp[[i]] <- 1
}

# We now loop through all individuals to estimate ranges:
# Note we don't specify isopleths until plotting

for( i in unique(ids) ){
  print(i)
  an <- rownames(summary( m.best[[i]][1]))
ctmm::plot( svf.t[[i]], m.best[[i]][[1]], 
            xlim =  xlimz ,
            main = paste( ids[i], an ) )
}
#check ranges
for( i in unique(ids) ){
  print(i)
  print(summary(akde.pfdp[[i]]), level.UD= 0.95)
}

par(mfrow = c(3,2))

##################################################################################
###Trouble shooting error

akde.pfdp[[i]] <- akde(data, ctmm.t[[i]])


akde.test <- akde(data, ctmm.t[[1]])
plot(akde.test)


for (i in unique(ids)) {
  print(i)
  plot(ctmm.t[[i]], UD = akde.pfdp[[i]], level.UD = 0.95, col = "blue")
  plot(ctmm.t[[i]], UD = akde.pfdp[[i]], level.UD = 0.50, add = TRUE, col = "red")
  title("Unweighted best model")
}
##################################################################################



##### Turning Home ranges into Polygons ##########
#extract mean HR estimates as sf polygon and combine 
pfdp.50_sf <- list()
pfdp.95_sf <- list()


for( i in unique(ids) ){
  #extract home range for each animal and turn into sf object
  sf.pfdp.50 <- as.sf( akde.pfdp[[i]], level.UD=0.50, level=0.50 )
  sf.pfdp.95 <- as.sf( akde.pfdp[[i]], level.UD=0.95, level=0.95 )
  sf.pfdp.t <- st_transform( sf.pfdp.50, crs = get_crs( trks.pfdp ) ) 
  sf.pfdp.t2 <- st_transform( sf.pfdp.95, crs = get_crs( trks.pfdp ) ) 
  pfdp.50_sf[[i]] <- sf.pfdp.t[2,]
  pfdp.95_sf[[i]] <- sf.pfdp.t2[2,]
}

pfdp.50_poly <- pfdp.50_sf %>%  dplyr::bind_rows()
pfdp.95_poly <- pfdp.95_sf %>%  dplyr::bind_rows()


### re-add attributes for each individual
#create dataframe of just the columns you need
iddf <- trks.pfdp %>% 
  group_by( id ) %>% 
  select( id, unit, territory, year ) %>% 
  slice(1 )
#view
iddf

iddf$unit <- as.factor(iddf$unit)
table(iddf$unit)
#add attributes as columns to the polygon
pfdp.50_poly$territory <- iddf$territory
pfdp.50_poly$id <- iddf$id
#pfdp.50_poly$sex <- iddf$sex
pfdp.50_poly$year <- iddf$year
pfdp.95_poly$territory <- iddf$territory
pfdp.95_poly$id <- iddf$id
#pfdp.95_poly$sex <- iddf$sex
pfdp.95_poly$year <- iddf$year
pfdp.50_poly$unit <- iddf$unit
pfdp.95_poly$unit <- iddf$unit



#import polygon of the NCA and other layers as sf spatial file:
NCA <- sf::st_read( "~/EEB 627/MoveEco//Data/BOPNCA_Boundary.shp")

#Plot range comparisons from the thined data for 50% and 95% isopleths
ggplot() +
  theme_bw( base_size = 15 ) +
  geom_sf( data = pfdp.50_poly,
           aes(fill = unit, color = unit, alpha = 0.7), linewidth = 0.8) +
  geom_sf( data = pfdp.95_poly, 
           fill= NA, aes(color = unit), linewidth = 0.8) +
  geom_sf( data = NCA, 
           fill = NA, col = "black", size = 1 )

# Facet wrap by year
ggplot() +
  theme_bw(base_size = 15) +
  geom_sf(data = pfdp.50_poly, 
          aes(fill = unit, color = unit, alpha = 0.7), linewidth = 0.8) +
  geom_sf(data = pfdp.95_poly, 
          fill = NA, aes(color = unit), linewidth = 0.8) +
  geom_sf(data = NCA, 
          fill = NA, col = "black", size = 1) +
  facet_wrap(~year) +
  labs(title = "Home Ranges by Year")

# Facet wrap by territory
ggplot() +
  theme_bw(base_size = 15) +
  geom_sf(data = pfdp.50_poly, 
          aes(fill = factor(year), color = factor(year), alpha = 0.7), linewidth = 0.8) +
  geom_sf(data = pfdp.95_poly, 
          fill = NA, aes(color = factor(year)), linewidth = 0.8) +
  geom_sf(data = NCA, 
          fill = NA, col = "black", size = 1) +
  facet_wrap(~unit) +
  labs(title = "Home Ranges by Territory")

####
### Save desired results #
# Create a directory to store the individual polygons
dir.create("individual_polygons", showWarnings = FALSE)

# Loop through each territory and year
for (territory in unique(pfdp.50_poly$territory)) {
  # Remove spaces and hyphens from territory name
  clean_territory <- gsub("[[:space:]-]", "", territory)
  
  for (year in unique(pfdp.50_poly$year)) {
    # Extract polygons for 50% and 95% isopleths
    poly_50 <- pfdp.50_poly[pfdp.50_poly$territory == territory & pfdp.50_poly$year == year, ]
    poly_95 <- pfdp.95_poly[pfdp.95_poly$territory == territory & pfdp.95_poly$year == year, ]
    
    # Save polygons for 50% and 95% isopleths with cleaned territory name
    sf::st_write(poly_50, paste0("individual_polygons/", clean_territory, "_", year, "_50.shp"))
    sf::st_write(poly_95, paste0("individual_polygons/", clean_territory, "_", year, "_95.shp"))
  }
}

##############################################################
# Estimating AKDE using atm package                          #
##############################################################
################
# We use amt package (talks to ctmm) to estimate AKDE #

##### What if you want to use atm for all individuals? #######
##### don't try this in class. It will take a long time #

#make sure individuals are in order so that they can be compared to ctmm results
nested.pfdp <- trks.pfdp %>% 
  arrange( id ) %>% 
  nest( data = -"id" ) #nest tibbles

#calculate home range using movement model you choose:
akde_all <- nested.pfdp %>% 
  mutate( hr_akde_all = map( data, ~hr_akde( ., 
                                             model = fit_ctmm(., "ou" ),
                                             levels = 0.95 ) ) )
# 
akde_all

amt::hr_area( akde_all$hr_akde_all[[1]] ) 

#Plot for all individuals against equivalent from ctmm
for( i in 1:length(ids) ){
  a <- ggplot() +
    theme_bw( base_size = 15 ) +
    #extract isopleths for ouf model using thinned data from amt
    geom_sf( data = hr_isopleths( akde_all$hr_akde_all[[i]] ),
             col = "black", size = 3 ) +
    # geom_sf( data = as.sf( akde.uw[[i]] ), fill = NA, 
    #          col = "purple", size = 3) +
    labs( title = ids[i] )
  print( a )
}
# Why did it not run all of the points? ########################################################################


# Save desired results
write_rds( akde_all, "akde.all" )

###################################################################################################
# Get area of home ranges. These steps are probably unnecessary. You can just look at the summary for akde.breed rather than extracting the area from the polygons. 

# Calculate the area in square kilometers for 50%
area_50 <- st_area(pfdp.50_poly) / 1e6  # Convert square meters to square kilometers

# Create a data frame with the associated territory and year for 50%
area_df <- st_set_geometry(pfdp.50_poly, NULL) %>%
  mutate(area_50 = round(area_50, 2))  # Round to 2 decimal places

# Remove the "[m^2]" suffix from the area_50 column
area_df$area_50 <- gsub("\\s*\\[m\\^2\\]", "", area_df$area_50)

# Print or use the updated data frame as needed
print(area_df)

# Calculate the area in square kilometers for 95%
area_95 <- st_area(pfdp.95_poly) / 1e6  # Convert square meters to square kilometers

# Create a data frame with the associated territory and year for 95%
area_95 <- st_set_geometry(pfdp.95_poly, NULL) %>%
  mutate(area_95 = round(area_95, 2))  # Round to 2 decimal places

# Remove the "[m^2]" suffix from the area_95 column
area_95$area_95 <- gsub("\\s*\\[m\\^2\\]", "", area_95$area_95)

# Print or use the data frame as needed
print(area_95)

# Bind the new column to the existing data frame
area_df <- left_join(area_df, area_95, by = c("id", "territory", "year"))

# Print or use the updated data frame as needed
print(area_df)

# Remove redundant columns
area_df <- area_df %>%
  select(-starts_with("name."))

# Print or use the updated data frame as needed
print(area_df)

#br_status <- c( "breeding", "unknown", "breeding", "breeding", "breeding", "breeding", "breeding", "non-breeding", "non-breeding", "breeding", "breeding")

#id <- 1:11

#add_stat <- data.frame(id = id, br_status = br_status)

#area_df <- left_join(area_df, add_stat, by = "id")

area_df$area_50 <- as.numeric(area_df$area_50)
area_df$area_95 <- as.numeric(area_df$area_95)

# Calculate average area_50 and area_95 for all individuals
overall_average <- area_df %>%
  summarise(avg_area_50_all = mean(area_50),
            avg_area_95_all = mean(area_95))

# Calculate average area_50 and area_95 by individual
individual_average <- area_df %>%
  group_by(id) %>%
  summarise(avg_area_50 = mean(area_50),
  avg_area_95 = mean(area_95))

# Round the values to two decimal places
individual_average$avg_area_50 <- format(round(individual_average$avg_area_50, 2), nsmall = 2)
individual_average$avg_area_95 <- format(round(individual_average$avg_area_95, 2), nsmall = 2)


print(overall_average)


#print(individual_average)


mean_home_ranges <- area_df %>%
  group_by(territory) %>%
  summarise(mean_area_50 = mean(area_50),
            mean_area_95 = mean(area_95))

# Round the values to two decimal places
mean_home_ranges$mean_area_50 <- format(round(mean_home_ranges$mean_area_50, 2), nsmall = 2)
mean_home_ranges$mean_area_95 <- format(round(mean_home_ranges$mean_area_95, 2), nsmall = 2)

# Print the results
print(mean_home_ranges)


save.image(file = "AKDEworkspace.RData")

########################### END OF SCRIPT ######################################
