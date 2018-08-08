
rm(list = ls())
# install required packages (if not already installed)
list.of.packages <- c("tidyverse", "rgeos", "gstat", "leaflet", "sp", "rgdal")
new.packages     <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)

# loading required packages
lapply(list.of.packages, require, character.only = T)
rm(list.of.packages, new.packages)

# 1. Import your data files
my_data <- read.csv("Site level univariate1.csv", stringsAsFactors = F) 

# Correct for NA values (make sure the variable is numeric)
my_data <- my_data %>% mutate(Exploit20 = ifelse(is.na(Exploit20) | Exploit20 == 0, NA, log(as.numeric(Exploit20))),
                              ExploitTot = as.numeric(ExploitTot)) # Warning message is fine 

# 2. Make sure your data is in the right format, e.g. each row is the scale you want the points to be at

my_data <-
  my_data %>%
  group_by(SiteCode) %>%
  summarise(Exploit20   = mean(Exploit20, na.rm = T),
            ExploitTot = mean(ExploitTot, na.rm = T),
            SiteLat  = mean(SiteLat),
            SiteLong = mean(SiteLong))

my_var    <- "Exploit20"
local_crs <- "+proj=utm +zone=55 +south +ellps=GRS80 +units=m" # The CRS format R is working in?
this_crs  <- sp::CRS("+init=epsg:4326") # The CRS format your data in is
map_res   <- c(300, 300) # Map resolution ( number of pixels/grid cells vertically and horizontally)
my_buffer <- 40 * 1000 # The size of the buffers around the points
map_output_width <- 4
asp_ratio <- 1 

# To map a variable you need three things.
# 1. A shape file of the region you're interested in
# 2. A

# Created a little function that converts a dataframe of lat and longs to a UTM
to_UTM <- function(x, crs = "+proj=utm +zone=55 +south +ellps=GRS80 +units=m") {
  sp::coordinates(x) = ~ SiteLong + SiteLat
  sp::proj4string(x) <- CRS("+init=epsg:4326")
  x <- spTransform(x, CRS(crs))
}


# Creating the hull shape (from lat and long values)
hull_shape <- 
  my_data %>% 
  dplyr::select(SiteLat, SiteLong) %>% 
  slice(c(chull(.), chull(.)[1])) %>% 
  to_UTM() %>% 
  SpatialPointsDataFrame(data.frame(id = 1:length(.))) %>% 
  spTransform(., raster::crs(local_crs))

# Import the reef and land shape files
reef_shape <- readOGR(dsn = getwd(), "gbrReef") %>% spTransform(raster::crs(local_crs))
land_shape <- readOGR(dsn = getwd(), "mainland") %>% spTransform(raster::crs(local_crs))

# Creating a UTM from ourdataframe
my_data_utm <- 
  my_data %>%
  filter(!is.na(Exploit20)) %>% 
  to_UTM() 
  
my_data_utm_df  <- 
  my_data_utm %>% 
  data.frame() %>% 
  mutate(x = SiteLong, 
         y = SiteLat)


# Adding a buffer to the hull - i.e. specifying outer bounds of map
map_coords <- raster::extent(hull_shape) 
map_xlim <- c(map_coords[1] - my_buffer, map_coords[2] + my_buffer)
map_ylim <- c(map_coords[3] - my_buffer, map_coords[4] + my_buffer)

empty_raster <- 
  raster::raster(map_coords,
                 ncols = map_res[1], nrows = map_res[2], 
                 crs = local_crs) %>%  
  raster::extend(c(20, 20)) # Blank mask

points_sites <- rgeos::gBuffer(my_data_utm, width = my_buffer) # Location of points and size of cookie cutter
mask         <- raster::rasterize(points_sites, empty_raster) # Applying cookie cutter to blank mask, Note: order is important

raster_vals <- 
  gstat::gstat(id = my_var, 
               formula = Exploit20 ~ 1, 
               locations = ~ x + y,
               data = my_data_utm_df) %>% 
  raster::interpolate(empty_raster, .) # painted canvas

masked_raster_vals <- raster::mask(raster_vals, mask) # Applying mask over painted canvas

# sorting the colours of the buffer of points
mycolpal <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))

# Prevents printing to screen, outputs as pdf
pdf('map_output.pdf', width = map_output_width *asp_ratio*1.5, height = map_output_width*1.5)

par(mar = c(1,1,1,1))
# plot(raster_vals, legend = T, xaxt = "n", yaxt = "n") # for illustrative 
plot(masked_raster_vals, col = (mycolpal(100)), legend = T, xaxt = "n", yaxt = "n") # what you'r seeing is the z-axis
plot(land_shape,  add = T, col = "darkolivegreen", border = F)
plot(reef_shape,  add = T, col = "darkgrey", border = F)
plot(my_data_utm, add = T, pch = 20, cex = 0.5, col = "deeppink") # add the locations of the sites


# pch = changes type of point (google 'pch in r')
# cex = size of points

# Adding a scale bar
raster::scalebar(200*1000, xy = c(10702,7333221),
                 label = "200 km")

# If you want to add a box
zm1  <- raster::extent(845589, 821417+290000, 7300289, 7300289+280000)
czm1 <- raster::crop(my_data_utm, zm1)
pzm1 <- as(zm1, "SpatialPolygons")
plot(pzm1, add = T, border = "grey", lty = 2)

# Make sure you turn the device off 
dev.off()


