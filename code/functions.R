## collection of functions to determine PDFs
# CREDITS: T. Kvalnes (NINA) 2024



# general PDF function for each raster cell (using matrix input)
pdf <- function(x, aorg, alost,z){
  (x*(1-(aorg-alost)/aorg)^z)/sum(values(x),na.rm = T)
}

## area lost for ES disturbances
#q min and q max can be either set by the user or determined from a questionnaire
A_lost_D<-function(tmpes,iar_rast,q_min,q_max, store_values, resolution){
  # takes the IAR raster, establishes a cumulative distribution function based on the distances between cells and es hot spots, 
  # based on the distribution and the given inputs (e.g. survey) of q_min and q_max (0-1) calculates a disturbance decay function
  # finally calc A lost based on decay
  
  #reproject and resample cv data
  iar_rast<-rast(iar_rast)
  iar_rast<-iar_rast%>%terra::project(crs("+init=epsg:3035"))
  bbox<-st_bbox(iar_rast)
  r1<-rast(res = resolution, xmin = bbox[1], xmax = bbox[3], ymin = bbox[2], ymax = bbox[4],crs="+init=epsg:4326")


  # cv_rast<-terra::project(cv_rast,st_crs(stud_site)$wkt)
  #resampling to standard raster
  iar_rast <- resample(iar_rast, r1)
  q_95<-quantile(values(iar_rast),0.95,na.rm=T)
  require(terra)
  require(stars)
  # transform rast to poly
  iar_rast<-as.polygons(iar_rast,aggregate = F)
  iar_rast<-st_as_sf(iar_rast)
  colnames(iar_rast)<-c("probability_stdDev","geometry")
  st_geometry(iar_rast)<-iar_rast$geometry
  
  ## here consensus should be relative to the ES. higher than Q95 -- hotspot
  
  # CV_poly_hot<-CV_poly%>%filter(probability_stdDev<median(CV_poly$probability_stdDev))
  IAR_poly_hot<-iar_rast%>%filter(probability_stdDev>q_95)
  # provide some space for stats for each cell
 
  
  ## inverse cumulative dist function
  inv_ecdf <- function(f){
    x <- environment(f)$x
    y <- environment(f)$y
    approxfun(y, x)
  }
  

  
  min_dist_list<-list()
  #### just use the min dist of each cell and make then the distribution
  for(i in 1: nrow(iar_rast)){
    print(paste0(round(i/nrow(iar_rast)*100,2)," %"))
    # print(i)
    a<-st_distance(iar_rast[i,],IAR_poly_hot, which = "Euclidean")
    a<-as.numeric(a)
    min_dist<-min(a)

    if(length(unique(a))<=1){
      next
    }
    #print(min_dist)
    min_dist_list[i]<-min_dist

    
  }
  
  a<-unlist(min_dist_list)

  b<-ecdf(a)
  #   ## inverse ecdf
  g <- inv_ecdf(b)
  

  if(q_min == 0){
    min_d<-0
    
    
  }else{
    min_d<-g(q_min)
    
  }
  
  max_d<-g(q_max)
  mean_d<-g((q_max-q_min)/2)
  
  ##2.2 area specific DK
  
  alpha<-0.1
  beta<-log((2-alpha)/alpha)/(min_d-mean_d)
  
  # function to integrate from d0 to dmax

  
  f <- function(x) {
    (1 - 1 / (1 + exp(beta * (x - max_d)))) / max_d
  }
  
  ## isn`t that min_d to max_d?
  #dK<-integrate(dk_fkt,0,max_d)
  dk<-integrate(f,0,max_d)
  
  ##2.3. area lost due to disturbance
  #A_lost_D<-pi*(max_d*dK$value)^2/10^6
  A_lost_D<-pi*(max_d*dk$value)^2/10^6

  #### write table to BQ just in case
  if(store_values == T){
    d_val<-as.data.frame(cbind(tmpes,min_d,max_d, dk$value, A_lost_D, resolution))
    colnames(d_val)<-c("esID","min_crit_d_m","max_crit_d_m","decay","A_lost_disturbance_km2","analysis_resolution_m")
    write.csv(d_val,here::here("output", paste0(tmpes,"_distPARAM_new.csv")))
  }

  
  return(A_lost_D)
  
}


## simpler calculation for area lost disturbance. Take values from literature distances in meters!
A_lost_D_simple<-function(min_d, max_d){
  mean_d<-mean(min_d,max_d)
  dk_fkt <- function(x) {    
    alpha<-0.1
    beta<-log((2-alpha)/alpha)/(min_d-mean_d)
    out <- (1-(1/(1+exp(beta*(x-mean_d)))))/max_d
  }
  
  ## isn`t that min_d to max_d?
  dK<-integrate(dk_fkt,min_d,max_d)
  
  ##2.3. area lost due to disturbance
  A_lost_D<-pi*(max_d*dK$value)^2/10^6
  return(A_lost_D)
}


# Given a set of turbine sites and an area lost per turbine - create the polygon for the wind farm by the union of turbine polygons
wind_farm_polygon <- function(turbine_sites, area_lost){
  # turbine_sites = coordinates for each wind turbine (points, sf object)
  # area_lost = area lost per turbine (a number, km^2)
  
  # Load required packages
  require(sf)
  
  # Buffer each turbine (meters). Assume that a circular area around each turbine is lost
  turbines_area_lost <- st_buffer(x = turbine_sites, dist = sqrt(area_lost/pi)*1000)
  
  # Wind farm polygon
  wind_farm_polygon <- st_union(turbines_area_lost)
  
  # Return
  wind_farm_polygon
  
}

# LCA: Normalize function ----
# Given a raster map, normalize the values in the first attribute (0-1 range)
# Normalize function
nor <- function(map, na.rm = TRUE){
  
  # Load required packages
  require(sf)
  
  # Normalize
  map <- (map-min(map[[1]], na.rm=na.rm))/(max(map[[1]], na.rm=na.rm) - min(map[[1]], na.rm=na.rm))
  
  # Return
  map
}


# LCA: Sliding window PDF (all possible locations of a given wind farm) ---- Reto: is the input of the sliding window function the normalized maps or is it normalized after?
pdf_sliding_window <- function(map, wind_farm_polygon, z, return_extremas = FALSE){
  
  # The function will find the pdf for each raster cell using a sliding window approach. The polygon for the wind farm is centered on the center of a raster cell in the map and the PDF calculated and summed for that cell, then the procedure is repeated for all cells in the map.
  
  # RETURN: A map with PDF values (default) or the extremas and central tendency (return_extremas = TRUE)
  
  # Load required packages
  require(terra)
  require(stars)
  
  # Calculate the mapsum
  mapsum <- sum(map[[1]], na.rm=TRUE)
  
  # Move wind farm polygon centroid to the center of a raster cell (to have the correct cell cover for the sliding window later)
  # Extract coordinates for the center of raster cells in the map
  coo <- st_coordinates(map)
  # Extract value of map for each coordinate and make coo a sf
  coo <- st_extract(x = map, at = st_as_sf(coo, coords = c("x", "y"), crs = st_crs(map)))
  names(coo)[1] <- "values"
  # Remove center coordinates for cells with missing values
  coo <- filter(coo, !is.na(values))
  # Choose coordinates in the center of the map
  coo <- coo[round(dim(coo)[1]/2),]
  # Extract only geometry
  coo <- st_geometry(coo)
  # Extract centroid for wind farm polygon
  wf_centroid <- st_geometry(st_centroid(wind_farm_polygon))
  # Move polygon for wind farm to the center of a raster cell given by "coo"
  wind_farm_polygon <- st_as_sf(st_geometry(wind_farm_polygon) - wf_centroid + coo, crs = st_crs(map))
  
  # Change to terra
  map <- as(map, "SpatRaster")
  wf <- as(wind_farm_polygon, "SpatVector")
  
  # Proportion of each raster cell covered by wind farm/turbine
  cellcover <- terra::rasterize(x = wind_farm_polygon, y = map, touches = T, cover=T)
  # Trim outer NA values
  cellcover <- terra::trim(cellcover)
  # Extract cover values as matrix
  cellcover <- as.matrix(cellcover, wide = TRUE)
  
  # Check for missing values and set these to 0.001
  cellcover[is.na(cellcover)] <- 0.001
  
  # Sliding window matrix ## 1 of covered 0 if not covered
  swmat <- ceiling(cellcover)
  
  # Check for length 1
  if(dim(swmat)[2] == 1) {
    swmat <- cbind(NA, swmat, NA)
  }
  if(dim(swmat)[1] == 1) {
    swmat <- rbind(NA, swmat, NA)
  }
  
  # Check for even number (need odd number)
  if((dim(swmat)[2] %% 2) == 0) {
    swmat <- cbind(swmat, NA)
  }
  if((dim(swmat)[1] %% 2) == 0) {
    swmat <- rbind(swmat, NA)
  }
  
  # PDF function (edge effect? - crop frame?)
  .pdf <- function(x, cc = c(cellcover), zz = z, ms = mapsum){
    sum((x-x*(1-cc)^zz)/ms,na.rm=T)
  }
  
  # Estimate the sliding window pdf
  sw_pdf <- terra::focal(x = map, w = swmat, fun = .pdf)
  # sw_pdf <- terra::focal(x = map, w = swmat, fun = mean, na.rm=T)
  # sw_pdf <- terra::focal(x = map, w = cellcover, fun = .pdf)
  
  
  # Change to stars
  sw_pdf <- st_as_stars(sw_pdf)
  names(sw_pdf) <- "pdf"
  
  # Return extremas
  if(return_extremas){
    # Calculate PDFs
    pdf_max <- max(sw_pdf[[1]], na.rm=TRUE)
    pdf_min <- min(sw_pdf[[1]], na.rm=TRUE)
    pdf_mean <- mean(sw_pdf[[1]], na.rm=TRUE)
    pdf_median <- median(sw_pdf[[1]], na.rm=TRUE)
    
    return(list(pdf_max = pdf_max, pdf_min = pdf_min, pdf_mean = pdf_mean, pdf_median = pdf_median))
  }
  
  # Return
  return(sw_pdf)
  
}

