## analysis helper functions

nor <- function(map, na.rm = TRUE){
  
  # Load required packages
  require(sf)
  
  # Normalize
  map <- (map-min(values(map), na.rm=na.rm))/(max(values(map), na.rm=na.rm) - min(values(map), na.rm=na.rm))
  
  # Return
  map
}


layer_stats <- function(x) {
  data.frame(
    mean = mean(x, na.rm = TRUE),
    q5   = quantile(x, probs = 0.05, na.rm = TRUE),
    q25   = quantile(x, probs = 0.25, na.rm = TRUE),
    q75  = quantile(x, probs = 0.75, na.rm = TRUE),
    q95  = quantile(x, probs = 0.95, na.rm = TRUE)
  )
}

get_jenks_breaks <- function(r, n = 5, name = "Raster", group = NA, wp) {
  library(classInt)
  library(terra)
  
  # get raster values, removing NAs
  vals <- values(r, na.rm = TRUE)
  vals <- vals[is.finite(vals)]
  
  # mean of raster values extracted by wp shapefile
  wt_pdf <- terra::extract(r, wp, fun = layer_stats)
  wt_pdf <- unlist(t(wt_pdf)[-1])
  wt_mean <- wt_pdf[1]
  wtQ25  <- wt_pdf[3]
  wtQ75  <- wt_pdf[4]
  
  # skip constant rasters
  if (length(unique(vals)) < 2) {
    return(NULL)
  }
  
  # Jenks breaks
  jenks <- classIntervals(vals, n = n, style = "fisher")
  brks  <- jenks$brks
  
  # Create dataframe with group + ID
  res <- data.frame(
    group   = group,
    ID      = name,
    wt_mean = wt_mean,
    wtQ25   = wtQ25,
    wtQ75   = wtQ75,
    matrix(NA, nrow = 1, ncol = n)
  )
  colnames(res)[6:(5+n)] <- paste0("Class", seq_len(n))
  
  # Fill each class column with "min - max"
  for (i in seq_len(n)) {
    res[[paste0("Class", i)]] <- paste0(brks[i], " - ", brks[i+1])
  }
  
  return(res)
}