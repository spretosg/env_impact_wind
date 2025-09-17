packages <- c(
  "terra",
  "dplyr",
  "readxl",
)

# Install missing packages
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org/")
}

# Load all packages
lapply(packages, library, character.only = TRUE)


### calculates the z coefficient (Invariability of ecosystem - area relationship) for each ecosystem service
### input: es_rating_participants.xlsx (just for es list), individual es maps (script 1_es_cap_cv)
### output: csv containing IAR and z for different areas and all ES (stored in data/4_z_coef)


es_vec<-read_excel("data/02_ecosystem_service_modelling/es_rating_participants.xlsx", sheet = "ID1",skip = 2)
es_vec<-es_vec %>%  filter(!row_number() %in% c(16:nrow(es_vec)))
es_vec<-es_vec[1:(length(es_vec)-5)]
lt_types<-unique(es_vec$Landskapstype)
es_vec<-colnames(es_vec)[-1]
es_vec<-es_vec[1:length(es_vec)-1]

z_list<-list()

for(i in 1: length(es_vec)){
  es_tmp<-es_vec[i]
  
  
  
  folder_path<-paste0("data/02_ecosystem_service_modelling/ind_maps/",es_tmp)
  raster_files <- list.files(here::here(folder_path), pattern = "\\.tif$", full.names = TRUE)
  
  es <- rast(raster_files)  
  
  ######## IAR0 = center cell of all cells
  
  center_row <- round(nrow(es)/2,0)+3
  center_col <- round(ncol(es)/2,0)-3
  # Get cell number of the center
  center_cell <- cellFromRowCol(es, center_row, center_col)
  
  # Get row/col of center
  cell_id <- rowColFromCell(es, center_cell)
  
  #cell_id <- cellFromRowCol(es, 28, 28)
  cent_vals <- values(es)[cell_id, ]
  mu <- mean(cent_vals, na.rm = TRUE)   # mean es for cent cell
  sigma <- sd(cent_vals,na.rm=TRUE)
  
  cv<-sigma/mu
  IAR_cen<-1/cv^2
  
  ### IAR1 3x3 window
  # Get surrounding 3x3 window cell IDs
  rows <- (center_row - 1):(center_row + 1)
  cols <- (center_col - 1):(center_col + 1)
  
  # Clip to raster bounds
  rows <- rows[rows >= 1 & rows <= nrow(es)]
  cols <- cols[cols >= 1 & cols <= ncol(es)]
  
  # Get all cell IDs in the 3x3 window
  cell_ids <- as.vector(outer(rows, cols, Vectorize(function(i, j) cellFromRowCol(es, i, j))))
  
  # Extract values for all 9 cells across all layers
  window_vals <- values(es)[cell_ids, ]
  #mean correl between all pixels in window
  p1<-mean(cor(t(window_vals)),na.rm=T)
  
  A1<-9
  I1 <- IAR_cen * (A1 / (((A1 - 1) * p1) + 1))
  
  
  
  ### IAR2 9x9 window
  # Get surrounding 3x3 window cell IDs
  rows <- (center_row - 4):(center_row + 4)
  cols <- (center_col - 4):(center_col + 4)
  
  # Clip to raster bounds
  rows <- rows[rows >= 1 & rows <= nrow(es)]
  cols <- cols[cols >= 1 & cols <= ncol(es)]
  
  # Get all cell IDs in the 3x3 window
  cell_ids <- as.vector(outer(rows, cols, Vectorize(function(i, j) cellFromRowCol(es, i, j))))
  
  # Extract values for all 9 cells across all layers
  window_vals <- values(es)[cell_ids, ]
  #mean correl between all pixels in window
  p2<-mean(cor(t(window_vals)),na.rm=T)
  
  A2<-81
  I2 <- I1 * (A2 / (((A2 - 1) * p2) + 1))
  
  
  ### IAR3 27x27 window
  # Get surrounding 3x3 window cell IDs
  rows <- (center_row - 13):(center_row + 13)
  cols <- (center_col - 13):(center_col + 13)
  
  # Clip to raster bounds
  rows <- rows[rows >= 1 & rows <= nrow(es)]
  cols <- cols[cols >= 1 & cols <= ncol(es)]
  
  # Get all cell IDs in the 3x3 window
  cell_ids <- as.vector(outer(rows, cols, Vectorize(function(i, j) cellFromRowCol(es, i, j))))
  
  # Extract values for all 9 cells across all layers
  window_vals <- values(es)[cell_ids, ]
  #mean correl between all pixels in window
  p3<-mean(cor(t(window_vals)),na.rm=T)
  
  A3<-729
  I3 <- I2 * (A3 / (((A3 - 1) * p3) + 1))
  
  
  

  ### IAR4 120*120 window
  # Get surrounding 3x3 window cell IDs
  rows <- (center_row - 60):(center_row + 60)
  cols <- (center_col - 60):(center_col + 60)
  
  # Clip to raster bounds
  rows <- rows[rows >= 1 & rows <= nrow(es)]
  cols <- cols[cols >= 1 & cols <= ncol(es)]
  
  # Get all cell IDs in the 3x3 window
  cell_ids <- as.vector(outer(rows, cols, Vectorize(function(i, j) cellFromRowCol(es, i, j))))
  
  # Extract values for all 9 cells across all layers
  window_vals <- values(es)[cell_ids, ]
  #mean correl between all pixels in window
  p4<-mean(cor(t(window_vals)),na.rm=T)
  
  A4<-14641
  I4 <- I3 * (A4 / (((A4 - 1) * p4) + 1))
  
  df<-data.frame(IAR = c(IAR_cen,I1,I2,I3,I4), area = c(3,9,81,729,14641), correl = c(NA,p1,p2,p3,p4))
  df$logI<-log(df$IAR)
  df$logA<-log(df$area)
  df$es<-c(es_tmp,es_tmp,es_tmp,es_tmp,es_tmp)
  
  lm_model <- lm(logI ~ logA, data = df)

  z<-lm_model$coefficients[2]
  
  df$z<-c(z,z,z,z,z)
  z_list[[i]]<-df
}


####For the sensitivity analysis, a random selection of n centre cells can performed.
set.seed(123) # For reproducibility


for (i in 1:length(es_vec)) {
  es_tmp <- es_vec[i]
  
  folder_path <- paste0("data/02_ecosystem_service_modelling/ind_maps/", es_tmp)
  raster_files <- list.files(here::here(folder_path), pattern = "\\.tif$", full.names = TRUE)
  es <- rast(raster_files)
  
  # Get all valid cells (non-NA across all layers)
  n_rows <- nrow(es)
  n_cols <- ncol(es)
  
  # Define central 10x10 window (centered in the middle of raster)
  center_row <- round(n_rows / 2)
  center_col <- round(n_cols / 2)
  
  rows_10 <- (center_row - 4):(center_row + 5)  # 10 rows
  cols_10 <- (center_col - 4):(center_col + 5)  # 10 columns
  
  # Ensure within bounds
  rows_10 <- rows_10[rows_10 >= 1 & rows_10 <= n_rows]
  cols_10 <- cols_10[cols_10 >= 1 & cols_10 <= n_cols]
  
  # Get cell indices of 10x10 grid
  central_cells <- as.vector(outer(rows_10, cols_10, Vectorize(function(i, j) cellFromRowCol(es, i, j))))
  
  # Filter to only valid (non-NA) cells
  valid_central_cells <- central_cells[which(rowSums(is.na(values(es)[central_cells, ])) == 0)]
  
  # Sample up to 10 from valid central cells
  if (length(valid_central_cells) < 40) {
    warning("Fewer than 10 valid central cells found.")
    sampled_cells <- valid_central_cells
  } else {
    sampled_cells <- sample(valid_central_cells, size = 40)
  }
  run_list <- list() 
  for (j in seq_along(sampled_cells)) {
    center_cell <- sampled_cells[j]
    cell_id <- rowColFromCell(es, center_cell)
    center_row <- cell_id[1]
    center_col <- cell_id[2]
    
    # ---- IAR0: Single center cell ----
    cent_vals <- values(es)[center_cell, ]
    mu <- mean(cent_vals, na.rm = TRUE)
    sigma <- sd(cent_vals, na.rm = TRUE)
    cv <- sigma / mu
    IAR_cen <- 1 / cv^2
    
    # Function to calculate IAR for a given window size
    compute_IAR <- function(es, center_row, center_col, window_radius, prev_IAR) {
      rows <- (center_row - window_radius):(center_row + window_radius)
      cols <- (center_col - window_radius):(center_col + window_radius)
      
      # Clip to raster boundaries
      rows <- rows[rows >= 1 & rows <= nrow(es)]
      cols <- cols[cols >= 1 & cols <= ncol(es)]
      
      cell_ids <- as.vector(outer(rows, cols, Vectorize(function(i, j) cellFromRowCol(es, i, j))))
      window_vals <- values(es)[cell_ids, ]
      
      A <- length(cell_ids)
      p <- mean(cor(t(window_vals)), na.rm = TRUE)
      
      IAR <- prev_IAR * (A / (((A - 1) * p) + 1))
      list(IAR = IAR, A = A, p = p)
    }
    
    # Calculate IARs for different window sizes
    I1_res <- compute_IAR(es, center_row, center_col, window_radius = 1, prev_IAR = IAR_cen)
    I2_res <- compute_IAR(es, center_row, center_col, window_radius = 4, prev_IAR = I1_res$IAR)
    I3_res <- compute_IAR(es, center_row, center_col, window_radius = 13, prev_IAR = I2_res$IAR)
    I4_res <- compute_IAR(es, center_row, center_col, window_radius = 60, prev_IAR = I3_res$IAR)
    
    df <- data.frame(
      IAR = c(IAR_cen, I1_res$IAR, I2_res$IAR, I3_res$IAR, I4_res$IAR),
      area = c(1, I1_res$A, I2_res$A, I3_res$A, I4_res$A),
      correl = c(NA, I1_res$p, I2_res$p, I3_res$p, I4_res$p),
      center_cell = center_cell,
      run = j
    )
    
    df$logI <- log(df$IAR)
    df$logA <- log(df$area)
    df$es <- es_tmp
    
    if (all(is.finite(df$logI)) && all(is.finite(df$logA)) && nrow(df) >= 2) {
      lm_model <- lm(logI ~ logA, data = df)
      z <- coef(lm_model)[2]
      df$z <- z
    } else {
      warning(paste("Skipping lm() for es =", es_tmp, "run =", j, ": insufficient data"))
      df$z <- NA
    }
    
    run_list[[j]] <- df 
  }
  df_all_runs <- do.call(rbind, run_list)
  z_list[[i]] <- df_all_runs 
}

# Combine results
results <- do.call(rbind, z_list)



# df_z <- bind_rows(z_list)
write.csv(results,"data/02_ecosystem_service_modelling/z_coef_sens.csv")
