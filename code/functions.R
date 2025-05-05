calculate_sex_specific_val <- function(pop_total, pop_se, sex_ratio, sex_ratio_se) {
  # Calculate population estimates
  pop_male <- pop_total * sex_ratio
  pop_female <- pop_total * (1 - sex_ratio)
  
  # Calculate standard errors using error propagation formula for products
  pop_male_rel_se <- sqrt((pop_se/pop_total)^2 + (sex_ratio_se/sex_ratio)^2)
  pop_male_se <- pop_male * pop_male_rel_se
  
  # For females, need to handle error propagation for (1-r)
  # The SE of (1-r) equals SE(r)
  # But the relative SE needs to use (1-r) in denominator
  pop_female_rel_se <- sqrt((pop_se/pop_total)^2 + (sex_ratio_se/(1-sex_ratio))^2)
  pop_female_se <- pop_female * pop_female_rel_se
  
  # Return results
  return(list(
    male = list(estimate = pop_male, se = pop_male_se),
    female = list(estimate = pop_female, se = pop_female_se)
  ))
}

update_sex_values <- function(df) {
  # Get unique IDs
  unique_ids <- unique(df$ID)
  
  # For each unique ID
  for (id in unique_ids) {
    # Get rows with this ID
    idx <- which(df$ID == id)
    
    # Get unique non-NA Sex values for this ID
    unique_sex <- unique(df$Sex[idx][!is.na(df$Sex[idx])])
    
    # Case 3: If there's exactly one non-NA Sex value and at least one NA
    if (length(unique_sex) == 1 && any(is.na(df$Sex[idx]))) {
      # Replace NAs with the non-NA value
      df$Sex[idx][is.na(df$Sex[idx])] <- unique_sex
    }
    # Cases 1 and 2 require no action
  }
  
  return(df)
}

convert_df_to_raster <- function(df, resolution = NULL, crs = NULL) {
  # Create a SpatVector from the data frame points
  points <- vect(df, geom = c("x", "y"), crs = crs)
  
  # Extract the extent
  ext <- ext(points)
  
  # If resolution is not provided, estimate it from the data
  if (is.null(resolution)) {
    # Estimate resolution based on average distance between points
    # or use a simple heuristic based on extent
    x_range <- ext[2] - ext[1]
    y_range <- ext[4] - ext[3]
    n_points <- nrow(df)
    resolution <- min(x_range, y_range) / sqrt(n_points) * 2  # Just an estimation
  }
  
  # Create an empty raster with the given resolution
  r <- rast(ext = ext, resolution = resolution, crs = crs)
  
  # Rasterize the points, using the D.0 values
  raster <- rasterize(points, r, field = "D.0", fun = mean)
  
  return(raster)
}