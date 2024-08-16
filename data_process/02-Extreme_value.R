# Extract the data and merge the data into one csv file 
# (if the data is separated into two files)

# load packages
library(ncdf4)
library(tidyverse)
library(dplyr)
library(lubridate)
library(stringr)
library(abind)

## function to extract the main variables ----- 
extract_variable = function(nc_name, two_month = FALSE){
  
  if (two_month) {
    path = paste0("combine/", nc_name)
  } else {
    path = paste0("data/", nc_name)
  }
  
  # open file
  nc = nc_open(path)
  
  # Extract dimensions
  time = ncvar_get(nc, "time")
  longitude = as.vector(ncvar_get(nc, "longitude"))
  latitude = as.vector(ncvar_get(nc, "latitude"))
  
  # extract variables
  tp = ncvar_get(nc, "tp")
  ws = ncvar_get(nc, "ws10")
  
  # Close the netCDF file
  nc_close(nc)
  
  # Create a list to store all data 
  data_list = list("time" = time, 
                   "latitude" = latitude, "longitude" = longitude, 
                   "ws" = ws, "tp" = tp)
  return(data_list)
}

## function to find the correlation matrix ----
replace_na_below_threshold = function(ws, tp, quantile_value) {
  
  ws_quantile = quantile(ws, probs = quantile_value, na.rm = TRUE)
  tp_quantile = quantile(tp, probs = quantile_value, na.rm = TRUE)
  
  # Loop through each location
  for (i in 1:dim(ws)[1]) {
    for (j in 1:dim(ws)[2]) {
      
      ws_series = ws[i, j, ]
      tp_series = tp[i, j, ]
      ws_time = rep(NA, length(ws_series))
      tp_time = rep(NA, length(tp_series))
      
      ws_index = ws_series >= ws_quantile
      tp_index = tp_series >= tp_quantile
      common_ind = as.integer(ws_index & tp_index)
      
      if (any(common_ind == 1)) {
        
        series_ind = which(common_ind == 1)
        change_ind = c(1, which(diff(series_ind) != 1) + 1, length(series_ind) + 1)
        max_diff = max(diff(change_ind))
        diff_ind = which(diff(change_ind) == max_diff)
        
        if (max_diff >= 5) {
          for (k in diff_ind) {
            first_ind = series_ind[change_ind[k]]
            last_ind = first_ind + max_diff - 1
            ws_time[first_ind:last_ind] = ws_series[first_ind:last_ind]
            tp_time[first_ind:last_ind] = tp_series[first_ind:last_ind]
          }
        }
      }
      ws[i, j, ] = ws_time
      tp[i, j, ] = tp_time
    } }
  
  return(list(ws = ws, tp = tp))
}

# Function to process each storm
process_storm = function(nc_name, base_time, quantile_value, two_month = FALSE) {
  if (two_month) {
    storm_first_half = extract_variable(nc_name[1])
    storm_second_half = extract_variable(nc_name[2])
    storm = list("time" = c(storm_first_half$time, storm_second_half$time),
                 "latitude" = storm_first_half$latitude,
                 "longitude" = storm_first_half$longitude,
                 "ws" = abind(storm_first_half$ws, storm_second_half$ws, along = 3),
                 "tp" = abind(storm_first_half$tp, storm_second_half$tp, along = 3))
  } else {
    storm = extract_variable(nc_name = nc_name)
  }
  
  storm_removed = replace_na_below_threshold(storm$ws, storm$tp, quantile_value)
  
  long_format = expand.grid(lon = storm$longitude,
                            lat = storm$latitude,
                            time = storm$time) %>%
    mutate(ws = as.vector(storm_removed$ws), 
           tp = as.vector(storm_removed$tp))
  
  storm_removed_na = na.omit(long_format)
  
  if (nrow(storm_removed_na) == 0) {
    return(NULL)
  }
  
  # Convert hours since 1900 to POSIXct datetime
  storm_removed_na = storm_removed_na %>%
    mutate(time_format = base_time + hours(time),
           time_format = as.POSIXct(time_format, tz = "UTC"),
           time_diff = as.numeric(difftime(time_format, min(time_format), units = "hours")),
           year = year(time_format))
  storm_removed_na$name = str_extract(nc_name[1], "^[^_]+_[^_]+")
  
  return(storm_removed_na)
}


