# Calculate the wind speed 
# Extract tp, longitude, latitude and time with units
# Saved it to a new nc file

# load packages
library(ncdf4)
library(tidyverse)
library(dplyr)


## function to calculate wind speed ----
calculate_wind_speed = function(u, v, 
                                scale_u, scale_v,
                                offset_u, offset_v) {
  if (!all(dim(u) == dim(v))) {
    stop("The dimensions of u and v must match")
  }
  u = u #* scale_u + offset_u
  v = v #* scale_v + offset_v
  wind_speed = sqrt(u^2 + v^2)
  return(wind_speed)
}

## extract wind speed and tp ----
extract_ws_tp = function(nc_name){
  # file path
  path = paste0("data_raw/", nc_name)
  new_path = paste0("data_extract/", nc_name)
  
  # open file
  nc = nc_open(path)
  
  # Extract dimensions
  time = ncvar_get(nc, "time")
  longitude = ncvar_get(nc, "longitude")
  latitude = ncvar_get(nc, "latitude")
  
  # extract variables
  u10 = ncvar_get(nc, "u10")
  v10 = ncvar_get(nc, "v10")
  tp = ncvar_get(nc, "tp")
  scale_u10 = ncatt_get(nc, 'u10', 'scale_factor')$value
  offset_u10 = ncatt_get(nc, 'u10', 'add_offset')$value
  scale_v10 = ncatt_get(nc, 'v10', 'scale_factor')$value
  offset_v10 = ncatt_get(nc, 'v10', 'add_offset')$value
  
  # find the wind speed
  ws = calculate_wind_speed(u = u10, v = v10, 
                            scale_u = scale_u10, scale_v = scale_v10,
                            offset_u = offset_u10, offset_v = offset_v10)
  
  # Extract dimensions from the original file
  lon_dim = ncdim_def("longitude", "degrees_east", longitude)
  lat_dim = ncdim_def("latitude", "degrees_north", latitude)
  time_units = ncatt_get(nc, "time", "units")$value
  time_dim = ncdim_def("time", time_units, time)
  
  # Extract units for v10 and tp from the original file
  v10_units = ncatt_get(nc, "v10", "units")$value
  tp_units = ncatt_get(nc, "tp", "units")$value
  
  # Define new variables ws and tp with the same dimensions
  ws_var = ncvar_def("ws", v10_units, list(lon_dim, lat_dim, time_dim), missval = -9999)
  tp_var = ncvar_def("tp", tp_units, list(lon_dim, lat_dim, time_dim), missval = -9999)
  
  # Create the new NetCDF file and add the new variables
  nc_new = nc_create(new_path, list(ws_var, tp_var))
  
  # Write data to the new file
  ncvar_put(nc_new, ws_var, ws)
  ncvar_put(nc_new, tp_var, tp)
  
  # Close the NetCDF files
  nc_close(nc)
  nc_close(nc_new)
  
}

# List all .nc files in the directory
nc_files = basename(list.files("data_raw", pattern = "\\.nc$", full.names = TRUE))

for (name in nc_files) {
  extract_ws_tp(name)
}

extract_ws_tp(nc_name = "IAN_2022_10.nc")

