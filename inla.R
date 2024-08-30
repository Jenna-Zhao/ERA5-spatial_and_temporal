### Read data ----
Wilma = extract_variable("WILMA_2005_10.nc")
Wilma_2005 = replace_na_below_threshold(ws = Wilma$ws[, , 220:240], 
                                        tp = Wilma$tp[, , 220:240], 
                                        quantile_value = .80)

### Generate data frame -----
base_time = as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
long_format = expand.grid(lon = Wilma$longitude,
                          lat = Wilma$latitude,
                          time = Wilma$time[220:240]) %>%
  mutate(ws = as.vector(Wilma_2005$ws), 
         tp = as.vector(Wilma_2005$tp))

storm_removed_na = na.omit(long_format)

# Convert hours since 1900 to POSIXct datetime
storm_data = storm_removed_na %>%
  arrange(lon, desc(lat)) %>%
  mutate(idx = as.numeric(factor(paste(lon, lat), levels = unique(paste(lon, lat))))) 

storm_data = storm_data %>%
  mutate(time_format = base_time + hours(time),
         time_format = as.POSIXct(time_format, tz = "UTC"),
         time_diff = as.numeric(difftime(time_format, min(time_format), units = "hours")) + 1,
         year = year(time_format))

storm_data$name = "Wilma_2005"

### matrix W ----
library(spdep)
library(sp)

storm_loc = storm_removed_na %>%
  group_by(lon, lat) %>%
  summarize(mean_ws = mean(ws, na.rm = TRUE),
            mean_tp = mean(tp, na.rm = TRUE),
            .groups = 'drop')

lat = sort(unique(storm_loc$lat), decreasing = TRUE)
lon = sort(unique(storm_loc$lon))

# save index
index_matrix = matrix(NA, nrow = length(lat), ncol = length(lon))

# fill index
for (i in 1:nrow(storm_loc)) {
  lat_idx = which(lat == storm_loc$lat[i])
  lon_idx = which(lon == storm_loc$lon[i])
  index_matrix[lat_idx, lon_idx] = i
}
# change to vector
idx.mapping = as.vector(index_matrix)

# remove na
idx.mapping = idx.mapping[!is.na(idx.mapping)]

# reorder coord_df
coord_df = storm_loc[idx.mapping, ]

# Convert the data frame to a spatial object
coordinates(coord_df) = ~ lon + lat
proj4string(coord_df) = CRS("+proj=longlat +datum=WGS84")

# Plot the spatial points
plot(coord_df, main = "Spatial Points")

# Create a distance-based neighbor object 
dthreshold = .36
nb = dnearneigh(coordinates(coord_df), 0, dthreshold)

# Convert the neighbor object to an adjacency matrix
W = nb2mat(nb, style = "B", zero.policy = TRUE)

D = matrix(0, nrow = nrow(W), ncol = nrow(W))
diag(D) = colSums(W)

### INLA data ---
df_inla = data.frame(OBS = c(storm_data$ws, storm_data$tp * 1000),
                     Intercept = rep(c("ws", "tp"), each = nrow(storm_data)),
                     spatial_idx = c(storm_data$idx, storm_data$idx + max(storm_data$idx)),
                     year = rep(storm_data$year, 2),
                     time_ws = c(storm_data$time_diff, rep(NA, nrow(storm_data))),
                     time_tp = c(rep(NA, nrow(storm_data)), storm_data$time_diff))
df_inla$idx = 1:nrow(df_inla)
df_inla$time_idx = c(storm_data$time_diff, storm_data$time_diff + max(storm_data$time_diff))
df_inla$spatial = df_inla$spatial_idx
df_inla$interaction = df_inla$spatial_idx * df_inla$time_idx
df_loc = data.frame(OBS = c(storm_loc$mean_ws, storm_loc$mean_tp),
                    Intercept = rep(c("ws", "tp"), each = nrow(storm_loc)),
                    spatial_idx = 1:(nrow(storm_loc) * 2))

### W_time ----
max_time = max(df_inla$time_idx)/2
W_time = matrix(0, nrow = max_time, ncol = max_time)
diag(W_time) = 2
for (i in 1:(max_time - 1)) {
  W_time[i, (i + 1)] = 1
  W_time[(i + 1), i] = 1
}

### location ----
INDPMCAR.loc = inla(OBS ~ 0 + Intercept + 
                      f(spatial_idx, model = model.indmcar),
                    data = df_loc, family = "gamma",  
                    #control.family = list(variant = 1),
                    control.predictor = list(compute = TRUE),
                    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))

summary(INDPMCAR.loc)

PMCAR.loc = inla(OBS ~ 0 + Intercept + 
                   f(spatial_idx, model = model.mcar),
                 data = df_loc, family = "gamma",  
                 control.predictor = list(compute = TRUE),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))
summary(PMCAR.loc)

### spatio-temporal ----
library(INLAMSM)
library(INLA)
alpha.min = 0
alpha.max = 1

shape_prior = list(prior = "loggamma",  
                   param = c(2, 0.01))
K = 2
u = 1
model.indmcar = inla.INDMCAR.model(k = K, W = W, alpha.min = alpha.min,
                                   alpha.max = alpha.max)
time.mcar = inla.MCAR.model(k = K, W = W_time, alpha.min = alpha.min,
                            alpha.max = alpha.max)
INDPMCAR.time = inla(OBS ~ 0 + Intercept + 
                       f(spatial_idx, model = model.indmcar) +
                       f(time_idx, model = time.mcar) + 
                       f(interaction, model = "iid"),
                     data = df_inla, family = "gamma",  
                     control.predictor = list(compute = TRUE),
                     control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))
summary(INDPMCAR.time)
indpmcar.pred.time = INDPMCAR.time$summary.fitted.values
indpmcar.prop = sum(df_inla$OBS < indpmcar.pred.time$`0.025quant` | df_inla$OBS > indpmcar.pred.time$`0.975quant`) / length(df_inla$OBS)

INDPMCAR.rw = inla(OBS ~ 0 + Intercept + 
                     f(spatial_idx, model = model.indmcar) +
                     f(interaction, model = "iid") +
                     f(time_ws, model = "rw1",
                       hyper = list(prec = list(prior = "loggamma", param = c(.01, 1)))) + 
                     f(time_tp, model = "rw1",
                       hyper = list(prec = list(prior = "loggamma", param = c(.01, 10)))),
                     data = df_inla, family = "gamma",  
                     control.predictor = list(compute = TRUE),
                     control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))
summary(INDPMCAR.time)
indpmcar.pred.rw = INDPMCAR.rw$summary.fitted.values
indpmcar.rw.prop = sum(df_inla$OBS < indpmcar.pred.rw$`0.025quant` | df_inla$OBS > indpmcar.pred.rw$`0.975quant`) / length(df_inla$OBS)

u = 1
model.mcar = inla.MCAR.model(k = K, W = W, alpha.min = alpha.min,
                             alpha.max = alpha.max)
PMCAR.time = inla(OBS ~ 0 + Intercept +
                    f(spatial_idx, model = model.mcar) +
                    f(time_idx, model = time.mcar) + 
                    f(interaction, model = "iid") , 
                  scale = c(runif(nrow(storm_data), 19, 30), runif(nrow(storm_data), 0.001, 0.01)),
                  #scale = c(rep(5.758812, nrow(storm_data)), rep(0.001265347, nrow(storm_data))),
                  data = df_inla, family = "gamma", 
                  control.predictor = list(compute = TRUE),
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))
summary(PMCAR.time)
pmcar.pred.time = PMCAR.time$summary.fitted.values
pmcar.prop = sum(df_inla$OBS < pmcar.pred.time$`0.025quant` | df_inla$OBS > pmcar.pred.time$`0.975quant`) / length(df_inla$OBS)

PMCAR.time = inla(OBS ~ 0 + Intercept +
                    f(spatial_idx, model = model.mcar) +
                    f(interaction, model = "iid") +
                    f(time_ws, model = "rw1",
                      hyper = list(prec = list(prior = "loggamma", param = c(.01, 1))),
                      scale.model = TRUE) + 
                    f(time_tp, model = "rw1",
                      hyper = list(prec = list(prior = "loggamma", param = c(.01, 1))),
                      scale.model = TRUE),
                  scale = c(runif(nrow(storm_data), 19, 20), runif(nrow(storm_data), 0.001, 0.01)),
                  data = df_inla, family = "gamma",  
                  control.predictor = list(compute = TRUE),
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))
summary(PMCAR.time)
pmcar.pred.time = PMCAR.time$summary.fitted.values
pmcar.prop = sum(df_inla$OBS < pmcar.pred.time$`0.025quant` | df_inla$OBS > pmcar.pred.time$`0.975quant`) / length(df_inla$OBS)

### mmmodel ------
library(Matrix)
M.model = inla.Mmodel.model(k = K, W = W,
                            alpha.min = alpha.min,
                            alpha.max = alpha.max)
Mmodel.loc = inla(OBS ~ 0 + Intercept + 
                   f(spatial_idx, model = M.model),
                 data = df_loc, family = "gamma",  
                 control.predictor = list(compute = TRUE),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))

M.time = inla.Mmodel.model(k = K, W = W_time,
                           alpha.min = alpha.min,
                           alpha.max = alpha.max)
Mmodel.time = inla(OBS ~ 0 + Intercept + 
                     f(spatial, model = M.model) +
                     f(time, model = M.time) + 
                     f(idx, model = "iid") + 
                     f(time_idx, model = "iid") + 
                     f(spatial_idx, model = "iid") + 
                     f(time_ws, model = "rw1",
                       hyper = list(prec = list(prior = "loggamma", param = c(u, 0.01)))) + 
                     f(time_tp, model = "rw1",
                       hyper = list(prec = list(prior = "loggamma", param = c(u, 0.01)))),
                   data = df_inla, family = "gamma",  
                   # control.family = list(variant = 1),
                   control.predictor = list(compute = TRUE),
                   control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))

