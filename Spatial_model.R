## Long Format Data ------------------------
library(tidyr)
library(dplyr)

long_format = expand.grid(lon = longitude,
                          lat = latitude,
                          time = time) %>%
  mutate(ws = as.vector(abby_remove$ws), 
         tp = as.vector(abby_remove$tp))

abby_Test = na.omit(long_format)

head(abby_Test)

# change time format
## Base time
base_time = as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
## Convert hours since 1900 to POSIXct datetime
converted_times = base_time + as.difftime(abby_Test$time, units = "hours")
## Format the datetime
abby_Test$time = format(converted_times, "%Y-%m-%d %H:%M")

## Model by mean value -----
### location mean value ------------------
abby_loc = abby_Test %>%
  group_by(lon, lat) %>%
  summarize(mean_ws = mean(ws, na.rm = TRUE),
            mean_tp = mean(tp, na.rm = TRUE),
            .groups = 'drop')

head(abby_loc)


### Prepare Data for INLA ------------
library(INLA)
library(INLAMSM)

d = data.frame(OBS = c(abby_loc$mean_ws, abby_loc$mean_tp),
               Intercept = rep(c("ws", "tp"), each = nrow(abby_loc)))
d$idx = 1:length(d$OBS)
d$Intercept = as.factor(d$Intercept)

# Check the structure of the transformed data
head(d)

###Create an adjacency matrix ---------------------
library(sp)
library(spdep)

# Convert the data frame to a spatial object
coordinates(abby_loc) = ~ lon + lat
proj4string(abby_loc) = CRS("+proj=longlat +datum=WGS84")

# Plot the spatial points
plot(abby_loc, main = "Spatial Points")

# Create a distance-based neighbor object (threshold set to 2 degrees)
dthreshold = .36
nb = dnearneigh(coordinates(abby_loc), 0, dthreshold)

# Convert the neighbor object to an adjacency matrix
W = nb2mat(nb, style = "B", zero.policy = TRUE)

# Define model parameters and constraints
k = 2 # ws and tp
alpha.min = 0
alpha.max = 1

# Constraint matrix
A = kronecker(Diagonal(k, 1), Matrix(1, ncol = nrow(W), nrow = 1))
e = rep(0, k)


### Fit the IMCAR model with intercept term ------------
model.imcar = inla.IMCAR.model(k = k, W = W)
IMCAR = inla(OBS ~ 0 + Intercept + 
               f(idx, model =  model.imcar,
                 extraconstr = list(A = as.matrix(A), e = e)),
             data = d, family = "gamma", 
            # control.family = list(hyper = list(shape = prior_phi)),  # Set the prior for phi here
             control.predictor = list(compute = TRUE, link = 1),
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE))

summary(IMCAR)

model.mcar = inla.MCAR.model(k = k, W = W, alpha.min = alpha.min,
                             alpha.max = alpha.max)
PMCAR = inla(OBS ~ 0 + Intercept + f(idx, model = model.mcar),
             data = d, family = "gamma",
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
             control.predictor = list(compute = TRUE, link = 1))

summary(PMCAR)

model.m = inla.Mmodel.model(k = k, W = W, alpha.min = alpha.min,
                            alpha.max = alpha.max)
Mmodel = inla(OBS ~ 0 + Intercept + f(idx, model = model.m), 
              data = d, family = "gamma",
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
              control.predictor = list(compute = TRUE, link = 1))
summary(Mmodel)


### plot -------
library(gridExtra)

n = nrow(abby_loc)
abby_loc$IMCAR.ws = IMCAR$summary.fitted[1:n, "mean"]
abby_loc$IMCAR.tp = IMCAR$summary.fitted[n + 1:n, "mean"]
abby_loc$PMCAR.ws = PMCAR$summary.fitted[1:n, "mean"]
abby_loc$PMCAR.tp = PMCAR$summary.fitted[n + 1:n, "mean"]
abby_loc$Mmodel.ws = Mmodel$summary.fitted[1:n, "mean"]
abby_loc$Mmodel.tp = Mmodel$summary.fitted[n + 1:n, "mean"]

abby_loc_df = as.data.frame(abby_loc)

Test_corr = Test_corr %>%
  left_join(abby_loc_df, by = c("Longitude" = "lon", "Latitude" = "lat"))

plot1 = plot_map(df = Test_corr, main = "IMCAR.ws", var = "IMCAR.ws")
plot2 = plot_map(df = Test_corr, main = "IMCAR.tp", var = "IMCAR.tp")
plot3 = plot_map(df = Test_corr, main = "PMCAR.ws", var = "PMCAR.ws")
plot4 = plot_map(df = Test_corr, main = "PMCAR.tp", var = "PMCAR.tp")
plot5 = plot_map(df = Test_corr, main = "Mmodel.ws", var = "Mmodel.ws")
plot6 = plot_map(df = Test_corr, main = "Mmodel.tp", var = "Mmodel.tp")
# Arrange the plots into a 3x2 grid
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 3, ncol = 2)

### predicted ----
imcar.pred = IMCAR$summary.fitted.values
pmcar.pred = PMCAR$summary.fitted.values
mmodel.pred = Mmodel$summary.fitted.values

# Proportion of data points outside the 95% prediction intervals
imcar.prop = sum(d$OBS < imcar.pred$`0.025quant` | d$OBS > imcar.pred$`0.975quant`) / length(d$OBS)
pmcar.prop = sum(d$OBS < pmcar.pred$`0.025quant` | d$OBS > pmcar.pred$`0.975quant`) / length(d$OBS)
mmodel.prop = sum(d$OBS < mmodel.pred$`0.025quant` | d$OBS > mmodel.pred$`0.975quant`) / length(d$OBS)


compare_prediction_intervals = function(model, upper_quantile = 0.975) {
  
  # Compute the prediction intervals using inla.posterior.sample
  num_samples = 1000
  posterior_samples = inla.posterior.sample(num_samples, model)
  
  N = nrow(model$.args$data)
  sampled_predictions = matrix(NA, nrow = N, ncol = num_samples)
  
  for (i in 1:num_samples) {
    sample = posterior_samples[[i]]
    linear_predictor = sample$latent[grep("Predictor", rownames(sample$latent))]
    sampled_predictions[, i] = exp(linear_predictor) # Apply the inverse link function 
  }
  
  mean_predictions_sample = apply(sampled_predictions, 1, mean)
  lower_bound_sample = apply(sampled_predictions, 1, quantile, 1 - upper_quantile)
  upper_bound_sample = apply(sampled_predictions, 1, quantile,  upper_quantile)
  
  posterior_intervals <- data.frame(
    mean = mean_predictions_sample,
    lower = lower_bound_sample,
    upper = upper_bound_sample,
    observed = model$.args$data$OBS
  )
  
  return(posterior_intervals)
}

imcar.post.pred = compare_prediction_intervals(model = IMCAR)
pmcar.post.pred = compare_prediction_intervals(model = PMCAR)
mmodel.post.pred = compare_prediction_intervals(model = Mmodel)

imcar.post.pred$inside_interval = with(imcar.post.pred, observed >= lower & observed <= upper)
pmcar.post.pred$inside_interval = with(pmcar.post.pred, observed >= lower & observed <= upper)
mmodel.post.pred$inside_interval = with(mmodel.post.pred, observed >= lower & observed <= upper)

sum(!imcar.post.pred$inside_interval)
sum(!pmcar.post.pred$inside_interval)
sum(!mmodel.post.pred$inside_interval)

### elpd -----
imcar.elpd_loo = sum(log(IMCAR$cpo$cpo), na.rm = TRUE)
imcar.elpd_waic = -IMCAR$waic$waic/2

pmcar.elpd_loo = sum(log(PMCAR$cpo$cpo), na.rm = TRUE)
pmcar.elpd_waic = -PMCAR$waic$waic/2

mmodel.elpd_loo = sum(log(Mmodel$cpo$cpo), na.rm = TRUE)
mmodel.elpd_waic = -Mmodel$waic$waic/2

