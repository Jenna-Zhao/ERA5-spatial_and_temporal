### Load data ----
Ian_first_half = extract_variable("newIAN_2022_09.nc", two_month = TRUE)
Ian_second_half = extract_variable("newIAN_2022_10.nc", two_month = TRUE)
Ian = list("time" = c(Ian_first_half$time, Ian_second_half$time),
             "latitude" = Ian_first_half$latitude,
             "longitude" = Ian_first_half$longitude,
             "ws" = abind(Ian_first_half$ws, Ian_second_half$ws, along = 3),
             "tp" = abind(Ian_first_half$tp, Ian_second_half$tp, along = 3))

Ian_2022 = replace_na_below_threshold(ws = Ian$ws[, , 120:180], 
                                      tp = Ian$tp[, , 120:180], 
                                      quantile_value = .90)

timepoints = c(16, 31, 56)  # 135 150 and 175
longitude = Ian$longitude
latitude = Ian$latitude

### function for map plotting ----
plot_map_matrix = function(var, t){
  long_format = expand.grid(longitude,
                            latitude) %>%
    mutate(var = as.vector(var[, , t]))
  colnames(long_format) = c("longitude", "latitude", "var")
  # Plotting the correlation on a map
  ggplot(long_format, 
           aes(x = longitude, y = latitude, fill = var)) +
      #geom_raster() +
      geom_tile() +
      borders("world", xlim = range(longitude), 
              ylim = range(latitude), colour = "black") +
      scale_fill_viridis_c(option = "D", na.value = "white") +
      labs(title = "",
           x = "Longitude",
           y = "Latitude",
           fill = "") +
      theme_minimal() + 
      coord_quickmap(xlim = c(-95, -70), ylim = c(22, 32), expand = FALSE) + 
      theme(axis.title = element_text(color='black', size = 10),
            legend.position = "none",
            plot.margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "mm"),
            plot.title = element_blank()) 
  }

### Visualization (before and after) -----
library(patchwork)
# ws
ws1 = plot_map_matrix(Ian$ws, timepoints[1] + 119)
ws2 = plot_map_matrix(Ian$ws, timepoints[2] + 119)
ws3 = plot_map_matrix(Ian$ws, timepoints[3] + 119)
# Combined plots
ws_combined <- ws1 / ws2 / ws3 +
  plot_layout(heights = c(1, 1, 1))
ws_combined 
ggsave("plot/ws_combined.pdf", ws_combined , width = 15, height = 20, units = "cm")

# tp
tp1 = plot_map_matrix(Ian$tp, timepoints[1] + 119)
tp2 = plot_map_matrix(Ian$tp, timepoints[2] + 119)
tp3 = plot_map_matrix(Ian$tp, timepoints[3] + 119)
tp_combined <- tp1 / tp2 / tp3 +
  plot_layout(heights = c(1, 1, 1))
tp_combined 
ggsave("plot/tp_combined.pdf", tp_combined , width = 15, height = 20, units = "cm")

# ws after
ws_2022_1 = plot_map_matrix(Ian_2022$ws, timepoints[1])
ws_2022_2 = plot_map_matrix(Ian_2022$ws, timepoints[2])
ws_2022_3 = plot_map_matrix(Ian_2022$ws, timepoints[3])
ws_2022_combined <- ws_2022_1 / ws_2022_2 / ws_2022_3  +
  plot_layout(heights = c(1, 1, 1))
ws_2022_combined 
ggsave("plot/ws_2022_combined.pdf", ws_2022_combined  , width = 15, height = 20, units = "cm")

# tp after
tp_2022_1 = plot_map_matrix(Ian_2022$tp, timepoints[1])
tp_2022_2 = plot_map_matrix(Ian_2022$tp, timepoints[2])
tp_2022_3 = plot_map_matrix(Ian_2022$tp, timepoints[3])
tp_2022_combined <- tp_2022_1 / tp_2022_2 / tp_2022_3  +
  plot_layout(heights = c(1, 1, 1))
tp_2022_combined 
ggsave("plot/tp_2022_combined.pdf", tp_2022_combined  , width = 15, height = 20, units = "cm")


### Density -----
par(mfrow = c(1, 2), 
    mar = c(5, 3, 2, 1))
hist(Ian$ws, breaks = 100, xlab = "Wind Speed", freq = FALSE, col = rgb(0.8, 0.8, 0.8, .7), border = NA, main = "")
lines(density(Ian$ws))
hist(Ian$tp, breaks = 100, xlab = "Total Precipitation", freq = FALSE, col = rgb(0.8, 0.8, 0.8, .7), border = NA, main = "")
lines(density(Ian$tp))
hist(Ian_2022$ws, breaks = 100, xlab = "Wind Speed", freq = FALSE, col = rgb(0.8, 0.8, 0.8, .7), border = NA,  main = "")
lines(density(Ian_2022$ws, na.rm = TRUE))
hist(Ian_2022$tp, breaks = 100, xlab = "Total Precipitation", freq = FALSE, col = rgb(0.8, 0.8, 0.8, .7), border = NA,  main = "")
lines(density(Ian_2022$tp, na.rm = TRUE))



### fitting ----
library(fitdistrplus)
library(mixtools)  # For mixture models
ws_90 = (as.vector(na.omit(as.vector(Ian_2022$ws))))
ws_90inv = 1/(as.vector(na.omit(as.vector(Ian_2022$ws))))
wsn = fitdist(ws_90, "norm")
wsw = fitdist(ws_90, "weibull")
wsg = fitdist(ws_90, "gamma")
wsninv = fitdist(ws_90inv, "norm")
wswinv = fitdist(ws_90inv, "weibull")
wsginv = fitdist(ws_90inv, "gamma")

# visulization ws
par(mfrow = c(1, 2), 
    mar = c(5, 3, 3, 1))
hist(ws_90, breaks = 50, freq = FALSE, main = "Empirical Density", xlab = "Wind Speed")
curve(dnorm(x, mean = wsn$estimate["mean"], sd = wsn$estimate["sd"]), add = TRUE, col = "red")
curve(dweibull(x, shape = wsw$estimate["shape"], scale = wsw$estimate["scale"]), add = TRUE, col = "green", lty = "dashed")
curve(dgamma(x, shape = wsg$estimate["shape"], rate = wsg$estimate["rate"]), add = TRUE, col = "blue", lty = "dotted")
legend("topright", legend = c("Normal", "Weibull", "Gamma"), col = c("red", "green", "blue"), 
       lty = c(1, 2, 3), bty = "n")
qqcomp(list(wsn, wsw, wsg), legendtext = c("Normal", "Weibull", "Gamma"))
# visulization ws inv
hist(ws_90inv, breaks = 50, freq = FALSE, main = "", xlab = "inverse of wind speed")
curve(dnorm(x, mean = wsninv$estimate["mean"], sd = wsninv$estimate["sd"]), add = TRUE, col = "red")
curve(dweibull(x, shape = wswinv$estimate["shape"], scale = wswinv$estimate["scale"]), add = TRUE, col = "green", lty = "dashed")
curve(dgamma(x, shape = wsginv$estimate["shape"], rate = wsginv$estimate["rate"]), add = TRUE, col = "blue", lty = "dotted")
legend("topleft", legend = c("Normal", "Weibull", "Gamma"), col = c("red", "green", "blue"), 
       lty = c(1, 2, 3), bty = "n")

# tp
tp_90 = as.vector(na.omit(as.vector(Ian_2022$tp)))
tp_90inv = 1/as.vector(na.omit(as.vector(Ian_2022$tp)))
tpn = fitdist(tp_90, "norm")
tpw = fitdist(tp_90, "weibull")
tpg = fitdist(tp_90, "gamma")
tpb = fitdist(tp_90, "beta")
tpninv = fitdist(tp_90inv, "norm")
tpwinv = fitdist(tp_90inv, "weibull")
tpginv = fitdist(tp_90inv, "gamma")
# visulization tp
par(mfrow = c(1, 2), 
    mar = c(5, 3, 3, 1))
hist(tp_90, breaks = 50, freq = FALSE, main = "Empirical Density", xlab = "Total Precipitation")
curve(dnorm(x, mean = tpn$estimate["mean"], sd = tpn$estimate["sd"]), add = TRUE, col = 2)
curve(dweibull(x, shape = tpw$estimate["shape"], scale = tpw$estimate["scale"]), add = TRUE, col = 3, lty = "dashed")
curve(dgamma(x, shape = tpg$estimate["shape"], rate = tpg$estimate["rate"]), add = TRUE, col = 4, lty = "dotted")
curve(dgamma(x, shape = tpb$estimate["shape1"], rate = tpb$estimate["shape2"]), add = TRUE, col = 5, lty = "dotted")
legend("topright", legend = c("Normal", "Weibull", "Gamma", "Beta"), col = 2:5, 
       lty = c(1, 2, 3, 4), bty = "n")
qqcomp(list(tpn, tpw, tpg, tpb), legendtext = c("Normal", "Weibull", "Gamma", "Beta"))
# visulization tp inv
par(mfrow = c(1, 2), 
    mar = c(5, 3, 3, 1))
hist(tp_90inv, breaks = 50, freq = FALSE, main = "Empirical Density", xlab = "Inverse of Total Precipitation")
curve(dnorm(x, mean = tpninv$estimate["mean"], sd = tpninv$estimate["sd"]), add = TRUE, col = "red")
curve(dweibull(x, shape = tpwinv$estimate["shape"], scale = tpwinv$estimate["scale"]), add = TRUE, col = "green", lty = "dashed")
curve(dgamma(x, shape = tpginv$estimate["shape"], rate = tpginv$estimate["rate"]), add = TRUE, col = "blue", lty = "dotted")
legend("topright", legend = c("Normal", "Weibull", "Gamma"), col = c("red", "green", "blue"), 
       lty = c(1, 2, 3), bty = "n")
qqcomp(list(tpninv, tpwinv, tpginv), legendtext = c("Normal", "Weibull", "Gamma"))



### correlation ----
ws = Ian_2022$ws
tp_inv = 1/Ian_2022$tp
tp = Ian_2022$tp
# Calculate correlation matrix
find_correlation_matrix = function(ws, tp) {
  # Initialize a matrix to store correlations
  correlation_matrix = matrix(NA, nrow = dim(ws)[1] * dim(ws)[2],
                              ncol = 6)
  colnames(correlation_matrix) = c("Longitude", "Latitude",
                                   "Corr", "Time",
                                   "ws_mean", "tp_mean")
  # Loop through each location
  for (i in 1:dim(ws)[1]) {
    for (j in 1:dim(ws)[2]) {
      ind = (j - 1) * dim(ws)[1] + i
      correlation_matrix[ind, 1:2] = c(longitude[i], latitude[j])
      # Extract the time series for the current location
      ws_series = ws[i, j, ]
      tp_series = tp[i, j, ]
      # Remove NA values
      valid_indices = !is.na(ws_series) & !is.na(tp_series)
      ws_series = ws_series[valid_indices]
      tp_series = tp_series[valid_indices]
      # Calculate the correlation if there are valid data points
      if (length(ws_series) > 1 && length(tp_series) > 1) {
        correlation_matrix[ind, 3] = cor(ws_series, tp_series)
        correlation_matrix[ind, 4] = length(ws_series)
        correlation_matrix[ind, 5] = mean(ws_series)
        correlation_matrix[ind, 6] = mean(tp_series)
      } else{
        correlation_matrix[ind, 3] = NA
        correlation_matrix[ind, 4] = 0
        correlation_matrix[ind, 5] = NA
        correlation_matrix[ind, 6] = NA
      }
    }
  }
  return(as.data.frame(correlation_matrix))
}

Ian_cor = find_correlation_matrix(ws, tp_inv)

cor_plot = ggplot(Ian_cor, aes(x = Longitude, y = Latitude, fill = Corr)) +
  geom_tile() +
  borders("world", xlim = range(longitude), 
          ylim = range(latitude), colour = "black") +
  scale_fill_viridis_c(option = "D", na.value = "white") +
  labs(title = "",
       x = "Longitude",
       y = "Latitude",
       fill = "") +
  theme_minimal() + 
  coord_quickmap(xlim = c(-95, -70), ylim = c(22, 32), expand = FALSE) + 
  theme(axis.title = element_text(color = 'black', size = 11),
        plot.margin = margin(t = 0, r = 1, b = 0, l = 5, unit = "mm"),
        plot.title = element_blank()) 
ggsave("plot/correlation.pdf", cor_plot , width = 15, height = 7, units = "cm")
