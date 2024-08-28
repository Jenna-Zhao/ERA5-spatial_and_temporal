### correlation ----
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

### Wilma -----
Wilma = extract_variable("newWILMA_2005_10.nc")

Wilma_2005 = replace_na_below_threshold(ws = Wilma$ws[, , 220:240], 
                                        tp = Wilma$tp[, , 220:240], 
                                        quantile_value = .80)
Wilma_ws = Wilma_2005$ws
Wilma_tp = Wilma_2005$tp
Wilma_cor = find_correlation_matrix(Wilma_ws, Wilma_tp)

Wilma_cor_plot = ggplot(Wilma_cor, aes(x = Longitude, y = Latitude, fill = Corr)) +
  geom_tile() +
  borders("world", xlim = range(longitude), 
          ylim = range(latitude), colour = "black") +
  scale_fill_viridis_c(option = "D", na.value = "white") +
  labs(title = "Wilma 2005",
       x = "Longitude",
       y = "Latitude",
       fill = "") +
  theme_minimal() + 
  coord_quickmap(xlim = c(-95, -70), ylim = c(22, 32), expand = FALSE) + 
  theme(axis.title = element_text(color = 'black', size = 11),
        plot.title = element_text(color = "black", size = 10),
        plot.margin = margin(t = 2, r = 3, b = 2, l = 3, unit = "mm"))
        #legend.position = "none")

### Ian -----
Ian_first_half = extract_variable("newIAN_2022_09.nc", two_month = TRUE)
Ian_second_half = extract_variable("newIAN_2022_10.nc", two_month = TRUE)
Ian = list("time" = c(Ian_first_half$time, Ian_second_half$time),
           "latitude" = Ian_first_half$latitude,
           "longitude" = Ian_first_half$longitude,
           "ws" = abind(Ian_first_half$ws, Ian_second_half$ws, along = 3),
           "tp" = abind(Ian_first_half$tp, Ian_second_half$tp, along = 3))

Ian_2022 = replace_na_below_threshold(ws = Ian$ws[, , 140:170], 
                                      tp = Ian$tp[, , 140:170], 
                                      quantile_value = .80)
Ian_ws = Ian_2022$ws
Ian_tp = Ian_2022$tp
Ian_cor = find_correlation_matrix(Ian_ws, Ian_tp)

Ian_cor_plot = ggplot(Ian_cor, aes(x = Longitude, y = Latitude, fill = Corr)) +
  geom_tile() +
  borders("world", xlim = range(longitude), 
          ylim = range(latitude), colour = "black") +
  scale_fill_viridis_c(option = "D", na.value = "white") +
  labs(title = "Ian 2022",
       x = "Longitude",
       y = "Latitude",
       fill = "") +
  theme_minimal() + 
  coord_quickmap(xlim = c(-95, -70), ylim = c(22, 32), expand = FALSE) + 
  theme(axis.title = element_text(color = 'black', size = 11),
        plot.title = element_text(color = "black", size = 10),
        plot.margin = margin(t = 2, r = 3, b = 2, l = 3, unit = "mm"))
        #legend.position = "none")

### IRMA -----
IRMA_first_half = extract_variable("newIRMA_2017_08.nc", two_month = TRUE)
IRMA_second_half = extract_variable("newIRMA_2017_09.nc", two_month = TRUE)
IRMA = list("time" = c(IRMA_first_half$time, IRMA_second_half$time),
           "latitude" = IRMA_first_half$latitude,
           "longitude" = IRMA_first_half$longitude,
           "ws" = abind(IRMA_first_half$ws, IRMA_second_half$ws, along = 3),
           "tp" = abind(IRMA_first_half$tp, IRMA_second_half$tp, along = 3))

IRMA_2017 = replace_na_below_threshold(ws = IRMA$ws[, , 280:310], 
                                      tp = IRMA$tp[, , 280:310], 
                                      quantile_value = .80)
IRMA_ws = IRMA_2017$ws
IRMA_tp = IRMA_2017$tp
IRMA_cor = find_correlation_matrix(IRMA_ws, IRMA_tp)

IRMA_cor_plot = ggplot(IRMA_cor, aes(x = Longitude, y = Latitude, fill = Corr)) +
  geom_tile() +
  borders("world", xlim = range(longitude), 
          ylim = range(latitude), colour = "black") +
  scale_fill_viridis_c(option = "D", na.value = "white") +
  labs(title = "Irma 2017",
       x = "Longitude",
       y = "Latitude",
       fill = "") +
  theme_minimal() + 
  coord_quickmap(xlim = c(-95, -70), ylim = c(22, 32), expand = FALSE) + 
  theme(axis.title = element_text(color = 'black', size = 11),
        plot.title = element_text(color = "black", size = 10),
        plot.margin = margin(t = 2, r = 3, b = 2, l = 3, unit = "mm"))
       # legend.position = "none")

### Charley -----
Charley = extract_variable("newCHARLEY_2004_08.nc")
Charley_2004 = replace_na_below_threshold(ws = Charley$ws[, , 100:140], 
                                          tp = Charley$tp[, , 100:140], 
                                          quantile_value = .80)
Charley_ws = Charley_2004$ws
Charley_tp = Charley_2004$tp
Charley_cor = find_correlation_matrix(Charley_ws, Charley_tp)

Charley_cor_plot = ggplot(Charley_cor, aes(x = Longitude, y = Latitude, fill = Corr)) +
  geom_tile() +
  borders("world", xlim = range(longitude), 
          ylim = range(latitude), colour = "black") +
  scale_fill_viridis_c(option = "D", na.value = "white") +
  labs(title = "Charley 2004",
       x = "Longitude",
       y = "Latitude",
       fill = "") +
  theme_minimal() + 
  coord_quickmap(xlim = c(-95, -70), ylim = c(22, 32), expand = FALSE) + 
  theme(axis.title = element_text(color = 'black', size = 11),
        plot.title = element_text(color = "black", size = 10),
        plot.margin = margin(t = 2, r = 3, b = 2, l = 3, unit = "mm"))
       # legend.position = "none")

### Correlation map -----
library(patchwork)
cor_combined = (Wilma_cor_plot + Ian_cor_plot) /
  (IRMA_cor_plot + Charley_cor_plot) +
  plot_layout(guides = "collect")
ggsave("plot/cor_combined.pdf", cor_combined , width = 20, height = 10, units = "cm")


## distribution
library(fitdistrplus)
library(sn)
wilmacor = as.vector(na.omit(as.vector(Wilma_cor$Corr)))
iancor = as.vector(na.omit(as.vector(Ian_cor$Corr)))
irmacor = as.vector(na.omit(as.vector(IRMA_cor$Corr)))
charleycor = as.vector(na.omit(as.vector(Charley_cor$Corr)))

wilman = fitdist(wilmacor, "norm")
iann = fitdist(iancor, "norm")
irman = fitdist(irmacor, "norm")
charleyn = fitdist(charleycor, "norm")

plot(wilman)
wsn = fitdist(ws_90, "norm")
wsw = fitdist(ws_90, "weibull")
wsg = fitdist(ws_90, "gamma")
wsninv = fitdist(ws_90inv, "norm")
wswinv = fitdist(ws_90inv, "weibull")
wsginv = fitdist(ws_90inv, "gamma")
hist(Wilma_cor$Corr, breaks = 50)

par(mfrow = c(2, 2))
hist(Wilma_cor$Corr, breaks = 50, freq = FALSE,
     xlab = "Correlation (N = 504)", main = "Wilma 2005")
#curve(dnorm(x, mean = wilman$estimate["mean"], sd = wilman$estimate["sd"]), add = TRUE, col = "red")
hist(Ian_cor$Corr, breaks = 50, freq = FALSE,
     xlab = "Correlation (N = 468)",  main = "Ian 2022")
#curve(dnorm(x, mean = iann$estimate["mean"], sd = iann$estimate["sd"]), add = TRUE, col = "red")
hist(IRMA_cor$Corr, breaks = 50, freq = FALSE,
     xlab = "Correlation (N = 889)", main = "Irma 2017")
#curve(dnorm(x, mean = irman$estimate["mean"], sd = irman$estimate["sd"]), add = TRUE, col = "red")
hist(Charley_cor$Corr, breaks = 50, freq = FALSE,
     xlab = "Correlation (N = 563)", main = "Charley 2004")
#curve(dnorm(x, mean = charleyn$estimate["mean"], sd = charleyn$estimate["sd"]), add = TRUE, col = "red")


