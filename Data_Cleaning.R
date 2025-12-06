# Load necessary libraries 
library(lubridate)
library(dplyr)
library(tsibble)

#Load Data
Energy = read.csv(file.choose(), header = TRUE)
str(Energy)

# Total number of missing values in the entire dataframe
sum(is.na(Energy))

# Count of duplicate rows
sum(duplicated(Energy))

#Combine Date and Hour
Energy <- Energy %>%
  mutate(time = ymd_h(paste(date, hour - 1)))  # hour 1 = 00:00

#Check if there are time gaps
Energy_ts <- Energy %>%
  as_tsibble(index = time)

Energy_ts %>% has_gaps()

#Drop irrelevant columns
Energy = subset(Energy, select = -c(date, hour))

#Energy monthly
Energy_daily <- Energy %>%
  mutate(day = as.Date(time)) %>%
  group_by(day) %>%
  summarise(
    demand_daily_avg = mean(hourly_demand, na.rm = TRUE),
    price_daily_avg  = mean(hourly_average_price, na.rm = TRUE),
    n_hours          = n(), .groups = "drop"
  )

Energy_monthly<- Energy_daily %>%
  mutate(month = floor_date(day, "month")) %>%
  group_by(month) %>%
  summarise(
    demand_monthly_avg = mean(demand_daily_avg, na.rm = TRUE),
    price_monthly_avg  = mean(price_daily_avg, na.rm = TRUE),
    days_obs           = n(), .groups = "drop"
  )

#Check for abnormal values/outliers 

numeric_columns = c("demand_monthly_avg", "price_monthly_avg")
plot_hist_with_npdf = function(x, col_name) {
  h_for_ylim = hist(Energy_monthly[[col_name]], plot = FALSE)   # new: get counts for ylim
  hist_data = hist(Energy_monthly[[col_name]],
                   main = paste("Histogram of", col_name),
                   xlab = col_name,
                   ylim = c(0, max(h_for_ylim$counts) * 1.2))   # new: apply ylim here
  
  xpt = seq(min(x), max(x), length.out = 50)
  n_den = dnorm(xpt, mean(x), sd(x))
  
  # Calculate the histogram scale factor
  bin_width = diff(hist_data$breaks[1:2])
  ypt = n_den * length(x) * bin_width
  lines(xpt, ypt, col = "blue")
}


par(mfrow = c(2, 1))
for (col in numeric_columns) {
  plot_hist_with_npdf(Energy_monthly[[col]], col)
  boxplot(Energy_monthly[[col]], xlab = col)
  print(summary(Energy_monthly[[col]]))
}
par(mfrow = c(1, 1))

start_y <- lubridate::year(min(Energy_monthly$month))
start_m <- lubridate::month(min(Energy_monthly$month))

price_ts <- ts(Energy_monthly$price_monthly_avg,
               start = c(start_y, start_m),
               frequency = 12)
# Plot the ts
plot(price_ts, ylab = "Avg HOEP (¢/kWh)", xlab = "Time",
     main = "Monthly Average Price (ts)")

#Required transform
model <- lm(Energy_monthly$price_monthly_avg ~ 1, data =Energy_monthly)
boxcox.model = MASS::boxcox(model, lambda = seq(-2, 2, 0.1), plotit = TRUE)
opt.lambda = boxcox.model$x[which.max(boxcox.model$y)]

Energy_monthly$price_monthly_avg_bc <- (Energy_monthly$price_monthly_avg^opt.lambda - 1) / opt.lambda

model <- lm(Energy_monthly$price_monthly_avg ~ 1, data =Energy_monthly)
boxcox.model = MASS::boxcox(model, lambda = seq(-2, 2, 0.1), plotit = TRUE)
opt.lambda = boxcox.model$x[which.max(boxcox.model$y)]

