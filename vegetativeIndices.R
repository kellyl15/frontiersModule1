# Vegetative Index Calculations
# NDVI, NDRE, EVI, CI, RTVIcore, SAVI, and MSAVI

rm(list=ls())

library("dplyr")
library("ggplot2")
library("mgcv") # for generalized additive model (gam)

# read in combined csv wide and long
long_data = read.csv("outputs/combined_csv_long_gdds.csv")

# filter the dataset for the MS red, nir, rededge, green, and blue bands - MEAN

    # using MS and mean - adjust if using another stat metric or diff bands
    # could make a function to filter and rename bands 

red_band = filter(long_data, grepl("MS_red", Variable) & 
                             grepl("mean", Variable) & 
                            !grepl("rededge", Variable))

nir_band = filter(long_data, grepl("MS_nir", Variable) & 
                             grepl("mean", Variable))

rededge_band = filter(long_data, grepl("MS_rededge", Variable) & 
                                 grepl("mean", Variable))

green_band = filter(long_data, grepl("MS_green", Variable) & 
                               grepl("mean", Variable))

blue_band = filter(long_data, grepl("MS_blue", Variable) & 
                              grepl("mean", Variable))


# rename columns - fixes duplication of rows
red_band = rename(red_band, Value_red = Value, 
                  Date_red = Date, Variable_red = Variable)

nir_band = rename(nir_band, Value_nir = Value, 
                  Date_nir = Date, Variable_nir = Variable)

rededge_band = rename(rededge_band, Value_rededge = Value, 
                      Date_rededge = Date, Variable_rededge = Variable)

green_band = rename(green_band, Value_green = Value, 
                    Date_green = Date, Variable_green = Variable)

blue_band = rename(blue_band, Value_blue = Value,
                   Date_blue = Date, Variable_blue = Variable)


# merge all necessary bands based on the date and ID
merged_data = Reduce(function(x, y) {
  merge(x, y, by = c("id", "Accumulated.GDD"), all = TRUE)
}, list(red_band, nir_band, rededge_band, green_band, blue_band))


##### calculate each vegetative index ###############################################

merged_data = transform(merged_data,
   NDVI = (Value_nir - Value_red) / (Value_nir + Value_red),
   
   NDRE = (Value_nir - Value_rededge) / (Value_nir + Value_rededge),
   
   EVI = 2.5 * (Value_nir - Value_red) / (Value_nir + 6 * Value_red - 7.5 * 
                                          Value_blue + 1),
   
   CI = (Value_nir / Value_green) - 1,
   
   RTVIcore = 100 * (Value_nir - Value_rededge) - 10 * (Value_nir - Value_green),
   
   SAVI = ((Value_nir - Value_red) / (Value_nir + Value_red + 0.5)) * (1 + 0.5),
   
   MSAVI = (2 * Value_nir + 1 - sqrt((2 * Value_nir + 1)^2 - 8 * 
                                    (Value_nir - Value_red))) / 2
)

    
    
###### Output results to CSV ########################################################

# Select columns to export
final_output = merged_data[, c("id", "Date_nir", "Accumulated.GDD", "NDVI", "NDRE",
                                "EVI", "CI", "RTVIcore", "SAVI", "MSAVI")]

# Rename columns
colnames(final_output) = c("Plot ID", "Date", "GDD", "NDVI", "NDRE", "EVI", "CI", 
                           "RTVIcore", "SAVI", "MSAVI")

# Write the results to a CSV file
write.csv(final_output, "outputs/vegetative_indices_results.csv", row.names = FALSE)


###### plot NDVI over time ##########################################################

# Create an empty list to store results
results_list = list()

# Create an empty data frame to store model comparison for each ID
model_comparison_all = data.frame(
  Plot_ID = character(),
  Model = character(),
  RMSE = numeric(),
  R_squared = numeric(),
  stringsAsFactors = FALSE
)

# plot ids vector to loop through
unique_ids = unique(merged_data$id)

# subset data for each plot, fit loess model, extract information, and plot ndvi over time
for (id in unique_ids) {
  df_subset = na.omit(merged_data[merged_data$id == id, ])
  
 ############################################################################ 
  # model fitting 
  
  # Fit the loess model and get smoothed values
  loess_model = loess(NDVI ~ Accumulated.GDD, data = df_subset,
                      span = 0.9, degree = 1) # may need to adjust span- was getting
                                              # warnings with default value
  smoothed_values = predict(loess_model, newdata = df_subset$Accumulated.GDD)
  
  # polynomial model
  poly_model = lm(NDVI ~ poly(Accumulated.GDD, 2), data = df_subset)
  smoothed_values_poly = predict(poly_model, newdata = df_subset)
  
  # Generalized Additive Model (GAM)
  gam_model = gam(NDVI ~ s(Accumulated.GDD, k = 4), data = df_subset)
  smoothed_values_gam = predict(gam_model, newdata = df_subset)

  # RMSE calculation function
  rmse = function(actual, predicted) {
    sqrt(mean((actual - predicted)^2))
  }
  
  # Calculate RMSE for each model
  rmse_loess = rmse(df_subset$NDVI, smoothed_values)
  
  rmse_poly = rmse(df_subset$NDVI, smoothed_values_poly)
  
  rmse_gam = rmse(df_subset$NDVI, smoothed_values_gam)
  
  # R-squared calculation
  r2_loess = 1 - sum((df_subset$NDVI - smoothed_values)^2) / 
                 sum((df_subset$NDVI - mean(df_subset$NDVI))^2)
  
  r2_poly = 1 - sum((df_subset$NDVI - smoothed_values_poly)^2) / 
                sum((df_subset$NDVI - mean(df_subset$NDVI))^2)
  
  r2_gam = 1 - sum((df_subset$NDVI - smoothed_values_gam)^2) / 
               sum((df_subset$NDVI - mean(df_subset$NDVI))^2)
  
  # Create model comparison for the current iteration
  model_comparison_iter = data.frame(
    Plot_ID = rep(id, 3),
    Model = c("Loess", "Polynomial", "GAM"),
    RMSE = c(rmse_loess, rmse_poly, rmse_gam),
    R_squared = c(r2_loess, r2_poly, r2_gam)
  )
  
  # Append each iteration's results to the overall data frame
  model_comparison_all = rbind(model_comparison_all, model_comparison_iter)
  
  df_subset$smoothed_loess = smoothed_values
  df_subset$smoothed_poly = smoothed_values_poly
  df_subset$smoothed_gam = smoothed_values_gam
  
  # Create the plot of model fits with a legend
  # model_plot = ggplot(df_subset, aes(x = Accumulated.GDD, y = NDVI)) +
  #   geom_point(color = "black") +  # Plot original data points
  #   geom_line(aes(y = smoothed_loess, color = "Loess", linetype = "Loess")) +  # Loess model
  #   geom_line(aes(y = smoothed_poly, color = "Polynomial", linetype = "Polynomial")) +  # Polynomial model
  #   geom_line(aes(y = smoothed_gam, color = "GAM", linetype = "GAM")) +  # GAM model
  #   scale_color_manual(values = c("Loess" = "blue", "Polynomial" = "red", "GAM" = "green")) +  # Custom colors
  #   scale_linetype_manual(values = c("Loess" = "dashed", "Polynomial" = "dashed", "GAM" = "dashed")) +  # Custom line types
  #   theme_minimal() +
  #   labs(title = "Model Comparison", y = "NDVI", x = "Accumulated GDD", color = "Model", linetype = "Model")  # Add legend

  # ggsave(filename = paste0("plots/model_fitting/model_comparison_", id, ".png"),
  #                 plot = model_plot,
  #                 width = 6, 
  #                 height = 3,   
  #                 units = "in",
  #                 dpi = 150)
  
  #############################################################################

  # Peak NDVI and GDD
  peak_NDVI = max(df_subset$NDVI)
  peak_GDD = df_subset$Accumulated.GDD[which.max(df_subset$NDVI)]

  # Trying to capture gdd at/around peak - not working rn
  # ndvi_threshold = 0.9 * peak_NDVI
  # gdd_near_peak = df_subset$Accumulated.GDD[df_subset$NDVI >= ndvi_threshold]
  # time_near_peak = range(gdd_near_peak)  # GDD range where NDVI is near peak


  # First derivative (NDVI rate of change)
  first_derivative = diff(smoothed_values) / diff(df_subset$Accumulated.GDD)

  # Calculate the mean rate of change (mean of first derivative)
  mean_rate_of_change = mean(first_derivative, na.rm = TRUE)
  
  # Calculate the maximum rate of change (max of first derivative)
  max_rate_of_change = max(first_derivative, na.rm = TRUE)

  # Second derivative 
  # second_derivative = diff(first_derivative) /
  #   diff(df_subset$Accumulated.GDD[-c(1, length(df_subset$Accumulated.GDD))][-1])

  # Identify inflection points where second derivative is zero (or maybe close to zero)
  # Could also maybe be used to find where NDVI is at/near peak. where/when curve is changing
  # not working
  # inflection_points = df_subset$Accumulated.GDD[which(abs(second_derivative) < 1e-6)]
  # inflection_points_str = paste(inflection_points, collapse = ", ")
  
  
  # Calculate the area under the curve
  loess_fun = function(x) predict(loess_model, 
                                  newdata = data.frame(Accumulated.GDD = x))
  auc = integrate(loess_fun, lower = min(df_subset$Accumulated.GDD), 
                  upper = max(df_subset$Accumulated.GDD))$value
  
  # Store results in list
  results_list[[as.character(id)]] = list(
    ID = as.character(id),
    Peak_NDVI = peak_NDVI,
    Peak_GDD = peak_GDD,
    AUC = auc,
    Mean_Rate_Of_Change = mean_rate_of_change,
    Max_Rate_Of_Change = max_rate_of_change
  )
  
  ##### plot NDVI over time #########################################################
   
  # # format date to add as secondary x axis 
  # df_subset$posix_dates = as.POSIXct(paste0("2023-", gsub("\\.", "-", 
  #                                    df_subset$Date_nir)), format = "%Y-%m-%d")
  #                                   
  # # Create the plot for the current ID with loess line
  # 
  # ndvi_plot = ggplot(data = df_subset) +
  #   geom_smooth(aes(x = Accumulated.GDD, y = NDVI), method = "loess", color = "mediumpurple3", size = 0.8) +
  #   geom_point(aes(x = Accumulated.GDD, y = NDVI)) +
  #   labs(title = paste("NDVI vs GDD for", id),
  #        x = "Accumulated GDD",
  #        y = "NDVI") +
  #   scale_x_continuous(
  #     name = "Accumulated GDD",
  #     breaks = seq(0, 2000, by = 200),
  #     sec.axis = sec_axis(
  #       trans = ~ .,  # Use the same scale as Accumulated.GDD
  #       breaks = df_subset$Accumulated.GDD,  # Set breaks corresponding to the actual GDD values
  #       labels = format(df_subset$posix_dates, "%b %d")  # Convert dates to a readable format
  #     )
  #   ) +
  #   theme_bw() +
  #   theme(
  #     axis.title.x.top = element_text(size = 10),  # Secondary axis title at the top
  #     axis.text.x.top = element_text(size = 8),   # Text for the secondary axis on top
  #     axis.title.x.bottom = element_text(size = 10),  # Bottom axis title
  #     axis.text.x.bottom = element_text(size = 10),  # Bottom axis text
  #     plot.title = element_text(hjust = 0.5)  # Center the main plot title
  #   )
  # 
  # # Show the plot for the current ID
  # print(ndvi_plot)
  # 
  # # Pause after the first iteration and ask for confirmation to continue
  # if (id == unique_ids[1]) {
  #   response = readline(prompt = "Does the plot look okay? (y/n): ")
  #   if (tolower(response) != "y") {
  #     cat("Adjust the plot and try again :) \n")
  #     break  # exit loop if plot is bad
  #   }
  # }
#   
#   # Save the plot to a file - commenting out while adjusting what we want to extract
#   ggsave(filename = paste0("plots/NDVI_", id, ".png"),
#          plot = plot,
#          width = 6,    # Moderate width (in inches)
#          height = 3,   # Moderate height (in inches)
#          units = "in", 
#          dpi = 150)   

}

# Summary of overall model fit (average RMSE and R-squared across all iterations)
model_summary = model_comparison_all %>%
  group_by(Model) %>%
  summarise(
    Mean_RMSE = mean(RMSE),
    Mean_R_squared = mean(R_squared),
    Min_RMSE = min(RMSE),
    Min_R_squared = min(R_squared),
    Max_RMSE = max(RMSE),
    Max_R_squared = max(R_squared),
    SD_RMSE = sd(RMSE),
    SD_R_squared = sd(R_squared)
  )

# Save to csv
write.csv(model_summary, "outputs/NDVI_model_summary.csv", row.names = FALSE)
write.csv(model_comparison_all, "outputs/NDVI_model_summary_allPlots.csv", row.names = FALSE)


# Convert the list of results into a data frame
results_df = do.call(rbind, lapply(results_list, function(x) {
  data.frame(ID = x$ID,
             Peak_NDVI = x$Peak_NDVI,
             Peak_GDD = x$Peak_GDD,
             Area_Under_Curve = x$AUC,
             Mean_Rate_Of_Change = x$Mean_Rate_Of_Change,
             Max_Rate_Of_Change = x$Max_Rate_Of_Change,
             stringsAsFactors = FALSE)
}))

# Remove duplicate rows and write the results to a CSV file
#results_df = unique(results_df)
write.csv(results_df, "outputs/NDVI_extract.csv", row.names = FALSE)



##### possible to dos or improvements ###############################################

# Quality control step - compare mean to median? Could do in data combine/clean script
# Try any additional models to fit curve?
    # Double logistic, cubic spline?
# Fix errors with calculating derivatives and inflection points - is it even useful?
# Plot any other vegetative indices?

