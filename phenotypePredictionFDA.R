# Phenotype prediction using NDVI over time as a predictor of Time to Pollen

# clear environment
rm(list=ls())

# load packages
library("fda")
library("dplyr")
library("ggplot2")
library("tibble")
library("randomForest")
library("caret")

# read in data
vegData = read.csv("data/MI23 and FL22 VI Accession Flowering StandCount.csv")

# remove last date - some plots NDVI plateaued or were trending upwards
vegData = vegData %>%
  filter(!near(Accumulated.GDD, 1868.45, tol = 0.0001))

# prepare matrices for FDA #########################################################
ndvi_wide = reshape(
  vegData[, c("Accession", "Accumulated.GDD", "NDVI")],
  idvar = "Accession",
  timevar = "Accumulated.GDD",
  direction = "wide"
)

# Set rownames as Accession and then remove the Accession column
rownames(ndvi_wide) = ndvi_wide$Accession
ndvi_wide$Accession = NULL

# Convert to matrix
ndvi_matrix = as.matrix(ndvi_wide)

time_points = as.numeric(gsub("NDVI_", "", colnames(ndvi_matrix)))  # time values
time_points = sort(unique(vegData$Accumulated.GDD))

# create a B-spline basis ##########################################################
nbasis = 4  # number of basis functions - only have 4 time points
time_range = range(time_points)
basis = create.bspline.basis(rangeval = time_range, nbasis = nbasis, norder = 3)

# smooth curves for each plot
fd_list = smooth.basis(argvals = time_points,
                       y = t(ndvi_matrix),
                       fdParobj = fdPar(basis))$fd

# Plot smoothed curves
plot(fd_list, main = "Smoothed NDVI curves for each plot")

# functional PCA to extract key functional modes ###################################
fpca = pca.fd(fd_list, nharm = 4)
fPCAplot = plot(fpca$harmonics, main = "Functional Principal Components")

# Add a legend to the plot
legend("top", legend = paste("PC", 1:4), 
       col = 1:4, lty = 1, 
       lwd = 2, cex = 0.8,
       ncol = 2)

# get scores and add to flowering data
scores = fpca$scores
colnames(scores) = paste0("fPC", 1:ncol(scores))

scores_df = data.frame(Accession = rownames(ndvi_matrix), scores, stringsAsFactors = FALSE)
veg_copy = vegData[!duplicated(vegData$Accession), c("Accession", "FL.GDD.to.50..Pollen")]
scores_df = merge(scores_df, veg_copy, by = "Accession", all.x = TRUE)

# Remove data where FL.GDD.to.50..Pollen is above 1750 - since we were not using last GDD as predictor
scores_df = scores_df %>% filter(FL.GDD.to.50..Pollen <= 1750)

# train random forest on these scores ##############################################
# Create training and test data
set.seed(5292013)
train_ids = sample(scores_df$Accession, size = ceiling(0.7 * nrow(scores_df)))
train_df = filter(scores_df, Accession %in% train_ids)
test_df = filter(scores_df, !Accession %in% train_ids)

# Remove the 'Accession' column before fitting the model
train_df_noAccession = dplyr::select(train_df, -Accession)
test_df_noAccession = dplyr::select(test_df, -Accession)

# Ensure FL.GDD.to.50..Pollen is numeric
train_df_noAccession$FL.GDD.to.50..Pollen = as.numeric(train_df_noAccession$FL.GDD.to.50..Pollen)
test_df_noAccession$FL.GDD.to.50..Pollen = as.numeric(test_df_noAccession$FL.GDD.to.50..Pollen)

# Train random forest model
rf_fd = randomForest(FL.GDD.to.50..Pollen ~ ., 
                      data = train_df_noAccession, 
                      ntree = 500, 
                      importance = TRUE)

# Predict and evaluate
pred_fd = predict(rf_fd, newdata = test_df_noAccession)

# RMSE and R² Evaluation
rmse_fd = sqrt(mean((test_df_noAccession$FL.GDD.to.50..Pollen - pred_fd)^2))
r2_fd = cor(test_df_noAccession$FL.GDD.to.50..Pollen, pred_fd)^2

cat("RMSE (FDA-based):", rmse_fd, "\n")
cat("R² (FDA-based):", r2_fd, "\n")

# Variable importance plot
varImpPlot(rf_fd)

# plot predicted vs. actual and residuals
plot_df = data.frame(
          Actual = test_df_noAccession$FL.GDD.to.50..Pollen,
          Predicted = pred_fd)
plot_df$Residuals = plot_df$Actual - plot_df$Predicted

actualVsPredPlot = ggplot(plot_df, aes(x = Actual, y = Predicted)) +
  geom_point(color = "purple") +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  labs(title = "Predicted vs. Actual (Random Forest - FDA)",
       x = "Actual FL GDD to 50% Pollen",
       y = "Predicted FL GDD to 50% Pollen") +
  theme_bw() +
  annotate("text", 
           x = min(plot_df$Actual) + 0.8 * diff(range(plot_df$Actual)), 
           y = max(plot_df$Predicted) - 0.1 * diff(range(plot_df$Predicted)), 
           label = paste("R² =", round(r2_fd, 2)),
           size = 5, color = "black")
plot(actualVsPredPlot)

# residual plot
ggplot(plot_df, aes(x = Predicted, y = Residuals)) +
  geom_point(color = "purple") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  labs(title = "Residuals vs. Predicted",
       x = "Predicted FL.GDD.to.50..Pollen",
       y = "Residuals (Actual - Predicted)") +
  theme_bw()
