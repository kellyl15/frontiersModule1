# Phenotype prediction using multiple vegetative indices to predict time to pollen. 
# Random forest model using caret package

# clear environment
rm(list=ls())

# load packages
library("randomForest")
library("ggplot2")
library("dplyr")
library("caret")

# read in data
vegData = read.csv("data/MI23 and FL22 VI Accession Flowering StandCount.csv")

# remove 1st and last dates
vegData = vegData %>%
  filter(!near(Accumulated.GDD, 498.4187, tol = 0.0001),
         !near(Accumulated.GDD, 1868.45, tol = 0.0001))

# divide data set into training and testing ########################################
# using 70% to train, can change

set.seed(5292013)

# list of unique ids - to group rows by plot
vegDataUnique = unique(vegData$id) 

trainingData = sample(x=vegDataUnique, size=ceiling(length(vegDataUnique)*0.7))

# create indicies of train & test rows
trainIndex = which(vegData$id %in% trainingData)
testIndex = setdiff(seq_len(nrow(vegData)), trainIndex)

# use k fold cross validation
trainControl = trainControl(method="cv", number=10)

# build model ######################################################################
rfModel = train(FL.GDD.to.50..Pollen ~ NDVI + NDRE + EVI + CI + RTVIcore + SAVI + 
                  MSAVI + Stand.Count, 
                data = vegData[trainIndex, ], method = "rf", trControl = trainControl)

# rfModel = randomForest(formula = FL.GDD.to.50..Pollen ~ NDVI + NDRE + EVI +  CI +
#                          RTVIcore + SAVI + MSAVI,
#                        data = vegData[trainIndex, ],
#                        ntree=500, importance=TRUE)

# predict using test data
#predictions = predict(rfModel, newdata=vegData[testIndex, ], type="response")

# predict using test data (using the caret with RF)
predictions = predict(rfModel, newdata = vegData[testIndex, ], type = "raw")


# evaluate model using rmse and r2 model ###########################################

# actual values
actualValues = vegData$FL.GDD.to.50..Pollen[testIndex]

# RMSE 
rmse = sqrt(mean((actualValues - predictions)^2))   

# R squared
r2 = cor(actualValues, predictions)^2  

# Create text annotation for plot
metricsText = paste0("RMSE: ", round(rmse, 2), "\nR²: ", round(r2, 3))

# plot actual vs predicted
actualVpredicted = ggplot(data = NULL, aes(x = actualValues, y = predictions)) +
  geom_point(color = "purple", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, max(actualValues), by = 200)) +
  scale_y_continuous(breaks = seq(0, max(predictions), by = 200)) +
  labs(x = "Actual FL GDD to 50% Pollen", y = "Predicted", 
       title = "Actual vs Predicted Values") +
  annotate("text", x = max(actualValues) * 0.95, y = max(predictions) * 0.85, 
           label = metricsText, hjust = 0.5, size = 3, color = "black") +
  theme_bw()
plot(actualVpredicted)
# ggsave("plots/actual_vs_predicted_initialRF_CV.png", actualVpredicted)
