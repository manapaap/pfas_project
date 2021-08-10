# To be run on ACCRE, to obtain XGBdart hyperparameters


library(doParallel)
library(caret)
library(xgboost)

extGrid <-  expand.grid(max_depth = c(1:10), 
                        eta = c(1:5),
                        rate_drop = c(1:5),
                        skip_drop = c(1:5),
                        subsample = c(1:5),
                        colsample_bytree = c(1:5),
                        nrounds = c(1:250),
                        gamma = c(1:5),
                        min_child_weight = c(1:10))

ky_esab_PFAS_normal <- read.csv('/home/manapaap/slurm_test/pfas_data.csv')

set.seed(420)
inTraining <- createDataPartition(ky_esab_PFAS_normal$PFAS_detect, p = .8, list = FALSE)
pfas_training <- ky_esab_PFAS_normal[ inTraining,]
pfas_testing  <- ky_esab_PFAS_normal[-inTraining,]

trainfolds <- trainControl(method = 'repeatedcv', number = 5, repeats = 10)


set.seed(2112)
# This tells caret to use 10-fold CV to validate model
# Hold back 1/10 of the data, 10 times

# Set a random seed so that results are reproducible.
set.seed(random_seed) 
# This fits the model on our data using the specified tuning parameters and 10-fold CV
fitModel <- train(PFAS_detect ~ ., data = pfas_testing, 
                  method = "xgbDART", 
                  trControl = fitControl,
                  tuneGrid = extGrid),
                  verbose = FALSE)





