# To be run on ACCRE, to obtain XGBdart hyperparameters


library(doParallel)
library(caret)
library(gbm)
library(klaR)
library(nnet)
library(caTools)
library(xgboost)

my.cluster <- parallel::makeCluster(14)

doParallel::registerDoParallel(cl = my.cluster)




ky_esab_PFAS_normal <- read.csv('pfas_data.csv')

ky_esab_PFAS_normal$PFAS_detect <- as.factor(ky_esab_PFAS_normal$PFAS_detect)


set.seed(420)
inTraining <- createDataPartition(ky_esab_PFAS_normal$PFAS_detect, p = .8, list = FALSE)
pfas_training <- ky_esab_PFAS_normal[ inTraining,]
pfas_testing  <- ky_esab_PFAS_normal[-inTraining,]

trainfolds <- trainControl(method = 'repeatedcv', number = 5, repeats = 10)

# GBM trials

extGrid <-  expand.grid(n.trees = c(1:20), # 9
                        interaction.depth = c(1:20), # 13 
                        shrinkage = c(1:10), # 1
                        n.minobsinnode = c(1:10)) # 8

set.seed(2112)
# This tells caret to use 10-fold CV to validate model
# Hold back 1/10 of the data, 10 times

# This fits the model on our data using the specified tuning parameters and 10-fold CV
fitModel_gbm <- train(PFAS_detect ~ ., data = pfas_training, 
                      method = "gbm", 
                      trControl = trainfolds,
                      tuneGrid = extGrid,
                      verbose = FALSE)

PFAS_predict_gb_norm <- predict(fitModel_gbm, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_gb_norm)


# Bayes trials

extGrid <-  expand.grid(fL = c(1:10), # 1
                        usekernel = c(1:10), # 1
                        adjust = c(1:10)) # 2


fitModel_tan <- train(PFAS_detect ~ ., data = pfas_training, 
                      method = "nb", 
                      trControl = trainfolds,
                      tuneGrid = extGrid,
                      verbose = FALSE)

PFAS_predict_by_norm <- predict(fitModel_tan, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_by_norm)


# Neural net trials


extGrid <-  expand.grid(size = c(1:10), # 3
                        decay = c(1:10)) # 1


fitModel_net <- train(PFAS_detect ~ ., data = pfas_training, 
                      method = "nnet", 
                      trControl = trainfolds,
                      tuneGrid = extGrid,
                      verbose = FALSE)

PFAS_predict_nn_norm <- predict(fitModel_net, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_nn_norm)


# Boosted logreg test

extGrid <-  expand.grid(nIter = c(1:20)) # 2


fitModel_lg <- train(PFAS_detect ~ ., data = pfas_training, 
                     method = "LogitBoost", 
                     trControl = trainfolds,
                     tuneGrid = extGrid,
                     verbose = FALSE)

PFAS_predict_lg_norm <- predict(fitModel_lg, newdata = pfas_testing)
PFAS_predict_lg_norm[is.na(PFAS_predict_lg_norm)] <- 0

table(pfas_testing$PFAS_detect, PFAS_predict_lg_norm)

# RFerns Test

extGrid <-  expand.grid(depth = c(1:20)) # 1

fitModel_Fe <- caret::train(PFAS_detect ~ .,
                            data = pfas_training,
                            method = 'rFerns',
                            tuneGrid = extGrid,
                            trControl = trainfolds)

PFAS_predict_Fe_norm <- predict(fitModel_Fe, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_Fe_norm)



# Random Forest test

extGrid <-  expand.grid(mtry = c(1:10), # 3
                        maxdepth = c(1:10)) # 1

fitModel_rf <- caret::train(PFAS_detect ~ . ,
                            data = pfas_training,
                            method = 'rfRules',
                            trControl = trainfolds,
                            tuneGrid = extGrid,
                            verbose = FALSE)

PFAS_predict_rf_norm <- predict(fitModel_rf, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_rf_norm)

# Exreme gradient boosting - some parameters

pfas_training$PFAS_detect <- pfas_training$PFAS_detect %>%
  as.factor() %>% as.numeric() %>% -1 %>%
  cut(breaks = c(-100, 0.5, 100), labels = c('no', 'yes'))

extGrid <-  expand.grid(max_depth = c(1:3), 
                        eta = 1,
                        rate_drop = 1,
                        skip_drop = 1,
                        subsample = c(1:3),
                        colsample_bytree = c(1:5),
                        nrounds = c(1:250),
                        gamma = c(1:5),
                        min_child_weight = c(1:5))

firModel_Xb <- caret::train(PFAS_detect ~ .,
                              data = pfas_training,
                              method = 'xgbDART',
                              tuneGrid = extGrid,
                              trControl = trainfolds)

PFAS_predict_Xb_norm <- predict(Xb_model_norm, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_Xb_norm)





























