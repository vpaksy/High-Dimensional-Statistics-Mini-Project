## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Install and load required packages

  if (!requireNamespace("mixOmics", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("mixOmics")}

  if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")}

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")}

  if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")}

  library(mixOmics)
  library(glmnet)
  library(ggplot2)
  library(reshape2)
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Load the liver.toxicity dataset from mixOmics package

  rm(list=ls())
  set.seed(1)
  data("liver.toxicity")
  X <- as.data.frame(liver.toxicity$gene) # creates covariates (genes) design matrix
  
  dim(X)  # X is a 64 x 3116 dataframe (n=64, p=3116)
  n <- nrow(liver.toxicity$gene)  # n=64
  Y <- as.data.frame(liver.toxicity$clinic$ALP) # creates column of the continuous response variable (measured ALP)
  

##' Center the covariates data to have zero mean and unit variance
##' There is no need toscale the response variable
  
  X <- scale(X, center=TRUE, scale=TRUE)
  Y <- scale(Y, center = FALSE, scale=FALSE)

  corrX <- cor(X)
  corrX <- round(corrX,4)
  hist(corrX) # histogram of the pairwise correlations
  # from the histogram, we can see the genes are correlated
  
  corrX <- corrX[1:50,1:50]
  melted_corrX <- melt(corrX)
  head(melted_corrX)


  ggplot(data = melted_corrX, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0,
                                     size= 8,hjust =0))+
    coord_fixed()
    # heat map of the pairwise correlations for a subset of 50 genes
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Create functions to return R^2

  r_squared <- function(y, yhat) {
    ybar <- mean(y)  
    ss_tot <- sum((y - ybar)^2)  # total sum of squares
    ss_res <- sum((y - yhat)^2)  # residual sum of squares
    1 - (ss_res / ss_tot)
  }
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



##############################################################################################################################################################################
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##                              Lasso Regression on full dataset                               
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################################################################################################################################################################

##' We will first investigate what happens when we perform Lasso regression on the full dataset, then
##' perform Lasso on the training set so we can make evaluate its performance using the test set. 

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Apply the lasso regression to the liver toxicity data, with λ being
##' chosen by cross validation (CV)

  set.seed(11)
  lasso <- glmnet(X, Y, alpha=1)  # applies lasso regression with X and Y
  lasso_cv <- cv.glmnet(X,Y,alpha=1) # applies the cross validation (CV)
  lambda_lasso <- lasso_cv$lambda.min # chooses optimal value of λ
  lambda_lasso  # returns optimal λ = 5.704953
  
  par(mfrow = c(1, 2))
  
  plot(lasso_cv) 
  title("Lasso on full dataset", line = 2) # plots the CV output
    
  plot(lasso, xvar = "lambda")
  title("Lasso on full dataset", line = 2) # plots a trace plot for λ
    
  
##' Extract the lasso estimates of parameters β
  
  betahat_lasso <- coef(lasso,s=lambda_lasso) # extracts the lasso estimates
  betahat_lasso_nointercept <- as.data.frame(betahat_lasso[-1,])  # removes the intercept term
  betahat_lasso_nointercept # returns estimates of the coefficients
    

##' Find how many parameters were "chosen"
  
  betahat_lasso_final <- betahat_lasso_nointercept[betahat_lasso_nointercept[,1]!=0,] # collects non-zero estimates
  length(betahat_lasso_final) # we find that 28 parameters are "chosen"
  which(betahat_lasso_nointercept[,1]!=0) # finds indices of non-zero estimates
    # the indices we find are: 22,65,66,501,577,745,761 .. etc
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



##############################################################################################################################################################################
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##                               Adaptive Lasso Regression on full dataset                            
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################################################################################################################################################################

##' We will first investigate what happens when we perform Adaptive Lasso regression on the full dataset, 
##' then perform Adaptive Lasso on the training set so we can make evaluate its performance using the test set. 

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Perform an initial lasso regression to the liver toxicity data, with λ being
##' chosen by cross validation (CV)

  set.seed(11)
  initial_lasso <- glmnet(x = X, y = Y,alpha = 1)
  initial_lasso_cv <- cv.glmnet(x = X, y = Y,type.measure = "mse",alpha = 1) # applies the cross validation (CV) (type.measure: loss to use for cv.)
  lambda_initial_lasso <- initial_lasso_cv$lambda.min # chooses optimal value of λ
  lambda_initial_lasso
    # Note: this produces same results as our lasso regression earlier, but I repeat the process for clarity reason
  
##' Extract the lasso estimates of parameters β
  betahat_initial_lasso <- coef(initial_lasso,s=lambda_initial_lasso) # The intercept estimate should be dropped.
  best_initial_lasso_coef <- as.numeric(betahat_initial_lasso)[-1] # removes the intercept term
    
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Perform adaptive LASSO

  set.seed(11)
  Adlasso <- glmnet(x = X, y = Y,alpha = 1,penalty.factor = 1 / abs(best_initial_lasso_coef)) # performs adaptive lasso with weights applied which we calculated in the previous step
  Adlasso_cv <- cv.glmnet(x = X, y = Y,type.measure = "mse",alpha = 1,penalty.factor = 1 / abs(best_initial_lasso_coef),keep = TRUE) # applies the cross validation (CV)
  lambda_Adlasso <- Adlasso_cv$lambda.min # chooses optimal value of λ
  lambda_Adlasso # we find the optimal value of λ to be 5.566721
    
  par(mfrow = c(1, 2))
  
  plot(Adlasso_cv)
  title("Adaptive Lasso on full dataset", line = 2) # plots the CV output
    
  plot(Adlasso, xvar = "lambda")
  title("Adaptive Lasso on full dataset", line = 2) # produces a trace plot for λ
    # produces a trace plot for λ
  
##' Extract the lasso estimates of parameters β
  betahat_Adlasso <- coef(Adlasso, s = lambda_Adlasso) # extracts the adaptive lasso estimates
  betahat_Adlasso_nointercept <- as.data.frame(betahat_Adlasso[-1,])  # removes the intercept term
  betahat_Adlasso_nointercept
  betahat_Adlasso_final <- betahat_Adlasso_nointercept[betahat_Adlasso_nointercept[,1]!=0,]
  
  length(betahat_Adlasso_final) # we find that now only 16 parameters are "chosen"
  which(betahat_Adlasso_nointercept[,1]!=0)  # finds indices of non-zero estimates
    # the indices we find are: 22,65,501,745,761,768... etc
    # these are a subset of the predictors chosen in lasso 
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# Now we will investigate what happens when we split the data into training and testing sets

##############################################################################################################################################################################
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##                              Lasso Regression on split dataset                               
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################################################################################################################################################################

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Randomly partition the data into training and test sets.

  set.seed(123)
  train_rows <- sample(1:n, 0.7*n) # samples indices from 1 to 64
    
  X.train <- X[train_rows, ]
  X.test <- X[-train_rows, ]
  Y.train <- Y[train_rows]
  Y.test <- Y[-train_rows]
    # creates training and test sets by partitioning the data with 70% of the data 
    # in the training set, 30% in the test set
  
  
##' Apply lasso regression to the training data
##' Use fitted model to predict the response values in the test data
  
  lasso_train <- glmnet(X.train, Y.train, alpha=1) # applies lasso regression with training sets
  lasso_train_cv <- cv.glmnet(X.train, Y.train, alpha=1) # applies the cross validation (CV) to the training set
  lambda_lasso_train <- lasso_train_cv$lambda.min # chooses optimal value of λ
  lambda_lasso_train  # returns optimal λ = 0.3636383
  
  par(mfrow = c(1, 2))
  plot(lasso_train_cv) 
  title("Lasso on training dataset", line = 2) # plots the CV output
    
  plot(lasso_train, xvar = "lambda")
  title("Lasso on training dataset", line = 2) # plots a trace plot for λ
    
  
  betahat_lasso_train <- coef(lasso_train,s=lambda_lasso_train)  # extracts the lasso estimates
  betahat_lasso_train_nointercept <- as.data.frame(betahat_lasso_train[-1,])  # removes the intercept term
  betahat_lasso_train_nointercept # returns estimates of the coefficients
    
  
##' Find how many parameters were "chosen"
  betahat_lasso_train_final <- betahat_lasso_train_nointercept[betahat_lasso_train_nointercept[,1]!=0,]  # collects non-zero estimates
    
  length(betahat_lasso_train_final)  # we find that 8 parameters are "chosen"
  which(betahat_lasso_train_nointercept[,1]!=0)  # finds indices of non-zero estimates
    # the indices we find are: 501,761,901,1817,1845,2384,2819,3032
    # note that there is an overlap in coefficients chosen when we compare with 
    # the previous run of lasso on the full dataset
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                            
##############################################################################################################################################################################





##############################################################################################################################################################################
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##                              Adaptive Lasso Regression on split dataset                              
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################################################################################################################################################################

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Perform initial lasso regression on the training set

  set.seed(123)
  initial_lasso_train <- glmnet(x = X.train, y = Y.train,alpha = 1) # applies lasso regression with X and Y
  initial_lasso_train_cv <- cv.glmnet(x = X.train, y = Y.train,type.measure = "mse",alpha = 1) # applies the cross validation (CV) (type.measure: loss to use for cv.)
  lambda_initial_lasso_train <- initial_lasso_train_cv$lambda.min  # chooses optimal value of λ
     
##' Extract the lasso estimates of parameters β
  
  betahat_initial_lasso_train <- coef(initial_lasso_train,s=lambda_initial_lasso_train) # The intercept estimate should be dropped.
  best_initial_lasso_train_coef <- as.numeric(betahat_initial_lasso_train)[-1] # removes the intercept term
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Perform Adaptive lasso regression on the training set

  set.seed(123)
  Adlasso_train <- glmnet(X.train, Y.train, alpha=1, penalty.factor = 1 / abs(best_initial_lasso_train_coef)) # applies adaptive lasso regression with training sets
  Adlasso_train_cv <- cv.glmnet(X.train, Y.train, alpha=1, penalty.factor = 1 / abs(best_initial_lasso_train_coef)) # applies the cross validation (CV) to the training set
  lambda_Adlasso_train <- Adlasso_train_cv$lambda.min # chooses optimal value of λ
  lambda_Adlasso_train
  
  par(mfrow = c(1, 2))
  
  plot(Adlasso_train_cv) 
  title("Adaptive Lasso on training dataset", line = 2) # plots the CV output
    
  plot(Adlasso_train, xvar = "lambda")
  title("Adaptive Lasso on training dataset", line = 2) # plots a trace plot for λ
    
  
  betahat_Adlasso_train <- coef(Adlasso_train,s=lambda_Adlasso_train)  # extracts the lasso estimates
  betahat_Adlasso_train_nointercept <- as.data.frame(betahat_Adlasso_train[-1,]) # removes the intercept term
  betahat_Adlasso_train_nointercept
  betahat_Adlasso_train_final <- betahat_Adlasso_train_nointercept[betahat_Adlasso_train_nointercept[,1]!=0,] # collects non-zero estimates
    
  length(betahat_Adlasso_train_final) # we find that 4 parameters are "chosen"
  which(betahat_Adlasso_train_nointercept[,1]!=0)  # finds indices of non-zero estimates
    # the indices we find are: 501,761,1845,2819
    # this is a subset of the indices chosen by lasso
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##############################################################################################################################################################################




##############################################################################################################################################################################
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##                                      Comparing results                                 
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################################################################################################################################################################



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Calculate the mean squared error of predictions (MSPE) and R^2 for Lasso

  set.seed(123)
  Yhat_lasso <- predict(lasso_train,X.test,s=lambda_lasso_train) # uses the fitted model to predict responses for the test data
  MSPE_lasso <- mean((Y.test- Yhat_lasso)^2) # calculates the mean square prediction error
  r_squared_lasso <- r_squared(as.vector(Y.test), as.vector(Yhat_lasso))
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Calculate the mean squared error of predictions (MSPE), R^2 for Adaptive Lasso

  set.seed(123)
  Yhat_Adlasso <- predict(Adlasso_train,X.test,s=lambda_Adlasso_train) # uses the fitted model to predict responses for the test data
  MSPE_Adlasso <- mean((Y.test- Yhat_Adlasso)^2) # calculates the MSPE
  r_squared_Adlasso <- r_squared(as.vector(Y.test), as.vector(Yhat_Adlasso))
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Compare the number of selected features on full dataset
  cat("Number of features selected by Lasso for full dataset:", length(betahat_lasso_final) )
  cat("Number of features selected by Adaptive Lasso for full dataset:", length(betahat_Adlasso_final) )
  
  
##' Plot coefficients for Lasso and Adaptive Lasso
  par(mfrow = c(1, 2))  # Plot side-by-side
  plot(betahat_lasso_final, main = "Lasso Coefficients for full dataset", xlab = "Features", ylab = "Coefficient Value", col = "blue", pch = 16)
  plot(betahat_Adlasso_final, main = "Adaptive Lasso Coefficients for full dataset", xlab = "Features", ylab = "Coefficient Value", col = "red", pch = 16)
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##' Second, present results from regression on training dataset

##' Coefficients from Lasso regression on training dataset
  lasso_coefficients <- betahat_lasso_train[c(502, 762, 902, 1818, 1846, 2385, 2820, 3033),]
  print("Lasso Coefficients for training dataset:")
  print(lasso_coefficients)

##' Coefficients from Adaptive Lasso regression on training dataset
  adaptive_lasso_coefficients <- betahat_Adlasso_train[c(502,762,2820),]
  print("Adaptive Lasso Coefficients for training dataset:")
  print(adaptive_lasso_coefficients)

##' Compare the number of selected features 
  cat("Number of features selected by Lasso for training dataset:", length(betahat_lasso_train_final) )
  cat("Number of features selected by Adaptive Lasso for training dataset:", length(betahat_Adlasso_train_final) )

##' Plot coefficients for Lasso and Adaptive Lasso
  par(mfrow = c(1, 2))  # Plot side-by-side
  plot(betahat_lasso_train_final, main = "Lasso Coefficients for training dataset", xlab = "Features", ylab = "Coefficient Value", col = "blue", pch = 16)
  plot(betahat_Adlasso_train_final, main = "Adaptive Lasso Coefficients for training dataset", xlab = "Features", ylab = "Coefficient Value", col = "red", pch = 16)


##' Compare MSPE and R^2
  cat("MSPE for Lasso:", MSPE_lasso)
  cat("R^2 for Lasso:", r_squared_lasso)
    
  cat("MSPE for adaptive Lasso:", MSPE_Adlasso)
  cat("R^2 for Lasso:", r_squared_Adlasso)


##' Plot predicted vs actual values
  par(mfrow = c(1, 2))
  plot(Y.test, Yhat_lasso, main = "Lasso: Predicted vs Actual", xlab = "Actual", ylab = "Predicted", col = "blue", pch = 16)
  abline(0, 1, col = "red")
  plot(Y.test, Yhat_Adlasso, main = "Adaptive Lasso: Predicted vs Actual", xlab = "Actual", ylab = "Predicted", col = "blue", pch = 16)
  abline(0, 1, col = "red")
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###############################################################################################################################################################################
