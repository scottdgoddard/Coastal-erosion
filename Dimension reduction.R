#Here we will address the high-dimensional issue caused by repeated measurements
#of the predictors over time by means of dimension reduction. Methods applied will
#be canonical correlation analysis (CCA), Sliced inverse regression (SIR), and 
#Sliced Average Variance Estimation (SAVE), and partial least squares. Also,
#potentially a Bayesian neural network, projective pursuit regression, and the 
#whiteboard dimension reduction method. Also, maybe something involving a data manifold.
library(tidyr)

#Created by the file Data loader.R
frame <- readRDS(file="/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Coastal erosion/Spring 2025/frame") #MacOS

### Apply cca
# library(mixOmics)
# subframe <- frame[,c(1:5,500:501,6,17)]#XWindv
# frame_wide <- pivot_wider(subframe,
#                           names_from="set_index",
#                           values_from="XWindv",
#                           id_cols="setyr",
#                           id_expand=TRUE,
#                           names_prefix="XWindv")
# .ind <- c(match(frame_wide[,"setyr"]$setyr,frame[,"setyr"]),nrow(frame))
# frame_wide <- cbind(frame_wide,frame[.ind[-1]-1,paste0("EPR_",1:161)])
# rcc_1 <- rcc(frame_wide[,1:30000],frame_wide[,124715:124875],ncomp=2,lambda1=10,lambda2=10)
# quantile(rcc_1$loadings$X,0.99)
#Persistent positive-definiteness problems; hyperparameter uncertainty

### Apply SIR to each response variable in sequence
# library(msir)
# msir(frame_wide[,1:3000],frame_wide[,124715])
#Takes a long time for a small subset of the data

### Apply SAVE to each response variable in sequence
#Standardize predictors, check for linearity and non-constant covariance, examine response vs. SAVE predictors for various choices of H

### Apply PLS to response and predictor matrices
library(pls)
subframe_X <- frame[,c(1:5,500:501,6,17)]#XWindv
frame_wide_X <- pivot_wider(subframe_X,
                          values_from="XWindv",
                          names_from="Tsec",
                          id_cols="Time",
                          id_expand=TRUE,
                          names_prefix="XWindv")
width_X <- grep("XWindv",colnames(frame_wide_X))
frame_wide_X <- as.matrix(frame_wide_X[,width_X])

subframe_Y <- frame[,c(17:177,500,501,1)] #Keep only a few columns #Should we also be pulling in the EPR uncertainty columns and the feature columns here?
#Check to see if there is data outside of the final setyr (FSy) rows 
FSy_ind <- sort(unique(c(which(!is.na(subframe_Y$EPR_1)),which(!is.na(subframe_Y$EPR_51)),which(!is.na(subframe_Y$EPR_101)))))
not_FSy_ind <- setdiff(1:nrow(subframe_Y),FSy_ind)
#which(!is.na(subframe_Y[not_FSy_ind,1:161]),arr.ind=TRUE) #indicates row #377267, col 34 has stray data
row_where_data_goes <- min(FSy_ind[which(FSy_ind>377267)]) #Put the stray data in the next FSy row
subframe_Y$EPR_34[row_where_data_goes] <- -2.7 #Put the stray data in the next FSy row
subframe_Y$EPR_34[377267] <- NA
frame_wide_Y <- subframe_Y[FSy_ind,1:(ncol(subframe_Y))] #Reduce responses to a single measurement per setyr--the FSy row
rownames(frame_wide_Y) <- frame_wide_Y$setyr

#Let's start small
small_frame_wide_X <- frame_wide_X[,1:50]
small_frame_wide_Y <- as.matrix(frame_wide_Y[,1:50])
set.seed(13)
pls_result <- plsr(small_frame_wide_Y~small_frame_wide_X,validation="CV",
                   segments=8,na.fail=na.omit,scale=TRUE)
summary(pls_result)
validationplot(pls_result)

#Scaling up
set.seed(13)
med_frame_wide_X <- frame_wide_X[,1:2970]
pls_result <- plsr(frame_wide_Y~med_frame_wide_X,validation="CV",
                   segments=8,scale=TRUE)
summary(pls_result)
#validationplot(pls_result)

