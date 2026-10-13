options(expressions = 50000)
#Here we will address the high-dimensional issue caused by repeated measurements
#of the predictors over time by means of dimension reduction. Methods applied will
#be canonical correlation analysis (CCA), Sliced inverse regression (SIR), and 
#Sliced Average Variance Estimation (SAVE), and partial least squares. Also,
#potentially a Bayesian neural network, projective pursuit regression, and the 
#whiteboard dimension reduction method. Also, maybe something involving a data manifold.
library(data.table)
library(here)
#Created by the file Data loader.R
frame <- readRDS(file=paste0(here(),"/frame"))

### Apply cca
#library(mixOmics)
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
# subframe_X <- frame[,c(1:5,500:501,6,17)]#XWindv
# frame_wide_X <- pivot_wider(subframe_X,
#                           values_from="XWindv",
#                           names_from="Tsec",
#                           id_cols="Time",
#                           id_expand=TRUE,
#                           names_prefix="XWindv")
# width_X <- grep("XWindv",colnames(frame_wide_X))
# frame_wide_X <- as.matrix(frame_wide_X[,width_X])
# 
# subframe_Y <- frame[,c(17:177,500,501,1)] #Keep only a few columns #Should we also be pulling in the EPR uncertainty columns and the feature columns here?
# #Check to see if there is data outside of the final setyr (FSy) rows 
# FSy_ind <- sort(unique(c(which(!is.na(subframe_Y$EPR_1)),which(!is.na(subframe_Y$EPR_51)),which(!is.na(subframe_Y$EPR_101)))))
# not_FSy_ind <- setdiff(1:nrow(subframe_Y),FSy_ind)
# #which(!is.na(subframe_Y[not_FSy_ind,1:161]),arr.ind=TRUE) #indicates row #377267, col 34 has stray data
# row_where_data_goes <- min(FSy_ind[which(FSy_ind>377267)]) #Put the stray data in the next FSy row
# subframe_Y$EPR_34[row_where_data_goes] <- -2.7 #Put the stray data in the next FSy row
# subframe_Y$EPR_34[377267] <- NA
# frame_wide_Y <- subframe_Y[FSy_ind,1:(ncol(subframe_Y))] #Reduce responses to a single measurement per setyr--the FSy row
# rownames(frame_wide_Y) <- frame_wide_Y$setyr

# ── 1. Prep ──────────────────────────────────────────────────────────────────
dt <- as.data.table(frame)
setorder(dt, stamp)

# Identify EPR occasion rows (assuming NA in EPR cols on wind-only rows)
dt[, is_epr := !is.na(ShorelineID_doy)]

# Each EPR row starts a new "period"; subsequent rows inherit that period ID.
# XWindv rows get the ID of the EPR occasion that just preceded them.
dt[, period_id := cumsum(is_epr)]

# ── 2. Extract EPR occasions ──────────────────────────────────────────────────
epr_cols <- c("period_id", "stamp", paste0("EPR_", 1:161))
epr_wide <- dt[is_epr == TRUE, ..epr_cols]
setnames(epr_wide, "stamp", "epr_stamp")

# ── 3. Extract XWindv measurements (drop rows before the first EPR occasion) ──
Xwind_long <- dt[is_epr == FALSE, .(period_id, stamp, XWindv)]
Ywind_long <- dt[is_epr == FALSE, .(period_id, stamp, YWindv)]
aice_long <- dt[is_epr == FALSE, .(period_id, stamp, aice)]
Hsig_long <- dt[is_epr == FALSE, .(period_id, stamp, Hsig)]
RTpeak_long <- dt[is_epr == FALSE, .(period_id, stamp, RTpeak)]
Tm01_long <- dt[is_epr == FALSE, .(period_id, stamp, Tm01)]
Dspr_long <- dt[is_epr == FALSE, .(period_id, stamp, Dspr)]
WavePower_long <- dt[is_epr == FALSE, .(period_id, stamp, WavePower)]
AirT_long <- dt[is_epr == FALSE, .(period_id, stamp, AirT)]
SeaT_long <- dt[is_epr == FALSE, .(period_id, stamp, SeaT)]

# Number each wind measurement within its period (1, 2, 3, ...)
Xwind_long[, xwind_idx := seq_len(.N), by = period_id]
Ywind_long[, ywind_idx := seq_len(.N), by = period_id]
aice_long[, aice_idx := seq_len(.N), by = period_id]
Hsig_long[, hsig_idx := seq_len(.N), by = period_id]
RTpeak_long[, rtpeak_idx := seq_len(.N), by = period_id]
Tm01_long[, tm01_idx := seq_len(.N), by = period_id]
Dspr_long[, dspr_idx := seq_len(.N), by = period_id]
WavePower_long[, wavepower_idx := seq_len(.N), by = period_id]
AirT_long[, airt_idx := seq_len(.N), by = period_id]
SeaT_long[, seat_idx := seq_len(.N), by = period_id]

# ── 4. Pivot wind measurements wide using data.table::dcast ──────────────────
#    dcast is far more memory-efficient than pivot_wider for large data
Xwind_wide <- dcast(Xwind_long, period_id ~ xwind_idx,
                    value.var = "XWindv", sep = "_")
Ywind_wide <- dcast(Ywind_long, period_id ~ ywind_idx,
                    value.var = "YWindv", sep = "_")
aice_wide <- dcast(aice_long, period_id ~ aice_idx,
                    value.var = "aice", sep = "_")
Hsig_wide <- dcast(Hsig_long, period_id ~ hsig_idx,
                    value.var = "Hsig", sep = "_")
RTpeak_wide <- dcast(RTpeak_long, period_id ~ rtpeak_idx,
                    value.var = "RTpeak", sep = "_")
Tm01_wide <- dcast(Tm01_long, period_id ~ tm01_idx,
                    value.var = "Tm01", sep = "_")
Dspr_wide <- dcast(Dspr_long, period_id ~ dspr_idx,
                    value.var = "Dspr", sep = "_")
WavePower_wide <- dcast(WavePower_long, period_id ~ wavepower_idx,
                    value.var = "WavePower", sep = "_")
AirT_wide <- dcast(AirT_long, period_id ~ airt_idx,
                    value.var = "AirT", sep = "_")
SeaT_wide <- dcast(SeaT_long, period_id ~ seat_idx,
                    value.var = "SeaT", sep = "_")

# Rename columns: "1" -> "XWindv_1", "2" -> "XWindv_2", etc.
Xwind_cols <- setdiff(names(Xwind_wide), "period_id")
setnames(Xwind_wide, Xwind_cols, paste0("XWindv_", Xwind_cols))
Ywind_cols <- setdiff(names(Ywind_wide), "period_id")
setnames(Ywind_wide, Ywind_cols, paste0("YWindv_", Ywind_cols))
aice_cols <- setdiff(names(aice_wide), "period_id")
setnames(aice_wide, aice_cols, paste0("aice_", aice_cols))
Hsig_cols <- setdiff(names(Hsig_wide), "period_id")
setnames(Hsig_wide, Hsig_cols, paste0("Hsig_", Hsig_cols))
RTpeak_cols <- setdiff(names(RTpeak_wide), "period_id")
setnames(RTpeak_wide, RTpeak_cols, paste0("RTpeak_", RTpeak_cols))
Tm01_cols <- setdiff(names(Tm01_wide), "period_id")
setnames(Tm01_wide, Tm01_cols, paste0("Tm01_", Tm01_cols))
Dspr_cols <- setdiff(names(Dspr_wide), "period_id")
setnames(Dspr_wide, Dspr_cols, paste0("Dspr_", Dspr_cols))
WavePower_cols <- setdiff(names(WavePower_wide), "period_id")
setnames(WavePower_wide, WavePower_cols, paste0("WavePower_", WavePower_cols))
AirT_cols <- setdiff(names(AirT_wide), "period_id")
setnames(AirT_wide, AirT_cols, paste0("AirT_", AirT_cols))
SeaT_cols <- setdiff(names(SeaT_wide), "period_id")
setnames(SeaT_wide, SeaT_cols, paste0("SeaT_", SeaT_cols))

# ── 5. Join EPR occasions with their trailing wind measurements ───────────────
result <- merge(epr_wide, Xwind_wide, by = "period_id", all.x = TRUE)
result <- merge(result, Ywind_wide, by = "period_id", all.x = TRUE)
result <- merge(result, aice_wide, by = "period_id", all.x = TRUE)
result <- merge(result, Hsig_wide, by = "period_id", all.x = TRUE)
result <- merge(result, RTpeak_wide, by = "period_id", all.x = TRUE)
result <- merge(result, Tm01_wide, by = "period_id", all.x = TRUE)
result <- merge(result, Dspr_wide, by = "period_id", all.x = TRUE)
result <- merge(result, WavePower_wide, by = "period_id", all.x = TRUE)
result <- merge(result, AirT_wide, by = "period_id", all.x = TRUE)
result <- merge(result, SeaT_wide, by = "period_id", all.x = TRUE)
rm(epr_wide)
rm(Xwind_wide)
rm(Ywind_wide)
rm(aice_wide)
rm(Hsig_wide)
rm(RTpeak_wide)
rm(Tm01_wide)
rm(Dspr_wide)
rm(WavePower_wide)
rm(AirT_wide)
rm(SeaT_wide)

# Optional: drop the bookkeeping column
result[, period_id := NULL]

#PLS
Y_cols <- grep("EPR_",colnames(result))
Y_mat <- as.matrix(result[,..Y_cols])
Y_mat <- scale(Y_mat,center=TRUE,scale=TRUE)
XWindv_cols <- grep("XWindv_",colnames(result))
XWindv_mat <- as.matrix(result[,..XWindv_cols])
XWindv_mat <- scale(XWindv_mat,center=TRUE,scale=TRUE)
YWindv_cols <- grep("YWindv_",colnames(result))
YWindv_mat <- as.matrix(result[,..YWindv_cols])
YWindv_mat <- scale(YWindv_mat,center=TRUE,scale=TRUE)
aice_cols <- grep("aice_",colnames(result))
aice_mat <- as.matrix(result[,..aice_cols])
aice_mat <- scale(aice_mat,center=TRUE,scale=TRUE)
Hsig_cols <- grep("Hsig_",colnames(result))
Hsig_mat <- as.matrix(result[,..Hsig_cols])
Hsig_mat <- scale(Hsig_mat,center=TRUE,scale=TRUE)
RTpeak_cols <- grep("RTpeak_",colnames(result))
RTpeak_mat <- as.matrix(result[,..RTpeak_cols])
RTpeak_mat <- scale(RTpeak_mat,center=TRUE,scale=TRUE)
Tm01_cols <- grep("Tm01_",colnames(result))
Tm01_mat <- as.matrix(result[,..Tm01_cols])
Tm01_mat <- scale(Tm01_mat,center=TRUE,scale=TRUE)
Dspr_cols <- grep("Dspr_",colnames(result))
Dspr_mat <- as.matrix(result[,..Dspr_cols])
Dspr_mat <- scale(Dspr_mat,center=TRUE,scale=TRUE)
WavePower_cols <- grep("WavePower_",colnames(result))
WavePower_mat <- as.matrix(result[,..WavePower_cols])
WavePower_mat <- scale(WavePower_mat,center=TRUE,scale=TRUE)
AirT_cols <- grep("AirT_",colnames(result))
AirT_mat <- as.matrix(result[,..AirT_cols])
AirT_mat <- scale(AirT_mat,center=TRUE,scale=TRUE)
SeaT_cols <- grep("SeaT_",colnames(result))
SeaT_mat <- as.matrix(result[,..SeaT_cols])
SeaT_mat <- scale(SeaT_mat,center=TRUE,scale=TRUE)
stamp <- result$epr_stamp
rm(result)

nco <- 1
XWindv_pls <- pls(XWindv_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(XWindv_pls$X) <- paste0("XWindv_",1:nco)
YWindv_pls <- pls(YWindv_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(YWindv_pls$X) <- paste0("YWindv_",1:nco)
aic_pls <- pls(aice_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(aic_pls$X) <- paste0("aice_",1:nco)
Hsig_pls <- pls(Hsig_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(Hsig_pls$X) <- paste0("Hsig_",1:nco)
RTpeak_pls <- pls(RTpeak_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(RTpeak_pls$X) <- paste0("RTpeak_",1:nco)
Tm01_pls <- pls(Tm01_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(Tm01_pls$X) <- paste0("Tm01_",1:nco)
Dspr_pls <- pls(Dspr_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(Dspr_pls$X) <- paste0("Dspr_",1:nco)
WavePower_pls <- pls(WavePower_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(WavePower_pls$X) <- paste0("WavePower_",1:nco)
AirT_pls <- pls(AirT_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(AirT_pls$X) <- paste0("AirT_",1:nco)
SeaT_pls <- pls(SeaT_mat,Y_mat,ncomp=nco,scale=FALSE)$variates
colnames(SeaT_pls$X) <- paste0("SeaT_",1:nco)

df_pls <- data.frame(cbind(stamp,Y_mat,XWindv_pls$X,YWindv_pls$X,aic_pls$X,Hsig_pls$X,
               RTpeak_pls$X,Tm01_pls$X,Dspr_pls$X,WavePower_pls$X,
               AirT_pls$X,SeaT_pls$X))

saveRDS(df_pls,file=paste0(here(),"/df_pls"))

#Is a linear model defensible?
X_pls <- data.frame(cbind(XWindv_pls$X,YWindv_pls$X,aic_pls$X,Hsig_pls$X,
                 RTpeak_pls$X,Tm01_pls$X,Dspr_pls$X,WavePower_pls$X,
                 AirT_pls$X,SeaT_pls$X))
png(paste0(here(),"/eda_plots/corrplot.png"))
corrplot::corrplot(cor(cbind(X_pls,Y=apply(Y_mat,MARGIN=1,FUN=mean,na.rm=TRUE))))
dev.off()
all_numbers <- 1:(nco * 10)
groups <- split(all_numbers, ceiling(seq_along(all_numbers) / 5))
for(i in 1:(10*nco/5)) {
  for(j in 1:161) {
    png(paste0(here(),"/eda_plots/scatterplot_i=",i,"_j=",j,".png"))
    pairs(cbind(X_pls[,groups[[i]]],Y=Y_mat[,j]),pch=16,cex=0.5)
    dev.off()
  }
}

#What probability model is defensible?
for(j in 1:161){
  png(paste0(here(),"/eda_plots/histogram_j=",j,".png"))
  hist(Y_mat[,j])
  dev.off()
}

#Do PCA also as a benchmark comparison #TODO finish this section
# library(pcaMethods)
# XWindv_mat[which(is.nan(XWindv_mat))] <- NA
# XWindv_nas <- apply(XWindv_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# XWindv_mat <- XWindv_mat[,-which(XWindv_nas)]
# XWindv_pca <- pcaMethods::pca(XWindv_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# YWindv_mat[which(is.nan(YWindv_mat))] <- NA
# YWindv_nas <- apply(YWindv_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# YWindv_mat <- YWindv_mat[,-which(YWindv_nas)]
# YWindv_pca <- pcaMethods::pca(YWindv_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# aice_mat[which(is.nan(aice_mat))] <- NA
# aice_nas <- apply(aice_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# aice_mat <- aice_mat[,-which(aice_nas)]
# aice_pca <- pcaMethods::pca(aice_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# Hsig_mat[which(is.nan(Hsig_mat))] <- NA
# Hsig_nas <- apply(Hsig_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# Hsig_mat <- Hsig_mat[,-which(Hsig_nas)]
# Hsig_pca <- pcaMethods::pca(Hsig_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# RTpeak_mat[which(is.nan(RTpeak_mat))] <- NA
# RTpeak_nas <- apply(RTpeak_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# RTpeak_mat <- RTpeak_mat[,-which(RTpeak_nas)]
# RTpeak_pca <- pcaMethods::pca(RTpeak_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# Tm01_mat[which(is.nan(Tm01_mat))] <- NA
# Tm01_nas <- apply(Tm01_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# Tm01_mat <- Tm01_mat[,-which(Tm01_nas)]
# Tm01_pca <- pcaMethods::pca(Tm01_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# Dspr_mat[which(is.nan(Dspr_mat))] <- NA
# Dspr_nas <- apply(Dspr_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# Dspr_mat <- Dspr_mat[,-which(Dspr_nas)]
# Dspr_pca <- pcaMethods::pca(Dspr_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# WavePower_mat[which(is.nan(WavePower_mat))] <- NA
# WavePower_nas <- apply(WavePower_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# WavePower_mat <- WavePower_mat[,-which(WavePower_nas)]
# WavePower_pca <- pcaMethods::pca(WavePower_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# AirT_mat[which(is.nan(AirT_mat))] <- NA
# AirT_nas <- apply(AirT_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# AirT_mat <- AirT_mat[,-which(AirT_nas)]
# AirT_pca <- pcaMethods::pca(AirT_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# SeaT_mat[which(is.nan(SeaT_mat))] <- NA
# SeaT_nas <- apply(SeaT_mat,MARGIN=2,FUN=\(x) all(is.na(x)))
# SeaT_mat <- SeaT_mat[,-which(SeaT_nas)]
# SeaT_pca <- pcaMethods::pca(SeaT_mat, method = "ppca", nPcs = nco,scale=c("none"))@scores
# 
# rm(XWindv_mat)
# rm(YWindv_mat)
# rm(aice_mat)
# rm(Hsig_mat)
# rm(RTpeak_mat)
# rm(Tm01_mat)
# rm(Dspr_mat)
# rm(WavePower_mat)
# rm(AirT_mat)
# rm(SeaT_mat)
# gc()