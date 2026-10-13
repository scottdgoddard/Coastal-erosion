library(mice)
library(brms)
library(here)
df_pls <- readRDS(file=paste0(here(),"/df_pls"))
df_pls$stamp <- as.numeric(format(as.POSIXct(df_pls$stamp), "%Y"))
nco <- 1

##########################
### Perform imputation ###
##########################

#Commented out in order to not always re-run when I run this script
df_pls_mice <- mice(df_pls, m=5, method="pmm", seed=123) #TODO check that the arguments here make sense (including defaults)
df_pls_imp1 <- complete(df_pls_mice,1)
saveRDS(df_pls_imp1,file=paste0(here(),"/df_pls_imp1"))
df_pls_imp1 <- readRDS(file=paste0(here(),"/df_pls_imp1"))

################################
### 0.0 "Basic linear model" ###
################################
#Fit model
formula_0.0 <- as.formula(paste0("cbind(",paste(paste0("EPR_",1:161),collapse=", "),") ~ ",
                                 paste(paste0("XWindv_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("YWindv_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("aice_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("Hsig_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("RTpeak_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("Tm01_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("Dspr_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("WavePower_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("AirT_",1:nco),collapse=" + ")," + ",
                                 paste(paste0("SeaT_",1:nco),collapse=" + ")," + stamp"))
(fit_0.0 <- lm(formula_0.0,data=df_pls_imp1))

#Examine diagnostics
png(paste0(here(),"/diagnostics/residuals_0-0.png"))
plot(rstandard(fit_0.0)~fitted(fit_0.0)) #Residual plot
dev.off()
png(paste0(here(),"/diagnostics/normprobplot_0-0.png"))
qqnorm(rstandard(fit_0.0)); qqline(rstandard(fit_0.0)) #Normal probability plot
dev.off()
sink(paste0(here(),"/diagnostics/tests.txt"))
shapiro.test(rstandard(fit_0.0)) #Normality test
lmtest::dwtest(fit_0.0,order.by=df_pls_imp1$stamp) #Durbin-Watson test by year
lmtest::dwtest(c(fit_0.0)) #Durbin-Watson test serially
png(paste0(here(),"/diagnostics/lag_residuals.png"))
plot(c(rstandard(fit_0.0)),data.table::shift(c(rstandard(fit_0.0)),n=1,type="lag",fill=NA),
     xlab="residuals",ylab="lagged residuals") #Lagged residual vs. residual scatterplot
abline(lm(data.table::shift(c(rstandard(fit_0.0)),n=1,type="lag",fill=NA)~c(rstandard(fit_0.0))),
          lwd=3,col="purple4") 
dev.off()

#Summarize model fit
temp <- 0
for(i in 1:161) {
  temp <- temp + summary(fit_0.0)[[i]]$adj.r.squared
}
cat("Adj. R2: ",temp/161)
sink()

#####################################################
### 1.0 "Random intercepts" for year and transect ###
#####################################################
#Fit model
formula_1.0 <- bf(mvbind(EPR_1,EPR_2,EPR_3,EPR_4,EPR_5,EPR_6,EPR_7,EPR_8,EPR_9,
                         EPR_10,EPR_11,EPR_12,EPR_13,EPR_14,EPR_15,EPR_16,EPR_17,
                         EPR_18,EPR_19,EPR_20,EPR_21,EPR_22,EPR_23,EPR_24,EPR_25,
                         EPR_26,EPR_27,EPR_28,EPR_29,EPR_30,EPR_31,EPR_32,EPR_33,
                         EPR_34,EPR_35,EPR_36,EPR_37,EPR_38,EPR_39,EPR_40,EPR_41,
                         EPR_42,EPR_43,EPR_44,EPR_45,EPR_46,EPR_47,EPR_48,EPR_49,
                         EPR_50,EPR_51,EPR_52,EPR_53,EPR_54,EPR_55,EPR_56,EPR_57,
                         EPR_58,EPR_59,EPR_60,EPR_61,EPR_62,EPR_63,EPR_64,EPR_65,
                         EPR_66,EPR_67,EPR_68,EPR_69,EPR_70,EPR_71,EPR_72,EPR_73,
                         EPR_74,EPR_75,EPR_76,EPR_77,EPR_78,EPR_79,EPR_80,EPR_81,
                         EPR_82,EPR_83,EPR_84,EPR_85,EPR_86,EPR_87,EPR_88,EPR_89,
                         EPR_90,EPR_91,EPR_92,EPR_93,EPR_94,EPR_95,EPR_96,EPR_97,
                         EPR_98,EPR_99,EPR_100,EPR_101,EPR_102,EPR_103,EPR_104,
                         EPR_105,EPR_106,EPR_107,EPR_108,EPR_109,EPR_110,EPR_111,
                         EPR_112,EPR_113,EPR_114,EPR_115,EPR_116,EPR_117,EPR_118,
                         EPR_119,EPR_120,EPR_121,EPR_122,EPR_123,EPR_124,EPR_125,
                         EPR_126,EPR_127,EPR_128,EPR_129,EPR_130,EPR_131,EPR_132,
                         EPR_133,EPR_134,EPR_135,EPR_136,EPR_137,EPR_138,EPR_139,
                         EPR_140,EPR_141,EPR_142,EPR_143,EPR_144,EPR_145,EPR_146,
                         EPR_147,EPR_148,EPR_149,EPR_150,EPR_151,EPR_152,EPR_153,
                         EPR_154,EPR_155,EPR_156,EPR_157,EPR_158,EPR_159,EPR_160,
                         EPR_161)~XWindv_1 + YWindv_1 + aice_1 + Hsig_1 + RTpeak_1 + 
                    Tm01_1 + Dspr_1 + WavePower_1 + AirT_1 + SeaT_1 + (1|stamp))
(fit_1.0 <- brm(formula_1.0,data=df_pls_imp1,family=gaussian(),set_rescor=FALSE))
