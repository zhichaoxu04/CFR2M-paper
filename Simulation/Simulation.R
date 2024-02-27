rm(list=ls())
# Load packages
library(foreach)
require(doSNOW)
library(parallel)
library(iterators)
library(snow)
library(tidyverse)
library(glmnet)
library(MASS)
library(Matrix)

# -----------------------------------
# Set the number of true mediators
args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])   # numcores
aID2 <- as.numeric(args[2])   # iter.max
aID3 <- as.numeric(args[3])  # Seed (default = 2022)
aID4 <- as.numeric(args[4])  # B
aID5 <- as.numeric(args[5])  # Large Sample for TrueR2
aID6 <- as.numeric(args[6])  # Settings
aID7 <- as.numeric(args[7])  # Whether M are correlated or independent given X # 1 2 3
aID8 <- as.numeric(args[8])  # Whether perform FDR
aID9 <- as.numeric(args[9])  # Correlation coef
aID10 <- as.numeric(args[10])  # Correlation coef

source("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/CoverP_Simulation/R2_Script/012924/R2_CP_012924_Cor.R")

if(aID7 == 1){
  CorM <- aID7
  cat("M are Independent", "\n")
}else if(aID7 == 2){
  CorM <- aID7
  cat("M are Independent given X", "\n")
}else{
  CorM <- aID7
  cat("M are Correlated", "\n")
}

if(aID8 == 1){
  FDR <- FALSE
  cat("Do NOT perform FDR", "\n")
}else if(aID8 == 2){
  FDR <- TRUE
  cat("Perfrom FDR", "\n")
}

if(aID5 == 1){
  LargeSam <- TRUE
  cat("Use Large Sample to compute TrueR2", "\n")
}else if(aID5 == 0){
  LargeSam <- FALSE
  cat("Not Use Large Sample to compute TrueR2", "\n")
}

numcores <- aID1
iter.max <- aID2
seed <- aID3
B2 <- aID4
mean <- 0
std <- 1.5
nsis <- NULL
tune <- "bic"
penalty <- "MCP"
corCoef1 <- aID9
corCoef2 <- aID10


cat("numcores= ", numcores, "\n")
cat("Iter.max = ", iter.max, "\n")
cat("seed = ", seed, "\n")
cat("tune = ", tune, "\n")
cat("penalty = ", penalty, "\n")
cat("B2 = ", B2, "\n")
cat("distribution = normal", "\n")
cat("mean = ", mean, "\n")
cat("std = ", std, "\n")
cat("corCoef1 = ", corCoef1, "\n")
cat("corCoef2 = ", corCoef2, "\n")

# -------------------------
summary <- data.frame()
summary <- rbind(summary, rep(0, 47))
colnames(summary) <- c("N", 
                       "CP_asym", "CI_asym", "STD_asym",
                       "TrueRsq", "MSE", "Bias", "STD_Est",
                       "SelectedM","FDRM",
                       "TrueM", "M1", "M2", "Noise", 
                       "SOS.true", "SOS", "SOS_bias", "SOS_MSE",
                       "AB", "AB_bias1", "AB_bias2", "AB2_bias1", "AB2_bias2", 
                       "AB_MSE1", "AB_MSE2", "AB2_MSE1", "AB2_MSE2", 
                       "prop","prop_bias1", "prop_bias2", "prop_MSE1", "prop_MSE2", 
                       "total","total_bias1", "total_bias2", "total_MSE1", "total_MSE2", 
                       "ratio", "ratio_bias1", "ratio_bias2","ratio_MSE1", "ratio_MSE2",
                       "TruePos", "FalsePos", "TimeUsed", "TimeSD",
                       "Scenario")

RunAnalysis <- function(){
  temResult <- R2_Sim_CoverProb_fixedAB(alpha = alpha, beta = beta, 
                                        r = 3, res.sd = 1, Perc = TrueM, p = p, p1 = p1, p2 = p2,
                                        N = SampleSize,
                                        B1 = 1, B2 = B2, numcores = numcores, 
                                        iter.max = iter.max, nsis = nsis, tune = tune, penalty = penalty, 
                                        CorM = CorM, corMat=corMat, FDR = FDR, FDRCutoff=0.2, LargeSam=LargeSam)
  # idx <- which(!is.na(as.numeric(temResult$R45)))
  # temResultCom <- temResult %>% slice(idx) %>% slice(1:200)
  temResultCom <- temResult %>% 
    filter(complete.cases(.)) %>%
    slice(1:200)
  return(temResultCom)
}

RunAnalysis2 <- function(){
  temResult <- R2_Sim_CoverProb_fixedAB(alpha = alpha, beta = beta, 
                                        r = 3, res.sd = 1, Perc = TrueM, p = p, p1 = p1, p2 = p2,
                                        N = SampleSize,
                                        B1 = 1, B2 = B2, numcores = numcores, 
                                        iter.max = iter.max, nsis = nsis, tune = tune, penalty = penalty, 
                                        CorM = CorM, corMat=corMat, FDR = FDR, FDRCutoff=0.2, LargeSam=LargeSam)
  # idx <- which(!is.na(as.numeric(temResult$R45)))
  # temResultCom <- temResult %>% slice(idx) %>% slice(1:200)
  numberB <- nrow(temResult)
  
  NArow <- temResult %>% 
    filter(!complete.cases(.)) %>% 
    nrow(.)
  
  if(NArow != 0){
    imputation <- temResult %>% 
      filter(!complete.cases(.)) %>% 
      mutate(R45=0, CI_width_asym=0, R45.true=0, indirect.true=0, diff.R45=0, diff_ori=0, SOS_bias=0, SOS=0, prop.true=0, total.true=0, 
             SOS.true=0, ratio.true=0, Rsq.YM.true=0, Rsq.YX.true=0, Rsq.YMX.true=0,
             pab1_FDR=0, pab1=0, pab1_0=0, pab1_1=0, pab1_2=0, pab1_noise=0, AB_bias1=0, prop_bias1=0, total_bias1=0, ratio_bias1=0,
             pab2_FDR=0, pab2=0, pab2_0=0, pab2_1=0, pab2_2=0, pab2_noise=0, AB_bias2=0, prop_bias2=0, total_bias2=0, ratio_bias2=0,
             std_asym=0, cover_asym=1, TruePos=-999, FalsePos=0, N=0, TimeUsed=0)
    
  }else{
    imputation <- temResult %>% 
      filter(!complete.cases(.))
  }
  
  
  temResultCom <- temResult %>% 
    filter(complete.cases(.)) %>% 
    bind_rows(imputation) %>% 
    slice(1:numberB)
    
  return(temResultCom)
}

if(aID6 == 1){
  TrueM <- 15
  p1 <- 0
  p2 <- 0
  p <- 1500

}else if(aID6 == 2){
  TrueM <- 150
  p1 <- 0
  p2 <- 0
  p <- 1500

}else if(aID6 == 3){
  TrueM <- 150
  p1 <- 1350
  p2 <- 0
  p <- 1500
  
}else if(aID6 == 4){
  TrueM <- 150
  p1 <- 0
  p2 <- 1350
  p <- 1500

}else if(aID6 == 5){
  TrueM <- 150
  p1 <- 150
  p2 <- 0
  p <- 1500

}else if(aID6 == 6){
  TrueM <- 150
  p1 <- 150
  p2 <- 150
  p <- 1500

}else if(aID6 == 7){
  TrueM <- 0
  p1 <- 20
  p2 <- 20
  p <- 1500

}else if(aID6 == 8){
  TrueM <- 5
  p1 <- 0
  p2 <- 0
  p <- 1500

}else if(aID6 == 9){
  TrueM <- 20
  p1 <- 0
  p2 <- 0
  p <- 1500

}else if(aID6 == 10){
  TrueM <- 20
  p1 <- 60
  p2 <- 0
  p <- 1500
  
}else if(aID6 == 11){
  TrueM <- 20
  p1 <- 0
  p2 <- 60
  p <- 1500

}else if(aID6 == 12){
  TrueM <- 20
  p1 <- 60
  p2 <- 60
  p <- 1500
  
}

set.seed(seed)
alpha <- c(rnorm(TrueM, mean = mean, sd = std),
           rep(0, p1),
           rnorm(p2, mean = mean, sd = std),
           rep(0, p - TrueM - p1 - p2))

beta <- c(rnorm(TrueM, mean = mean, sd = std),
          rnorm(p1, mean = mean, sd = std),
          rep(0, p - TrueM - p1))

corMat <- diag(p)
# corMatDat <- matrix(runif((TrueM+p1)^2, corCoef1, corCoef2), nrow = (TrueM+p1), ncol = (TrueM+p1))
corMatDat <- matrix(rnorm((TrueM+p1)^2, corCoef1, corCoef2), nrow = (TrueM+p1), ncol = (TrueM+p1))
corMatDat[upper.tri(corMatDat)] <- t(corMatDat)[upper.tri(corMatDat)]
diag(corMatDat) <- 1
corMatDat <- Matrix::nearPD(corMatDat, corr = TRUE, base.matrix = TRUE, maxit = 300)$mat
summary(corMatDat[upper.tri(corMatDat, diag = F)])
corMat[1:(TrueM+p1), 1:(TrueM+p1)] <- corMatDat
diag(corMat) <- 1


SampleSize <- 750
Scenario <- paste0("N_", SampleSize, "_Perc_", TrueM, "_", p1, "_", p2, "_", p - TrueM - p1 - p2, "_M_", mean, "_S_", std)
tem <- RunAnalysis2()
result_tem <- c(summary_tem(tem), Scenario)
summary <- rbind(summary, result_tem)
cat("Current time is: ", as.character(Sys.time()), "\n")

SampleSize <- 1500
Scenario <- paste0("N_", SampleSize, "_Perc_", TrueM, "_", p1, "_", p2, "_", p - TrueM - p1 - p2, "_M_", mean, "_S_", std)
tem <- RunAnalysis2()
result_tem <- c(summary_tem(tem), Scenario)
summary <- rbind(summary, result_tem)
cat("Current time is: ", as.character(Sys.time()), "\n")

SampleSize <- 3000
Scenario <- paste0("N_", SampleSize, "_Perc_", TrueM, "_", p1, "_", p2, "_", p - TrueM - p1 - p2, "_M_", mean, "_S_", std)
tem <- RunAnalysis2()
result_tem <- c(summary_tem(tem), Scenario)
summary <- rbind(summary, result_tem)
cat("Current time is: ", as.character(Sys.time()), "\n")

cat("Job Done!", "\n")

# ---------------------
# Remove Unnecessary Variabels
rm(TrueM, SampleSize, p1, p2, p, numcores, iter.max, nsis, B2,
   alpha, beta, tem, tune, penalty, mean, std,
   result_tem)

summary <- summary %>% dplyr::slice(-1)

outputDir <- "/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/CoverP_Simulation/R2_Script/012924/Result/"
setwd(outputDir)

save(list = ls(.GlobalEnv),
     file = paste0(outputDir, "Sim_", aID1, "_", aID2, "_", aID3, "_",aID4,"_",aID5,"_",aID6,"_",aID7,"_",aID8,"_",aID9,"_",aID10, ".RData"))   #### 

save(list = "summary",
     file = paste0(outputDir, "Sim_", aID1, "_", aID2, "_", aID3, "_",aID4,"_",aID5,"_",aID6,"_",aID7,"_",aID8,"_",aID9,"_",aID10, ".txt"))   #### 

