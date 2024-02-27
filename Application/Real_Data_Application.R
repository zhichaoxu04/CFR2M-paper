# ----- Load the packages
# .libPaths("/home/zxu7/R/ubuntu/4.2.0")
library(GMMAT)
library(SIS)
library(tidyverse)
library(HDMT)
library(tidyr)
cat("Load Packages Completed", "\n")

cat("This result is using the first half of the dataset.", "\n")

### Use x.res.std for regression and estimation

# ------ Load the data from local or HPC

# load("S:/Rotation/PW/Cohort_Cleaning/RDA_033122.rData")
# load("S:/Rotation/PW/Cohort_Cleaning/RDA_033122_meta.rData")
load("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/Cohort_Cleaning/RDA_033122.rData")  # Change Rdata
load("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/Cohort_Cleaning/Top10PCs.rData")
source("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/Cohort_Cleaning/RDA_script/022624/RDA_Infer_101023.R") # Change Function 
cat("Data Loaded", "\n")


args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])   # Random Seed
aID2 <- as.numeric(args[2])   # iter.max
aID3 <- as.numeric(args[3])   # Outcome and exposure and FDR

seed <- aID1
iter.max <- aID2
OutExpDat <- tidyr::expand_grid(data.frame(outcome = c("SBP_adj", "HDL"),
                                           exposure = c("Age", "Sex")),
                                method = c("iSIS", "HDMT"),
                                FDR = c(TRUE,FALSE))
outcome <- OutExpDat$outcome[aID3]
exposure <- OutExpDat$exposure[aID3]
FDR <- OutExpDat$FDR[aID3]
method <- OutExpDat$method[aID3]

if(outcome == "SBP_adj"){
  PCsNum <- 1
}else{
  PCsNum <- 2
  }

PCs <- Top10PCs[[PCsNum]]
colnames(PCs) <- paste("PC", 1:10, sep = "")

if(FDR == TRUE){
  FDR_out <- "FDR"
}else if (FDR == FALSE){
  FDR_out <- "NOFDR"
}

cat("Random seed= ", seed, "\n")
cat("Iter.max = ", iter.max, "\n")

CovarName <- c("Age", "Sex", "Currsmk", "BMI", "Alcohol2")

cat("outcome is ", outcome, "\n")
cat("exposure is ", exposure, "\n")
cat("FDR is ", FDR, "\n")
cat("Method is ", method, "\n")

set.seed(seed)
data <- Save_033122 %>%  
  # dplyr::filter(cohort == "offspring") %>% 
  mutate(cohort = as.numeric(as.factor(cohort))) %>%
  dplyr::select(shareid, 
                all_of(CovarName),
                all_of(outcome), 
                starts_with("ID_")) %>% 
  filter(complete.cases(.))  %>% 
  sample_frac(size = 1, replace = F) 

# data[1:10, 1:10]

# cat("N Row of data = ", dim(data)[1], "\n") 
# cat("N Col of data = ", dim(data)[2], "\n") 
# cat("Last colnames of dataset:", tail(colnames(data)), "\n")
# cat("First 9 colnames of dataset:", head(colnames(data), 9), "\n")

# ----- 
n <- nrow(data)   # Sample Size
cat("Sample Size = ", n, "\n") 
X <- as.numeric(data %>% dplyr::select(all_of(exposure)) %>% pull())   # X is the exposure
M <- as.matrix(data %>% dplyr::select(starts_with("ID_")) )     # M is mediators, n * d
cat("First 3 meditors name is:", colnames(M)[1:3], "\n")

d <- ncol(M)   
cat("number of M is:", d, "\n")

Covar <- data %>% dplyr::select(all_of(CovarName), 
                                -all_of(exposure))  # Other 5 variables
Covar <- Covar %>% bind_cols(PCs)
cat("Covariates names are:", colnames(Covar), "\n")

Y <- as.numeric(data %>% dplyr::select(all_of(outcome)) %>% pull()) # Y is the outcome

cat("Started. Current time is: ", as.character(Sys.time()), "\n")

result <- R2(Y = Y, M = M, Covar = Covar, X = X, d = d, n = n, 
             iter.max = iter.max, nsis = NULL, first.half = TRUE, seed = seed,
             FDR = FDR, FDRCutoff = 0.2, method=method)
result$output

cat("Completed. Current time is: ", as.character(Sys.time()), "\n")


outputDir <- "/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/Cohort_Cleaning/RDA_script/022624/Result/"
setwd(outputDir)
save(result, file = paste0(outputDir, "CFOLS2_", aID1, "_", aID2, "_",method, "_",outcome, "_",exposure,"_",FDR_out, ".RData"))

