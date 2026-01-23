## ODDESSI - Marginal Hazard Ratio for primary Cox model (time to relapse)
## Author: ORIGINAL CODE FROM RHIAN DANIEL 2020 PAPER (https://onlinelibrary.wiley.com/doi/10.1002/bimj.201900297)
    ## EDIT/REVISIONS MADE FOR ODDESSI TRIAL BY: Petrina Chu
    ## please excuse my code and annotations...
## Date/version: 18 Mar 2025

library("haven")
library("tidyverse")
library("survival")
library("survminer")
library("boot")
library("simsurv")

setwd("C:/Users/k2150120/OneDrive - King's College London/Documents/ODDESSI/out")
## read in rr_analysis .dta file
data1 <- read_dta("rr_set.dta")
  table(data1$arm)
  colnames(data1) <- c("pid", "cid", "catchment", "list_size", "imd", "arm", "DEM_17", "PH_03", "relap", "time1", "cen_lsize", "cen_imd", "catch1", "catch2", "catch3", "catch4", "catch6", "catch7", "catch8", "ph1", "house1", "house3", "house4")
    ## predictors of full medical record status: housing status, previous contact with secondary care mental health services
    ## catch5 (largest catchment area - design factor) and house2 (living with relatives/others) are reference groups. having previous contact is reference group
    ## contains original factor variables but i'm using indicators here ('catchX', etc.) to make sure i'm calling the right thing
  

########## USE RHIAN DANIEL ET AL EXAMPLE CODE TO FIND MARGINAL HAZARD RATIO FROM FULLY ADJUSTED MODEL ########
  ## Define function 'HR_new'. It's long!!
  HR_new <- function(data, indices, N_new_half){#N_new_half, survival times be simulated, i.e. m in the paper
    
    odd <- data[indices, ]
    cox.fit <- coxph(Surv(odd$time1, odd$relap==1) ~ arm + catch1 + catch2 + catch3 + catch4 + catch6 + catch7 + catch8 + cen_lsize + cen_imd + house1 + house3 + house4 + ph1, data=odd) #Conditional model here
    bh <- basehaz(cox.fit, centered = FALSE) #extract baseline hazard from model above
    
    ###eventtime of relapse in original dataset
    N <- nrow(odd)
    s1 <- matrix(0, nrow = nrow(bh), ncol = N)#all exposure, i,e, when all x==1 
    s0 <- matrix(0, nrow = nrow(bh), ncol = N)#all non exposure, i,e, when all x==0
    for (i in 1:N) {
      s1[, i] <- exp(-bh$hazard) ^ (exp(cox.fit$coef[1] + (cox.fit$coef[2] * odd$catch1[i]) + (cox.fit$coef[3] * odd$catch2[i]) + (cox.fit$coef[4] * odd$catch3[i]) + (cox.fit$coef[5] * odd$catch4[i]) + (cox.fit$coef[6] * odd$catch6[i]) + (cox.fit$coef[7] * odd$catch7[i]) + (cox.fit$coef[8] * odd$catch8[i]) + (cox.fit$coef[9] * odd$cen_lsize[i]) + (cox.fit$coef[10] * odd$cen_imd[i]) + (cox.fit$coef[11] * odd$house1[i]) + (cox.fit$coef[12] * odd$house3[i]) + (cox.fit$coef[13] * odd$house4[i]) + (cox.fit$coef[14] * odd$ph1[i]) ))
      s0[, i] <- exp(-bh$hazard) ^ (exp((cox.fit$coef[2] * odd$catch1[i]) + (cox.fit$coef[3] * odd$catch2[i]) + (cox.fit$coef[4] * odd$catch3[i]) + (cox.fit$coef[5] * odd$catch4[i]) + (cox.fit$coef[6] * odd$catch6[i]) + (cox.fit$coef[7] * odd$catch7[i]) + (cox.fit$coef[8] * odd$catch8[i]) + (cox.fit$coef[9] * odd$cen_lsize[i]) + (cox.fit$coef[10] * odd$cen_imd[i]) + (cox.fit$coef[11] * odd$house1[i]) + (cox.fit$coef[12] * odd$house3[i]) + (cox.fit$coef[13] * odd$house4[i]) + (cox.fit$coef[14] * odd$ph1[i]) ))
    }
    s1_cond <- vector(mode="numeric", length = nrow(bh))
    s0_cond <- vector(mode="numeric", length = nrow(bh))
    s1_Mean <- rowMeans(s1)
    s0_Mean <- rowMeans(s0)
    
    ## step 3 in section 4.2?
    s1_cond[1] <- s1_Mean[1]
    s0_cond[1] <- s0_Mean[1]
    for (j in 2:nrow(bh)) {
      s1_cond[j] <- s1_Mean[j] / s1_Mean[j - 1]
      s0_cond[j] <- s0_Mean[j] / s0_Mean[j - 1]
    }
    
    #### IN PRESENCE OF CENSORING, SECTION 4.2.4 IN PAPER
    ## STEP 4' IN SECTION 4.2.4
    odd_c <- odd
    odd_c$relap <- 1 - odd$relap #being censored is event of interest in 'pbc_c'
    ## STEP 5' IN SECTION 4.2.4
    cox.fit_c <- coxph(Surv(odd_c$time1, odd_c$relap ==1 ) ~ arm + catch1 + catch2 + catch3 + catch4 + catch6 + catch7 + catch8 + cen_lsize + cen_imd + house1 + house3 + house4 + ph1, odd_c)
       #conditional Cox model using same event time but new censoring status
    ## STEP 6' IN SECTION 4.2.4
    ##find baseline hazards, then find marginal causal survivor function for censoring
    bh_c <- basehaz(cox.fit_c, centered = FALSE)
    N_c <- nrow(odd_c)
    s1_c <- matrix(0, nrow = nrow(bh_c), ncol = N_c)#all exposure, i,e, when all x==1 
    s0_c <- matrix(0, nrow = nrow(bh_c), ncol = N)#all non exposure, i,e, when all x==0
    for (i in 1:N_c) {
      s1_c[, i] <- exp(-bh_c$hazard) ^ (exp(cox.fit_c$coef[1] + (cox.fit_c$coef[2] * odd_c$catch1[i]) + (cox.fit_c$coef[3] * odd_c$catch2[i]) + (cox.fit_c$coef[4] * odd_c$catch3[i]) + (cox.fit_c$coef[5] * odd_c$catch4[i]) + (cox.fit_c$coef[6] * odd_c$catch6[i]) + (cox.fit_c$coef[7] * odd_c$catch7[i]) + (cox.fit_c$coef[8] * odd_c$catch8[i]) + (cox.fit_c$coef[9] * odd_c$cen_lsize[i]) + (cox.fit_c$coef[10] * odd_c$cen_imd[i]) + (cox.fit_c$coef[11] * odd_c$house1[i]) + (cox.fit_c$coef[12] * odd_c$house3[i]) + (cox.fit_c$coef[13] * odd_c$house4[i]) +  (cox.fit_c$coef[14] * odd_c$ph1[i])))
      s0_c[, i] <- exp(-bh_c$hazard) ^ (exp((cox.fit_c$coef[2] * odd_c$catch1[i]) + (cox.fit_c$coef[3] * odd_c$catch2[i]) + (cox.fit_c$coef[4] * odd_c$catch3[i]) + (cox.fit_c$coef[5] * odd_c$catch4[i]) + (cox.fit_c$coef[6] * odd_c$catch6[i]) + (cox.fit_c$coef[7] * odd_c$catch7[i]) + (cox.fit_c$coef[8] * odd_c$catch8[i]) + (cox.fit_c$coef[9] * odd_c$cen_lsize[i]) + (cox.fit_c$coef[10] * odd_c$cen_imd[i]) + (cox.fit_c$coef[11] * odd_c$house1[i]) + (cox.fit_c$coef[12] * odd_c$house3[i]) + (cox.fit_c$coef[13] * odd_c$house4[i]) + (cox.fit_c$coef[14] * odd_c$ph1[i]) ))
    }
    s1_cond_c <- vector(mode = "numeric", length = nrow(bh_c))
    s0_cond_c <- vector(mode = "numeric", length = nrow(bh_c))
    s1_Mean_c <- rowMeans(s1_c)
    s0_Mean_c <- rowMeans(s0_c)
    s1_cond_c[1] <- s1_Mean_c[1]
    s0_cond_c[1] <- s0_Mean_c[1]
    for (j in 2:nrow(bh_c)) {
      s1_cond_c[j] <- s1_Mean_c[j] / s1_Mean_c[j-1]
      s0_cond_c[j] <- s0_Mean_c[j] / s0_Mean_c[j-1]
    }
    
    #############################Simulated Data 
    N_1 <- N_new_half#100#round(N/2),100,500,1000,2500,5000,10000,25000
    N_0 <- N_new_half#100#N-N_1,100,500,1000,2500,5000,10000
    ######Exposure
    #eventtime of dead
    newd1 <- matrix(0, nrow = nrow(bh), ncol = N_1)
    for (i in 1:N_1) {
      newd1[,i] <- rbinom(nrow(bh), 1, s1_cond)
    }
    
    index_event1_d <- apply(newd1, 2, function(x) which(x == 0)[1])
    eventtime1_d <- bh$time[index_event1_d]
    eventtime1_d[which(is.na(index_event1_d))] <- max(bh$time)
    status1_d <- rep(1, N_1)
    status1_d[which(is.na(index_event1_d))] <- 0
    #####censoring eventime as death time
    newd1_c <- matrix(0, nrow = nrow(bh_c), ncol = N_1)
    for (i in 1:N_1) {
      newd1_c[,i] <- rbinom(nrow(bh_c), 1, s1_cond_c)
    }
    
    index_event1_c <- apply(newd1_c, 2, function(x) which(x == 0)[1])
    eventtime1_c <- bh_c$time[index_event1_c]
    eventtime1_c[which(is.na(index_event1_c))] <- max(bh_c$time)
    status1_c <- rep(1, N_1)
    status1_c[which(is.na(index_event1_c))] <- 0
    ###############choose event time for simulated data (step 10' ?)
    eventtime1 <- vector(mode="numeric", length = length(eventtime1_c))
    status1 <- vector(mode = "numeric", length = nrow(bh))# 0=alive, 1=dead
    ## see 4.2.5 in paper. Choose the earlier of 1) simulated event time, or 2) simulated censoring time
    for (k in 1:length(eventtime1_c))
    {
      if (eventtime1_c[k] < eventtime1_d[k]){
        
        eventtime1[k] <- eventtime1_c[k]
        status1[k] <- 0#censoring
        
      }else if (eventtime1_c[k] > eventtime1_d[k] | eventtime1_c[k] == eventtime1_d[k]){
        
        eventtime1[k] <- eventtime1_d[k]
        status1[k] <- 1#relap
        
      }
    }
    ########UnExposure
    #eventtime of dead
    newd0 <- matrix(0, nrow = nrow(bh), ncol = N_0)
    for (i in 1:N_0) {
      newd0[,i] <- rbinom(nrow(bh), 1, s0_cond)
    }
    
    index_event0_d <- apply(newd0, 2, function(x) which(x == 0)[1])
    
    eventtime0_d <- bh$time[index_event0_d]
    eventtime0_d[which(is.na(index_event0_d))] <- max(bh$time)
    status0_d <- rep(1, N_0)
    status0_d[which(is.na(index_event0_d))] <- 1
    #######censoring eventime as death time
    newd0_c <- matrix(0, nrow = nrow(bh_c), ncol = N_0)
    for (i in 1:N_0) {
      newd0_c[,i] <- rbinom(nrow(bh_c), 1, s0_cond_c)
    }
    
    index_event0_c <- apply(newd0_c, 2, function(x) which(x == 0)[1])
    
    eventtime0_c <- bh_c$time[index_event0_c]
    eventtime0_c[which(is.na(index_event0_c))] <- max(bh_c$time)
    status0_c <- rep(1, N_0)
    status0_c[which(is.na(index_event0_c))] <- 1
    #####choose event time for simulated data
    eventtime0 <- vector(mode = "numeric", length = length(eventtime0_c))
    status0 <- vector(mode = "numeric", length = nrow(bh))# 0=alive, 1=dead
    for (k in 1:length(eventtime0_c))
    {
      if (eventtime0_c[k] < eventtime0_d[k]){
        
        eventtime0[k] <- eventtime0_c[k]
        status0[k] <- 0#censoring
        
      }else if (eventtime0_c[k] > eventtime0_d[k] | eventtime0_c[k] == eventtime0_d[k]){
        
        eventtime0[k] <- eventtime0_d[k]
        status0[k] <- 1#relap
        
      }
    }
    
    #####New data set and final cox model (Steps 5 and 6 of 4.2.1)
    newdata <- data.frame(rbind(cbind(eventtime1, status1, rep(1, N_1)), 
                                cbind(eventtime0, status0, rep(0, N_0))))
    names(newdata) <- c("eventtime", "status", "x")
    cox.fit.new <- coxph(Surv(eventtime, status) ~ x, newdata) 
    return(coef(cox.fit.new))
  }
  ## end of function!!
  
## create results matrix to save output later
results <- matrix(0, 3, 4)
  colnames(results) <- c("reps", "mhr", "ll", "ul")

sink("relapse_mhr_log.txt")    
######## USE FUNCTION 'HR_new' with varying '2m' simulations ######
set.seed(11022025)
  
## m=100, so 2m=200 replications. 
HR_new.boot <- boot(data1, HR_new, R = 1000, N_new_half = 100)
  ## calc 95% CIs using bootstrap percentile intervals
  HR_new.boot.ci <- boot.ci(HR_new.boot, conf = c(0.95), type = c("perc"))
print(HR_new.boot)
  results[1,1] <- 2*HR_new.boot[["call"]][["N_new_half"]]
  results[1,2] <- exp(HR_new.boot[["t0"]][["x"]])
print(HR_new.boot.ci)
  ci_vector <- HR_new.boot.ci[["percent"]]
  results[1,3] <- exp(ci_vector[4])
  results[1,4] <- exp(ci_vector[5])
  
## m=500, so 2m=1000 replications. 
HR_new.boot <- boot(data1, HR_new, R = 1000, N_new_half = 500)
  ## calc 95% CIs using bootstrap percentile intervals
  HR_new.boot.ci <- boot.ci(HR_new.boot, conf = c(0.95), type = c("perc"))
print(HR_new.boot)
  results[2,1] <- 2*HR_new.boot[["call"]][["N_new_half"]]
  results[2,2] <- exp(HR_new.boot[["t0"]][["x"]])
print(HR_new.boot.ci)
  ci_vector <- HR_new.boot.ci[["percent"]]
  results[2,3] <- exp(ci_vector[4])
  results[2,4] <- exp(ci_vector[5])
  
## m=1000, so 2m=2000 replications. 
HR_new.boot <- boot(data1, HR_new, R=1000, N_new_half = 1000)
  ## calc 95% CIs using bootstrap percentile intervals
  HR_new.boot.ci <- boot.ci(HR_new.boot, conf = c(0.95), type = c("perc"))
print(HR_new.boot)
  results[3,1] <- 2*HR_new.boot[["call"]][["N_new_half"]]
  results[3,2] <- exp(HR_new.boot[["t0"]][["x"]])
print(HR_new.boot.ci)
  ci_vector <- HR_new.boot.ci[["percent"]]
  results[3,3] <- exp(ci_vector[4])
  results[3,4] <- exp(ci_vector[5])
  
## export results matrix into dataframe and then save as .dta file to read back in Stata report code
est_all <- as.data.frame(results)
  library(foreign)
  write.dta(est_all,"est_primary.dta")
  ## LAST LINE
  sink()