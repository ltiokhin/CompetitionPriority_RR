
#########
library(tidyverse)
library(rethinking)
library(lme4)
library(lmerTest)
library(effects)
library(sjPlot)
library(stringr)

rm(list=ls()) #remove previous objects
load("Pilot_Analyses_RR_without_pwr.RData") # if not already loaded in working memory

#####
#Effect of competition and effort on time until guessing 
#use model m.1.b
####

df_HPDI_time <- data.frame(c_low = rep(0, 200), c_high = rep(0, 200), ce_low = rep(0, 200), ce_high = rep(0, 200),
                           c_plusce_low = rep(0, 200), c_plusce_high = rep(0, 200))
#extract samples from posterior for all parameters
post <- extract.samples(m.1.b)

m_list <- list()

#creating 100 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(runs in 1:100) {
  
  #generate empty data frame to store simulated data
  df_guess <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), 
                         effort = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)), n_major.s = NA, ElapsedTime_Guess = NA)
  
  #extract number of observations per participant, for effort and no effort conditions
  g_num_e <- aggregate(Guess_Number ~ ID_Player, data=data.c[data.c$Effort==1,], FUN=max)
  g_num_e <- g_num_e$Guess_Number
  g_num_noe <- aggregate(Guess_Number ~ ID_Player, data=data.c[data.c$Effort==0,], FUN=max)
  g_num_noe <- g_num_noe$Guess_Number
  
  #generate multiple observations per participant
  df_guess$nobs <- NA
  #separate data frames into effort and no effort
  df_guess_e <- df_guess[df_guess$effort==0,]
  df_guess_noe <- df_guess[df_guess$effort==1,]
  
  #repeat observations in effort and no effort conditions based on number of observations in original data set
  df_guess_e$nobs <- sample(g_num_e, size = nrow(df_guess_e), replace = TRUE)
  df_guess_noe$nobs <- sample(g_num_noe, size = nrow(df_guess_noe), replace = TRUE)
  df_guess <- rbind(df_guess_e, df_guess_noe)
  n.times <- df_guess$nobs
  df_guess_full <- df_guess[rep(seq(1, nrow(df_guess)), n.times),] #full empty data frame, with varying numbers of obs per participant
  
  #generate n_major.s values for each problem that a participant solves and store in the full data
  df_guess_full$n_major.s <- sample(data.c$n_major.s, size = nrow(df_guess_full), replace = TRUE)
  
  ##generate simulated ElapsedTime_Guess values for each participant
  #generate global intercept values
  for(i in 1:nrow(df_guess_full)){
    df_guess_full$a[i] <- post$a[i]
  }
  #generate intercept values for each player
  df_guess_list <- list() 
  
  for(i in 1:200) {
    #subset data set for each unique participant
    a <- df_guess_full[df_guess_full$id==i,]
    #generate unique intercept values for that participant
    unique_int <- sample(1:ncol(post$a_player), size = 1)
    intercepts <- post$a_player[1:nrow(a),unique_int]
    a$a_player <- intercepts
    df_guess_list[[i]] <- a
  }
  
  #combine into full data set
  df_guess_full <- do.call(rbind, df_guess_list)
  
  #generate simulated mean elapsedtime for 200 players
  df_guess_full$mean_time <- NA
  
  for(i in 1:nrow(df_guess_full)) {
    
    df_guess_full$mean_time[i] <- df_guess_full$a[i] + df_guess_full$a_player[i] + post$bC[i]*df_guess_full$competition[i] + 
      post$bE[i]*df_guess_full$effort[i] + post$bCE[i]*df_guess_full$competition[i]*df_guess_full$effort[i] + 
      post$bNs[i]*df_guess_full$n_major.s[i]
  }
  
  #generate individual-level data for each participant
  df_guess_full$ElapsedTime_Guess <- rnorm(nrow(df_guess_full), mean = df_guess_full$mean_time, sd = post$sigma)
  
  #round elapsed time less than 0 to be 0
  df_guess_full$ElapsedTime_Guess[df_guess_full$ElapsedTime_Guess < 0] <- 0
  
  #subset data frame for analysis
  df <- df_guess_full[,c("id", "competition", "effort", "n_major.s", "ElapsedTime_Guess")]
  
  if(runs == 1) {
    
    m.1.power <- map2stan(
      alist(
        ElapsedTime_Guess ~ dnorm(mu, sigma), 
        mu <- a + a_player[id] + bC*competition + bE*effort + bCE*competition*effort + bNs*n_major.s, 
        bC ~ dnorm(0, 10), 
        bE ~ dnorm(0, 30), 
        bCE ~ dnorm(0, 10),
        bNs ~ dnorm(0, 10),
        a ~ dgamma(1.5, 0.05), 
        a_player[id] ~ dnorm(0, sigma_player),
        sigma_player ~ dgamma(1.5, 0.05),
        sigma ~ dgamma(2, 0.5)
      ), data=df, iter=5000, warmup=1000, chains = 4, cores = 4)
    
    samples <- extract.samples(m.1.power)
  } else if(runs > 1) {
    m_list[[runs]] <- map2stan(m.1.power, data=df, iter=5000, chains=4, cores=4, warmup=1000)
    samples <- extract.samples(m_list[[runs]])
  }
  
  df_HPDI_time$c_low[runs] <- HPDI(samples$bC, prob = 0.95)[1]
  df_HPDI_time$c_high[runs] <- HPDI(samples$bC, prob = 0.95)[2]
  df_HPDI_time$ce_low[runs] <- HPDI(samples$bCE, prob = 0.95)[1]
  df_HPDI_time$ce_high[runs] <- HPDI(samples$bCE, prob = 0.95)[2]
  df_HPDI_time$c_plusce_low[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[1]
  df_HPDI_time$c_plusce_high[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[2]
  
  #saving 
  save(df_HPDI_time, file = "df_HPDI_time.RData")
  
  ##shuffling samples from m.1.b
  post_a_player <- post["a_player"] # just samples for intercept values for each player
  post["a_player"] <- NULL #updates post to not include intercept values
  
  #shuffle samples in post
  ran_pos <- sample(1:nrow(post$a), nrow(post$a), replace = FALSE)
  post <- lapply(post, function(a) a[ran_pos])
  
  #shuffle samples in post_a_player
  for(i in 1:ncol(post_a_player$a_player)) {
    
    post_a_player$a_player[,i] <- post_a_player$a_player[ran_pos,i]
  }
  
  #add post_a_player back into post
  open_spot <- 1+length(post)
  post[[open_spot]] <- post_a_player$a_player
  names(post)[open_spot] <- "a_player"
  
}

