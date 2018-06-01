
rm(list=ls())

#########
library(tidyverse)
library(rethinking)
library(lme4)
library(lmerTest)
library(effects)
library(sjPlot)
library(stringr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#load workspace from main analyses
load("Pilot_Analyses_RR_without_pwr.RData") 


#####
#Effect of competition and effort on accuracy
#use model m.3.b
####

df_HPDI_accuracy <- data.frame(c_low = rep(0, 200), c_high = rep(0, 200), ce_low = rep(0, 200), ce_high = rep(0, 200),
                               c_plusce_low = rep(0, 200), c_plusce_high = rep(0, 200))

#extract samples from posterior for all parameters
post <- extract.samples(m.3.b)

m_list <- vector("list", length = 100)

#creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(runs in 1:200) {
  
  #generate empty data frame to store simulated data
  df_guess <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), 
                         effort = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)), 
                         a_player = NA, a = NA, n_major.s = NA, Correct_Guess = NA)
  
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
  
  #generate simulated data for all 200 participants, for the probability that they get each problem correct
  df_guess_full$probc <- NA
  
  for(i in 1:nrow(df_guess_full)) {
    df_guess_full$probc[i] <- 
      logistic(
        df_guess_full$a[i] + df_guess_full$a_player[i] + post$bC[i]*df_guess_full$competition[i] + 
          post$bE[i]*df_guess_full$effort[i] + post$bCE[i]*df_guess_full$competition[i]*df_guess_full$effort[i] + 
          post$bNs[i]*df_guess_full$n_major.s[i])
  }
  
  #generate binary 1/0 values for each trial, based on the probability of getting the problem correct
  df_guess_full$Correct_Guess <- rbinom(nrow(df_guess_full), size = 1, prob = df_guess_full$probc)
  
  #analyze data
  df <- df_guess_full[,c("id", "competition", "effort", "n_major.s", "Correct_Guess")]
  
  if(runs == 1) {
    
    m.3.power <- map2stan(
      alist(
        Correct_Guess ~ dbinom(1, theta), 
        logit(theta) <- a + a_player[id] + bC*competition + bE*effort + bCE*competition*effort + bNs*n_major.s, 
        bC ~ dnorm(0, 10), 
        bE ~ dnorm(0, 10), 
        bCE ~ dnorm(0, 10),
        bNs ~ dnorm(0, 10),
        a ~ dnorm(0, 10), 
        a_player[id] ~ dnorm(0, sigma_player),
        sigma_player ~ dgamma(1.5, 0.05)
      ), data=df, iter=6000, warmup=1000, chains=3, cores = 4)
    samples <- extract.samples(m.3.power)
    
  } else if(runs > 1) {
    m_list[[runs]] <- map2stan(m.3.power, data=df, iter=6000, warmup=1000, chains=3, cores = 4)
    samples <- extract.samples(m_list[[runs]])
  }
  
  df_HPDI_accuracy$c_low[runs] <- HPDI(samples$bC, prob = 0.95)[1]
  df_HPDI_accuracy$c_high[runs] <- HPDI(samples$bC, prob = 0.95)[2]
  df_HPDI_accuracy$ce_low[runs] <- HPDI(samples$bCE, prob = 0.95)[1]
  df_HPDI_accuracy$ce_high[runs] <- HPDI(samples$bCE, prob = 0.95)[2]
  df_HPDI_accuracy$c_plusce_low[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[1]
  df_HPDI_accuracy$c_plusce_high[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[2]
  
  #saving 
  save(df_HPDI_accuracy, file = "df_HPDI_accuracy.RData")
  
  ##shuffling samples from m.3.b
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
#########
#ROPES for Accuracy
#########
HPDI_Accuracy <- df_HPDI_accuracy

#visualizing
ggplot(HPDI_Accuracy, aes(x = 1:200)) + 
  geom_errorbar(aes(ymax = c_high, ymin = c_low))

#ROPE for C
sum(HPDI_Accuracy$c_high > -0.12)

#visualizing
ggplot(HPDI_Accuracy, aes(x = 1:200)) + 
  geom_errorbar(aes(ymax = ce_high, ymin = ce_low))

#ROPE for CE
sum(HPDI_Accuracy$ce_low > 0.05)


#####
#Effect of competition and effort on tiles
####

#original model
# m.2.b <- map2stan(
#   alist(
#     TilesRevealed ~ dnorm(mu, sigma), 
#     mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*n_major.s, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     bNs ~ dnorm(0, 10),
#     a ~ dunif(0, 25), 
#     a_player[ID_Player] ~ dnorm(0, sigma_player),
#     sigma_player ~ dgamma(1.5, 0.05),
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.c, iter=10000, chains=3, cores=4, warmup=1000)


#modified model where data  generated from a gamma distribution
m.2.b.gamma <- map2stan(
  alist(
    TilesRevealed ~ dgamma(shape, rate),
    shape <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*n_major.s,
    bC ~ dnorm(0, 10),
    bE ~ dnorm(0, 10),
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05),
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    rate ~ dgamma(1, 0.05)
  ), data=data.c, start=list(rate = 1),  iter=2000, chains=1, cores=4, warmup=1000)

precis(m.2.b.gamma)
pairs(m.2.b.gamma, pars=c("bC", "bCE", "bNs", "a"))
plot(m.2.b.gamma)

save(m.2.b.gamma, file = "m.2.b.gamma.RData")
load("m.2.b.gamma.RData")

df_HPDI_tiles <- data.frame(c_low = rep(0, 200), c_high = rep(0, 200), ce_low = rep(0, 200), ce_high = rep(0, 200),
                            c_plusce_low = rep(0, 200), c_plusce_high = rep(0, 200))

#stores each run
m_list <- vector("list", length = 200)

#extract samples from posterior for all parameters
post <- extract.samples(m.2.b.gamma)

#creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(runs in 1:200) {
  
  #generate empty data frame to store simulated data
  df_guess <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), 
                         effort = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)), n_major.s = NA, TilesRevealed = NA)
  
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
  
  ##generate simulated TilesRevealed values for each participant
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
  
  #generate simulated mean tilesrevealed for 200 players
  df_guess_full$mean_tiles <- NA
  
  for(i in 1:nrow(df_guess_full)) {
    
    df_guess_full$tiles_shape.m[i] <- df_guess_full$a[i] + df_guess_full$a_player[i] + post$bC[i]*df_guess_full$competition[i] + 
      post$bE[i]*df_guess_full$effort[i] + post$bCE[i]*df_guess_full$competition[i]*df_guess_full$effort[i] + 
      post$bNs[i]*df_guess_full$n_major.s[i]
  }
  
  #generate individual-level data for each participant
  df_guess_full$TilesRevealed <- rgamma(nrow(df_guess_full), shape = df_guess_full$tiles_shape.m, rate = post$rate)
  
  #round tiles less than 0 to be 0
  df_guess_full$TilesRevealed[df_guess_full$TilesRevealed < 0] <- 0
  
  #subset data frame for analysis
  df <- df_guess_full[,c("id", "competition", "effort", "n_major.s", "TilesRevealed")]
  
  if(runs == 1){
    
    m.2.power <- map2stan(
      alist(
        TilesRevealed ~ dnorm(mu, sigma), 
        mu <- a + a_player[id] + bC*competition + bE*effort + bCE*competition*effort + bNs*n_major.s, 
        bC ~ dnorm(0, 10), 
        bE ~ dnorm(0, 10), 
        bCE ~ dnorm(0, 10),
        bNs ~ dnorm(0, 10),
        a ~ dunif(0, 25), 
        a_player[id] ~ dnorm(0, sigma_player),
        sigma_player ~ dgamma(1.5, 0.05),
        sigma ~ dgamma(2, 0.5)
      ), data=df, iter=3000, chains=3, cores = 4, warmup=500)
    
    samples <- extract.samples(m.2.power)
  } else if(runs > 1) {
    m_list[[runs]] <- map2stan(m.2.power, data=df, iter=3000, chains=3, cores=4, warmup=500)
    samples <- extract.samples(m_list[[runs]])
  }
  
  df_HPDI_tiles$c_low[runs] <- HPDI(samples$bC, prob = 0.95)[1]
  df_HPDI_tiles$c_high[runs] <- HPDI(samples$bC, prob = 0.95)[2]
  df_HPDI_tiles$ce_low[runs] <- HPDI(samples$bCE, prob = 0.95)[1]
  df_HPDI_tiles$ce_high[runs] <- HPDI(samples$bCE, prob = 0.95)[2]
  df_HPDI_tiles$c_plusce_low[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[1]
  df_HPDI_tiles$c_plusce_high[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[2]
  
  #saving 
  save(df_HPDI_tiles, file = "df_HPDI_tiles_gamma1.RData")
  
  ##shuffling samples from m.2.b
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

#####
#Effect of competition on time to guess for a single math problem
#####

#model for below power analysis
# m.mathtime <- map2stan(
#   alist(
#     ElapsedTime_Math ~ dnorm(mu, sigma), 
#     mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s, 
#     bC ~ dnorm(0, 10), 
#     bNs ~ dnorm(0, 10),
#     a ~ dgamma(1, 0.05), 
#     a_player[ID_Player] ~ dnorm(0, sigma_player),
#     sigma_player ~ dgamma(1, 0.05),
#     sigma ~ dgamma(2, 0.5)
#   ), data=data_math.c, iter=18000, chains=3, cores=4, warmup=1000)

load("m.mathtime") 

#extract samples from posterior for all parameters
post <- extract.samples(m.mathtime)

df_HPDI_mathtime <- data.frame(c_low = rep(0, 200), c_high = rep(0, 200))

#stores each run
m_list <- vector("list", length = 200)

#creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(runs in 1:200) {
  
  #generate empty data frame to store simulated data
  df_guess <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), n_major.s = NA, ElapsedTime_Math = NA)
  
  #extract number of observations per participant
  g_num_e <- aggregate(ElapsedTime_Math ~ ID_Player, data=data_math.c, FUN=length)
  g_num_e <- g_num_e$ElapsedTime_Math
  
  #extract number of guesses per participant 
  guesses <- aggregate(Guess_Number ~ ID_Player, data=data_math.c, FUN=max)
  guesses <- guesses$Guess_Number
  
  #generate multiple observations per participant
  df_guess$nobs <- NA
  df_guess$nguess <- NA
  
  #repeat observations based on number of observations in original data set
  df_guess$nobs <- sample(g_num_e, size = nrow(df_guess), replace = TRUE)
  n.times <- df_guess$nobs
  df_guess_full <- df_guess[rep(seq(1, nrow(df_guess)), n.times),] #full empty data frame, with varying numbers of obs per participant
  
  df_guess$nguess <- sample(guesses, size = nrow(df_guess), replace = TRUE)
  n.times <- df_guess$nguess
  
  #number of math problems per guess for each participant
  for(i in unique(df_guess_full$id)){
        df_guess_full$nguess[df_guess_full$id == i] <- n.times[i]
  }
  
  df_guess_full$mathperguess <- round(df_guess_full$nobs / df_guess_full$nguess)
  
  #generate n_major.s values for each problem that a participant solves and store in the full data
  for(i in unique(df_guess_full$id)){
    maths <- unique(df_guess_full$mathperguess[df_guess_full$id==i])
    rows <- nrow(df_guess_full[df_guess_full$id == i,])
    effs <- sample(data_math.c$n_major.s, size = 1000, replace = TRUE)
    effs <- rep(effs, each = maths )
    effs <- effs[1:rows]
    df_guess_full$n_major.s[df_guess_full$id==i] <- effs
  }
  
  ##generate simulated ElapsedTime_Math values for each participant
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
  
  #generate simulated mean ElapsedTime_Math for 200 players
  df_guess_full$meanmath <- NA
  
  for(i in 1:nrow(df_guess_full)) {
    df_guess_full$meanmath[i] <- df_guess_full$a[i] + df_guess_full$a_player[i] + post$bC[i]*df_guess_full$competition[i] + 
      post$bNs[i]*df_guess_full$n_major.s[i]
  }
  
  #generate individual-level data for each participant
  df_guess_full$ElapsedTime_Math <- rnorm(nrow(df_guess_full), mean = df_guess_full$meanmath, sd = post$sigma)
  
  #round tiles less than 0 to be 0
  df_guess_full$ElapsedTime_Math[df_guess_full$ElapsedTime_Math < 0] <- 0
  
  #subset data frame for analysis
  df <- df_guess_full[,c("id", "competition", "n_major.s", "ElapsedTime_Math")]
  
  if(runs == 1){
    
    m.mathpower <- map2stan(
      alist(
        ElapsedTime_Math ~ dnorm(mu, sigma), 
        mu <- a + a_player[id] + bC*competition + bNs*n_major.s, 
        bC ~ dnorm(0, 10), 
        bNs ~ dnorm(0, 10),
        a ~ dgamma(1, 0.05), 
        a_player[id] ~ dnorm(0, sigma_player),
        sigma_player ~ dgamma(1, 0.05),
        sigma ~ dgamma(2, 0.5)
      ), data=df, iter=6000, chains=3, cores=4, warmup=1000)
    
    samples <- extract.samples(m.mathpower)
    
  } else if(runs > 1) {
    m_list[[runs]] <- map2stan(m.mathpower, data=df, iter=6000, chains=3, cores=4, warmup=1000)
    samples <- extract.samples(m_list[[runs]])
  }
  
  df_HPDI_mathtime$c_low[runs] <- HPDI(samples$bC, prob = 0.95)[1]
  df_HPDI_mathtime$c_high[runs] <- HPDI(samples$bC, prob = 0.95)[2]

  #saving 
  save(df_HPDI_mathtime, file = "df_HPDI_mathtime1.RData")
  
  ##shuffling samples from m.2.b
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

#####
#Effect of competition and effort on time until guessing 
#use model m.1.b
# ####
# 
# df_HPDI_time <- data.frame(c_low = rep(0, 200), c_high = rep(0, 200), ce_low = rep(0, 200), ce_high = rep(0, 200),
#                            c_plusce_low = rep(0, 200), c_plusce_high = rep(0, 200))
# #extract samples from posterior for all parameters
# post <- extract.samples(m.1.b)
# 
# m_list <- vector("list", length = 200)
# 
# #creating 100 simulated data sets, analyzing, and comparing their CI's to the ROPE
# for(runs in 1:200) {
#   
#   #generate empty data frame to store simulated data
#   df_guess <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), 
#                          effort = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)), n_major.s = NA, ElapsedTime_Guess = NA)
#   
#   #extract number of observations per participant, for effort and no effort conditions
#   g_num_e <- aggregate(Guess_Number ~ ID_Player, data=data.c[data.c$Effort==1,], FUN=max)
#   g_num_e <- g_num_e$Guess_Number
#   g_num_noe <- aggregate(Guess_Number ~ ID_Player, data=data.c[data.c$Effort==0,], FUN=max)
#   g_num_noe <- g_num_noe$Guess_Number
#   
#   #generate multiple observations per participant
#   df_guess$nobs <- NA
#   #separate data frames into effort and no effort
#   df_guess_e <- df_guess[df_guess$effort==0,]
#   df_guess_noe <- df_guess[df_guess$effort==1,]
#   
#   #repeat observations in effort and no effort conditions based on number of observations in original data set
#   df_guess_e$nobs <- sample(g_num_e, size = nrow(df_guess_e), replace = TRUE)
#   df_guess_noe$nobs <- sample(g_num_noe, size = nrow(df_guess_noe), replace = TRUE)
#   df_guess <- rbind(df_guess_e, df_guess_noe)
#   n.times <- df_guess$nobs
#   df_guess_full <- df_guess[rep(seq(1, nrow(df_guess)), n.times),] #full empty data frame, with varying numbers of obs per participant
#   
#   #generate n_major.s values for each problem that a participant solves and store in the full data
#   df_guess_full$n_major.s <- sample(data.c$n_major.s, size = nrow(df_guess_full), replace = TRUE)
#   
#   ##generate simulated ElapsedTime_Guess values for each participant
#   #generate global intercept values
#   for(i in 1:nrow(df_guess_full)){
#     df_guess_full$a[i] <- post$a[i]
#   }
#   #generate intercept values for each player
#   df_guess_list <- list() 
#   
#   for(i in 1:200) {
#     #subset data set for each unique participant
#     a <- df_guess_full[df_guess_full$id==i,]
#     #generate unique intercept values for that participant
#     unique_int <- sample(1:ncol(post$a_player), size = 1)
#     intercepts <- post$a_player[1:nrow(a),unique_int]
#     a$a_player <- intercepts
#     df_guess_list[[i]] <- a
#   }
#   
#   #combine into full data set
#   df_guess_full <- do.call(rbind, df_guess_list)
#   
#   #generate simulated mean elapsedtime for 200 players
#   df_guess_full$mean_time <- NA
#   
#   for(i in 1:nrow(df_guess_full)) {
#     
#     df_guess_full$mean_time[i] <- df_guess_full$a[i] + df_guess_full$a_player[i] + post$bC[i]*df_guess_full$competition[i] + 
#       post$bE[i]*df_guess_full$effort[i] + post$bCE[i]*df_guess_full$competition[i]*df_guess_full$effort[i] + 
#       post$bNs[i]*df_guess_full$n_major.s[i]
#   }
#   
#   #generate individual-level data for each participant
#   df_guess_full$ElapsedTime_Guess <- rnorm(nrow(df_guess_full), mean = df_guess_full$mean_time, sd = post$sigma)
#   
#   #round elapsed time less than 0 to be 0
#   df_guess_full$ElapsedTime_Guess[df_guess_full$ElapsedTime_Guess < 0] <- 0
#   
#   #subset data frame for analysis
#   df <- df_guess_full[,c("id", "competition", "effort", "n_major.s", "ElapsedTime_Guess")]
#   
#   if(runs == 1) {
#     
#     m.1.power <- map2stan(
#       alist(
#         ElapsedTime_Guess ~ dnorm(mu, sigma), 
#         mu <- a + a_player[id] + bC*competition + bE*effort + bCE*competition*effort + bNs*n_major.s, 
#         bC ~ dnorm(0, 10), 
#         bE ~ dnorm(0, 30), 
#         bCE ~ dnorm(0, 10),
#         bNs ~ dnorm(0, 10),
#         a ~ dgamma(1.5, 0.05), 
#         a_player[id] ~ dnorm(0, sigma_player),
#         sigma_player ~ dgamma(1.5, 0.05),
#         sigma ~ dgamma(2, 0.5)
#       ), data=df, iter=5000, warmup=1000, chains = 4, cores = 4)
#     
#     samples <- extract.samples(m.1.power)
#   } else if(runs > 1) {
#     m_list[[runs]] <- map2stan(m.1.power, data=df, iter=5000, chains=4, cores=4, warmup=1000)
#     samples <- extract.samples(m_list[[runs]])
#   }
#   
#   df_HPDI_time$c_low[runs] <- HPDI(samples$bC, prob = 0.95)[1]
#   df_HPDI_time$c_high[runs] <- HPDI(samples$bC, prob = 0.95)[2]
#   df_HPDI_time$ce_low[runs] <- HPDI(samples$bCE, prob = 0.95)[1]
#   df_HPDI_time$ce_high[runs] <- HPDI(samples$bCE, prob = 0.95)[2]
#   df_HPDI_time$c_plusce_low[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[1]
#   df_HPDI_time$c_plusce_high[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[2]
#   
#   #saving 
#   save(df_HPDI_time, file = "df_HPDI_time.RData")
#   
#   ##shuffling samples from m.1.b
#   post_a_player <- post["a_player"] # just samples for intercept values for each player
#   post["a_player"] <- NULL #updates post to not include intercept values
#   
#   #shuffle samples in post
#   ran_pos <- sample(1:nrow(post$a), nrow(post$a), replace = FALSE)
#   post <- lapply(post, function(a) a[ran_pos])
#   
#   #shuffle samples in post_a_player
#   for(i in 1:ncol(post_a_player$a_player)) {
#     
#     post_a_player$a_player[,i] <- post_a_player$a_player[ran_pos,i]
#   }
#   
#   #add post_a_player back into post
#   open_spot <- 1+length(post)
#   post[[open_spot]] <- post_a_player$a_player
#   names(post)[open_spot] <- "a_player"
#   
# }

#####
#Effect of competition and effort on guess rate
#use model m.4.b
####
# 
# df_HPDI_guessrate <- data.frame(c_low = rep(0, 200), c_high = rep(0, 200), ce_low = rep(0, 200), ce_high = rep(0, 200),
#                                 c_plusce_low = rep(0, 200), c_plusce_high = rep(0, 200))
# 
# m_list <- vector("list", length = 100)
# 
# #extract samples from posterior for all parameters
# post <- extract.samples(m.4.b)
# 
# #creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
# for(runs in 1:200) {
#   
#   #generate empty data frame to store simulated data
#   df_guess_full <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), 
#                               effort = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)), mean_n_major.s = NA, guessrate = NA)
#   
#   #generate n_major.s values for each problem that a participant solves and store in the full data
#   df_guess_full$mean_n_major.s <- sample(data.rate$mean_n_major.s, size = nrow(df_guess_full), replace = TRUE)
#   
#   #generate simulated mean guessrate for 200 players
#   df_guess_full$mean_guessrate <- NA
#   
#   for(i in 1:nrow(df_guess_full)) {
#     
#     df_guess_full$mean_guessrate[i] <- post$a[i] + post$bC[i]*df_guess_full$competition[i] + 
#       post$bE[i]*df_guess_full$effort[i] + post$bCE[i]*df_guess_full$competition[i]*df_guess_full$effort[i] + 
#       post$bNs[i]*df_guess_full$mean_n_major.s[i]
#   }
#   
#   #generate individual-level data for each participant
#   df_guess_full$guessrate <- rnorm(nrow(df_guess_full), mean = df_guess_full$mean_guessrate, sd = post$sigma)
#   
#   #round guessrates less than 0 to be 0
#   df_guess_full$guessrate[df_guess_full$guessrate < 0] <- 0
#   
#   #subset data frame for analysis
#   df <- df_guess_full[,c("id", "competition", "effort", "mean_n_major.s", "guessrate")]
#   
#   if(runs == 1) {
#     
#     m.4.power <- map2stan(
#       alist(
#         guessrate ~ dnorm(mu, sigma), 
#         mu <- a + bC*competition + bE*effort + bCE*competition*effort + bNs*mean_n_major.s, 
#         bC ~ dnorm(0, 10), 
#         bE ~ dnorm(0, 10), 
#         bCE ~ dnorm(0, 10),
#         bNs ~ dnorm(0, 10),
#         a ~ dgamma(2, 0.05), 
#         sigma ~ dgamma(2, 0.5)
#       ), data=df, iter=5000, chains=4, cores = 4, warmup=1000)
#     
#     samples <- extract.samples(m.4.power) 
#   } else if(runs > 1) {
#     m_list[[runs]] <- map2stan(m.4.power, data=df, iter=5000, chains=4, cores=4, warmup=1000)
#     samples <- extract.samples(m_list[[runs]])
#   }
#   
#   df_HPDI_guessrate$c_low[runs] <- HPDI(samples$bC, prob = 0.95)[1]
#   df_HPDI_guessrate$c_high[runs] <- HPDI(samples$bC, prob = 0.95)[2]
#   df_HPDI_guessrate$ce_low[runs] <- HPDI(samples$bCE, prob = 0.95)[1]
#   df_HPDI_guessrate$ce_high[runs] <- HPDI(samples$bCE, prob = 0.95)[2]
#   df_HPDI_guessrate$c_plusce_low[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[1]
#   df_HPDI_guessrate$c_plusce_high[runs] <- HPDI(samples$bC + samples$bCE, prob = 0.95)[2]
#   
#   #save 
#   save(df_HPDI_guessrate, file = "df_HPDI_guessrate.RData")
#   
#   #shuffle samples in post
#   ran_pos <- sample(1:nrow(post$a), nrow(post$a), replace = FALSE)
#   post <- lapply(post, function(a) a[ran_pos])
# }
# 
# 
# 
# 



