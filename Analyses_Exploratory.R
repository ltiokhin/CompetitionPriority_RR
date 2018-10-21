
rm(list=ls())

#########
library(tidyverse)
library(knitr)
library(kableExtra)
library(dplyr)
library(rethinking)
library(sjPlot)
library(stringr)
library(RColorBrewer)

###Load Data###

#load("XXXXX")
#load("XXXXX")
sequences <- read.table("Sequences1.txt", header=FALSE) #the sequence of blue and yellow tiles for each grid

###Cleaning and Data Prep###

###Store original data in separate objects###
d.conf_1_2.b <- d.conf_1_2
d.conf_3_math.b <- d.conf_3_math

#Exclude participants who indicated technical difficulties
d.conf_1_2.b <- d.conf_1_2.b[d.conf_1_2.b$ID_Player != 26,]
d.conf_3_math.b <- d.conf_3_math.b[d.conf_3_math.b$ID_Player != 26,]
d.conf_1_2.b <- d.conf_1_2.b[d.conf_1_2.b$ID_Player != 61,]
d.conf_3_math.b <- d.conf_3_math.b[d.conf_3_math.b$ID_Player != 61,]

#Exclude rows with NA's, but keep participants who did not complete study
d.conf_1_2.b <- d.conf_1_2.b[complete.cases(d.conf_1_2.b[,c(1:4, 6:12)]),]
d.conf_3_math.b <- d.conf_3_math.b[complete.cases(d.conf_3_math.b[,c(1:4, 6:17)]),]

###for each of the 600 possible sequences, the number of majority tiles (i.e. effect size)###
seq.df <- data.frame(seq = sequences)
t <- str_split(seq.df$V1, ";")
t <- do.call(rbind, t)
t <- as.data.frame(t)

for(i in 1:ncol(t)){
  t[,i] <- as.numeric(t[,i]) - 1
}

#number of majority tiles for each grid
n_majority <- apply(t, 1, sum) 

#identifies which grids correspond to which effect size
g_13 <- which(n_majority == 13)
g_15 <- which(n_majority == 15)
g_17 <- which(n_majority == 17)

#adds this to relevant data sets
d.conf_1_2.b$n_major <- rep(NA, nrow(d.conf_1_2.b))
d.conf_3_math.b$n_major <- rep(NA, nrow(d.conf_3_math.b))

for(i in 1:nrow(d.conf_1_2.b)){
  if(d.conf_1_2.b$Guess_Number[i] %in% g_13){d.conf_1_2.b$n_major[i] <- 13}
  if(d.conf_1_2.b$Guess_Number[i] %in% g_15){d.conf_1_2.b$n_major[i] <- 15}
  if(d.conf_1_2.b$Guess_Number[i] %in% g_17){d.conf_1_2.b$n_major[i] <- 17}
}

for(i in 1:nrow(d.conf_3_math.b)){
  if(d.conf_3_math.b$Guess_Number[i] %in% g_13){d.conf_3_math.b$n_major[i] <- 13}
  if(d.conf_3_math.b$Guess_Number[i] %in% g_15){d.conf_3_math.b$n_major[i] <- 15}
  if(d.conf_3_math.b$Guess_Number[i] %in% g_17){d.conf_3_math.b$n_major[i] <- 17}
}

#Standardize number of majority tiles
d.conf_1_2.b$n_major.s <- (d.conf_1_2.b$n_major - mean(d.conf_1_2.b$n_major)) / sd(d.conf_1_2.b$n_major)
d.conf_3_math.b$n_major.s <- (d.conf_3_math.b$n_major - mean(d.conf_3_math.b$n_major)) / sd(d.conf_3_math.b$n_major)

#Aggregate data for each participant for confirmatory analyses 1 and 2
d.conf.agg <- aggregate(cbind(TilesRevealed, Correct_Guess) ~ Effort + Competition +
                          Sex + Guess_Number + n_major.s + ID_Player, data=d.conf_1_2.b, FUN = mean)

d.conf.agg$ID_Player <- coerce_index(d.conf.agg$ID_Player)
d.conf.agg$Sex <- as.integer(d.conf.agg$Sex)

###########################
###Exploratory Analyses###
##########################

###Repeating confirmatory analyses without excluding 1)  outliers or 2) participants who did not complete the study

######
###Effect of competition and effort on number of tiles revealed###
######

m.tiles.noexclusion <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + gamma*Competition + bE*Effort + bNs*n_major.s, 
    gamma <- bC + bCE*Effort,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=20000, chains=4, cores=4, warmup=500)

plot(m.tiles.noexclusion)
par(mfrow=c(1,1))
precis(m.tiles.noexclusion, prob = 0.95)
plot(precis(m.tiles.noexclusion, prob = 0.95), xlab = "Tiles Revealed")

######
###Effect of competition and effort on accuracy###
######

m.accuracy.noexclusion <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=6000, chains = 2, cores = 2, warmup=500)

plot(m.accuracy)
par(mfrow=c(1,1))
precis(m.accuracy, prob=0.95)
plot(precis(m.accuracy, prob=0.95), xlab = "Log Odds of Correct Guess")

##########
###Effect of competition on effort (i.e. time to accurately solve an arithmetic problem)###
#########

d.conf_3_math.b$ID_Player <- coerce_index(d.conf_3_math.b$ID_Player)
d.conf_3_math.b$Sex <- as.integer(d.conf_3_math.b$Sex)

d.math.agg <- aggregate(cbind(ElapsedTime_MathSolved) ~ Effort + Competition + Sex + Tile_Number +
                          Guess_Number + n_major.s + ID_Player, data=d.conf_3_math.b, FUN = mean)

### Model 3 ###
m.effort.noexclusion <- map2stan(
  alist(
    ElapsedTime_MathSolved ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.math.agg, iter=20000, chains = 3, cores = 3, warmup=500)

plot(m.effort.noexclusions)
par(mfrow=c(1,1))
precis(m.effort.noexclusions, prob=0.95)
plot(precis(m.effort.noexclusions, prob=0.95), xlab = "Time (seconds) to Accurately Solve One Arithmetic Problem")


##############
###Table with parameter estimates from models without exclusions
##############

precis_tiles <- precis(m.tiles.noexclusion, depth=1, prob = 0.95)
precis_tiles <- precis_tiles@output
precis_tiles <- round(precis_tiles, 2)
precis_tiles <- precis_tiles[1:5, c(1, 3, 4)]
precis_tiles <- rename(precis_tiles, `Lower 0.95` = `lower 0.95`)
precis_tiles <- rename(precis_tiles, `Upper 0.95` = `upper 0.95`)

row.names(precis_tiles) <- c("1. Competition", "1. Effort", "1. Competition x Effort Interaction", 
                             "1. Effect Size", "1. Intercept")

precis_accuracy <- precis(m.accuracy.noexclusion, depth=1, prob = 0.95)
precis_accuracy <- precis_accuracy@output
precis_accuracy <- round(precis_accuracy, 2)
precis_accuracy <- precis_accuracy[1:5, c(1, 3, 4)]
precis_accuracy <- rename(precis_accuracy, `Lower 0.95` = `lower 0.95`)
precis_accuracy <- rename(precis_accuracy, `Upper 0.95` = `upper 0.95`)

row.names(precis_accuracy) <- c("2. Competition", "2. Effort", "2. Competition x Effort Interaction", 
                                "2. Effect Size", "2. Intercept")

precis_e <- precis(m.effort.noexclusion, depth=1, prob = 0.95)
precis_e <- precis_e@output
precis_e <- round(precis_e, 2)
precis_e <- precis_e[1:3, c(1, 3, 4)]
precis_e <- rename(precis_e, `Lower 0.95` = `lower 0.95`)
precis_e <- rename(precis_e, `Upper 0.95` = `upper 0.95`)

row.names(precis_e) <- c("3. Competition",  
                         "3. Effect Size", "3. Intercept")

mcomb_df <- rbind(precis_tiles, precis_accuracy, precis_e)

mcomb_df %>% kable(caption = "CONFIRMATORY MODELS WITHOUT EXCLUSIONS")  %>% 
  kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), 
                full_width = TRUE) %>%
  column_spec(1, width = "5cm") %>%
  column_spec(2:4, width = "2cm") %>%
  group_rows("Model 1: Tiles Revealed", 1, 5) %>%
  group_rows("Model 2: Accuracy", 6, 10) %>%
  group_rows("Model 3: Time (seconds) to Accurately Solve One Arithmetic Problem", 11, 13)

#############
###Repeating confirmatory analyses with extra exclusions: in addition to excluding participants who did###
#not complete the study, and excluding time-until-guess and arithmetic-solving-times more than 5 standard deviations
#away from the mean, this analysis also excludes participants who removed 0 tiles when guessing.
############

#Exclude participants who did not complete study
d.conf_1_2.complete <- d.conf_1_2.b[complete.cases(d.conf_1_2.b),]
d.conf_3_math.complete <- d.conf_3_math.b[complete.cases(d.conf_3_math.b),]

#remove rows with time-until-guess values that are more than 5 standard deviations away from the mean
agg <- aggregate(ElapsedTime_Guess ~ Guess_Number + ID_Player + Effort + Competition, data=d.conf_1_2.complete, FUN = mean)
sd5.times <- mean(agg$ElapsedTime_Guess) + (5 * sd(agg$ElapsedTime_Guess))
nrow(agg[agg$ElapsedTime_Guess > sd5.times,]) 

#excluding observations from all datasets
d.conf_1_2.complete <- d.conf_1_2.complete[d.conf_1_2.complete$ElapsedTime_Guess < sd5.times,]
d.conf_3_math.complete <- d.conf_3_math.complete[d.conf_3_math.complete$ElapsedTime_Guess < sd5.times,]

#remove rows with arithmetic-solving-times that are more than 5 standard deviations away from the mean
sd5.times <- mean(d.conf_3_math.complete$ElapsedTime_Math) + (5 * sd(d.conf_3_math.complete$ElapsedTime_Math))
nrow(d.conf_3_math.complete[d.conf_3_math.complete$ElapsedTime_Math > sd5.times,]) 

#excluding observations from math dataset
d.conf_3_math.complete <- d.conf_3_math.complete[d.conf_3_math.complete$ElapsedTime_Math < sd5.times,]

##############
#exclude participants who revealed 0 tiles
##############
d.conf_1_2.complete_no0 <- d.conf_1_2.complete[d.conf_1_2.complete$TilesRevealed > 0,]
d.conf_3_math.complete_no0 <- d.conf_3_math.complete[d.conf_3_math.complete$TilesRevealed > 0,]

#Re-standardize number of majority tiles
d.conf_1_2.complete_no0$n_major.s <- (d.conf_1_2.complete_no0$n_major - mean(d.conf_1_2.complete_no0$n_major)) 
                                      / sd(d.conf_1_2.complete_no0$n_major)
d.conf_3_math.complete_no0$n_major.s <- (d.conf_3_math.complete_no0$n_major - mean(d.conf_3_math.complete_no0$n_major)) 
                                      / sd(d.conf_3_math.complete_no0$n_major)

#Aggregate data for each participant for confirmatory analyses 1 and 2
d.conf.agg <- aggregate(cbind(TilesRevealed, Correct_Guess) ~ Effort + Competition +
                          Sex + Guess_Number + n_major.s + ID_Player, data=d.conf_1_2.complete_no0, FUN = mean)

d.conf.agg$ID_Player <- coerce_index(d.conf.agg$ID_Player)
d.conf.agg$Sex <- as.integer(d.conf.agg$Sex)

######
###Effect of competition and effort on number of tiles revealed###
######

m.tiles.superexclusion <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + gamma*Competition + bE*Effort + bNs*n_major.s, 
    gamma <- bC + bCE*Effort,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=15000, chains=4, cores=4, warmup=500)

plot(m.tiles.superexclusion)
par(mfrow=c(1,1))
precis(m.tiles.superexclusion, prob = 0.95)
plot(precis(m.tiles.superexclusion, prob = 0.95), xlab = "Tiles Revealed")

######
###Effect of competition and effort on accuracy###
######

m.accuracy.superexclusion <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=4500, chains = 2, cores = 2, warmup=500)

plot(m.accuracy.superexclusion)
par(mfrow=c(1,1))
precis(m.accuracy.superexclusion, prob=0.95)
plot(precis(m.accuracy.superexclusion, prob=0.95), xlab = "Log Odds of Correct Guess")

#########
###Effect of competition on effort (i.e. time to accurately solve an arithmetic problem)###
#########

d.conf_3_math.complete_no0$ID_Player <- coerce_index(d.conf_3_math.complete_no0$ID_Player)
d.conf_3_math.complete_no0$Sex <- as.integer(d.conf_3_math.complete_no0$Sex)

d.math.agg <- aggregate(cbind(ElapsedTime_MathSolved) ~ Effort + Competition + Sex + Tile_Number +
                          Guess_Number + n_major.s + ID_Player, data=d.conf_3_math.complete_no0, FUN = mean)

### Model 3 ###
m.effort.superexclusion <- map2stan(
  alist(
    ElapsedTime_MathSolved ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.math.agg, iter=20000, chains = 3, cores = 3, warmup=500)

plot(m.effort.superexclusion)
par(mfrow=c(1,1))
precis(m.effort.superexclusion, prob=0.95)
plot(precis(m.effort.superexclusion, prob=0.95), xlab = "Time (seconds) to Accurately Solve One Arithmetic Problem")



##############
###Table with parameter estimates from models with super exclusions
##############

precis_tiles <- precis(m.tiles.superexclusion, depth=1, prob = 0.95)
precis_tiles <- precis_tiles@output
precis_tiles <- round(precis_tiles, 2)
precis_tiles <- precis_tiles[1:5, c(1, 3, 4)]
precis_tiles <- rename(precis_tiles, `Lower 0.95` = `lower 0.95`)
precis_tiles <- rename(precis_tiles, `Upper 0.95` = `upper 0.95`)

row.names(precis_tiles) <- c("1. Competition", "1. Effort", "1. Competition x Effort Interaction", 
                             "1. Effect Size", "1. Intercept")

precis_accuracy <- precis(m.accuracy.superexclusion, depth=1, prob = 0.95)
precis_accuracy <- precis_accuracy@output
precis_accuracy <- round(precis_accuracy, 2)
precis_accuracy <- precis_accuracy[1:5, c(1, 3, 4)]
precis_accuracy <- rename(precis_accuracy, `Lower 0.95` = `lower 0.95`)
precis_accuracy <- rename(precis_accuracy, `Upper 0.95` = `upper 0.95`)

row.names(precis_accuracy) <- c("2. Competition", "2. Effort", "2. Competition x Effort Interaction", 
                                "2. Effect Size", "2. Intercept")

precis_e <- precis(m.effort.superexclusion, depth=1, prob = 0.95)
precis_e <- precis_e@output
precis_e <- round(precis_e, 2)
precis_e <- precis_e[1:3, c(1, 3, 4)]
precis_e <- rename(precis_e, `Lower 0.95` = `lower 0.95`)
precis_e <- rename(precis_e, `Upper 0.95` = `upper 0.95`)

row.names(precis_e) <- c("3. Competition",  
                         "3. Effect Size", "3. Intercept")

mcomb_df <- rbind(precis_tiles, precis_accuracy, precis_e)

mcomb_df %>% kable(caption = "CONFIRMATORY MODELS EXCLUDING CASES WHEN PARTICIPANTS GUESSED RANDOMLY (i.e. EXCLUDING TILES REVEALED = 0)")  %>% 
  kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), 
                full_width = TRUE) %>%
  column_spec(1, width = "5cm") %>%
  column_spec(2:4, width = "2cm") %>%
  group_rows("Model 1: Tiles Revealed", 1, 5) %>%
  group_rows("Model 2: Accuracy", 6, 10) %>%
  group_rows("Model 3: Time (seconds) to Accurately Solve One Arithmetic Problem", 11, 13)

##########
###Sensitivity checks using the exclusion criteria in the main paper
##########

#Re-standardize number of majority tiles
d.conf_1_2.complete$n_major.s <- (d.conf_1_2.complete$n_major - mean(d.conf_1_2.complete$n_major)) /
                                          sd(d.conf_1_2.complete$n_major)
d.conf_3_math.complete$n_major.s <- (d.conf_3_math.complete$n_major - mean(d.conf_3_math.complete$n_major)) / 
                                          sd(d.conf_3_math.complete$n_major)

#Aggregate data for each participant for confirmatory analyses 1 and 2
d.conf.agg <- aggregate(cbind(TilesRevealed, Correct_Guess) ~ Effort + Competition +
                          Sex + Guess_Number + n_major.s + ID_Player, data=d.conf_1_2.complete, FUN = mean)

d.conf.agg$ID_Player <- coerce_index(d.conf.agg$ID_Player)
d.conf.agg$Sex <- as.integer(d.conf.agg$Sex)
d.conf.agg$Sex[d.conf.agg$Sex==1] <- 0 #female
d.conf.agg$Sex[d.conf.agg$Sex==2] <- 1 #male

###Tiles: No effort 
m.tiles.noE <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

precis(m.tiles.noE, prob = 0.95)

###Tiles: No effort interaction

m.tiles.noE_interaction <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

###Tiles: original confirmatory analyses 

m.tiles.orig <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

#compare alternative models
compare(m.tiles.noE, m.tiles.noE_interaction, m.tiles.orig)
#compare just interaction vs no-interaction models
compare(m.tiles.noE_interaction, m.tiles.orig)

multi_models <- coeftab(m.tiles.noE, m.tiles.noE_interaction, m.tiles.orig)
coeftab_plot(multi_models, pars=c("bC", "bE", "bNs","bCE"), 
             prob = 0.95, col.ci = rangi2)

###Tiles: No effort interaction + guess number

d.conf.agg$Guess_Number.s <- (d.conf.agg$Guess_Number - mean(d.conf.agg$Guess_Number)) / sd(d.conf.agg$Guess_Number)

m.tiles.guessnum <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s + bGs*Guess_Number.s, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

###Tiles: main effect of sex####
m.tiles.sex.nointer <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s + bGs*Guess_Number.s + 
      bS*Sex,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10), 
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

###Tiles: main effect of sex and interaction####
m.tiles.sex.inter <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + 
      bNs*n_major.s + bGs*Guess_Number.s + 
      bS*Sex + bCS*Competition*Sex + bES*Effort*Sex,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10), 
    bCS ~ dnorm(0, 10),
    bES ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

m.tiles.sex.inter2 <- map2stan(
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort + 
      bNs*n_major.s + bGs*Guess_Number.s + 
      bS*Sex + bCS*Competition*Sex + bES*Effort*Sex,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10),
    bCE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10), 
    bCS ~ dnorm(0, 10),
    bES ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)


#Compare alternative models using information criteria and plot parameter estimates
compare(m.tiles.guessnum, m.tiles.sex.nointer, m.tiles.sex.inter)

###Interactions with Guess Number and n_major.s
m.tiles.guess.nmajor.inter <- map2stan( 
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s + bGs*Guess_Number.s + 
      bEGs*Effort*Guess_Number.s + bCGs*Competition*Guess_Number.s +
      bENs*Effort*n_major.s + bCNs*Competition*n_major.s,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bEGs ~ dnorm(0, 10),
    bCGs ~ dnorm(0, 10),
    bENs ~ dnorm(0, 10),
    bCNs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

plot(precis(m.tiles.guess.nmajor.inter))

m.tiles.guess.nmajor.inter2only <- map2stan( 
  alist(
    TilesRevealed ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s + bGs*Guess_Number.s + 
      bENs*Effort*n_major.s + bCNs*Competition*n_major.s,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bENs ~ dnorm(0, 10),
    bCNs ~ dnorm(0, 10),
    a ~ dunif(0, 25), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.conf.agg, iter=5000, chains=3, cores=3, warmup=500)

precis(m.tiles.guess.nmajor.inter2only)

#Final comparison of alternative models for tiles revealed (see if sex should go in here)
compare(m.tiles.orig, m.tiles.guessnum, m.tiles.sex.inter, m.tiles.sex.inter2,
        m.tiles.guess.nmajor.inter, m.tiles.guess.nmajor.inter2only, 
        m.tiles.noE)

####RUUNNING SEX INTERACTION MODEL - REDO THE MULTIMODELS AND COEFPLOT FOR TILES AND ACCURACY
multi_models <- coeftab( 
                        m.tiles.guess.nmajor.inter, 
                        m.tiles.guess.nmajor.inter2only, 
                        m.tiles.guessnum,
                        m.tiles.sex.inter,
                        m.tiles.sex.inter2,
                        m.tiles.orig)

coeftab_plot(multi_models, pars=c("bC", "bE", "bNs", "bGs", "bCE", "bES", "bCS",
                                  "bENs", "bCNs", "bCGs", "bEGs"), 
                            prob = 0.95, col.ci = rangi2, main = "Tiles Revealed")

#concise
coeftab_plot(multi_models, pars=c("bC", "bE", "bCE", "bNs"), 
             prob = 0.95, col.ci = rangi2, main = "Tiles Revealed")

###Visualize raw data: guess number against tiles revealed###
plot(d.conf.agg$TilesRevealed[d.conf.agg$Competition==0 &
                                 d.conf.agg$Effort==0] ~ 
       d.conf.agg$Guess_Number[d.conf.agg$Competition==0 &
                                 d.conf.agg$Effort==0], 
     col=col.alpha(rangi2, alpha = 0.4))

###Accuracy###
m_accuracy_orig.guess <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bE*Effort + 
      bCE*Competition*Effort + bNs*n_major.s + bGs*Guess_Number.s,  
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10),
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=4500, chains = 2, cores = 2, warmup=500)

m_accuracy_componly <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s + bGs*Guess_Number.s,
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=4500, chains = 2, cores = 2, warmup=500)

m_accuracy_sex <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s + bGs*Guess_Number.s + bS*Sex, 
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=4500, chains = 2, cores = 2, warmup=500)

m_accuracy_sex_inter <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s + 
      bGs*Guess_Number.s + bS*Sex + bNsS*n_major.s*Sex + bCS*Competition*Sex + bES*Effort*Sex,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10),
    bNsS ~ dnorm(0, 10), 
    bCS ~ dnorm(0, 10), 
    bES ~ dnorm(0, 10), 
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=4500, chains = 2, cores = 2, warmup=500)

m_accuracy_nmajor_sex_inters <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s + 
      bGs*Guess_Number.s + bS*Sex + bNsS*n_major.s*Sex + 
      bENs*Effort*n_major.s + bCNs*Competition*n_major.s,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10), 
    bNsS ~ dnorm(0, 10), 
    bENs ~ dnorm(0, 10),
    bCNs ~ dnorm(0, 10),
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=4500, chains = 2, cores = 2, warmup=500)

m_accuracy_compeff_guess_nmajor_inter <- map2stan(
  alist(
    Correct_Guess ~ dbinom(1, theta), 
    logit(theta) <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bNs*n_major.s + 
      bGs*Guess_Number.s + bCE*Competition*Effort + 
      bENs*Effort*n_major.s + bCNs*Competition*n_major.s,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bCE ~ dnorm(0, 10), 
    bENs ~ dnorm(0, 10),
    bCNs ~ dnorm(0, 10),
    a ~ dnorm(0, 10), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05)
  ), data=d.conf.agg, iter=4500, chains = 2, cores = 2, warmup=500)

plot(precis(m_accuracy_compeff_guess_nmajor_inter))

#compare alternative models using information criteria
compare(
  m_accuracy_nmajor_sex_inters,
  m_accuracy_compeff_guess_nmajor_inter,
  m_accuracy_sex_inter, 
  m_accuracy_componly, 
  m_accuracy_sex, 
  m_accuracy_orig.guess)

multi_models_accuracy <- coeftab(
                        m_accuracy_nmajor_sex_inters,
                        m_accuracy_compeff_guess_nmajor_inter,
                        m_accuracy_sex_inter, 
                        m_accuracy_sex, 
                        m_accuracy_orig.guess)

#full plot
coeftab_plot(multi_models_accuracy, 
             pars=c("bC", "bE", "bNs", "bGs", "bS", "bCE", "bENs", "bCNs",  "bCS", "bES", "bNsS"), 
             prob = 0.95, col.ci = rangi2, 
             main = "Accuracy: Log Odds of Correct Guess")
#concise
coeftab_plot(multi_models_accuracy, 
             pars=c("bC", "bE", "bNs", "bCE"),
             prob = 0.95, col.ci = rangi2, main = "Accuracy: Log Odds of Correct Guess")

#Arithmetic problems
d.conf_3_math.complete$ID_Player <- coerce_index(d.conf_3_math.complete$ID_Player)
d.conf_3_math.complete$Sex <- as.integer(d.conf_3_math.complete$Sex)

d.math.agg <- aggregate(cbind(ElapsedTime_MathSolved, ElapsedTime_Math) ~ Effort + Competition + 
                          Sex + Tile_Number + Guess_Number + n_major.s + ID_Player, 
                        data=d.conf_3_math.complete, FUN = mean)

d.math.agg$Guess_Number.s <- (d.math.agg$Guess_Number - mean(d.math.agg$Guess_Number)) / sd(d.math.agg$Guess_Number)
d.math.agg$Sex[d.math.agg$Sex==1] <- 0
d.math.agg$Sex[d.math.agg$Sex==2] <- 1

m_effort_orig <- map2stan(
  alist(
    ElapsedTime_MathSolved ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.math.agg, iter=5000, chains = 3, cores = 3, warmup=500)

m_effort_guessnum <- map2stan(
  alist(
    ElapsedTime_MathSolved ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s + bGs*Guess_Number.s,
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.math.agg, iter=5000, chains = 3, cores = 3, warmup=500)

compare(m_effort_orig, m_effort_guessnum)

m_effort_sex_inter <- map2stan(
  alist(
    ElapsedTime_MathSolved ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s + bGs*Guess_Number.s +
      bS*Sex + bCS*Competition*Sex, 
    bC ~ dnorm(0, 10), 
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10), 
    bCS ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.math.agg, iter=5000, chains = 3, cores = 3, warmup=500)

m_effort_inters <- map2stan(
  alist(
    ElapsedTime_MathSolved ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s + bGs*Guess_Number.s +
      bS*Sex + bCS*Competition*Sex + bGsS*Guess_Number.s*Sex + bCGs*Competition*Guess_Number.s,
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10), 
    bCS ~ dnorm(0, 10), 
    bGsS ~ dnorm(0, 10), 
    bCGs ~ dnorm(0, 10), 
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.math.agg, iter=5000, chains = 3, cores = 3, warmup=500)

plot(precis(m_effort_inters))

#compare alternative models
compare(m_effort_inters, m_effort_sex_inter, m_effort_guessnum, m_effort_orig)

multi_models_effort <- coeftab(m_effort_inters, m_effort_sex_inter, m_effort_guessnum, m_effort_orig)

#full plot
coeftab_plot(multi_models_effort, 
             pars=c("bC", "bNs", "bGs", "bS",  "bCS", "bGsS", "bCGs"), 
             prob = 0.95, col.ci = rangi2)
#part of plot
coeftab_plot(multi_models_effort, 
             pars=c("bC", "bNs", "bGs", "bS",  "bCS"), 
             prob = 0.95, col.ci = rangi2, main = "Time (seconds) to Accurately Solve One Arithmetic Problem")

#time to produce any answer (i.e correct or incorrect)
m_math_any_sex_inter <- map2stan(
  alist(
    ElapsedTime_Math ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s + bGs*Guess_Number.s +
      bS*Sex + bCS*Competition*Sex, 
    bC ~ dnorm(0, 10), 
    bGs ~ dnorm(0, 10),
    bS ~ dnorm(0, 10), 
    bCS ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.math.agg, iter=5000, chains = 3, cores = 3, warmup=500)

plot(precis(m_math_any_sex_inter, prob=0.95), 
     main = "Time (seconds) to Produce Any Answer to One Arithmetic Problem", 
     xlab = "Estimate", 
     col.ci = rangi2)

###########
###Reward per unit time, without outlier removal###
###########

df.reward <- d.conf_1_2.b[complete.cases(d.conf_1_2.b),]

d.reward.agg <- aggregate(cbind(Reward) ~ Effort + Competition + Sex +
                          Correct_Guess + ElapsedTime_Guess + ID_Player, 
                        data=df.reward, FUN = mean)

df_sum_times <- aggregate(ElapsedTime_Guess ~ ID_Player, data=d.reward.agg, FUN = sum)

#total time for each player
d.reward.agg$totaltime <- NA
for (i in unique(df_sum_times$ID_Player)){
              d.reward.agg$totaltime[d.reward.agg$ID_Player == i] <- 
                    df_sum_times$ElapsedTime_Guess[df_sum_times$ID_Player==i]
                }

#correct guesses per unit time
d.reward.agg.acc <- aggregate(cbind(Correct_Guess) ~ Effort + Competition + Sex + totaltime + ID_Player, 
                                                    data=d.reward.agg, FUN = sum)

#total number of problems that each player attempted
d.reward.agg.acc$totalgrids <- NA
for(i in unique(d.reward.agg$ID_Player)){
  ROWS <- nrow(d.reward.agg[d.reward.agg$ID_Player == i,])
  d.reward.agg.acc$totalgrids[d.reward.agg.acc$ID_Player==i] <- ROWS
}

#proportion of true results for each player
d.reward.agg.acc$trueresult_prop <- d.reward.agg.acc$Correct_Guess / d.reward.agg.acc$totalgrids

#######
###proportion of true results as a function of competition and effort
#######

m_trueresult_prop <- map2stan(
  alist(
    trueresult_prop ~ dnorm(mu, sigma), 
    mu <- a + bC*Competition + bE*Effort + 
      bCE*Competition*Effort, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10), 
    a ~ dgamma(1, 0.05), 
    sigma ~ dgamma(2, 0.5)
  ), data=d.reward.agg.acc, iter=5000, chains = 1, cores = 1, warmup=500)

plot(precis(m_trueresult_prop))

#######
###Reward per unit time as a function of competition and effort
#######

###Store original data in separate objects###
d.reward <- d.conf_1_2

#Exclude participants who indicated technical difficulties
d.reward <- d.reward[d.reward$ID_Player != 26,]
d.reward <- d.reward[d.reward$ID_Player != 61,]

#Exclude participants who did not complete study and rows with NA (i.e. grids that were never finished)
d.reward <- d.reward[complete.cases(d.reward),]

#aggregate data for time elapsed per guess, by first removing duplicate times, then summing up total times per player
d.reward.agg <- aggregate(ElapsedTime_Guess ~ Effort + Competition + Sex + Guess_Number + ID_Player, data = d.reward, FUN = mean)
d.reward.agg <- aggregate(ElapsedTime_Guess ~ Effort + Competition + Sex + ID_Player, data = d.reward.agg, FUN = sum)

#find total number of guesses for each participant and add to d.rewrd.agg
d.reward.agg$Reward <- NA
for(i in unique(d.reward.agg$ID_Player)){
  payout <- d.reward$Reward[d.reward$ID_Player == i]
  d.reward.agg$Reward[d.reward.agg$ID_Player == i] <- unique(payout)
}

#find total number of grids for each participant and add to d.rewrd.agg
d.reward.agg$totalgrids <- NA
for(i in unique(d.reward.agg$ID_Player)){
  grids <- d.reward$Guess_Number[d.reward$ID_Player == i]
  d.reward.agg$totalgrids[d.reward.agg$ID_Player == i] <- max(grids)
}

#calculate reward per unit time, for aggregated data (add 5 seconds for each grid complated)
d.reward.agg$totaltime <- d.reward.agg$ElapsedTime_Guess + d.reward.agg$totalgrids * 5
#calculate reward per unit time
d.reward.agg$reward_per_time <- d.reward.agg$Reward / d.reward.agg$totaltime

#model
m_reward_per_second <- map2stan(
  alist(
    reward_per_time ~ dnorm(mu, sigma), 
    mu <- a + bC*Competition + bE*Effort + bCE*Competition*Effort, 
    bC ~ dnorm(0, 0.1), 
    bE ~ dnorm(0, 0.1), 
    bCE ~ dnorm(0, 0.1), 
    sigma ~ dgamma(2, 1), 
    a ~ dgamma(1, 1)
  ), data=d.reward.agg, iter=5000, chains = 1, cores = 1, warmup=500)

#Plot Predictions
d.pred.a <- list(
  Competition = c(0, 0, 1, 1), 
  Effort = c(0, 1, 0, 1))

#replace varying intercept samples with zeros
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(d.reward.agg$ID_Player)))
m.reward.link <- link(m_reward_per_second, n=2000, data=d.pred.a, 
                        replace = list(a_player=a_player_zeros))

#summarize and plot with plot function
pred.p.mean.reward <- apply( m.reward.link , 2 , mean )
pred.p.PI.reward <- apply( m.reward.link , 2 , HPDI , prob=0.95)

d.gg <- data.frame(mean = pred.p.mean.reward, low_ci = pred.p.PI.reward[1,], high_ci = pred.p.PI.reward[2,], 
                   competition = as.factor(c(0, 0, 1, 1)), effort = as.factor(c(0, 1, 0, 1)))

#plot
d.reward.agg.plot <- d.reward.agg
d.reward.agg.plot$Effort <- as.factor(d.reward.agg.plot$Effort)

gg.reward <- ggplot(d.reward.agg.plot, aes(x = as.factor(Competition), y = reward_per_time, 
                                           color = Effort)) +
   geom_jitter(width = 0.25) + facet_grid(. ~ Effort) + 
   ylab("Payoff per Second") +
  scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes")) +
  theme_classic(base_size = 14) + theme(strip.text.x = element_blank()) + 
  scale_color_brewer(name="Effort", labels = c("No", "Yes"), palette = "Set1")

names(d.gg) <- c("reward_per_time", "low_ci", "high_ci", "Competition", "Effort")

gg.reward + geom_pointrange(data = d.gg, aes(ymin=low_ci, ymax=high_ci), lwd=0.8, 
                  colour = "black", alpha = 0.9) 
  
###compare alternative models###
m_reward_per_second_nointer <- map2stan(
  alist(
    reward_per_time ~ dnorm(mu, sigma), 
    mu <- a + bC*Competition + bE*Effort,  
    bC ~ dnorm(0, 0.1), 
    bE ~ dnorm(0, 0.1), 
    a ~ dgamma(1, 1), 
    sigma ~ dgamma(2, 1)
  ), data=d.reward.agg, iter=5000, chains = 1, cores = 1, warmup=500)

m_reward_per_second_noC <- map2stan(
  alist(
    reward_per_time ~ dnorm(mu, sigma), 
    mu <- a + bE*Effort,  
    bE ~ dnorm(0, 0.1), 
    a ~ dgamma(1, 1), 
    sigma ~ dgamma(2, 1)
  ), data=d.reward.agg, iter=5000, chains = 1, cores = 1, warmup=500)

reward_df <- compare(m_reward_per_second, m_reward_per_second_nointer, m_reward_per_second_noC)
reward_df <- round(reward_df@output, 2)
reward_df <- reward_df[,1:4]

row.names(reward_df) <- c("Reward_E",  
                           "Reward_E_C", "Reward_E_C_EC")

reward_df %>% kable(caption = "REWARD PER SECOND")  %>% 
  kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), 
                full_width = TRUE) %>%
  column_spec(1:5, width = "2cm")


##################
###Time to reveal 1 tile
#################

#Load data for quality checks
#load("XXXXX")

#remove missing observations
d.quality.effortcheck <- d.quality[complete.cases(d.quality),]

#exclude participants who indicated technical difficulties
d.quality.effortcheck <- d.quality.effortcheck[d.quality.effortcheck$ID_Player != 26,]
d.quality.effortcheck <- d.quality.effortcheck[d.quality.effortcheck$ID_Player != 61,]

#remove rows with time-to-remove 1 tile values that are more than 5 standard deviations away from the mean
agg <- aggregate(ElapsedTime_Tile ~ Guess_Number + ID_Player + Effort + Competition, data=d.quality.effortcheck, FUN = mean)
sd5.times <- mean(agg$ElapsedTime_Tile) + (5 * sd(agg$ElapsedTime_Tile))
nrow(agg[agg$ElapsedTime_Tile > sd5.times,]) #101 observations excluded

#excluding observations 
d.quality.effortcheck <- d.quality.effortcheck[d.quality.effortcheck$ElapsedTime_Tile < sd5.times,]
#making index
d.quality.effortcheck$ID_Player <- coerce_index(d.quality.effortcheck$ID_Player)


######
###Effect of Effort Treatment on Time to Click 1 Tile
######

m.quality.original.nooutliers <- map2stan(
  alist(
    ElapsedTime_Tile ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10), 
    a ~ dgamma(1.5, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.quality.effortcheck, iter=15000, chains=4, cores=4, warmup=500)

plot(precis(m.quality.original.nooutliers, prob = 0.95), 
     main = "Time (seconds) to Reveal 1 Tile")

m.quality.original.nointeraction <- map2stan(
  alist(
    ElapsedTime_Tile ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    a ~ dgamma(1.5, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.quality.effortcheck, iter=15000, chains=4, cores=4, warmup=500)

m.quality.original.nocomp <- map2stan(
  alist(
    ElapsedTime_Tile ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bE*Effort,
    bE ~ dnorm(0, 10), 
    a ~ dgamma(1.5, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=d.quality.effortcheck, iter=15000, chains=4, cores=4, warmup=500)

quality_df <- compare(m.quality.original.nocomp, m.quality.original.nointeraction, m.quality.original.nooutliers)
quality_df <- round(quality_df@output, 2)
quality_df <- quality_df[,1:4]

row.names(quality_df) <- c("Time_Tile_E",  
                         "Time_Tile_E_C_EC", "Time_Tile_E_C")

quality_df %>% kable(caption = "TIME TO REVEAL 1 TILE")  %>% 
  kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), 
                full_width = TRUE) %>%
  column_spec(1:5, width = "2cm")

###########
#plot tiles revealed as function of competition, effort, and effect size -----------
############

precis(m.tiles.guess.nmajor.inter, prob = 0.95)

EE <- unique(d.conf.agg$n_major.s)
median.guess.s <- median(d.conf.agg$Guess_Number.s)

d.pred <- list(
  Competition = c(rep(0, 6), rep(1, 6)),
  Effort = rep(c(0, 1), 6),
  n_major_s = c(rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2), rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2)),
  Guess_Number_s = rep(median.guess.s, 12),
  ID_Player = rep(2, 12) #placeholder
)

#replace varying intercept samples with zeros
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(d.conf.agg$ID_Player)))

m.quality1.link <- link(m.tiles.guess.nmajor.inter, n=2000, data=d.pred, 
                        replace = list(a_player=a_player_zeros))

#summarize#
pred.p.mean <- apply(m.quality1.link , 2 , mean)
pred.p.PI <- apply( m.quality1.link , 2 , HPDI , prob=0.95)

##plot
d.gg <- data.frame(mean = pred.p.mean, low_ci = pred.p.PI[1,], high_ci = pred.p.PI[2,], 
                   competition = c(rep(0, 6), rep(1, 6)),
                   effort = rep(c(0, 1), 6),
                   n_major_s = c(rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2), rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2)),
                   Guess_Number_s = rep(median.guess.s, 12)
                   )
#plot
d.gg$competition <- as.factor(d.gg$competition)
d.gg$competition <- c(rep("No Competition", 6), rep("Competition", 6))
d.gg$effort <- as.factor(d.gg$effort)
d.gg$effort <- rep(c("No Effort", "Effort"), 6)

ggplot(data=d.gg, aes(x = as.factor(n_major_s), y = mean, 
                      colour = as.factor(n_major_s), group=as.factor(n_major_s))) + 
  ylim(0, 15) + ylab("Tiles Revealed") +
  facet_grid(cols = vars(fct_rev(effort), fct_rev(competition))) +
  geom_pointrange(aes(ymin = low_ci, ymax=high_ci), lwd=1) +
  theme_bw(base_size = 14) +
  scale_x_discrete(name ="Effect Size", labels=c("Small","Medium","Large", "Small", "Medium", "Large")) +
  scale_color_brewer(name="Effect Size", palette = "Reds", 
                     labels = c("Small", "Medium", "Large"))

###########
#plot accuracy as function of competition, effort, and effect size -----------
############

precis(m_accuracy_compeff_guess_nmajor_inter, prob = 0.95)

EE <- unique(d.conf.agg$n_major.s)
median.guess.s <- median(d.conf.agg$Guess_Number.s)

d.pred <- list(
  Competition = c(rep(0, 6), rep(1, 6)),
  Effort = rep(c(0, 1), 6),
  n_major_s = c(rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2), rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2)),
  Guess_Number_s = rep(median.guess.s, 12),
  ID_Player = rep(2, 12) #placeholder
)

#replace varying intercept samples with zeros
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(d.conf.agg$ID_Player)))

m.acc.effect.link <- link(m_accuracy_compeff_guess_nmajor_inter, n=2000, data=d.pred, 
                        replace = list(a_player=a_player_zeros))

#summarize#
pred.p.mean <- apply(m.acc.effect.link , 2 , mean)
pred.p.PI <- apply( m.acc.effect.link , 2 , HPDI , prob=0.95)

##plot
d.gg <- data.frame(mean = pred.p.mean, low_ci = pred.p.PI[1,], high_ci = pred.p.PI[2,], 
                   competition = c(rep(0, 6), rep(1, 6)),
                   effort = rep(c(0, 1), 6),
                   n_major_s = c(rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2), rep(EE[1], 2), rep(EE[2], 2), rep(EE[3], 2)),
                   Guess_Number_s = rep(median.guess.s, 12)
)
#plot
d.gg$competition <- as.factor(d.gg$competition)
d.gg$competition <- c(rep("No Competition", 6), rep("Competition", 6))
d.gg$effort <- as.factor(d.gg$effort)
d.gg$effort <- rep(c("No Effort", "Effort"), 6)


##plot with ggplot
ggplot(data=d.gg, aes(x = as.factor(n_major_s), y = mean, 
                      colour = as.factor(n_major_s), group=as.factor(n_major_s))) + 
  ylim(0, 1) + ylab("Probability of Correct Guess") +
  facet_grid(cols = vars(fct_rev(effort), fct_rev(competition))) +
  geom_pointrange(aes(ymin = low_ci, ymax=high_ci), lwd=1) +
  theme_bw(base_size = 14) +
  scale_x_discrete(name ="Effect Size", labels=c("Small","Medium","Large", "Small", "Medium", "Large")) +
  scale_color_brewer(name="Effect Size", palette = "Reds", 
                     labels = c("Small", "Medium", "Large"))



