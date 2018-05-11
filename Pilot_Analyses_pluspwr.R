
rm(list=ls())

#########
library(tidyverse)
library(rethinking)
library(lme4)
library(lmerTest)
library(effects)
library(sjPlot)
library(stringr)

load("PilotData1.RData")
load("PilotData2_MathProbs.RData")
sequences <- read.table("Sequences1.txt", header=FALSE) #the sequence of blue and yellow tiles for each grid

#removing rows with missing observations (also removes participants who did not complete the study
#(i.e. had NA values for compcheck)
data.c <- (data[complete.cases(data),])
data_math.c <- (data_math[complete.cases(data_math),])

#simply histogram of time until guess
simplehist(data.c$ElapsedTime_Guess, xlab = "Time Until Guess (seconds), All Data.")

#remove times until guess that are more than 5 standard deviations away from the mean
sd5.times <- mean(data.c$ElapsedTime_Guess) + ( 5 * sd(data.c$ElapsedTime_Guess))
nrow(data.c[data.c$ElapsedTime_Guess > sd5.times,]) #number of observations excluded
data.c <- data.c[data.c$ElapsedTime_Guess < sd5.times,]
simplehist(data.c$ElapsedTime_Guess, xlab = "Time Until Guess (seconds), Outliers Exluded (5 SD)")

#simply histogram of arithmetic problem solving times
simplehist(data_math.c$ElapsedTime_Math, xlab = "Arithmetic Problem Solving Times (seconds), All Data.")

#remove math-problem times that are more than 3 standard deviations away from the mean
sd5 <- mean(data_math.c$ElapsedTime_Math) + ( 5 * sd(data_math.c$ElapsedTime_Math))
nrow(data_math.c[data_math.c$ElapsedTime_Math > sd5,]) #number of observations excluded
data_math.c <- data_math.c[data_math.c$ElapsedTime_Math < sd5,]
simplehist(data_math.c$ElapsedTime_Math, xlab = "Arithmetic Problem Solving Times (seconds), Outliers Exluded (5 SD)")


#copying data into separate objects for frequentist analyses
data.cf <- data.c
data_math.cf <- data_math.c

#changing sex and player ID to numeric
data.c$Sex <- as.numeric(data.c$Sex)
data_math.c$Sex <- as.numeric(data_math.c$Sex)
data.c$ID_Player <- as.numeric(data.c$ID_Player)
data_math.c$ID_Player <- as.numeric(data_math.c$ID_Player)

#changing variables to factors for frequentist analyses
factorfxn <- function(x) as.factor(x)
cols <- c("Sex", "Effort", "Competition", "CompCheck", "Correct_Guess", "Faster", "ID_Player")
data.cf[,cols] <- lapply(data.cf[,cols], factorfxn)
cols.m <- c("Sex", "Effort", "Competition", "CompCheck", "Correct_Guess", "Correct_Math", "Treatment", "ID_Player")
data_math.cf[,cols.m] <- lapply(data_math.cf[,cols.m], factorfxn)

#for each of the 600 possible sequences, the number of majority tiles (i.e. effect size)
seq.df <- data.frame(seq = sequences)
t <- str_split(seq.df$V1, ";")
t <- do.call(rbind, t)
t <- as.data.frame(t)

for(i in 1:ncol(t)){
  t[,i] <- as.numeric(t[,i]) - 1
}
n_majority <- apply(t, 1, sum)
#the number of majority tiles, only for those grids that had at least 1 player attempt them
n_majority <- n_majority[1:length(unique(data.c$Guess_Number))]

#identifies which grids correspond to which effect size
g_13 <- which(n_majority == 13)
g_15 <- which(n_majority == 15)
g_17 <- which(n_majority == 17)

#adds this to all data sets
data.c$n_major <- rep(NA, nrow(data.c))
data.cf$n_major <- rep(NA, nrow(data.cf))
data_math.c$n_major <- rep(NA, nrow(data_math.c))
data_math.cf$n_major <- rep(NA, nrow(data_math.cf))

for(i in 1:nrow(data.c)){
  if(data.c$Guess_Number[i] %in% g_13){data.c$n_major[i] <- 13}
  if(data.c$Guess_Number[i] %in% g_15){data.c$n_major[i] <- 15}
  if(data.c$Guess_Number[i] %in% g_17){data.c$n_major[i] <- 17}
}
for(i in 1:nrow(data.cf)){
  if(data.cf$Guess_Number[i] %in% g_13){data.cf$n_major[i] <- 13}
  if(data.cf$Guess_Number[i] %in% g_15){data.cf$n_major[i] <- 15}
  if(data.cf$Guess_Number[i] %in% g_17){data.cf$n_major[i] <- 17}
}
for(i in 1:nrow(data_math.c)){
  if(data_math.c$Guess_Number[i] %in% g_13){data_math.c$n_major[i] <- 13}
  if(data_math.c$Guess_Number[i] %in% g_15){data_math.c$n_major[i] <- 15}
  if(data_math.c$Guess_Number[i] %in% g_17){data_math.c$n_major[i] <- 17}
}
for(i in 1:nrow(data_math.cf)){
  if(data_math.cf$Guess_Number[i] %in% g_13){data_math.cf$n_major[i] <- 13}
  if(data_math.cf$Guess_Number[i] %in% g_15){data_math.cf$n_major[i] <- 15}
  if(data_math.cf$Guess_Number[i] %in% g_17){data_math.cf$n_major[i] <- 17}
}

#standardizing n_major
data.c$n_major.s <- (data.c$n_major - mean(data.c$n_major)) / sd(data.c$n_major)
data.cf$n_major.s <- (data.cf$n_major - mean(data.cf$n_major)) / sd(data.cf$n_major)
data_math.c$n_major.s <- (data_math.c$n_major - mean(data_math.c$n_major)) / sd(data_math.c$n_major)
data_math.cf$n_major.s <- (data_math.cf$n_major - mean(data_math.cf$n_major)) / sd(data_math.cf$n_major)

#change to factor for frequentist analyses
data.cf$n_major <- as.factor(data.cf$n_major)
data_math.cf$n_major.s <- as.factor(data_math.cf$n_major)
data.cf$n_major.s <- as.factor(data.cf$n_major.s)
data_math.cf$n_major.s <- as.factor(data_math.cf$n_major.s)

###################################
#Bayesian analyses (Rethinking) followed by frequentist implementations of the same model (lmer)
###################################

#Quality Checks#
df_comp <- aggregate(cbind(CompCheck, Competition) ~ ID_Player,
                                  data=data.c, FUN=unique)
q2 <- map2stan(
  alist(
    CompCheck ~ dbinom(1, theta), 
    logit(theta) <- a + bC*Competition, 
    bC ~ dnorm(0, 10), 
    a ~ dnorm(0, 10)
  ), data=df_comp, iter=10000, chains=3, warmup=1000)

plot(q2)
pairs(q2)
par(mfrow=c(1,1))
precis(q2, prob=0.95)
plot(precis(q2, prob=0.95), xlab = "Log Odds of Answering 'Yes'")
logistic(5)
logistic(-2.37)

#Frequentist
#Model fails to converge because the predictor (competition) almost perfectly
#separates the outcome (CompCheck1). All individuals in the competition condition 
#responded "yes" and the model breaks down as it considers log-odds values up to
#infinity as plausible. 

df_comp.f <- df_comp
df_comp.f$Competition <- as.factor(df_comp.f$Competition)
df_comp.f$CompCheck <- as.factor(df_comp.f$CompCheck)

#For illustrative purposes only: estimation problem will 
# disappear if you replace one of the "yes" answers with a "no" in  competition treatment
# as follows: df_comp.f$CompCheck[33] <- as.factor(0)

q2.f <- glm(CompCheck ~ Competition, family=binomial, data=df_comp.f)
summary(q2.f)

sjt.glm(q2.f,
          show.icc=FALSE, show.col.header = TRUE,
          string.est = "Estimate",
          string.ci = "CI", 
          string.p = "P", 
          separate.ci.col = FALSE, 
          group.pred = TRUE)
plot(allEffects(q2.f))


#######
#Effect of competition and effort on accuracy
#######

m.3.b <- map2stan(
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
  ), data=data.c, iter=10000, chains=3, cores = 4, warmup=1000)

plot(m.3.b)
pairs(m.3.b)
par(mfrow=c(1,1))
precis(m.3.b, prob=0.95)
plot(precis(m.3.b, prob=0.95), xlab = "Log odds of Correct Guess")

#Frequentist
m.3.f <- glmer(Correct_Guess ~ Competition*Effort + n_major.s + (1|ID_Player), 
               data=data.cf, family=binomial)
sjt.glmer(m.3.f,
          show.icc=FALSE, show.col.header = TRUE,
          string.est = "Estimate",
          string.ci = "CI", 
          string.p = "P", 
          separate.ci.col = FALSE, 
          group.pred = TRUE, 
          exp.coef = FALSE, 
          file = "CorrectGuess_LogOdds.doc")
plot(allEffects(m.3.f))
plot(effect("Competition*Effort", m.3.f))
summary(m.3.f)

logistic(.99)

#Effect of competition and effort on time to guess 
m.1.b <- map2stan(
  alist(
    ElapsedTime_Guess ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 30), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1.5, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=data.c, iter=10000, chains=3, warmup=1000)

plot(m.1.b)
pairs(m.1.b, pars=c("bC", "bE", "bCE", "bNs", "a", "sigma"))
par(mfrow=c(1,1))
precis(m.1.b, prob = 0.95)
plot(precis(m.1.b, prob = 0.95), xlab = "Time Until Guess")

#Frequentist
m.1.f <- lmer(ElapsedTime_Guess ~ Competition*Effort + n_major.s + (1|ID_Player),
              data=data.cf, REML=FALSE)
sjt.lmer(m.1.f,
         show.icc=FALSE, show.col.header = TRUE,
         string.est = "Estimate",
         string.ci = "CI", 
         string.p = "P", 
         separate.ci.col = FALSE, 
         group.pred = TRUE, 
         file = "ElapsedTime.doc")
plot(effect("Competition*Effort", m.1.f))

#######
#Effect of competition and effort on number of tiles revealed
######

m.2.b <- map2stan(
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
  ), data=data.c, iter=10000, chains=3, warmup=1000)

plot(m.2.b)
pairs(m.2.b)
par(mfrow=c(1,1))
precis(m.2.b, prob = 0.95)
plot(precis(m.2.b, prob = 0.95), xlab = "Tiles Revealed")

#Frequentist
m.2.f <- lmer(TilesRevealed ~ Competition*Effort + n_major.s + (1|ID_Player), 
              data=data.cf, REML=FALSE)
sjt.lmer(m.2.f,
         show.icc=FALSE, show.col.header = TRUE,
         string.est = "Estimate",
         string.ci = "CI", 
         string.p = "P", 
         separate.ci.col = FALSE, 
         group.pred = TRUE, 
         file = "TilesRevealed.doc")
plot(allEffects(m.2.f))
plot(effect("Competition*Effort", m.2.f))
summary(m.2.f)

#######
#Effect of competition and effort on quantity (i.e. rate of guessing the majority color)
#######
data.rate <- aggregate(cbind(ElapsedTime_Guess, Correct_Guess) ~ ID_Player,
          data=data.c, FUN=sum)

#mean problem effect size (i.e. n_major)
mean_n_major <- aggregate(n_major ~ ID_Player, data=data.c, FUN=mean)
data.rate$mean_n_major <- mean_n_major$n_major

a <- aggregate(Correct_Guess ~ ID_Player,data=data.c, FUN=length)
data.rate$total_guesses <- a$Correct_Guess
data.rate$totalminutes <- data.rate$ElapsedTime_Guess / 60

#adjusting time of participation to account for 5 second delay between grids
data.rate$totalminutes <- data.rate$totalminutes + ((4/60) * data.rate$total_guesses)
############################

data.rate$guessrate <- data.rate$total_guesses / data.rate$totalminutes
data.rate$correct_guessrate <- data.rate$Correct_Guess / data.rate$totalminutes

#extract which condition each player was in
data.rate$Competition <- rep(NA, nrow(data.rate))
data.rate$Effort <- rep(NA, nrow(data.rate))
for(i in unique(data.c$ID_Player)) {
  dat <- data.c[data.c$ID_Player==i,]
    data.rate$Competition[data.rate$ID_Player==i] <- unique(dat$Competition)
    data.rate$Effort[data.rate$ID_Player==i] <- unique(dat$Effort)
    
}

#rate of acquiring points
data.rate$wrong_guess <- data.rate$total_guesses - data.rate$Correct_Guess
data.rate$points <- data.rate$Correct_Guess - data.rate$wrong_guess
data.rate$points_rate <- data.rate$points / data.rate$totalminutes

#standardizing mean_n_major
data.rate$mean_n_major.s <- (data.rate$mean_n_major - mean(data.rate$mean_n_major)) / sd(data.rate$mean_n_major)

#Model
m.4.b <- map2stan(
  alist(
    guessrate ~ dnorm(mu, sigma), 
    mu <- a + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*mean_n_major.s, 
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dgamma(2, 0.05), 
    sigma ~ dgamma(2, 0.5)
  ), data=data.rate, iter=10000, chains=3, warmup=1000)

plot(m.4.b)
pairs(m.4.b)
precis(m.4.b)
par(mfrow=c(1,1))
plot(precis(m.4.b, prob=0.95), xlab = "Guesses Per Minute")

#Frequentist
data.rate.f <- data.rate
data.rate.f$ID_Player <- as.factor(data.rate.f$ID_Player)
data.rate.f$Correct_Guess <- as.factor(data.rate.f$Correct_Guess)
data.rate.f$Competition <- as.factor(data.rate.f$Competition)
data.rate.f$Effort <- as.factor(data.rate.f$Effort)

m.4.f <- lm(guessrate ~ Competition*Effort + mean_n_major.s, 
            data=data.rate.f)
sjt.lm(m.4.f,
       show.icc=FALSE, show.col.header = TRUE,
       string.est = "Estimate",
       string.ci = "CI", 
       string.p = "P", 
       separate.ci.col = FALSE, 
       group.pred = TRUE, 
       file = "Guessrate.doc")
plot(allEffects(m.4.f))
plot(effect("Competition*Effort", m.4.f))


#identical model without mean_n_major
m.4.f.no_nmajor <- lm(guessrate ~ Competition*Effort, 
                      data=data.rate.f)
summary(m.4.f.no_nmajor)
plot(allEffects(m.4.f.no_nmajor))

#######
#Effect of competition on effort
#######

#Rate#

#percentage instances when math problems answered correctly
nrow(data_math.c[data_math.c$Correct_Math==1,]) / 
  nrow(data_math.c) #96.5 percent of math answers were correct in pilot study

#rate of solving arithmetic problems
data.mathrate <- aggregate(cbind(ElapsedTime_Math, Correct_Math) ~ ID_Player,
                       data=data_math.c, FUN=sum)

#mean problem effect size (i.e. n_major)
mean_n_major_math <- aggregate(n_major ~ ID_Player, data=data_math.c, FUN=mean)
data.mathrate$mean_n_major <- mean_n_major_math$n_major

data.mathrate$ElapsedTime_Math_Minutes <- data.mathrate$ElapsedTime_Math / 60
data.mathrate$mathprobs_perminute <- data.mathrate$Correct_Math / data.mathrate$ElapsedTime_Math_Minutes

#extract which condition each player was in
data.mathrate$Competition <- rep(NA, nrow(data.mathrate))
data.mathrate$Effort <- rep(NA, nrow(data.mathrate))
for(i in unique(data.c$ID_Player)) {
  dat <- data.c[data.c$ID_Player==i,]
  data.mathrate$Competition[data.mathrate$ID_Player==i] <- unique(dat$Competition)
  data.mathrate$Effort[data.mathrate$ID_Player==i] <- unique(dat$Effort)
}
#standardizing mean_n_major
data.mathrate$mean_n_major.s <- (data.mathrate$mean_n_major - mean(data.mathrate$mean_n_major)) / sd(data.mathrate$mean_n_major)

#Model
m.5.b <- map2stan(
  alist(
    mathprobs_perminute ~ dnorm(mu, sigma), 
    mu <- a + bC*Competition + bNs*mean_n_major.s, 
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1.5, 0.05), 
    sigma ~ dgamma(2, 0.5)
  ), data=data.mathrate, iter=10000, chains=3, warmup=1000)

plot(m.5.b)
pairs(m.5.b)
par(mfrow=c(1,1))
precis(m.5.b, prob=0.95)
plot(precis(m.5.b, prob=0.95), xlab = "Arithmetic Problems Solved Per Minute")

#Frequentist
data.mathrate.f <- data.mathrate
data.mathrate.f$Competition <- as.factor(data.mathrate.f$Competition)
data.mathrate.f$ID_Player <- as.factor(data.mathrate.f$ID_Player)

m.5.f <- lm(mathprobs_perminute ~ Competition + mean_n_major.s, 
              data=data.mathrate.f)
sjt.lm(m.5.f,
         show.icc=FALSE, show.col.header = TRUE,
         string.est = "Estimate",
         string.ci = "CI", 
         string.p = "P", 
         separate.ci.col = FALSE, 
         group.pred = TRUE, 
         file = "math_per_min.doc")

plot(allEffects(m.5.f))
plot(effect("Competition", m.5.f))

#save workspace
save.image(file = "Pilot_Analysis_UpdatedRR.RData")

###################
#Power Analysis
###################

rm(list=ls()) #remove previous objects
load("Pilot_Analysis_UpdatedRR.RData") # if not already loaded in working memory

#####
#Effect of competition and effort on accuracy
#use model m.3.b
####

ropelow <- -0.2
ropehigh <- 0.2

df_rope <- data.frame(c_outside = rep(0, 200), ce_outside = rep(0, 200), c_plusce_outside = rep(0, 200), 
                      c_inside = rep(0, 200), ce_inside = rep(0, 200), c_plusce_inside = rep(0, 200))

#extract samples from posterior for all parameters
post <- extract.samples(m.3.b)

#creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(z in 1:200) {

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
  ), data=df, iter=5000, warmup=1000, chains = 4, cores = 4)

samples <- extract.samples(m.3.power)
beta_c <- HPDI(samples$bC, prob = 0.95)
beta_c_plusce <- HPDI(samples$bC + samples$bCE, prob = 0.95)
beta_ce <- HPDI(samples$bCE, prob = 0.95)

#competition beta
if(all(beta_c < ropelow) | all(beta_c > ropehigh)) {
  df_rope_tiles$betac_outside[runs] <- 1
}
if(all(beta_c > ropelow) & all(beta_c < ropehigh)) {
  df_rope_tiles$betac_inside[runs] <- 1
}

#difference between no comp x effort and comp x effort
if(all(beta_c_plusce < ropelow) | all(beta_c_plusce > ropehigh)) {
  df_rope_tiles$c_plusce_outside[runs] <- 1
}
if(all(beta_c_plusce > ropelow) & all(beta_c_plusce < ropehigh)) {
  df_rope_tiles$c_plusce_outside[runs] <- 1
}

#interaction between competition and effort
if(all(beta_ce < ropelow) | all(beta_ce > ropehigh)) {
  df_rope_tiles$ce_outside[runs] <- 1
}
if(all(beta_ce > ropelow) & all(beta_ce < ropehigh)) {
  df_rope_tiles$ce_inside[runs] <- 1
}

#saving rope tracker, then removing model for power
save(df_rope, "df_rope_accuracy.RData")

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

#####
#Effect of competition and effort on time until guessing 
#use model m.1.b
####

ropelow <- -1
ropehigh <- 1

df_rope_time <- data.frame(c_outside = rep(0, 200), e_outside = rep(0, 200), ce_outside = rep(0, 200), 
                      c_inside = rep(0, 200), e_inside = rep(0, 200), ce_inside = rep(0, 200))

#extract samples from posterior for all parameters
post <- extract.samples(m.1.b)

#creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(runs in 1:200) {
  
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
  beta_c <- HPDI(samples$bC, prob = 0.95)
  beta_c_plusce <- HPDI(samples$bC + samples$bCE, prob = 0.95)
  
  #competition beta
  if(all(beta_c < ropelow) | all(beta_c > ropehigh)) {
    df_rope_time$c_outside[1] <- 1
  }
  if(all(beta_c > ropelow) & all(beta_c < ropehigh)) {
    df_rope_time$c_inside[1] <- 1
  }
  #competition*effort beta
  if(all(beta_c_plusce < ropelow) | all(beta_c_plusce > ropehigh)) {
    df_rope_time$ce_outside[1] <- 1
  }
  if(all(beta_c_plusce > ropelow) & all(beta_c_plusce < ropehigh)) {
    df_rope_time$ce_inside[1] <- 1
  }
  
  #save rope cound, then remove already compiled models
  save(df_rope_time, "df_rope_time.RData")
  
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

#####
#Effect of competition and effort on tiles
#use model m.2.b
####

ropelow <- -0.5
ropehigh <- 0.5

df_rope_tiles <- data.frame(c_outside = rep(0, 200), e_outside = rep(0, 200), ce_outside = rep(0, 200), 
                            c_inside = rep(0, 200), e_inside = rep(0, 200), ce_inside = rep(0, 200))

#extract samples from posterior for all parameters
post <- extract.samples(m.2.b)

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
    
    df_guess_full$mean_tiles[i] <- df_guess_full$a[i] + df_guess_full$a_player[i] + post$bC[i]*df_guess_full$competition[i] + 
      post$bE[i]*df_guess_full$effort[i] + post$bCE[i]*df_guess_full$competition[i]*df_guess_full$effort[i] + 
      post$bNs[i]*df_guess_full$n_major.s[i]
  }
  
  #generate individual-level data for each participant
  df_guess_full$TilesRevealed <- rnorm(nrow(df_guess_full), mean = df_guess_full$mean_tiles, sd = post$sigma)
  
  #round tiles less than 0 to be 0
  df_guess_full$TilesRevealed[df_guess_full$TilesRevealed < 0] <- 0
  
  #subset data frame for analysis
  df <- df_guess_full[,c("id", "competition", "effort", "n_major.s", "TilesRevealed")]
  
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
    ), data=df, iter=5000, chains=4, cores = 4, warmup=1000)
  
  samples <- extract.samples(m.2.power)
  beta_c <- HPDI(samples$bC, prob = 0.95)
  beta_c_plusce <- HPDI(samples$bC + samples$bCE, prob = 0.95)
  
  #competition beta
  if(all(beta_c < ropelow) | all(beta_c > ropehigh)) {
    df_rope_tiles$c_outside[runs] <- 1
  }
  if(all(beta_c > ropelow) & all(beta_c < ropehigh)) {
    df_rope_tiles$c_inside[runs] <- 1
  }
  #competition*effort beta
  if(all(beta_c_plusce < ropelow) | all(beta_c_plusce > ropehigh)) {
    df_rope_tiles$ce_outside[runs] <- 1
  }
  if(all(beta_c_plusce > ropelow) & all(beta_c_plusce < ropehigh)) {
    df_rope_tiles$ce_inside[runs] <- 1
  }
  
  #save rope tracker, then remove already compiled model
  save(df_rope_tiles, "df_rope_tiles.RData")
  
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
#Effect of competition and effort on guess rate
#use model m.4.b
####

ropelow <- -0.5
ropehigh <- 0.5

df_rope_guessrate <- data.frame(c_outside = rep(0, 200), c_plusce_outside = rep(0, 200), 
                            c_inside = rep(0, 200), c_plusce_inside = rep(0, 200))

#extract samples from posterior for all parameters
post <- extract.samples(m.4.b)

#creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(runs in 1:2) {
  
  #generate empty data frame to store simulated data
  df_guess_full <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), 
                         effort = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)), mean_n_major.s = NA, guessrate = NA)
  
  #generate n_major.s values for each problem that a participant solves and store in the full data
  df_guess_full$mean_n_major.s <- sample(data.rate$mean_n_major.s, size = nrow(df_guess_full), replace = TRUE)
  
  
  #generate simulated mean guessrate for 200 players
  df_guess_full$mean_guessrate <- NA
  
  for(i in 1:nrow(df_guess_full)) {
    
    df_guess_full$mean_guessrate[i] <- post$a[i] + post$bC[i]*df_guess_full$competition[i] + 
      post$bE[i]*df_guess_full$effort[i] + post$bCE[i]*df_guess_full$competition[i]*df_guess_full$effort[i] + 
      post$bNs[i]*df_guess_full$mean_n_major.s[i]
  }

  #generate individual-level data for each participant
  df_guess_full$guessrate <- rnorm(nrow(df_guess_full), mean = df_guess_full$mean_guessrate, sd = post$sigma)
  
  #round guessrates less than 0 to be 0
  df_guess_full$guessrate[df_guess_full$guessrate < 0] <- 0
  
  #subset data frame for analysis
  df <- df_guess_full[,c("id", "competition", "effort", "mean_n_major.s", "guessrate")]
  
  m.4.power <- map2stan(
    alist(
      guessrate ~ dnorm(mu, sigma), 
      mu <- a + bC*competition + bE*effort + bCE*competition*effort + bNs*mean_n_major.s, 
      bC ~ dnorm(0, 10), 
      bE ~ dnorm(0, 10), 
      bCE ~ dnorm(0, 10),
      bNs ~ dnorm(0, 10),
      a ~ dgamma(2, 0.05), 
      sigma ~ dgamma(2, 0.5)
    ), data=df, iter=5000, chains=4, cores = 4, warmup=1000)
  
  samples <- extract.samples(m.4.power)
  beta_c <- HPDI(samples$bC, prob = 0.95)
  beta_c_plusce <- HPDI(samples$bC + samples$bCE, prob = 0.95)
  
  #competition beta
  if(all(beta_c < ropelow) | all(beta_c > ropehigh)) {
    df_rope_guessrate$c_outside[runs] <- 1
  }
  if(all(beta_c > ropelow) & all(beta_c < ropehigh)) {
    df_rope_guessrate$c_inside[runs] <- 1
  }
  #difference between competition x no effort and competition x effort
  if(all(beta_c_plusce < ropelow) | all(beta_c_plusce > ropehigh)) {
    df_rope_guessrate$c_plusce_outside[runs] <- 1
  }
  if(all(beta_c_plusce > ropelow) & all(beta_c_plusce < ropehigh)) {
    df_rope_guessrate$c_plusce_inside[runs] <- 1
  }
  
  #save rope tracker, then remove already compiled model
  save(df_rope_guessrate, file = "df_rope_guessrate.RData")
  
  #shuffle samples in post
  ran_pos <- sample(1:nrow(post$a), nrow(post$a), replace = FALSE)
  post <- lapply(post, function(a) a[ran_pos])
}

########
#Effect of competition on rate of solving arithmetic problems
#use model m.5.b

ropelow <- -0.2
ropehigh <- 0.2

df_rope_mathrate <- data.frame(c_outside = rep(0, 200), c_plusce_outside = rep(0, 200), 
                                c_inside = rep(0, 200), c_plusce_inside = rep(0, 200))

#extract samples from posterior for all parameters
post <- extract.samples(m.5.b)

#creating 200 simulated data sets, analyzing, and comparing their CI's to the ROPE
for(runs in 1:2) {
  
  #generate empty data frame to store simulated data
  df_guess_full <- data.frame(id = 1:200, competition = c(rep(0, 100), rep(1, 100)), 
                              effort = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)), mean_n_major.s = NA, mathrate = NA)
  
  #generate n_major.s values for each problem that a participant solves and store in the full data
  df_guess_full$mean_n_major.s <- sample(data.mathrate$mean_n_major.s, size = nrow(df_guess_full), replace = TRUE)
  
  #generate simulated mean guessrate for 200 players
  df_guess_full$mean_mathrate <- NA
  
  for(i in 1:nrow(df_guess_full)) {
    
    df_guess_full$mean_mathrate[i] <- post$a[i] + post$bC[i]*df_guess_full$competition[i] + 
                                                          post$bNs[i]*df_guess_full$mean_n_major.s[i]
  }
  #generate individual-level data for each participant
  df_guess_full$mathrate <- rnorm(nrow(df_guess_full), mean = df_guess_full$mean_mathrate, sd = post$sigma)
  
  #round guessrates less than 0 to be 0
  df_guess_full$mathrate[df_guess_full$mathrate < 0] <- 0
  
  #subset data frame for analysis
  df <- df_guess_full[,c("id", "competition", "effort", "mean_n_major.s", "mathrate")]
  
  m.5.power <- map2stan(
    alist(
      mathrate ~ dnorm(mu, sigma), 
      mu <- a + bC*competition + bNs*mean_n_major.s, 
      bC ~ dnorm(0, 10), 
      bNs ~ dnorm(0, 10),
      a ~ dgamma(1.5, 0.05), 
      sigma ~ dgamma(2, 0.5)
    ), data=df, iter=5000, chains=4, cores = 4, warmup=1000)
  
######STOPPED HERE
  ###to do is to 1) finish the rest of this section and 2) find a way to recompile without crashing R session
  
  samples <- extract.samples(m.4.power)
  beta_c <- HPDI(samples$bC, prob = 0.95)
  beta_c_plusce <- HPDI(samples$bC + samples$bCE, prob = 0.95)
  
  #competition beta
  if(all(beta_c < ropelow) | all(beta_c > ropehigh)) {
    df_rope_guessrate$c_outside[runs] <- 1
  }
  if(all(beta_c > ropelow) & all(beta_c < ropehigh)) {
    df_rope_guessrate$c_inside[runs] <- 1
  }
  #difference between competition x no effort and competition x effort
  if(all(beta_c_plusce < ropelow) | all(beta_c_plusce > ropehigh)) {
    df_rope_guessrate$c_plusce_outside[runs] <- 1
  }
  if(all(beta_c_plusce > ropelow) & all(beta_c_plusce < ropehigh)) {
    df_rope_guessrate$c_plusce_inside[runs] <- 1
  }
  
  #save rope tracker, then remove already compiled model
  save(df_rope_guessrate, file = "df_rope_guessrate.RData")
  rm(m.4.power)
  
  #shuffle samples in post
  ran_pos <- sample(1:nrow(post$a), nrow(post$a), replace = FALSE)
  post <- lapply(post, function(a) a[ran_pos])
}
























##############
#Supplemental
##############
# 
# ####Exploratory######
# 
# #limiting analysis to instances when math problems answered correctly
# data_math.c.correct <- data_math.c[data_math.c$Correct_Math==1,]
# 
# m.6.b <- map2stan(
#   alist(
#     ElapsedTime_Math ~ dnorm(mu, sigma), 
#     mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s, 
#     bC ~ dnorm(0, 10), 
#     bNs ~ dnorm(0, 10),
#     a ~ dgamma(1.5, 0.05), 
#     a_player[ID_Player] ~ dnorm(0, sigma_player),
#     sigma_player ~ dgamma(1.5, 0.05), 
#     sigma ~ dgamma(2, 0.5)
#   ), data=data_math.c.correct, iter=10000, chains=3, warmup=1000)
# 
# plot(m.6.b)
# pairs(m.6.b, pars=c("bC", "bNs", "a", "sigma"))
# par(mfrow=c(1,1))
# precis(m.6.b, prob=0.95)
# plot(precis(m.6.b, prob=0.95), xlab = "Time (seconds) Spent on Accurately-Solved Arithmetic Problems")
# 
# #Frequentist
# data_math.c.correct.f <- data_math.c.correct
# data_math.c.correct.f$Competition <- as.factor(data_math.c.correct.f$Competition)
# data_math.c.correct.f$ID_Player <- as.factor(data_math.c.correct.f$ID_Player)
# 
# m.6.f <- lmer(ElapsedTime_Math ~ Competition + n_major.s + (1|ID_Player), 
#               data=data_math.c.correct.f, REML=FALSE)
# sjt.lmer(m.6.f,
#          show.icc=FALSE, show.col.header = TRUE,
#          string.est = "Estimate",
#          string.ci = "CI", 
#          string.p = "P", 
#          separate.ci.col = FALSE, 
#          group.pred = TRUE, 
#          file = "ElapsedTime_Math.doc")
# 
# plot(allEffects(m.6.f))
# plot(effect("Competition", m.6.f))
# 
# #Rate of Accurate Guesses
# m.7.b <- map2stan(
#   alist(
#     correct_guessrate ~ dnorm(mu, sigma), 
#     mu <- a + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*mean_n_major.s, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     bNs ~ dnorm(0, 10),
#     a ~ dgamma(1.5, 0.05), 
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.rate, iter=10000, chains=3, warmup=1000)
# 
# plot(m.7.b)
# pairs(m.7.b)
# precis(m.7.b)
# par(mfrow=c(1,1))
# plot(precis(m.7.b, prob=0.95), xlab = "Accurate Guesses per Minute")
# 
# #Frequentist
# data.rate.f <- data.rate
# data.rate.f$Competition <- as.factor(data.rate.f$Competition)
# data.rate.f$ID_Player <- as.factor(data.rate.f$ID_Player)
# data.rate.f$Effort <- as.factor(data.rate.f$Effort)
# 
# m.7.f <- lm(correct_guessrate ~ Competition*Effort + mean_n_major.s, 
#               data=data.rate.f)
# sjt.lm(m.7.f,
#          show.icc=FALSE, show.col.header = TRUE,
#          string.est = "Estimate",
#          string.ci = "CI", 
#          string.p = "P", 
#          separate.ci.col = FALSE, 
#          group.pred = TRUE, 
#        file = "CorrectGuessRate.doc")
# 
# plot(allEffects(m.7.f))
# plot(effect("Competition*Effort", m.7.f))
# 
# #Rate of Acquiring Points
# m.8.b <- map2stan(
#   alist(
#     points_rate ~ dnorm(mu, sigma), 
#     mu <- a + bC*Competition + bE*Effort + bCE*Competition*Effort + bNs*mean_n_major.s, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     bNs ~ dnorm(0, 10),
#     a ~ dgamma(1.5, 0.05), 
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.rate, iter=10000, chains=3, warmup=1000)
# 
# plot(m.8.b)
# pairs(m.8.b)
# precis(m.8.b)
# par(mfrow=c(1,1))
# plot(precis(m.8.b, prob=0.95), xlab = "Points per Minute")
# 
# #Frequentist
# m.8.f <- lm(points_rate ~ Competition*Effort + mean_n_major.s, 
#             data=data.rate.f)
# sjt.lm(m.8.f,
#        show.icc=FALSE, show.col.header = TRUE,
#        string.est = "Estimate",
#        string.ci = "CI", 
#        string.p = "P", 
#        separate.ci.col = FALSE, 
#        group.pred = TRUE, 
#        file = "PointsRate.doc")
# 
# plot(allEffects(m.8.f))
# plot(effect("Competition*Effort", m.8.f))
# 
# #################
# #Alternate (bad) model specifications 
# #result in high levels of autocorelaton between samples for at least 1 of the parameters
#####################################
# 
# m <- map2stan(
#   alist(
#     ElapsedTime_Guess ~ dnorm(mu, sigma), 
#     mu <- a + a_player[ID_Player] + a_effect[n_major] + bC*Competition + bE*Effort + bCE*Competition*Effort, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     a ~ dgamma(1.5, 0.05), 
#     a_player[ID_Player] ~ dnorm(0, sigma_player),
#     sigma_player ~ dgamma(1.5, 0.05),
#     a_effect[n_major] ~ dnorm(0, sigma_n_major), 
#     sigma_n_major ~ dgamma(1.5, 0.05),
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.c, iter=3000, chains=1, warmup=1000)
# 
# plot(m)
# pairs(m, pars=c("bC", "bE", "bCE", "a", "sigma_n_major", "a_effect"))
# 
# g <- map2stan(
#   alist(
#     ElapsedTime_Guess ~ dnorm(mu, sigma), 
#     mu <- a + a_player[ID_Player] + a_guess[Guess_Number] + bC*Competition + bE*Effort + bCE*Competition*Effort, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     a ~ dgamma(1.5, 0.05), 
#     a_player[ID_Player] ~ dnorm(0, sigma_player),
#     sigma_player ~ dgamma(1.5, 0.05),
#     a_guess[Guess_Number] ~ dnorm(0, sigma_guess), 
#     sigma_guess ~ dgamma(1.5, 0.05),
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.c, iter=3000, chains=1, warmup=1000)
# 
# plot(g)
# pairs(g, pars=c("bC", "bE", "bCE", "a", "sigma_n_major", "a_effect"))
# 
# z <- map2stan(
#   alist(
#     ElapsedTime_Guess ~ dnorm(mu, sigma), 
#     mu <- a + a_player[ID_Player] + a_effect[n_major] + bC*Competition + bE*Effort + bCE*Competition*Effort, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     a ~ dgamma(1.5, 0.05), 
#     a_player[ID_Player] ~ dnorm(0, sigma_player),
#     sigma_player ~ dgamma(1.5, 0.05),
#     a_effect[n_major] ~ dnorm(0, 10), 
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.c, iter=3000, chains=1, warmup=1000)
# 
# plot(z)
# pairs(z, pars=c("bC", "bE", "bCE", "a", "a_effect", "sigma"))
# precis(z)
# 
# zz <- map2stan(
#   alist(
#     ElapsedTime_Guess ~ dnorm(mu, sigma), 
#     mu <- a_player[ID_Player] + a_effect[n_major] + bC*Competition + bE*Effort + bCE*Competition*Effort, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     a_player[ID_Player] ~ dnorm(0, sigma_player),
#     sigma_player ~ dgamma(1.5, 0.05),
#     a_effect[n_major] ~ dnorm(0, 10), 
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.c, iter=3000, chains=1, warmup=1000)
# 
# plot(zz)
# pairs(zz, pars=c("bC", "bE", "bCE", "a_effect", "sigma"))
# precis(zz)
# 
# guessrate_bad <- map2stan(
#   alist(
#     guessrate ~ dnorm(mu, sigma), 
#     mu <- a + bC*Competition + bE*Effort + bCE*Competition*Effort + bN*mean_n_major, 
#     bC ~ dnorm(0, 10), 
#     bE ~ dnorm(0, 10), 
#     bCE ~ dnorm(0, 10),
#     bN ~ dnorm(0, 10),
#     a ~ dnorm(0, 10), 
#     sigma ~ dgamma(2, 0.5)
#   ), data=data.rate, iter=4000, chains=1, warmup=1000)
# 
# plot(guessrate_bad)
# pairs(guessrate_bad)
# precis(guessrate_bad)
# 
# plot(coeftab(guessrate_bad))
# 




