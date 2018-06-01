
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

###

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
  ), data=df_comp, iter=10000, chains=3, cores = 3, warmup=1000)

plot(q2)
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
  ), data=data.c, iter=10000, chains=3, cores=4, warmup=1000)

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
  ), data=data.c, iter=10000, chains=3, cores=4, warmup=1000)

plot(m.2.b)
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
  ), data=data.rate, iter=10000, chains=3, cores=4, warmup=1000)

plot(m.4.b)
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
  ), data=data.mathrate, iter=10000, chains=3, cores=4, warmup=1000)

plot(m.5.b)
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
save.image(file = "Pilot_Analyses_RR_without_pwr.RData")


##############
#Supplemental Exploratory
##############
# 
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
