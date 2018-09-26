
rm(list=ls())

#########
library(tidyverse)
library(rethinking)
library(lme4)
library(lmerTest)
library(effects)
library(sjPlot)
library(stringr)

###Load Data###

#load("PilotData1.RData")
#load("PilotData2_MathProbs.RData")
sequences <- read.table("Sequences1.txt", header=FALSE) #the sequence of blue and yellow tiles for each grid

###Cleaning and Data Prep###
data <- Data
data_math <- Data_MathProblem

#removing redundant data
rm(Data)
rm(Data_MathProblem)

#removing rows with missing observations (also removes participants who did not complete the study
#(i.e. had NA values for compcheck)
data.c <- data[complete.cases(data),] #remove 5631 by this (ignore colummn 14 if want to keep 0 tile guesses)
data_math.c <- data_math[complete.cases(data_math),] #also make this not do the last column

#alternatively
# data <- data[complete.cases(data_math[, 1:13]), ]
# data_math.c <- data_math[complete.cases(data_math[, c(1:14, 16)]), ]

##################
###check numbers in each condition
data_f <- data.c[data.c$Sex == "Female",]
data_m <- data.c[data.c$Sex == "Male",]

#females
length(unique(data_f$ID_Player[data_f$Effort == 0 & data_f$Competition == 0])) # 10
length(unique(data_f$ID_Player[data_f$Effort == 0 & data_f$Competition == 1])) # 10
length(unique(data_f$ID_Player[data_f$Effort == 1 & data_f$Competition == 0])) # 10
length(unique(data_f$ID_Player[data_f$Effort == 1 & data_f$Competition == 1])) # 9

#males
length(unique(data_m$ID_Player[data_m$Effort == 0 & data_m$Competition == 0])) # 12
length(unique(data_m$ID_Player[data_m$Effort == 0 & data_m$Competition == 1])) # 4
length(unique(data_m$ID_Player[data_m$Effort == 1 & data_m$Competition == 0])) # 9
length(unique(data_m$ID_Player[data_m$Effort == 1 & data_m$Competition == 1])) # 3

rm(data_f)
rm(data_m)
########################

#simple histogram of time until guess
simplehist(data.c$ElapsedTime_Guess, xlab = "Time Until Guess (seconds), All Data.")

#remove times until guess that are more than 5 standard deviations away from the mean
sd5.times <- mean(data.c$ElapsedTime_Guess) + ( 5 * sd(data.c$ElapsedTime_Guess))
nrow(data.c[data.c$ElapsedTime_Guess > sd5.times,]) #number of observations excluded
data.c <- data.c[data.c$ElapsedTime_Guess < sd5.times,]
simplehist(data.c$ElapsedTime_Guess, xlab = "Time Until Guess (seconds), Outliers Exluded (5 SD)")

#simple histogram of arithmetic problem solving times
simplehist(data_math.c$ElapsedTime_Math, xlab = "Arithmetic Problem Solving Times (seconds), All Data.")

#remove math-problem times that are more than 3 standard deviations away from the mean
sd5 <- mean(data_math.c$ElapsedTime_Math) + ( 5 * sd(data_math.c$ElapsedTime_Math))
nrow(data_math.c[data_math.c$ElapsedTime_Math > sd5,]) #number of observations excluded
data_math.c <- data_math.c[data_math.c$ElapsedTime_Math < sd5,]
simplehist(data_math.c$ElapsedTime_Math, xlab = "Arithmetic Problem Solving Times (seconds), Outliers Exluded (5 SD)")

#copying data into separate objects for frequentist analyses
# data.cf <- data.c
# data_math.cf <- data_math.c

#changing sex and player ID to numeric
data.c$Sex <- as.numeric(data.c$Sex)
data_math.c$Sex <- as.numeric(data_math.c$Sex)
data.c$ID_Player <- as.numeric(data.c$ID_Player)
data_math.c$ID_Player <- as.numeric(data_math.c$ID_Player)

###changing variables to factors for frequentist analyses##
# factorfxn <- function(x) as.factor(x)
# cols <- c("Sex", "Effort", "Competition", "CompCheck", "Correct_Guess", "Faster", "ID_Player", "ID_Group")
# data.cf[,cols] <- lapply(data.cf[,cols], factorfxn)
# cols.m <- c("Sex", "Effort", "Competition", "CompCheck", "Correct_Guess", "Correct_Math", "ID_Player", 
#             "ID_Group")
# data_math.cf[,cols.m] <- lapply(data_math.cf[,cols.m], factorfxn)

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

#adds this to all data sets
data.c$n_major <- rep(NA, nrow(data.c))
#data.cf$n_major <- rep(NA, nrow(data.cf))
data_math.c$n_major <- rep(NA, nrow(data_math.c))
#data_math.cf$n_major <- rep(NA, nrow(data_math.cf))

for(i in 1:nrow(data.c)){
  if(data.c$Guess_Number[i] %in% g_13){data.c$n_major[i] <- 13}
  if(data.c$Guess_Number[i] %in% g_15){data.c$n_major[i] <- 15}
  if(data.c$Guess_Number[i] %in% g_17){data.c$n_major[i] <- 17}
}
# for(i in 1:nrow(data.cf)){
#   if(data.cf$Guess_Number[i] %in% g_13){data.cf$n_major[i] <- 13}
#   if(data.cf$Guess_Number[i] %in% g_15){data.cf$n_major[i] <- 15}
#   if(data.cf$Guess_Number[i] %in% g_17){data.cf$n_major[i] <- 17}
# }
for(i in 1:nrow(data_math.c)){
  if(data_math.c$Guess_Number[i] %in% g_13){data_math.c$n_major[i] <- 13}
  if(data_math.c$Guess_Number[i] %in% g_15){data_math.c$n_major[i] <- 15}
  if(data_math.c$Guess_Number[i] %in% g_17){data_math.c$n_major[i] <- 17}
}
# for(i in 1:nrow(data_math.cf)){
#   if(data_math.cf$Guess_Number[i] %in% g_13){data_math.cf$n_major[i] <- 13}
#   if(data_math.cf$Guess_Number[i] %in% g_15){data_math.cf$n_major[i] <- 15}
#   if(data_math.cf$Guess_Number[i] %in% g_17){data_math.cf$n_major[i] <- 17}
# }

#standardizing n_major
data.c$n_major.s <- (data.c$n_major - mean(data.c$n_major)) / sd(data.c$n_major)
#data.cf$n_major.s <- (data.cf$n_major - mean(data.cf$n_major)) / sd(data.cf$n_major)
data_math.c$n_major.s <- (data_math.c$n_major - mean(data_math.c$n_major)) / sd(data_math.c$n_major)
#data_math.cf$n_major.s <- (data_math.cf$n_major - mean(data_math.cf$n_major)) / sd(data_math.cf$n_major)

# #change to factor for frequentist analyses
# data.cf$n_major <- as.factor(data.cf$n_major)
# data_math.cf$n_major.s <- as.factor(data_math.cf$n_major)
# data.cf$n_major.s <- as.factor(data.cf$n_major.s)
# data_math.cf$n_major.s <- as.factor(data_math.cf$n_major.s)

###########################
###Quality Checks###
##########################

######
###Effect of Effort Treatment on Time to Click 1 Tile
######

qq1 <- map2stan(
  alist(
    ElapsedTime_Tile ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bE*Effort + bCE*Competition*Effort,
    bC ~ dnorm(0, 10), 
    bE ~ dnorm(0, 10), 
    bCE ~ dnorm(0, 10),
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1.5, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1.5, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=data.c, iter=3000, chains=3, cores=4, warmup=500)

plot(qq1)
par(mfrow=c(1,1))
precis(qq1, prob=0.95)
plot(precis(qq1, prob=0.95), xlab = "Time (seconds) to Click 1 Tile")

######
###Competition Attention Check###
######

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

###########################
###Confirmatory Analyses###
##########################

######
###Effect of competition and effort on number of tiles revealed###
######

d.tiles <- data.c
d.tiles <- aggregate(TilesRevealed ~ ID_Player + Guess_Number + Effort + Competition + n_major.s, 
                     data=d.tiles, FUN = mean)

m.tiles <- map2stan(
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
  ), data=d.tiles, iter=3500, chains=1, cores=1, warmup=500)

plot(m.tiles)
par(mfrow=c(1,1))
precis(m.tiles, prob = 0.95)
plot(precis(m.tiles, prob = 0.95), xlab = "Tiles Revealed")

###Plot: Tiles Revealed###

#Plot Predictions
d.pred <- list(
  Competition = c(0, 0, 1, 1), 
  Effort = c(0, 1, 0, 1), 
  n_major_s = c(0, 0, 0, 0),
  ID_Player = rep(2, 4) #placeholder
)

#replace varying intercept samples with zeros
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(data.c$ID_Player)))

m.tiles.link <- link(m.tiles, n=2000, data=d.pred, 
                 replace = list(a_player=a_player_zeros))

#summarize#
pred.p.mean <- apply(m.tiles.link , 2 , mean)
pred.p.PI <- apply( m.tiles.link , 2 , HPDI , prob=0.95)

### plot the raw data and the 95% HPDI from m.tiles ###
d.gg.f <- data.c
d.gg.f$mean <- NA
d.gg.f$low_ci <- NA
d.gg.f$high_ci <- NA

#means
d.gg.f$mean[d.gg.f$Competition == 0 & d.gg.f$Effort == 0] <- d.gg$mean[d.gg$competition == 0 & d.gg$effort == 0]
d.gg.f$mean[d.gg.f$Competition == 0 & d.gg.f$Effort == 1] <- d.gg$mean[d.gg$competition == 0 & d.gg$effort == 1]
d.gg.f$mean[d.gg.f$Competition == 1 & d.gg.f$Effort == 0] <- d.gg$mean[d.gg$competition == 1 & d.gg$effort == 0]
d.gg.f$mean[d.gg.f$Competition == 1 & d.gg.f$Effort == 1] <- d.gg$mean[d.gg$competition == 1 & d.gg$effort == 1]

#CI's
d.gg.f$low_ci[d.gg.f$Competition == 0 & d.gg.f$Effort == 0] <- d.gg$low_ci[d.gg$competition == 0 & d.gg$effort == 0]
d.gg.f$low_ci[d.gg.f$Competition == 0 & d.gg.f$Effort == 1] <- d.gg$low_ci[d.gg$competition == 0 & d.gg$effort == 1]
d.gg.f$low_ci[d.gg.f$Competition == 1 & d.gg.f$Effort == 0] <- d.gg$low_ci[d.gg$competition == 1 & d.gg$effort == 0]
d.gg.f$low_ci[d.gg.f$Competition == 1 & d.gg.f$Effort == 1] <- d.gg$low_ci[d.gg$competition == 1 & d.gg$effort == 1]

d.gg.f$high_ci[d.gg.f$Competition == 0 & d.gg.f$Effort == 0] <- d.gg$high_ci[d.gg$competition == 0 & d.gg$effort == 0]
d.gg.f$high_ci[d.gg.f$Competition == 0 & d.gg.f$Effort == 1] <- d.gg$high_ci[d.gg$competition == 0 & d.gg$effort == 1]
d.gg.f$high_ci[d.gg.f$Competition == 1 & d.gg.f$Effort == 0] <- d.gg$high_ci[d.gg$competition == 1 & d.gg$effort == 0]
d.gg.f$high_ci[d.gg.f$Competition == 1 & d.gg.f$Effort == 1] <- d.gg$high_ci[d.gg$competition == 1 & d.gg$effort == 1]

#plot
ggplot(d.gg.f, aes(x = as.factor(Competition), y = TilesRevealed, 
                   fill = as.factor(Effort))) + 
  geom_violin(trim = FALSE) + 
  facet_grid(. ~ Effort) + theme_bw(base_size = 12) +
  theme(strip.text.x = element_blank()) + 
  geom_pointrange(aes(ymin = low_ci, ymax=high_ci), lwd=0.8) +
  ylim(0, 25) + xlab("Competition") + ylab("Tiles Revealed") +
  scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes")) +
  scale_fill_brewer(name="Effort", palette = "Set1", 
                    labels = c("No Effort", "Effort"))

######
###Effect of competition and effort on accuracy###
######

d.acc <- data.c
d.acc <- aggregate(Correct_Guess ~ ID_Player + Guess_Number + Effort + Competition + n_major.s, 
                     data=d.acc, FUN = mean)

m.accuracy <- map2stan(
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
  ), data=d.acc, iter=4000, chains = 2, cores = 2, warmup=500)

plot(m.accuracy)
par(mfrow=c(1,1))
precis(m.accuracy, prob=0.95)
plot(precis(m.accuracy, prob=0.95), xlab = "Log Odds of Correct Guess")

#Plot Predictions
d.pred.a <- list(
  Competition = c(0, 0, 1, 1), 
  Effort = c(0, 1, 0, 1), 
  n_major_s = c(0, 0, 0, 0),
  ID_Player = rep(2, 4) #placeholder
)

#replace varying intercept samples with zeros
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(data.c$ID_Player)))

m.accuracy.link <- link(m.accuracy, n=2000, data=d.pred.a, 
                        replace = list(a_player=a_player_zeros))

#summarize and plot with plot function
pred.p.mean.accuracy <- apply( m.accuracy.link , 2 , mean )
pred.p.PI.accuracy <- apply( m.accuracy.link , 2 , HPDI , prob=0.95)

##plot with ggplot
d.gg <- data.frame(mean = pred.p.mean.accuracy, low_ci = pred.p.PI.accuracy[1,], high_ci = pred.p.PI.accuracy[2,], 
                   competition = as.factor(c(0, 0, 1, 1)), effort = as.factor(c(0, 1, 0, 1)))

n <- ggplot(data=d.gg, aes(x=competition, y=mean, color=effort, group=effort)) +
  geom_pointrange(aes(ymin=low_ci, ymax=high_ci), lwd=0.8) +
  facet_grid(. ~ effort) + theme_bw(base_size=12) +
  theme(strip.text.x = element_blank()) + 
  geom_line(lwd=0.8) +
  ylim(0.5, 1) + xlab("Competition") + ylab("Accuracy")

n <- n + scale_colour_brewer(name="Effort",
                             labels = c("No Effort", "Effort"),
                             palette = "Set1") 

n <- n + scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes"))

n

##########
###Effect of competition on effort (i.e. time to accurately solve an arithmetic problem)###
#########

m.effort <- map2stan(
  alist(
    ElapsedTime_MathSolved ~ dnorm(mu, sigma), 
    mu <- a + a_player[ID_Player] + bC*Competition + bNs*n_major.s, 
    bC ~ dnorm(0, 10), 
    bNs ~ dnorm(0, 10),
    a ~ dgamma(1, 0.05), 
    a_player[ID_Player] ~ dnorm(0, sigma_player),
    sigma_player ~ dgamma(1, 0.05),
    sigma ~ dgamma(2, 0.5)
  ), data=data_math.c, iter=4000, chains = 2, cores = 2, warmup=500)

plot(m.effort)
par(mfrow=c(1,1))
precis(m.effort, prob=0.95)
plot(precis(m.effort, prob=0.95), xlab = "Time (seconds) to Accurately Solve One Arithmetic Problem")

#Plot Predictions
d.pred.math <- list(
  Competition = c(0, 1), 
  n_major_s = c(0, 0),
  ID_Player = rep(2, 2) #placeholder
)

#replace varying intercept samples with zeros
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(data_math.c$ID_Player)))

m.math.link <- link(m.effort, n=2000, data=d.pred.math, 
                    replace = list(a_player=a_player_zeros))

#summarize#
pred.p.mean.math <- apply(m.math.link , 2 , mean)
pred.p.PI.math <- apply( m.math.link , 2 , HPDI , prob=0.95)

##plot with ggplot
d.gg <- data.frame(mean = pred.p.mean.math, low_ci = pred.p.PI.math[1,], high_ci = pred.p.PI.math[2,], 
                   competition = as.factor(c(0, 1)))

### plot the raw data and the 95% HPDI from m.arithmetic ###
d.gg.math <- data_math.c
d.gg.math$mean <- NA
d.gg.math$low_ci <- NA
d.gg.math$high_ci <- NA

#means
d.gg.math$mean[d.gg.math$Competition == 0] <- d.gg$mean[d.gg$competition == 0]
d.gg.math$mean[d.gg.math$Competition == 1] <- d.gg$mean[d.gg$competition == 1]

#CI's
d.gg.math$low_ci[d.gg.math$Competition == 0] <- d.gg$low_ci[d.gg$competition == 0]
d.gg.math$low_ci[d.gg.math$Competition == 1] <- d.gg$low_ci[d.gg$competition == 1]
d.gg.math$high_ci[d.gg.math$Competition == 0] <- d.gg$high_ci[d.gg$competition == 0]
d.gg.math$high_ci[d.gg.math$Competition == 1] <- d.gg$high_ci[d.gg$competition == 1]

#plot
d.gg$ElapsedTime_MathSolved <- d.gg$mean
d.gg$Competition <- d.gg$competition

ggplot(d.gg.math, aes(x = as.factor(Competition), y = ElapsedTime_MathSolved)) +
  geom_violin(trim = FALSE, fill = "#377eb8") + theme_bw(base_size = 12) +
  stat_summary(fun.y=mean, data = d.gg, geom="point", size=2, color = "black") +
  geom_pointrange(data = d.gg, aes(ymin = low_ci, ymax=high_ci), lwd=0.8) +
  ylim(0, 25) + 
  xlab("Competition") + ylab("Time (seconds) to Accurately Solve an Arithmetic Problem") +
  scale_x_discrete(name ="Competition", labels=c("No","Yes")) +
  scale_fill_brewer(palette="Set1")

#################
###Exploratory###
################

####frequentist
data_math.cF <- data_math.c
data_math.cF$ID_Player <- as.factor(data_math.cF$ID_Player)
data_math.cF$Competition <- as.factor(data_math.cF$Competition)
data_math.cF$Sex <- as.factor(data_math.cF$Sex)

m.math.freq <- lmer(ElapsedTime_MathSolved ~ Competition*Sex + n_major.s + (1|ID_Player), data=data_math.cF)
summary(m.math.freq)
plot(allEffects(m.math.freq))

#save workspace
save.image(file = "Confirmatory_plus_Exploratory_Analyses.RData")


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



############

agg <- aggregate(TilesRevealed ~ Guess_Number + Competition + Effort, data=data.c, FUN = mean)

plot(TilesRevealed ~ Guess_Number, data = agg[agg$Competition==0 & agg$Effort==0,], main = "no comp no effort")
plot(TilesRevealed ~ Guess_Number, data = agg[agg$Competition==0 & agg$Effort==1,], main = "no comp,yes effort")

plot(TilesRevealed ~ Guess_Number, data = agg[agg$Competition==1 & agg$Effort==0,], main = "yes comp no effort")
plot(TilesRevealed ~ Guess_Number, data = agg[agg$Competition==1 & agg$Effort==1,], main = "yes comp,yes effort")

cor(agg$Guess_Number[agg$Competition==1], agg$TilesRevealed[agg$Competition==1])
cor(agg$Guess_Number[agg$Competition==0], agg$TilesRevealed[agg$Competition==0])

#ggplot of just the points
ggplot(data=agg, aes(x=Guess_Number, y=TilesRevealed, color=as.factor(Competition), 
                     group=as.factor(Competition))) +
  geom_point(colour=rangi2, alpha = 0.4)+ 
  theme_classic(base_size=12) +
  ylim(0, 15) + 
  scale_colour_brewer(name="Competition",
                      labels = c("No Competition", "Competition"),
                      palette = "Dark2") 
# 
