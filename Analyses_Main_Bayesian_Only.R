
rm(list=ls())

#########
library(tidyverse)
library(rethinking)
library(effects)
library(sjPlot)
library(stringr)
library(RColorBrewer)

###Load Data###

#load("PilotData1.RData")
#load("PilotData2_MathProbs.RData")
sequences <- read.table("Sequences1.txt", header=FALSE) #the sequence of blue and yellow tiles for each grid

###Cleaning and Data Prep###
d.quality.c <- d.quality[complete.cases(d.quality[c(1:8, 10)]),]
d.conf_1_2.c <- d.conf_1_2[complete.cases(d.conf_1_2),]
d.conf_3_math.c <- d.conf_3_math[complete.cases(d.conf_3_math),]

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
d.conf_1_2.c$n_major <- rep(NA, nrow(d.conf_1_2.c))
d.conf_3_math.c$n_major <- rep(NA, nrow(d.conf_3_math.c))

for(i in 1:nrow(d.conf_1_2.c)){
  if(d.conf_1_2.c$Guess_Number[i] %in% g_13){d.conf_1_2.c$n_major[i] <- 13}
  if(d.conf_1_2.c$Guess_Number[i] %in% g_15){d.conf_1_2.c$n_major[i] <- 15}
  if(d.conf_1_2.c$Guess_Number[i] %in% g_17){d.conf_1_2.c$n_major[i] <- 17}
}

for(i in 1:nrow(d.conf_3_math.c)){
  if(d.conf_3_math.c$Guess_Number[i] %in% g_13){d.conf_3_math.c$n_major[i] <- 13}
  if(d.conf_3_math.c$Guess_Number[i] %in% g_15){d.conf_3_math.c$n_major[i] <- 15}
  if(d.conf_3_math.c$Guess_Number[i] %in% g_17){d.conf_3_math.c$n_major[i] <- 17}
}

#remove rows with time-until-guess values that are more than 5 standard deviations away from the mean
agg <- aggregate(ElapsedTime_Guess ~ Guess_Number + ID_Player + Effort + Competition, data=d.quality.c, FUN = mean)
sd5.times <- mean(agg$ElapsedTime_Guess) + (5 * sd(agg$ElapsedTime_Guess))
nrow(agg[agg$ElapsedTime_Guess > sd5.times,]) #104 observations excluded

#excluding observations from all datasets
d.quality.c <- d.quality.c[d.quality.c$ElapsedTime_Guess < sd5.times,]
d.conf_1_2.c <- d.conf_1_2.c[d.conf_1_2.c$ElapsedTime_Guess < sd5.times,]
d.conf_3_math.c <- d.conf_3_math.c[d.conf_3_math.c$ElapsedTime_Guess < sd5.times,]

#remove rows with arithmetic-solving-times that are more than 5 standard deviations away from the mean
sd5.times <- mean(d.conf_3_math.c$ElapsedTime_Math) + (5 * sd(d.conf_3_math.c$ElapsedTime_Math))
nrow(d.conf_3_math.c[d.conf_3_math.c$ElapsedTime_Math > sd5.times,]) #93 observations excluded

#excluding observations from math dataset
d.conf_3_math.c <- d.conf_3_math.c[d.conf_3_math.c$ElapsedTime_Math < sd5.times,]

###########################
###Quality Checks###
##########################

#Remove remaining NA's
d.quality.final <- d.quality.c[complete.cases(d.quality.c),]

######
###Effect of Effort Treatment on Time to Click 1 Tile
######
d.quality.final$ID_Player <- coerce_index(d.quality.final$ID_Player)

m.quality.2 <- map2stan(
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
  ), data=d.quality.final, iter=30000, chains=5, cores=5, warmup=500)

plot(m.quality.2)
par(mfrow=c(1,1))
precis(m.quality.2, prob=0.95)
plot(precis(m.quality.2, prob=0.95), xlab = "Time (seconds) to Click 1 Tile")

######
###Competition Attention Check###
######

df_comp <- aggregate(cbind(CompCheck, Competition) ~ ID_Player,
                     data=d.quality.final, FUN=unique)
q2 <- map2stan(
  alist(
    CompCheck ~ dbinom(1, theta), 
    logit(theta) <- a + bC*Competition, 
    bC ~ dnorm(0, 10), 
    a ~ dnorm(0, 10)
  ), data=df_comp, iter=10000, chains=2, cores = 2, warmup=500)

plot(q2)
par(mfrow=c(1,1))
precis(q2, prob=0.95)
plot(precis(q2, prob=0.95), xlab = "Log Odds of Answering 'Yes'")
logistic(5)

###########################
###Confirmatory Analyses###
##########################

#Standardize number of majority tiles
d.conf_1_2.c$n_major.s <- (d.conf_1_2.c$n_major - mean(d.conf_1_2.c$n_major)) / sd(d.conf_1_2.c$n_major)
d.conf_3_math.c$n_major.s <- (d.conf_3_math.c$n_major - mean(d.conf_3_math.c$n_major)) / sd(d.conf_3_math.c$n_major)

#Aggregate data for each participant for confirmatory analyses 1 and 2
d.conf.agg <- aggregate(cbind(TilesRevealed, Correct_Guess) ~ Effort + Competition + Faster +
                       Sex + Guess_Number + n_major.s + ID_Player, data=d.conf_1_2.c, FUN = mean)

d.conf.agg$ID_Player <- coerce_index(d.conf.agg$ID_Player)
d.conf.agg$Sex <- as.integer(d.conf.agg$Sex)

######
###Effect of competition and effort on number of tiles revealed###
######

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
  ), data=d.conf.agg, iter=8000, chains=2, cores=2, warmup=500)

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
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(d.conf.agg$ID_Player)))

m.tiles.link <- link(m.tiles, n=2000, data=d.pred, 
                 replace = list(a_player=a_player_zeros))

#summarize#
pred.p.mean <- apply(m.tiles.link$mu , 2 , mean)
pred.p.PI <- apply( m.tiles.link$mu , 2 , HPDI , prob=0.95)

### plot the raw data and the 95% HPDI from m.tiles ###
d.gg.f <- d.conf.agg
d.gg.f$mean <- NA
d.gg.f$low_ci <- NA
d.gg.f$high_ci <- NA

d.gg <- data.frame(mean = pred.p.mean, low_ci = pred.p.PI[1,], high_ci = pred.p.PI[2,], 
                   competition = as.factor(c(0, 0, 1, 1)), effort = as.factor(c(0, 1, 0, 1)))

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
d.gg2 <- d.gg
d.gg2$Competition <- d.gg2$competition
d.gg2$Effort <- d.gg2$effort
d.gg2$TilesRevealed <- d.gg2$mean
d.gg2$competition <- NULL
d.gg2$effort <- NULL
d.gg.f$Effort <- as.factor(d.gg.f$Effort)

ggplot(d.gg.f, aes(x = as.factor(Competition), y = TilesRevealed, 
                   fill = Effort)) + 
  geom_violin(trim = FALSE) + 
  facet_grid(. ~ Effort) + theme_bw(base_size = 12) + theme(strip.text.x = element_blank()) +
  ylim(0, 25) + xlab("Competition") + ylab("Tiles Revealed") +
  scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes")) +
  geom_pointrange(aes(ymin = low_ci, ymax=high_ci), data=d.gg2, lwd=0.8) +
  scale_fill_brewer(name="Effort", palette = "Set1", 
                    labels = c("No Effort", "Effort"))

###Alternative plots for tiles revealed###
pred.p.mean_sorted <- c(pred.p.mean[1], pred.p.mean[3], pred.p.mean[2], 
                        pred.p.mean[4])
pred.p.PI_sorted <- matrix(NA, 2, 4)
pred.p.PI_sorted[,1] <- pred.p.PI[,1]
pred.p.PI_sorted[,2] <- pred.p.PI[,3]
pred.p.PI_sorted[,3] <- pred.p.PI[,2]
pred.p.PI_sorted[,4] <- pred.p.PI[,4]

plot( 0 , 0 , type="n" , xlab="Competition/Effort Treatment" ,
      ylab="Tiles Revealed" , ylim=c(0,25) , xaxt="n" , main = "Mean Tiles Revealed in 4 Treatments",
      xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("0/0","1/0","0/1","1/1") )
lines( 1:4 , pred.p.mean_sorted, type = "b")
shade( pred.p.PI_sorted , 1:4 , col=col.alpha(rangi2, 0.15))

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
