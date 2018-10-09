
rm(list=ls())

#########
library(tidyverse)
library(knitr)
library(kableExtra)
library(rethinking)
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

#exclude participants who indicated technical difficulties
d.quality.c <- d.quality.c[d.quality.c$ID_Player != 26,]
d.conf_1_2.c <- d.conf_1_2.c[d.conf_1_2.c$ID_Player != 26,]
d.conf_3_math.c <- d.conf_3_math.c[d.conf_3_math.c$ID_Player != 26,]

d.quality.c <- d.quality.c[d.quality.c$ID_Player != 61,]
d.conf_1_2.c <- d.conf_1_2.c[d.conf_1_2.c$ID_Player != 61,]
d.conf_3_math.c <- d.conf_3_math.c[d.conf_3_math.c$ID_Player != 61,]

#Number of participants in each treatment, before excluding individual observations
#within participants
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Female" & d.conf_1_2.c$Effort==1])) 
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Female" & d.conf_1_2.c$Effort==0])) 
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Male" & d.conf_1_2.c$Effort==1])) 
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Male" & d.conf_1_2.c$Effort==0])) 

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
nrow(agg[agg$ElapsedTime_Guess > sd5.times,]) #101 observations excluded

#excluding observations from all datasets
d.quality.c <- d.quality.c[d.quality.c$ElapsedTime_Guess < sd5.times,]
d.conf_1_2.c <- d.conf_1_2.c[d.conf_1_2.c$ElapsedTime_Guess < sd5.times,]
d.conf_3_math.c <- d.conf_3_math.c[d.conf_3_math.c$ElapsedTime_Guess < sd5.times,]

#remove rows with arithmetic-solving-times that are more than 5 standard deviations away from the mean
sd5.times <- mean(d.conf_3_math.c$ElapsedTime_Math) + (5 * sd(d.conf_3_math.c$ElapsedTime_Math))
nrow(d.conf_3_math.c[d.conf_3_math.c$ElapsedTime_Math > sd5.times,]) #93 observations excluded

#excluding observations from math dataset
d.conf_3_math.c <- d.conf_3_math.c[d.conf_3_math.c$ElapsedTime_Math < sd5.times,]

#Number of participants in each treatment, after excluding individual observations
#within participants
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Female" & d.conf_1_2.c$Effort==1])) 
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Female" & d.conf_1_2.c$Effort==0])) 
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Male" & d.conf_1_2.c$Effort==1])) 
length(unique(d.conf_1_2.c$ID_Player[d.conf_1_2.c$Sex=="Male" & d.conf_1_2.c$Effort==0])) 



###########################
###Quality Checks###
##########################

#Remove remaining NA's
d.quality.effortcheck <- d.quality.c[complete.cases(d.quality.c),]

######
###Effect of Effort Treatment on Time to Click 1 Tile
######

d.quality.effortcheck$ID_Player <- coerce_index(d.quality.effortcheck$ID_Player)

m.quality1 <- map2stan(
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
  ), data=d.quality.effortcheck, iter=30000, chains=5, cores=5, warmup=500)

plot(m.quality1)
par(mfrow=c(1,1))
precis(m.quality1, prob=0.95)
plot(precis(m.quality1, prob=0.95), xlab = "Time (seconds) to Click 1 Tile")

### Means and CI in each treatment ###

d.pred <- list(
  Competition = c(0, 0, 1, 1), 
  Effort = c(0, 1, 0, 1), 
  ID_Player = rep(2, 4) #placeholder
)

#replace varying intercept samples with zeros
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(d.quality.effortcheck$ID_Player)))

m.quality1.link <- link(m.quality1, n=2000, data=d.pred, 
                     replace = list(a_player=a_player_zeros))

#summarize#
pred.p.mean <- apply(m.quality1.link , 2 , mean)
pred.p.PI <- apply( m.quality1.link , 2 , HPDI , prob=0.95)



######
###Competition Attention Check###
######

df_comp <- aggregate(cbind(CompCheck, Competition) ~ ID_Player,
                     data=d.quality.c, FUN=unique)
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

#number of participants in each treatment
length(unique(d.conf.agg$ID_Player[d.conf.agg$Competition==0 & d.conf.agg$Effort==0])) 
length(unique(d.conf.agg$ID_Player[d.conf.agg$Competition==1 & d.conf.agg$Effort==0]))
length(unique(d.conf.agg$ID_Player[d.conf.agg$Competition==0 & d.conf.agg$Effort==1])) 
length(unique(d.conf.agg$ID_Player[d.conf.agg$Competition==1 & d.conf.agg$Effort==1]))

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

#compare Effort, No-Competition treatment to Effort, Competition treatment
mean_diff_Enoc_Ec <- mean(m.tiles.link$mu[,4] - m.tiles.link$mu[,2])
HPDI_Enoc_Ec <- HPDI(m.tiles.link$mu[,4] - m.tiles.link$mu[,2], prob = 0.95)

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
  geom_violin(trim = FALSE, bw = 0.8) + 
  facet_grid(. ~ Effort) + theme_bw(base_size = 14) + theme(strip.text.x = element_blank()) +
  ylim(0, 25) + xlab("Competition") + ylab("Tiles Revealed") +
  scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes")) +
  geom_pointrange(aes(ymin = low_ci, ymax=high_ci), data=d.gg2, lwd=0.8) +
  scale_fill_brewer(name="Effort", palette = "Set1", 
                    labels = c("No", "Yes"))

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
  ), data=d.conf.agg, iter=6000, chains = 3, cores = 3, warmup=500)

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
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(d.conf.agg$ID_Player)))

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
  facet_grid(. ~ effort) + theme_bw(base_size=14) +
  theme(strip.text.x = element_blank()) + 
  geom_line(lwd=0.8) +
  ylim(0.5, 1) + xlab("Competition") + ylab("Accuracy")

n <- n + scale_colour_brewer(name="Effort",
                             labels = c("No", "Yes"),
                             palette = "Set1") 

n <- n + scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes"))

n



##########
###Effect of competition on effort (i.e. time to accurately solve an arithmetic problem)###
#########

d.conf_3_math.c$ID_Player <- coerce_index(d.conf_3_math.c$ID_Player)
d.conf_3_math.c$Sex <- as.integer(d.conf_3_math.c$Sex)

d.math.agg <- aggregate(cbind(ElapsedTime_MathSolved) ~ Effort + Competition + Sex + Tile_Number +
                          Guess_Number + n_major.s + ID_Player, data=d.conf_3_math.c, FUN = mean)

### Model 3 ###
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
  ), data=d.math.agg, iter=20000, chains = 3, cores = 3, warmup=500)

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
a_player_zeros <- matrix(0, nrow=2000, ncol = length(unique(d.math.agg$ID_Player)))

m.math.link <- link(m.effort, n=2000, data=d.pred.math, 
                    replace = list(a_player=a_player_zeros))

#summarize#
pred.p.mean.math <- apply(m.math.link , 2 , mean)
pred.p.PI.math <- apply( m.math.link , 2 , HPDI , prob=0.95)

##plot with ggplot
d.gg <- data.frame(mean = pred.p.mean.math, low_ci = pred.p.PI.math[1,], high_ci = pred.p.PI.math[2,], 
                   competition = as.factor(c(0, 1)))

### plot the raw data and the 95% HPDI from m.arithmetic ###
d.gg.math <- d.math.agg
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
  geom_violin(trim = FALSE, fill = "#377eb8", bw = 0.5) + theme_bw(base_size = 14) +
  stat_summary(fun.y=mean, data = d.gg, geom="point", size=2, color = "black") +
  geom_pointrange(data = d.gg, aes(ymin = low_ci, ymax=high_ci), lwd=0.8) +
  ylim(0, 25) + 
  xlab("Competition") + ylab("Time (seconds) to Accurately Solve an Arithmetic Problem") +
  scale_x_discrete(name ="Competition", labels=c("No","Yes")) +
  scale_fill_brewer(palette="Set1")

##############
###Table with parameter estimates from all models
##############

precis_tiles <- precis(m.tiles, depth=1, prob = 0.95)
precis_tiles <- precis_tiles@output
precis_tiles <- round(precis_tiles, 2)
precis_tiles <- precis_tiles[1:5, c(1, 3, 4)]
precis_tiles <- rename(precis_tiles, `Lower 0.95` = `lower 0.95`)
precis_tiles <- rename(precis_tiles, `Upper 0.95` = `upper 0.95`)

row.names(precis_tiles) <- c("Competition", "Effort", "Competition x Effort Interaction", 
                             "Effect Size", "Intercept")

precis_tiles %>% kable %>% 
  kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), 
                full_width = TRUE) %>%
  column_spec(1, width = "5cm") %>%
  column_spec(2:4, width = "2cm") %>%
  group_rows("Model 1: Tiles Revealed", 1, 5) 

precis_accuracy <- precis(m.accuracy, depth=1, prob = 0.95)
precis_accuracy <- precis_accuracy@output
precis_accuracy <- round(precis_accuracy, 2)
precis_accuracy <- precis_accuracy[1:5, c(1, 3, 4)]
precis_accuracy <- rename(precis_accuracy, `Lower 0.95` = `lower 0.95`)
precis_accuracy <- rename(precis_accuracy, `Upper 0.95` = `upper 0.95`)

precis_tiles %>% kable %>% 
  kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), 
                full_width = TRUE) %>%
  column_spec(1, width = "5cm") %>%
  column_spec(2:4, width = "2cm") %>%
  group_rows("Model 2: Accuracy", 1, 5) 
    
  

