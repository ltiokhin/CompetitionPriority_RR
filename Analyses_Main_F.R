
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
data.c <- data[complete.cases(data),] #remove 5631 by this; by ignoring coliumn 14, remove 1731
data_math.c <- data_math[complete.cases(data_math),] #also make this not do the last column

#alternatively
# data <- data[complete.cases(data_math[, 1:13]), ]
# data_math.c <- data_math[complete.cases(data_math[, c(1:14, 16)]), ]

#temporary#####participant who told me that there was a bug#####
data.c <- data.c[data.c$ID_Player != 26,]
data_math.c <- data_math.c[data_math.c$ID_Player != 26,]

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
data.cf <- data.c
data_math.cf <- data_math.c

#changing sex and player ID to numeric
data.c$Sex <- as.numeric(data.c$Sex)
data_math.c$Sex <- as.numeric(data_math.c$Sex)
data.c$ID_Player <- as.numeric(data.c$ID_Player)
data_math.c$ID_Player <- as.numeric(data_math.c$ID_Player)

###changing variables to factors for frequentist analyses##
factorfxn <- function(x) as.factor(x)
cols <- c("Sex", "Effort", "Competition", "CompCheck", "Correct_Guess", "Faster", "ID_Player", "ID_Group")
data.cf[,cols] <- lapply(data.cf[,cols], factorfxn)
cols.m <- c("Sex", "Effort", "Competition", "CompCheck", "Correct_Guess", "Correct_Math", "ID_Player", 
            "ID_Group")
data_math.cf[,cols.m] <- lapply(data_math.cf[,cols.m], factorfxn)

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
  ), data=data.c, iter=3000, chains=1, cores=1, warmup=500)

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

m.tiles <- map2stan(
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
  ), data=data.c, iter=2500, chains=3, cores=3, warmup=500)

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
  geom_errorbar(aes(ymin = low_ci, ymax=high_ci), width=0.1, lwd=0.7) +
  ylim(0, 25) + xlab("Competition") + ylab("Tiles Revealed") +
scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes")) +
  scale_fill_brewer(name="Effort", palette = "Set1", 
                    labels = c("No Effort", "Effort"))


######
###Effect of competition and effort on accuracy##
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
  ), data=data.c, iter=6000, chains=2, cores = 3, warmup=1000)

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

plot( 0 , 0 , type="n" , xlab="comp/effort condition" ,
      ylab="probability correct" , ylim=c(0,1) , xaxt="n" ,
      xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("0/0","0/1","1/0","1/1") )
lines( 1:4 , pred.p.mean )
shade( pred.p.PI , 1:4 )

##plot with ggplot
d.gg <- data.frame(mean = pred.p.mean, low_ci = pred.p.PI[1,], high_ci = pred.p.PI[2,], 
                   competition = as.factor(c(0, 0, 1, 1)), effort = as.factor(c(0, 1, 0, 1)))
pd <- position_dodge(0.1) #move each point .05 to right

n <- ggplot(data=d.gg, aes(x=competition, y=mean, color=effort, group=effort)) +
  geom_errorbar(aes(ymin=low_ci, ymax=high_ci), width=0.1, lwd=0.8, position=pd) + 
  geom_line(position=pd, lwd=0.8) +
  theme_classic(base_size=12) +
  ylim(0, 1) 

n <- n + scale_colour_brewer(name="effort",
                             labels = c("No Effort", "Effort"),
                             palette = "Dark2") 

n <- n + scale_x_discrete(labels=c("0" = "No Competition",
                                   "1" = "Competition"))
n

##violin plot with ggplot
ggplot(data.c, aes(x = as.factor(Competition), y = TilesRevealed, fill = as.factor(Effort))) + 
  geom_violin(trim = FALSE) + 
  facet_grid(. ~ Effort) + guides(fill=FALSE) +
  theme_bw(base_size = 12) +
  geom_point(stat = "summary", fun.y = "mean", size = 2, color = "black") +
   xlab("Competition") + 
  scale_fill_brewer(name="Effort", palette = "Dark2", guides(fill=FALSE)) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
               width=0.1, lwd=0.5)

##aggregated for each player
d.c.agg <- aggregate(TilesRevealed ~ Competition + Effort + ID_Player, data=data.c, 
                     FUN = mean)
d.c.agg$Effort <- as.factor(d.c.agg$Effort)
levels(d.c.agg$Effort) <- c("No Effort", "Effort")
##violin plot with ggplot
ggplot(d.c.agg, aes(x = as.factor(Competition), y = TilesRevealed, fill = Effort)) + 
  geom_violin(trim = FALSE) +
  # #optional
  # geom_dotplot(binaxis="y", stackdir='center', 
  #                                          stackratio = 1, dotsize = 0.5) +
  facet_grid(. ~ Effort) + guides(fill=FALSE) +
  theme_bw(base_size = 12) + ylim(0, 25) +
  geom_point(stat = "summary", fun.y = "mean", size = 2, color = "black") +
  xlab("Competition") + ylab("Tiles Revealed") +
  scale_fill_brewer(name="Effort", palette = "Set1", guides(fill=FALSE)) +
  scale_x_discrete(name ="Competition", labels=c("No","Yes","No", "Yes")) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
               width=0.1, lwd=0.5)

#standardize guess number and run the same analyses - get the same results
data.cf$Guess_Number.s <- (data.cf$Guess_Number - mean(data.cf$Guess_Number)) / sd(data.cf$Guess_Number)

#Frequentist
m.3.f <- glmer(Correct_Guess ~ Competition*Effort + n_major.s + Guess_Number.s*Effort + (1|ID_Player), 
               data=data.cf, family=binomial)
sjt.glmer(m.3.f,
          show.icc=FALSE, show.col.header = TRUE,
          string.est = "Estimate",
          string.ci = "CI", 
          string.p = "P", 
          separate.ci.col = FALSE, 
          group.pred = TRUE, 
          exp.coef = FALSE)

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
