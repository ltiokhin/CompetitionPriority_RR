library(ggplot2)
library(reshape2)
library(viridisLite)
library(viridis)
library(RColorBrewer)

setwd() #set to location of data file
rm(list=ls())
load('ScriptInformation2_25tiles_500sims.RData')

N_Tiles <- ncol(InformationFinal)
plot(NULL, xlab = 'Number of tiles revealed', ylab ='Average Information', xlim =c(0,N_Tiles),  ylim = c(0.5,1))
Color <- c("blue","red", "black")
Information_Subset <- InformationFinal[c(1, 3, 5),]

for (i in 1:length(Color)){
  lines(Information_Subset[i,], type = "l", col = Color[i])
}
############
#Expected Payoffs for Individual Players
############
plot(NULL, xlab = 'Number of Tiles Revealed', ylab ='Reward', xlim =c(0,N_Tiles), ylim=c(0, 50))
Curve <- c(1, 2, 3) #Curve 1 is 8 yellow; Curve 2 is 10 yellow; Curve 3 is 12 yellow
Color <- c("#8856a7","#9ebcda", "#e0ecf4")
Reward_wait <- NULL
Optimal_Tiles_Revealed <- NULL
reward_list <- list()
TimeStep <- 60 * 20 #length of experiment
TilesWait <- 5 # 5 second waiting time between each problem

for(curv in 1:length(Curve)) {
  for(i in 1 : N_Tiles) {
    Reward_wait[i] <- 
      ((Information_Subset[curv,i] * (TimeStep/i)) - (1 * ((1-Information_Subset[curv,i])*(TimeStep/i)))) * (i/ (i + TilesWait)) 
    
  }
  lines(Reward_wait, type = 'l', col = Color[curv], lwd = 3 )
  Optimal_Tiles_Revealed[curv] <- which.max(Reward_wait)
  reward_list[[curv]] <- Reward_wait
}

Color <- c("#e0ecf4","#9ebcda", "#8856a7")
legend(0.1, 50, title="Effect Size", legend=c("Small", "Medium", "Large"), 
       col=c(Color), lty=1, lwd=3)

####average payoff across different tile ratios###
reward_list_f <- rbind(reward_list[[1]], reward_list[[2]],reward_list[[3]])
reward_list_f <- as.matrix(reward_list_f)
reward_f <- colMeans(reward_list_f) #average reward to opponent across all tiles
average_opp <- which.max(reward_f)

#####################
#How competition affects payoff
#############################
Competitor_Tile <- 1:25
Player_Tile <- 1:25
Reward_player <- NULL
Rewards_list <- list()
Rewards_Curves <- list()

for(curv in 1:length(Curve)) {
for(comp in 1:length(Competitor_Tile)) {
  for (player in 1:length(Player_Tile)) {
    
    if(comp < player) {          #when competitor guesses first
      Reward_player[player] <- (1 * (1 - Information_Subset[curv, comp] ) * Information_Subset[curv, player]) - 
                               (1 * (1 - Information_Subset[curv, comp] ) * (1 - Information_Subset[curv, player]))
    }
    
    if(comp >= player) {          #when player guesses first
      Reward_player[player] <- (1 * Information_Subset[curv, player]) - (1 * (1 - Information_Subset[curv, player]))
    }
  }
  Rewards_list[[comp]] <- Reward_player
}
  Rewards_Curves[[curv]] <- melt(Rewards_list)
  Rewards_Curves[[curv]]$PlayerTiles <- rep(1:25, 25)
}
#plotting
Curve_Plots <- list()
ytiles <- c(8, 10, 12)
ratio <- c("8:17", "10:15", "12:13")
for(i in 1:length(Rewards_Curves)) {
  
  pp <- ggplot(Rewards_Curves[[i]], aes(L1, PlayerTiles, value)) +
    geom_raster(aes(fill = value), interpolate = FALSE) +
    scale_x_continuous(name = "Number of Tiles Revealed by Opponent", expand = c(0, 0)) +
    scale_y_continuous(name = "Number of Tiles Revealed by Player in Competition Treatment", expand = c(0, 0)) +
    scale_fill_viridis(name = "Player \nPayoff\n ", limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1))  + 
    ggtitle(paste(ratio[i], 'Tile Ratio')) +
    theme_bw()  +
    #eliminates background, gridlines, and chart border
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank() , 
      plot.title = element_text(hjust = 0.5))
  Curve_Plots[[i]] <- pp
}
Curve_Plots[[1]]
Curve_Plots[[2]]
Curve_Plots[[3]]

#########
#assuming that the opponent (i.e. player in the no-competition condition) plays optimally 
#(i.e. get the value from Optimal_Tiles_Revealed); then for this value, calculate the expected payoff \
#to a participant who guesses after 1:25 tiles, and store this in a dataframe
############

df <- as.data.frame(cbind(Curve, Optimal_Tiles_Revealed))
df$Optimal_Tiles_Revealed <- average_opp
Reward_wait <- NULL
Reward_list <- list()

for(z in 1:nrow(df)) {
  
  opt_tiles <- df$Optimal_Tiles_Revealed[z]
  crv <- df$Curve[z]
  
  for(i in 1:opt_tiles) {
    Reward_wait[i] <- 
      ((Information_Subset[crv,i]  - (1 * ((1-Information_Subset[crv,i])))))
  }
  
  if(opt_tiles < 25) {
    for(i in (opt_tiles+1):25) {
      Reward_wait[i] <- (1 * (1 - Information_Subset[crv,opt_tiles] ) * (Information_Subset[crv,i]))
      - (1 * (1 - Information_Subset[crv,opt_tiles]) * (1 - Information_Subset[curv, i]))
    } }
  
  Reward_list[[z]] <- Reward_wait 
}

Curve_8yellow <- data.frame(Payoff = Reward_list[[1]], Tiles = 1:25)
Curve_10yellow <- data.frame(Payoff = Reward_list[[2]], Tiles = 1:25)
Curve_12yellow <- data.frame(Payoff = Reward_list[[3]], Tiles = 1:25)  

plot(Curve_8yellow$Tiles, Curve_8yellow$Payoff)
plot(Curve_10yellow$Tiles, Curve_10yellow$Payoff)
plot(Curve_12yellow$Tiles, Curve_12yellow$Payoff)

df_payoff <- rbind(Curve_8yellow, Curve_10yellow, Curve_12yellow)
df_payoff$NYellow <- c(rep(17, 25), rep(15, 25), rep(13, 25))

#plotting
pp <- ggplot(df_payoff, aes(Tiles, NYellow, Payoff)) +
  geom_raster(aes(fill = Payoff), interpolate = FALSE) +
  scale_y_continuous(name = "Tile Ratio", expand = c(0, 0), breaks = c(17, 15, 13), 
                     labels = c("8:17", "10:15", "12:13")) +
  scale_x_continuous(name = "Number of Tiles Revealed by Player", expand = c(0, 0)) +
  scale_fill_viridis(name = "Expected \nPayoff")  + 
  theme_bw()  +
  #eliminates background, gridlines, and chart border
  ggtitle("SD = 0") +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank() , 
    plot.title = element_text(hjust = 0.5))

pp

#############################
#assuming that an opponent chooses how
#many tiles to reveal by drawing a random number from a normal distribution centered
#around the optimal number of tiles, what is the expected payoff to players who
#guess after anywhere from 1:25 tiles?
#########################

df <- as.data.frame(cbind(Curve, Optimal_Tiles_Revealed))
df$Optimal_Tiles_Revealed <- average_opp
Reward_list <- list()
sd <- 8 # change to increase or decrease uncertainty in when an opponent guesses

for(z in 1:nrow(df)) {
  
  Reward_wait <- NULL
  all_rewards <- NULL
  
  opt_tiles <- round(rnorm(1e4, mean = df$Optimal_Tiles_Revealed[z], sd = sd))
  opt_tiles <- pmax(opt_tiles, 1)
  opt_tiles <- pmin(opt_tiles, 25)
  
  crv <- df$Curve[z]
  
  for(opp_guess in 1:length(opt_tiles)){
    
    for(i in 1:opt_tiles[opp_guess]) {
      Reward_wait[i] <- 
        ((Information_Subset[crv,i] )) - (1 * ((1-Information_Subset[crv,i])))
    }
    
    if(opt_tiles[opp_guess] < 25) {
      for(i in (opt_tiles[opp_guess]+1):25) {
        Reward_wait[i] <- (1 * (1 - Information_Subset[crv,opt_tiles[opp_guess]] ) * (Information_Subset[crv,i]))
        - (1 * (1 - Information_Subset[crv,opt_tiles[opp_guess]] * (1 - Information_Subset[curv, i])))
      } }
    
    all_rewards <- c(all_rewards, Reward_wait)
    Reward_wait <- NULL
  } 
  
  payoffs_to_player <- matrix(all_rewards, nrow = length(opt_tiles), byrow = TRUE)
  
  # average reward to player for guessing anywhere from 1:25 tiles when there is uncertainty about when opponent guesses
  Reward_list[[z]] <-colMeans(payoffs_to_player)
  
}

Curve_8yellow <- data.frame(Payoff = Reward_list[[1]], Tiles = 1:25)
Curve_10yellow <- data.frame(Payoff = Reward_list[[2]], Tiles = 1:25)
Curve_12yellow <- data.frame(Payoff = Reward_list[[3]], Tiles = 1:25)  

plot(Curve_8yellow$Tiles, Curve_8yellow$Payoff)
plot(Curve_10yellow$Tiles, Curve_10yellow$Payoff)
plot(Curve_12yellow$Tiles, Curve_12yellow$Payoff)

df_payoff <- rbind(Curve_8yellow, Curve_10yellow, Curve_12yellow)
df_payoff$NYellow <- c(rep(17, 25), rep(15, 25), rep(13, 25))

#plotting
pp <- ggplot(df_payoff, aes(Tiles, NYellow, Payoff)) +
  geom_raster(aes(fill = Payoff), interpolate = FALSE) +
  scale_y_continuous(name = "Tile Ratio", expand = c(0, 0), breaks = c(17, 15, 13), 
                     labels = c("8:17", "10:15", "12:13")) +
  scale_x_continuous(name = "Number of Tiles Revealed by Player", expand = c(0, 0)) +
  scale_fill_viridis(name = "Expected \nPayoff")  + 
  theme_bw()  +
  ggtitle("SD = 8") +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank() , 
    plot.title = element_text(hjust = 0.5))

pp




