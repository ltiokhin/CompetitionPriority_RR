library(ggplot2)
library(reshape2)
library(viridisLite)
library(viridis)
library(RColorBrewer)
library(rethinking)

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
    
    if(comp > player) {          #when player guesses first
      Reward_player[player] <- (1 * Information_Subset[curv, player]) - (1 * (1 - Information_Subset[curv, player]))
    }
    
    if(comp == player) {          #when players guess at same time,
      Reward_player[player] <- 0.5 * Information_Subset[curv, player]^2 - ((1 - Information_Subset[curv, player])^2) 
            }
  }
  Rewards_list[[comp]] <- Reward_player
}
  Rewards_Curves[[curv]] <- melt(Rewards_list)
  Rewards_Curves[[curv]]$PlayerTiles <- rep(1:25, 25)
}
#plotting
# Curve_Plots <- list()
# ytiles <- c(8, 10, 12)
# ratio <- c("8:17", "10:15", "12:13")
# for(i in 1:length(Rewards_Curves)) {
#   
#   pp <- ggplot(Rewards_Curves[[i]], aes(L1, PlayerTiles, value)) +
#     geom_raster(aes(fill = value), interpolate = FALSE) +
#     scale_x_continuous(name = "Number of Tiles Revealed by Opponent", expand = c(0, 0)) +
#     scale_y_continuous(name = "Number of Tiles Revealed by Player in Competition Treatment", expand = c(0, 0)) +
#     scale_fill_viridis(name = "Player \nPayoff\n ", limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1))  + 
#     ggtitle(paste(ratio[i], 'Tile Ratio')) +
#     theme_bw()  +
#     #eliminates background, gridlines, and chart border
#     theme(
#       plot.background = element_blank()
#       ,panel.grid.major = element_blank()
#       ,panel.grid.minor = element_blank()
#       ,panel.border = element_blank() , 
#       plot.title = element_text(hjust = 0.5))
#   Curve_Plots[[i]] <- pp
# }
# Curve_Plots[[1]]
# Curve_Plots[[2]]
# Curve_Plots[[3]]

##########################
#evolution towards equilibrium
######################
payoffs_m_1 <- matrix(Rewards_Curves[[1]]$value, nrow = 25)
payoffs_m_1[,16]

###################
# function to play a vector of scientists against each other #####
#######################
lifespan <- 100

#probability of removing anywhere from 1 to 25 tiles
play_ess <- function(playera, playerb) {
  
  payoffs <- rep(0, 2)

  for(i in 1:lifespan) {
    
    tile_a <- sample(1:25, 1, prob = playera)
    tile_b <- sample(1:25, 1, prob = playerb)
    
    payoffs[1] <- payoffs[1] + payoffs_m_1[tile_a, tile_b]
    payoffs[2] <- payoffs[2] + payoffs_m_1[tile_b, tile_a]
    
}
  return(payoffs) 
}

##################
#1 species; mixed strategies; 1:25
################
run_simulation <- function(N, RR, G, R) {
  
  popsize <- N # sets population size
  rounds <- RR # how many other scientists each scientist faces per generation
  gens <- G # number of generations
  repeats <- R # number of simulations for every unique combo of effect and startup cost
  
    eq.tileprob <- matrix(rep(0, 25), nrow = gens, ncol = 25)
    eq.totalfitness <- vector()
            
            mean_tile_prob <- matrix(rep(0, 25), nrow = gens, ncol = 25)
            mean_total_fitness <- rep(0, gens)
            all_players <- matrix(rep(0, 25), nrow = popsize, ncol = 25)
            
            for(rep in 1:repeats) {
              
              # initialize the population, for each repeat. Each player is a row in the matrix. 
              for(i in 1:popsize){
                playr <- runif(25, 0, 1)
                all_players[i,] <- round(playr / sum(playr), 3)
              }
              # start looping
              for (gen in 1:gens) {
                
                #make sure fitness values are 0 at the start of each generation
                fitness <- rep(0.0000001, popsize)
                
                for (round in 1:rounds) {
                  # for each round randomize the partners
                  indexes <- sample(1:popsize, size=popsize)
                  all_players <- all_players[indexes,]
                  fitness <- fitness[indexes]
                  
                  for (i in 1:(popsize/2)) {
                    # pairs play against each other
                    competitor_a <- all_players[(i*2 - 1),]
                    competitor_b <- all_players[(2*i),]
                    dum <- play_ess(competitor_a, competitor_b)
                    fitness[(i*2 - 1):(2*i)] <- fitness[(i*2 - 1):(2*i)] + dum
                  } # end of round
                }
                #save state of the population sample sizes and fitness
                mean_tile_prob[gen,] <- mean_tile_prob[gen,] + apply(all_players, 2, mean)
                mean_total_fitness[gen] <- mean_total_fitness[gen] + sum(fitness)
                
                # calculate fitness and manage reproduction
                fitness2 <- fitness/sum(fitness)
                ss <- sample(1:nrow(all_players), size = popsize, replace=TRUE, prob=fitness2)
                all_players <- all_players[ss,]

                for(i in 1:popsize){
                  all_players[i,] <- all_players[i,] + rnorm(25, 0, 0.001)
                  all_players[i,] <- pmin(pmax(all_players[i,], 0), 1)
                  all_players[i,] <- all_players[i,] / sum(all_players[i,])
                }
                
              } # end of all generations
              eq.tileprob <- c(eq.tileprob, mean_tile_prob[gens,])
              eq.totalfitness <- c(eq.totalfitness, sum(fitness))
            } 
            # end of repeats
              mean_tile_prob <- mean_tile_prob/repeats
              mean_total_fitness <- mean_total_fitness/repeats

  return(list(all_players, mean_tile_prob))
}

run_simulation(10,     5,     10,        1)

r_1species_mixed_03mutation_5runs <- run_simulation(100,     5,     1000,        20)
r_1species_mixed_03mutation_40runs <- run_simulation(100,     10,     1000,        40)


# #                               popsize, rounds, generations, repeats
# results_1species_1 <- run_simulation(100,     2 ,     1000,        4)
# results_1species_2 <- run_simulation(100,     5,      1000,        4)
# results_1species_3 <- run_simulation(100,     5,      1000,        4)
# results_1species_4 <- run_simulation(100,     5,      1000,        4)

for(i in c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)){
  
plot(1:25, r_1species_mixed_03mutation_5runs[[1]][i,] / sum(r_1species_mixed_03mutation_5runs[[1]][i,]), 
     ylim = c(0, 0.1), type = "b", main = paste("5 runs, mixed species probabalistic, generation", i))
  
}

plot(1:25, results_1species_1[[1]][1000,], ylim = c(0, 0.1))
plot(1:25, results_1species_2[[1]][1000,], ylim = c(0, 0.1))
plot(1:25, results_1species_3[[1]][1000,], ylim = c(0, 0.1))


for(i in c(1, 100, 200, 300, 400)){
  
  plot(1:25, aaa[[1]][i,], ylim = c(0, 0.2))
  
}

###take home message when populations are mixed - - - nothing happens
######################

#################
#1 species; mixed strategies; 15:25
#This shit works. Everyone goes to guessing the smallest number of allowable tiles, 15 or 16. 
################
run_simulation_15_25 <- function(N, RR, G, R) {
  
  popsize <- N # sets population size
  rounds <- RR # how many other scientists each scientist faces per generation
  gens <- G # number of generations
  repeats <- R # number of simulations for every unique combo of effect and startup cost
  
  eq.tileprob <- matrix(rep(0, 25), nrow = gens, ncol = 25)
  eq.totalfitness <- vector()
  
  mean_tile_prob <- matrix(rep(0, 25), nrow = gens, ncol = 25)
  mean_total_fitness <- rep(0, gens)
  all_players <- matrix(rep(0, 25), nrow = popsize, ncol = 25)
  
  for(rep in 1:repeats) {
    
    # initialize the population, for each repeat. Each player is a row in the matrix. 
    for(i in 1:popsize){
      playr <- runif(25, 0, 1)
      all_players[i,] <- round(playr / sum(playr), 3)
    }
    
    #sets anything below tile of 14 to 0. 
    all_players[,1:14] <- 0
    
    # start looping
    for (gen in 1:gens) {
      
      #make sure fitness values are 0 at the start of each generation
      fitness <- rep(0.0000001, popsize)
      
      for (round in 1:rounds) {
        # for each round randomize the partners
        indexes <- sample(1:popsize, size=popsize)
        all_players <- all_players[indexes,]
        fitness <- fitness[indexes]
        
        for (i in 1:(popsize/2)) {
          # pairs play against each other
          competitor_a <- all_players[(i*2 - 1),]
          competitor_b <- all_players[(2*i),]
          dum <- play_ess(competitor_a, competitor_b)
          fitness[(i*2 - 1):(2*i)] <- fitness[(i*2 - 1):(2*i)] + dum
        } # end of round
      }
      
      #save state of the population sample sizes and fitness
      mean_tile_prob[gen,] <- mean_tile_prob[gen,] + apply(all_players, 2, mean)
      mean_total_fitness[gen] <- mean_total_fitness[gen] + sum(fitness)
      
      # calculate fitness and manage reproduction
      fitness2 <- fitness/sum(fitness)
      ss <- sample(1:nrow(all_players), size = popsize, replace=TRUE, prob=fitness2)
      all_players <- all_players[ss,]
      
      for(i in 1:popsize){
        all_players[i,] <- all_players[i,] + rnorm(25, 0, 0.01)
        all_players[i,] <- pmin(pmax(all_players[i,], 0), 1)
      }
      #sets all probabilities for tiles 14 or lower to 0. 
      all_players[,1:14] <- 0
      
    } # end of all generations
    eq.tileprob <- c(eq.tileprob, mean_tile_prob[gens,])
    eq.totalfitness <- c(eq.totalfitness, sum(fitness))
  } 
  # end of repeats
  mean_tile_prob <- mean_tile_prob/repeats
  mean_total_fitness <- mean_total_fitness/repeats
  
  return(list(mean_tile_prob))
}

aaa <- run_simulation_15_25(100,     5,     1000,        5)

for(i in seq(1, 1000, by = 50)){
  
  plot(1:25, aaa[[1]][i,] / sum(aaa[[1]][i,]), ylim = c(0, 0.1), main = i, type = "b")
  
}


###fixed strategy play ess###

play_ess_fixed <- function(tile_a, tile_b) {
  
  payoffs <- rep(0, 2)
  
  for(i in 1:2) {
    
    payoffs[1] <- payoffs[1] + payoffs_m_1[tile_a, tile_b]
    payoffs[2] <- payoffs[2] + payoffs_m_1[tile_b, tile_a]
    
  }
  return(payoffs) 
}


#################
#1 species; fixed strategies
################
run_simulation_fixed <- function(N, RR, G, R) {
  
  popsize <- N # sets population size
  rounds <- RR # how many other scientists each scientist faces per generation
  gens <- G # number of generations
  repeats <- R # number of simulations for every unique combo of effect and startup cost

  eq.tiles <- matrix(rep(0, popsize), nrow = repeats, ncol = popsize)
  fitness_acrossgens <- matrix(rep(0, popsize), nrow = gens, ncol = popsize)
  eq.totalfitness <- matrix(rep(0, popsize), nrow = repeats, ncol = popsize)
  mean_tiles <- vector()
  all_players <- vector()

  for(rep in 1:repeats) {
    
    everyone_tiles <- matrix(rep(0, popsize), nrow = gens, ncol = popsize)
    
    # initialize the population, for each repeat. Each player is a row in the matrix. 
    all_players <- round(runif(popsize, 1, 25), 0)
    
    # start looping
    for (gen in 1:gens) {
      
      #make sure fitness values are 0 at the start of each generation
      fitness <- rep(0.0000001, popsize)
      
      for (round in 1:rounds) {
        # for each round randomize the partners
        indexes <- sample(1:popsize, size=popsize)
        all_players <- all_players[indexes]
        fitness <- fitness[indexes]
        
        for (i in 1:(popsize/2)) {
          # pairs play against each other
          competitor_a <- all_players[(i*2 - 1)]
          competitor_b <- all_players[(2*i)]
          dum <- play_ess_fixed(competitor_a, competitor_b)
          fitness[(i*2 - 1):(2*i)] <- fitness[(i*2 - 1):(2*i)] + dum
        } # end of round
      }
      #save state of the population sample sizes and fitness
      everyone_tiles[gen,] <- everyone_tiles[gen,] + all_players
      fitness_acrossgens[gen,] <- fitness/sum(fitness)

      # calculate fitness and manage reproduction
      fitness2 <- fitness/sum(fitness)
      all_players <- sample(all_players, size = popsize, replace=TRUE, prob=fitness2)
        all_players <- all_players + rnorm(length(all_players), 0, 0.4)
        all_players <- round(all_players, 0)
        all_players <- pmin(pmax(all_players, 1), 25)

    } # end of all generations
    eq.tiles[rep,] <- everyone_tiles[gen,]
    eq.totalfitness[rep,] <- fitness/sum(fitness)
  }
  # end of repeats
  mean_tiles <- apply(eq.tiles, 2, mean)

  return(list(everyone_tiles, eq.tiles, mean_tiles, fitness_acrossgens, eq.totalfitness))
}

bbb <- run_simulation_fixed(200,     100,     300,        3)

#plotting 1 repetition, across the time of the simulation
###aaa has results for when mutation is 2
for(i in seq(1, 300, by = 50)){
  
  hist(bbb[[4]][i,], main = paste("fitness at across gens, gen", i))
}

hist(bbb[[3]], main = "averaging across runs", breaks = c(1:25))

hist[]


##plotting the distribution of sample sizes at the end of 2000 generations of evolution
for(i in 1:10){
  
  hist(bbb[[2]][i,], main = paste("mutation 3, distribution at the end of repeat", i), 
       breaks = c(1:25))
  
}
















































####ess###

payoffs_big_e <- Rewards_Curves[[3]]

best_strat <- vector()
for(i in 1:25){
  subb <- payoffs_big_e[payoffs_big_e$L1 == i,]
  best_strat[i] <- which.max(subb$value)
}

######
###OR###
payoffs_av <- cbind(Rewards_Curves[[1]]$L1, Rewards_Curves[[1]]$PlayerTiles, Rewards_Curves[[1]]$value, 
                    Rewards_Curves[[2]]$value, Rewards_Curves[[3]]$value)
best_strat <- vector()

for(i in 1:25){
  subb <- payoffs_av[payoffs_av[,1] == i,]
  mean_pay <- apply(subb[,c(3,4,5)], 1, mean)
  best_strat[i] <- which.max(mean_pay)
}
########

#plotting###

plot(1:25, best_strat, xlim = c(0, 25), ylim=c(0, 25),
     ylab = "Best Strategy Player", 
     xlab = "Number Tiles Revealed by Competitor", 
     pch = 17, 
     main = "Effect Size 12:13")
lines(1:25, best_strat, col=col.alpha(rangi2,0.7), lwd=3)
grid(5, 5, lwd = 2) # grid only in y-direction

plot(1:25, best_strat, xlim = c(0, 25), ylim=c(0, 25),
     ylab = "Best Strategy Player", 
     xlab = "Number Tiles Revealed by Competitor", 
     pch = 17, 
     main = "Averaging across effect sizes")
lines(1:25, best_strat, col=col.alpha(rangi2,0.7), lwd=3)
grid(5, 5, lwd = 2) # grid only in y-direction


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

##averaging across different tile-ratios; including this assumes that a player has no information whatsoever about the
#underlying tile ratio and they just try to maximize their expected payoff across the 3 different possible tile ratios that they
#encounter
list_f <- rbind(Reward_list[[1]], Reward_list[[2]],Reward_list[[3]])
list_f <- as.matrix(list_f)
reward_comp_f <- colMeans(list_f) #average reward to opponent across all tiles
average_against_opp <- which.max(reward_comp_f)

Curve_8yellow <- data.frame(Payoff = Reward_list[[1]], Tiles = 1:25)
Curve_10yellow <- data.frame(Payoff = Reward_list[[2]], Tiles = 1:25)
Curve_12yellow <- data.frame(Payoff = Reward_list[[3]], Tiles = 1:25)  

plot(Curve_8yellow$Tiles, Curve_8yellow$Payoff)
plot(Curve_10yellow$Tiles, Curve_10yellow$Payoff)
plot(Curve_12yellow$Tiles, Curve_12yellow$Payoff)
#plotting average; obvious that, if you have no information at all about which effect you are encountering, then you should
#guess as late as possible. 
plot(1:25, reward_comp_f)

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




