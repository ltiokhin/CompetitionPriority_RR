#This is the 2nd version of Information_25Tiles.R . The first version uploaded to OSF had the following
#bug, which leads to NA information values:
## Remove mirrored sequences (Keep only those that start with 0)
## SpaceOfPossibilities <- SpaceOfPossibilities[-which(SpaceOfPossibilities[,1] == 0),]
##############################
library(rethinking)
library(tidyverse)
library(knitr)
library(kableExtra)
library(dplyr)
library(sjPlot)
library(stringr)
library(RColorBrewer)
library(reshape2)

#############
###Data from RR study (no excluding individual observations)
############

###Load Data###
sequences <- read.table("Sequences1.txt", header=FALSE) #the sequence of blue and yellow tiles for each grid

###Cleaning and Data Prep###
d.conf_1_2.c <- d.conf_1_2[complete.cases(d.conf_1_2),]

#exclude participants who indicated technical difficulties
d.conf_1_2.c <- d.conf_1_2.c[d.conf_1_2.c$ID_Player != 26,]
d.conf_1_2.c <- d.conf_1_2.c[d.conf_1_2.c$ID_Player != 61,]

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

for(i in 1:nrow(d.conf_1_2.c)){
  if(d.conf_1_2.c$Guess_Number[i] %in% g_13){d.conf_1_2.c$n_major[i] <- 13}
  if(d.conf_1_2.c$Guess_Number[i] %in% g_15){d.conf_1_2.c$n_major[i] <- 15}
  if(d.conf_1_2.c$Guess_Number[i] %in% g_17){d.conf_1_2.c$n_major[i] <- 17}
}

#Aggregate data for each participant
d.conf.agg <- aggregate(TilesRevealed ~ Effort + Competition+
                          Guess_Number + n_major + ID_Player, data=d.conf_1_2.c, FUN = mean)

################
#Plot tiles revealed in experiment as function of effect size -------------------
###############

#separate objects just based on competition (not segregating on effort)
esize <- c(13, 15, 17)
e_label <- c("Small Effect.", "Medium Effect.", "Large Effect.")

for(e in 1:length(esize)) {

df.comp <- d.conf.agg$TilesRevealed[d.conf.agg$Competition==1 & d.conf.agg$n_major==esize[e]]
df.comp <- as.data.frame(df.comp)
mean_comp <- round(mean(df.comp$df.comp), 2)

df.nocomp <- d.conf.agg$TilesRevealed[d.conf.agg$Competition==0 & d.conf.agg$n_major==esize[e]]
df.nocomp <- as.data.frame(df.nocomp)
mean_nocomp <- round(mean(df.nocomp$df.nocomp), 2)

a <- ggplot(df.nocomp, aes(x = df.nocomp))
a <- a + geom_histogram(fill = "#3182bd", alpha = 0.8, bins = 26) + guides(fill=FALSE) +
  scale_x_continuous(name = "Tiles Revealed", breaks = seq(0, 25, by = 5), limits=c(-1, 26)) +
  ylab("Frequency") + 
  ggtitle(paste("No Competition,",  e_label[e], "Mean = ", mean_nocomp)) + 
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(data=df.nocomp, aes(xintercept=mean_nocomp),
             linetype="dashed", size=1, colour="#e41a1c")

plot(a)

a <- ggplot(df.comp, aes(x = df.comp))
a <- a + geom_histogram(fill = "#3182bd", alpha = 0.8, bins = 26) + guides(fill=FALSE) +
  scale_x_continuous(name = "Tiles Revealed", breaks = seq(0, 25, by = 5), limits=c(-1, 26)) +
  ylab("Frequency") + 
  ggtitle(paste("Competition,",  e_label[e], "Mean = ", mean_comp)) + 
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(data=df.comp, aes(xintercept=mean_comp),
             linetype="dashed", size=1, colour="#e41a1c")

plot(a)

}

################
#Comparing optimal bayesian to actual data -------------------
###############

#Medium effect size (change this to small, large, and all 3 effects and see what happens)
#also change the belief threshold as wanted

#ptm <- proc.time()
NumberOfTiles <- 25
replicats <- 1e4 # Number of replicates
tiles_revealed <- rep(0, replicats)
correct <- rep(NA, replicats)

#generate matrix to hold all tile sequences
tile_seq_matrix <- matrix(ncol = 25, nrow = 1e4)
N_Yellow <- sample(c(10, 10, 10), 1e4, replace=TRUE) # Number of Yellow tiles (has to be lower than NumberOfTiles/2)
N_Blue <- NumberOfTiles - N_Yellow # Number of Blue tiles

#sampling sequences, with all possible effects
for(i in 1:replicats){

  Tiles <- c(rep(0, N_Yellow[i]), rep(1, N_Blue[i]))
  Observation <- sample(Tiles, NumberOfTiles, replace = F) # Shuffle the grid to generate a sequence of observation
  tile_seq_matrix[i,] <- Observation
}

#how confident do you want the opponent (no-comp condition) to be before guessing (e.g. 80% = 0.2)
threshold_l <- 0.18
threshold_u <- 1 - threshold_l

## Simulations
theta<-seq(0,1,0.001) #create theta range from 0 to 1

# Start                    
for (j in 1 :replicats){
  
  aposterior<-1 #Set the alpha for the Beta distribution for the prior
  bposterior<-1 #Set the beta for the Beta distribution for the prior
  LL <- 0
  
  for(i in 1:length(Tiles)){
    n <- i #total trials
    if(tile_seq_matrix[j,i] == 1){
      aposterior <- aposterior + 1
    } else{
      bposterior <- bposterior + 1
    } 
    posterior <- dbeta(theta, aposterior, bposterior) #determine posterior distribution
    LL <- qbeta(threshold_l, aposterior, bposterior) #calculate lower limit credible interval
    UL <- qbeta(threshold_u ,aposterior, bposterior) #calculate upper limit credible interval
    
    #if there is at least an XX% probability that there are more yellow than blue, then guess at that number of tiles
    if(LL > 0.5) {
      tiles_revealed[j] <- n
      correct[j] <- 1
      break
      
    } else if(UL < 0.5){
      tiles_revealed[j] <- n
      correct[j] <- 0
      break
      
    } else if(aposterior == 14){
      tiles_revealed[j] <- n
      correct[j] <- 1
      break
  }
  } } 

#opponent guesses and whether they were correct in a data frame
df_opponent <- data.frame(tiles_opp = tiles_revealed, 
                          correct = correct, 
                          esize = N_Yellow)

#kernal density plot of tiles revealed by bayesian opponent ----------
#* rerun above simulated data to get whatever data you wnat for plots
#just simulated data
a <- ggplot(df_opponent, aes(x = tiles_opp))
a <- a + geom_density(fill = "#984ea3", alpha = 0.5, bw = 1) + guides(fill=FALSE) +
  scale_x_continuous(name = "Tiles Revealed", breaks = seq(0, 25, by = 5), limits=c(-1, 26)) +
  labs(y = "Density", title = "Medium Effect. Optimal Opponent = 82%") + 
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(data=df_opponent, aes(xintercept=round(mean(df_opponent$tiles_opp)), 2),
             linetype="dashed", size=1, colour="black")

a

#adding actual data
df.nocomp <- d.conf.agg$TilesRevealed[d.conf.agg$Competition==0 & d.conf.agg$n_major==13]
df.nocomp <- as.data.frame(df.nocomp)
df.nocomp$tiles_opp <- df.nocomp$df.nocomp
mean_nocomp <- round(mean(df.nocomp$df.nocomp), 2)

a + geom_density(data = df.nocomp, fill = "#ff7f00", alpha = 0.5, bw = 1) + 
  geom_vline(data=df.nocomp, aes(xintercept=mean_nocomp, 2),
             linetype="dashed", size=1, colour="#ff7f00")

#histogram
# a <- ggplot(df_opponent, aes(x = tiles_opp))
# a + geom_histogram(fill = "#984ea3", alpha = 0.8, bins = 25) + guides(fill=FALSE) +
#   scale_x_continuous(name = "Tiles Revealed", breaks = seq(0, 25, by = 5), limits=c(-1, 26)) +
#   ylab("Frequency") + 
#   ggtitle(paste("Medium Effect. Guessing with 85% Confidence. Mean = ", 
#                 round(mean(df_opponent$tiles_opp), 2))) + 
#   theme_bw(base_size = 14) +
#   theme(plot.title = element_text(hjust = 0.5), 
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   geom_vline(data=df_opponent, aes(xintercept=round(mean(df_opponent$tiles_opp)), 2),
#              linetype="dashed", size=1, colour="#e41a1c")


#############
#payoff for a bayesian who reveals tiles, varying the level of confidence ------------
############ 

#make sure this is the same object as above
tile_seq_matrix <- tile_seq_matrix
confidence_levels <- seq(0.6, 1, by = 0.01)
tiles_revealed <- rep(0, replicats)
correct <- rep(NA, replicats)

##matrixes to store player tiles and whether or not they were correct
matrix_playertiles <- matrix(ncol = length(confidence_levels), nrow = 1e4)
matrix_playercorrect <- matrix(ncol = length(confidence_levels), nrow = 1e4)

for(z in 1:length(confidence_levels)){

#how confident do you want to be before guessing (e.g. 80% = 0.2)
threshold_u <- confidence_levels[z]
threshold_l <- 1 - threshold_u

## Simulations
theta<-seq(0,1,0.001) #create theta range from 0 to 1

# Start                    
for (j in 1 :replicats){
  
  aposterior<-1 #Set the alpha for the Beta distribution for the prior
  bposterior<-1 #Set the beta for the Beta distribution for the prior
  LL <- 0
  
  for(i in 1:length(Tiles)){
    n <- i #total trials
    if(tile_seq_matrix[j,i] == 1){
      aposterior <- aposterior + 1
    } else{
      bposterior <- bposterior + 1
    } 
    posterior <- dbeta(theta, aposterior, bposterior) #determine posterior distribution
    LL <- qbeta(threshold_l, aposterior, bposterior) #calculate lower limit credible interval
    UL <- qbeta(threshold_u ,aposterior, bposterior) #calculate upper limit credible interval
    
    #if there is at least an XX% probability that there are more yellow than blue, then guess at that number of tiles
    if(LL > 0.5) {
      tiles_revealed[j] <- n
      correct[j] <- 1
      break
      
    } else if(UL < 0.5){
      tiles_revealed[j] <- n
      correct[j] <- 0
      break
      
    } else if(aposterior == 14){
      tiles_revealed[j] <- n
      correct[j] <- 1
      break
    }
  } }

matrix_playertiles[,z] <- tiles_revealed
matrix_playercorrect[,z] <- correct
}

#######
#calculate payoffs to the player, when playing against the opponent from df_opponent--------
#######

#data frame to store results
df_payoff_confidence <- matrix(ncol=length(confidence_levels), 
                               nrow = 1e4)

#currently this is correct, but assumes that the player only gets scooped if they guess after opponent
for(z in 1:length(confidence_levels)){
  for(i in 1:nrow(df_payoff_confidence)){
    
    #if opponent scoops you by guessing earlier
    if(df_opponent$tiles_opp[i] < matrix_playertiles[i,z] & df_opponent$correct[i] == 1){
      df_payoff_confidence[i,z] <- 0
    } else if(matrix_playercorrect[i,z] == 1){
      df_payoff_confidence[i,z] <- 1 #else if opponent doesn't scoop you and you get it right, you get a point
    } else {
      df_payoff_confidence[i,z] <- -1 #else if you don't get it right, it means you got it wrong
    }
  }
}

df_final_player_payoffs <- data.frame(mean_payoff = colMeans(df_payoff_confidence), 
                                      confidence = confidence_levels, 
                                      mean_tiles = colMeans(matrix_playertiles))

plot(df_final_player_payoffs$mean_payoff ~ df_final_player_payoffs$confidence, type = "b", lwd = 2,
     col="#cb181d", ylab = "Mean Payoff", xlab = "Player's Confidence Level When Guessing", 
     main = "65% Confident Opponent. 3 Effects.")
abline(v = 0.65, col="black", lwd=2, lty=2)

######plot of overlapping distributions of simulated datasets ---------

#select the matrix player tiles confidence level that you want
df_overlapping <- data.frame(Tiles = c(df_opponent$tiles_opp, matrix_playertiles[,1]), 
                             Who = c(rep("Opponent", 1e4), rep("Player", 1e4)))

####kernal density plot
a <- ggplot(df_overlapping, aes(x = Tiles, fill = Who))
a + geom_density(alpha = 0.4, bw = 0.7) +
  scale_x_continuous(name = "Tiles Revealed", breaks = seq(0, 25, by = 5), limits=c(-1, 26)) +
  ylab("Density") + 
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

############
#Individually optimal number of tiles to reveal in the game, disocunting by number of seconds ------
############ 

rm(list=ls())  

#ptm <- proc.time()
NumberOfTiles <- 25
replicats <- 1e3 # Number of replicates

yellow_tile_matrix <- matrix(ncol = 4, nrow = 1e4)
yellow_tile_matrix[,1] <- sample(c(8, 10, 12), 1e4, replace=TRUE)
yellow_tile_matrix[,2]  <- rep(8, 1e4)
yellow_tile_matrix[,3]  <- rep(10, 1e4)
yellow_tile_matrix[,4]  <- rep(12, 1e4)

names <- c("All Effects,", "Large Effect,", "Medium Effect,", "Small Effect,")
openlist <- list()

for(E in 4){

N_Yellow <- yellow_tile_matrix[,E]
N_Blue <- NumberOfTiles - N_Yellow # Number of Blue tiles

#confidence levels
confidence_levels <- seq(0.6, 1, by = 0.01)

#keep track of score and runs
runs <- 1e4 #change to 1e4
matrix_playerscore <- matrix(ncol = length(confidence_levels), nrow = runs)

for(z in 1:length(confidence_levels)){
  
  #how confident do you want to be before guessing (e.g. 80% = 0.2)
  threshold_u <- confidence_levels[z]
  threshold_l <- 1 - threshold_u
  theta<-seq(0,1,0.001) #create theta range from 0 to 1
  
  #start an individual run
  for(run in 1:runs){
    
    #reset score and time
    time_remaining <- 1200
    score <- 0
    
    #generate matrix to hold all tile sequences
    tile_seq_matrix <- matrix(ncol = 25, nrow = 1e3)
    
  #sampling sequences, with all possible effects, to start one individual replicate
  for(i in 1:replicats){
    
    Tiles <- c(rep(0, N_Yellow[i]), rep(1, N_Blue[i]))
    Observation <- sample(Tiles, NumberOfTiles, replace = F) # Shuffle the grid to generate a sequence of observation
    tile_seq_matrix[i,] <- Observation
  }
  
  while(time_remaining > 0){
    
    aposterior<-1 #Set the alpha for the Beta distribution for the prior
    bposterior<-1 #Set the beta for the Beta distribution for the prior
    
    for(i in 1:length(Tiles)){
      
      n <- i #total trials
      
      if(tile_seq_matrix[1,i] == 1){
        aposterior <- aposterior + 1
      } else{
        bposterior <- bposterior + 1
      } 
      
      posterior <- dbeta(theta, aposterior, bposterior) #determine posterior distribution
      LL <- qbeta(threshold_l, aposterior, bposterior) #calculate lower limit credible interval
      UL <- qbeta(threshold_u ,aposterior, bposterior) #calculate upper limit credible interval
      
      #if there is at least an XX% probability that there are more yellow than blue, then guess at that number of tiles
      if(LL > 0.5) {
        time_remaining <- time_remaining - (n+5) # 5 seconds for delay between grids
        score <- score + 1
        break
        
      } else if(UL < 0.5){
        time_remaining <- time_remaining - (n+5)
        score <- score - 1
        break
        
      } else if(aposterior == 14){
        time_remaining <- time_remaining - (n+5)
        score <- score + 1
        break
      }
    }

        #remove first row of tile sequence matrix to prevent redundant tiles
    tile_seq_matrix <- tile_seq_matrix[-1,]
    
  } # end while loop
    matrix_playerscore[run,z] <- score
    
  } #end of run 
    } # end loop for different confidence levels
  
df.indvpay <- data.frame(Reward = colSums(matrix_playerscore) / runs, 
                         Confidence = confidence_levels)

plot(df.indvpay$Reward ~ df.indvpay$Confidence, type = "b", main = paste(names[E], "10k Repeats."), 
 lwd = 2, col="#cb181d", ylab = "Mean Payoff", xlab = "Player's Confidence Level When Guessing")

openlist[[4]] <- df.indvpay
}

save(openlist, file="payoff_optimal_bayesian_all.RData")
lapply(openlist, max)

############
#Best response of a bayesian to another bayesian who plays at any given level of confidence -------
###########

list_bestresponseplayer <- list()

#ptm <- proc.time()
NumberOfTiles <- 25
replicats <- 1e4 # Number of replicates
tiles_revealed <- rep(0, replicats)
correct <- rep(NA, replicats)

#generate matrix to hold all tile sequences
tile_seq_matrix <- matrix(ncol = 25, nrow = 1e4)
N_Yellow <- sample(c(8, 8, 8), 1e4, replace=TRUE) # Number of Yellow tiles (has to be lower than NumberOfTiles/2)
N_Blue <- NumberOfTiles - N_Yellow # Number of Blue tiles

#sampling sequences, with all possible effects
for(i in 1:replicats){
  
  Tiles <- c(rep(0, N_Yellow[i]), rep(1, N_Blue[i]))
  Observation <- sample(Tiles, NumberOfTiles, replace = F) # Shuffle the grid to generate a sequence of observation
  tile_seq_matrix[i,] <- Observation
}

#run everything from here. make sure to customize effect size
opp_conf_levels <- seq(0.6, 1, by = 0.01)

for(RUN in 1:length(opp_conf_levels)){

#how confident do you want the opponent (no-comp condition) to be before guessing (e.g. 80% = 0.2)
threshold_u_opp <- opp_conf_levels[RUN]
threshold_l_opp <- 1 - threshold_u_opp

## Simulations
theta<-seq(0,1,0.001) #create theta range from 0 to 1

# Start                    
for (j in 1 :replicats){
  
  aposterior<-1 #Set the alpha for the Beta distribution for the prior
  bposterior<-1 #Set the beta for the Beta distribution for the prior
  LL <- 0
  
  for(i in 1:length(Tiles)){
    n <- i #total trials
    if(tile_seq_matrix[j,i] == 1){
      aposterior <- aposterior + 1
    } else{
      bposterior <- bposterior + 1
    } 
    posterior <- dbeta(theta, aposterior, bposterior) #determine posterior distribution
    LL <- qbeta(threshold_l_opp, aposterior, bposterior) #calculate lower limit credible interval
    UL <- qbeta(threshold_u_opp ,aposterior, bposterior) #calculate upper limit credible interval
    
    #if there is at least an XX% probability that there are more yellow than blue, then guess at that number of tiles
    if(LL > 0.5) {
      tiles_revealed[j] <- n
      correct[j] <- 1
      break
      
    } else if(UL < 0.5){
      tiles_revealed[j] <- n
      correct[j] <- 0
      break
      
    } else if(aposterior == 14){
      tiles_revealed[j] <- n
      correct[j] <- 1
      break
    }
  } } 

#opponent guesses and whether they were correct in a data frame
df_opponent <- data.frame(tiles_opp = tiles_revealed, 
                          correct = correct, 
                          esize = N_Yellow)

##########
#now the bayesian player who varies level of confidence --------
#########

#make sure this is the same object as above
tile_seq_matrix <- tile_seq_matrix
confidence_levels <- seq(0.6, 1, by = 0.01)
tiles_revealed <- rep(0, replicats)
correct <- rep(NA, replicats)

##matrixes to store player tiles and whether or not they were correct
matrix_playertiles <- matrix(ncol = length(confidence_levels), nrow = 1e4)
matrix_playercorrect <- matrix(ncol = length(confidence_levels), nrow = 1e4)

for(z in 1:length(confidence_levels)){
  
  #how confident do you want to be before guessing (e.g. 80% = 0.2)
  threshold_u <- confidence_levels[z]
  threshold_l <- 1 - threshold_u
  
  ## Simulations
  theta<-seq(0,1,0.001) #create theta range from 0 to 1
  
  # Start                    
  for (j in 1 :replicats){
    
    aposterior<-1 #Set the alpha for the Beta distribution for the prior
    bposterior<-1 #Set the beta for the Beta distribution for the prior
    LL <- 0
    
    for(i in 1:length(Tiles)){
      n <- i #total trials
      if(tile_seq_matrix[j,i] == 1){
        aposterior <- aposterior + 1
      } else{
        bposterior <- bposterior + 1
      } 
      posterior <- dbeta(theta, aposterior, bposterior) #determine posterior distribution
      LL <- qbeta(threshold_l, aposterior, bposterior) #calculate lower limit credible interval
      UL <- qbeta(threshold_u ,aposterior, bposterior) #calculate upper limit credible interval
      
      #if there is at least an XX% probability that there are more yellow than blue, then guess at that number of tiles
      if(LL > 0.5) {
        tiles_revealed[j] <- n
        correct[j] <- 1
        break
        
      } else if(UL < 0.5){
        tiles_revealed[j] <- n
        correct[j] <- 0
        break
        
      } else if(aposterior == 14){
        tiles_revealed[j] <- n
        correct[j] <- 1
        break
      }
    } }
  
  matrix_playertiles[,z] <- tiles_revealed
  matrix_playercorrect[,z] <- correct
}

#######
#calculate payoffs to the player, when playing against the opponent from df_opponent--------
#######

#data frame to store results
df_payoff_confidence <- matrix(ncol=length(confidence_levels), 
                               nrow = 1e4)

#currently this is correct, but assumes that the player only gets scooped if they guess after opponent
for(z in 1:length(confidence_levels)){
  for(i in 1:nrow(df_payoff_confidence)){
    
    #if opponent scoops you by guessing earlier
    if(df_opponent$tiles_opp[i] < matrix_playertiles[i,z] & df_opponent$correct[i] == 1){
      df_payoff_confidence[i,z] <- 0
    } else if(matrix_playercorrect[i,z] == 1){
      df_payoff_confidence[i,z] <- 1 #else if opponent doesn't scoop you and you get it right, you get a point
    } else {
      df_payoff_confidence[i,z] <- -1 #else if you don't get it right, it means you got it wrong
    }
  }
}

list_bestresponseplayer[[RUN]] <- data.frame(mean_payoff = colMeans(df_payoff_confidence), 
                                      confidence = confidence_levels, 
                                      mean_tiles = colMeans(matrix_playertiles), 
                                      confidence_opp = threshold_u_opp)

  } #end for loop for all confidence levels of the opponent
 
#plot optimal confidence level of bayesian player as a function of the confidence level of opponent
data_df <- df_bestresponseplayer[complete.cases(df_bestresponseplayer),]
data_df$confidence[1:15] <- seq(0.6, 0.74, by = 0.01)
  
ggplot(data=data_df, aes(x = confidence_opp, y = confidence, group = 1)) +
  geom_line(aes(y = confidence), size = 1.2, color = "black") + 
  theme_classic(base_size=14) +
  labs(title = "Large Effect", y = "Bayesian Player Confidence Level") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggplot(data=zzz, aes(x = Generation, y = SampleSize)) +
  geom_line(size = 1, color = "#cb181d") + 
  theme_bw(base_size=14) +
  scale_colour_brewer(name="Startup Cost", palette = "Reds")





df_final_player_payoffs[df_final_player_payoffs$mean_payoff == max_pay_player,]






















####To redo the same analysis, but assuming that you get scooped if you guess at the exact same time
###as opponent, run code below
# 
# df_payoff_confidence <- matrix(ncol=length(confidence_levels), 
#                                nrow = 1e4)
# 
# for(z in 1:length(confidence_levels)){
#   for(i in 1:nrow(df_payoff_confidence)){
#     
#     #if opponent scoops you by guessing earlier or at the same time
#     if(df_opponent$tiles_opp[i] <= matrix_playertiles[i,z] & df_opponent$correct[i] == 1){
#       df_payoff_confidence[i,z] <- 0
#     } else if(matrix_playercorrect[i,z] == 1){
#       df_payoff_confidence[i,z] <- 1 #else if opponent doesn't scoop you and you get it right, you get a point
#     } else {
#       df_payoff_confidence[i,z] <- -1 #else if you don't get it right, it means you got it wrong
#     }
#   }
# }
# 
# df_final_player_payoffs <- data.frame(mean_payoff = colMeans(df_payoff_confidence), 
#                                       confidence = confidence_levels)
# 
# plot(df_final_player_payoffs$mean_payoff ~ df_final_player_payoffs$confidence, type = "b", lwd = 2,
#      col="#cb181d", ylab = "Payoff", xlab = "Player's Confidence Level When Guessing", 
#      main = "Payoff against 85% Confident Opponent (Easy to be Scooped)")
