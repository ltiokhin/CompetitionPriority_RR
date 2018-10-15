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

#separate plots of tiles for each condition
hist(d.conf.agg$TilesRevealed[d.conf.agg$Effort==0 & d.conf.agg$Competition==0])
hist(d.conf.agg$TilesRevealed[d.conf.agg$Effort==0 & d.conf.agg$Competition==1])
hist(d.conf.agg$TilesRevealed[d.conf.agg$Effort==1 & d.conf.agg$Competition==0])
hist(d.conf.agg$TilesRevealed[d.conf.agg$Effort==1 & d.conf.agg$Competition==1])

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

#####Histogram of guesses as a function of belief#####
####medium effect size##### (change this to small and large and see what happens)
####shuffle the tile_seq_matrix and mix in a bunch of effects

#ptm <- proc.time()
NumberOfTiles <- 25
replicats <- 1e4 # Number of replicates
tiles_revealed <- rep(0, replicats)
correct <- rep(NA, replicats)

#generate matrix to hold all tile sequences
tile_seq_matrix <- matrix(ncol = 25, nrow = 1e4)
N_Yellow <- sample(c(8, 10, 12), 1e4, replace=TRUE) # Number of Yellow tiles (has to be lower than NumberOfTiles/2)
N_Blue <- NumberOfTiles - N_Yellow # Number of Blue tiles

#sampling sequences, with all possible effects
for(i in 1:replicats){

  Tiles <- c(rep(0, N_Yellow[i]), rep(1, N_Blue[i]))
  Observation <- sample(Tiles, NumberOfTiles, replace = F) # Shuffle the grid to generate a sequence of observation
  tile_seq_matrix[i,] <- Observation
}

#how confident do you want the opponent (no-comp condition) to be before guessing (e.g. 80% = 0.2)
threshold_l <- 0.05
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

#plot tiles revealed by a bayesian opponent (player in the no-competition condition)
a <- ggplot(df_opponent, aes(x = tiles_opp))
a + geom_histogram(fill = "#984ea3", alpha = 0.8, bins = 25) + guides(fill=FALSE) +
  scale_x_continuous(name = "Tiles Revealed", breaks = seq(0, 25, by = 5), limits=c(-1, 26)) +
  ylab("Frequency") + 
  ggtitle(paste("Medium Effect. Guessing with 85% Confidence. Mean = ", 
                round(mean(df_opponent$tiles_opp), 2))) + 
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(data=df_opponent, aes(xintercept=round(mean(df_opponent$tiles_opp)), 2),
             linetype="dashed", size=1, colour="#e41a1c")

#############
####player, revealing tiles in the same way. varying the level of confidence
############ 

#make sure this is the same object as above
tile_seq_matrix <- tile_seq_matrix
confidence_levels <- seq(0.6, 1, by = 0.01)

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
######calculate payoffs to the player, when playing against the opponent from df_opponent#####
#########

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
     main = "95% Confident Opponent. 3 Effects.")
abline(v = 0.95, col="black", lwd=2, lty=2)


###############
plot(df_final_player_payoffs$mean_tiles ~ df_final_player_payoffs$confidence, type = "b", lwd = 2,
     col="#cb181d", ylab = "Mean Tiles", xlab = "Player's Confidence Level When Guessing", 
     main = "Mean Tiles Revealed. 3 Effects.")



save.image(file = "bayesiancompetitors_payoffs_killinit.RData")


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
