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


#####Histogram of guesses as a function of belief

rm(list=ls())
#ptm <- proc.time()
NumberOfTiles <- 25

replicats <- 1000 # Number of replicates
tiles_revealed <- rep(0, replicats)

## Generation the underlying state of the grid 
N_Yellow <- 12 # Number of Yellow tiles (has to be lower than NumberOfTiles/2)
N_Blue <- NumberOfTiles - N_Yellow # Number of Blue tiles
Tiles <- c(rep(0, N_Yellow), rep(1, N_Blue))

#how confident do you want to be before guessing (e.g. 80% = 0.2)
threshold_l <- 0.2
threshold_u <- 1 - threshold_l

## Simulations
theta<-seq(0,1,0.001) #create theta range from 0 to 1

# Start                    
for (j in 1 :replicats){
  
  Observation <- sample(Tiles, NumberOfTiles, replace = F) # Shuffle the grid to generate a sequence of observation
  aposterior<-1 #Set the alpha for the Beta distribution for the prior
  bposterior<-1 #Set the beta for the Beta distribution for the prior
  LL <- 0
  
  for(i in 1:length(Observation)){
    n <- i #total trials
    if(Observation[i] == 1){
      aposterior <- aposterior + 1
    } else{
      bposterior <- bposterior + 1
    } 
    posterior <- dbeta(theta, aposterior, bposterior) #determine posterior distribution
    LL <- qbeta(threshold_l, aposterior, bposterior) #calculate lower limit credible interval
    UL <- qbeta(threshold_u ,aposterior, bposterior) #calculate upper limit credible interval
    
    #if there is at least an XX% probability that there are more yellow than blue, then guess at that number of tiles
    if(LL > 0.5 | UL < 0.5 | aposterior == 14){
      tiles_revealed[j] <- n
      break
    } }
  
} 

simplehist(tiles_revealed, 
           main = paste("Number Minority = 12, 1000 repeats; 
                        Players guess with", threshold_u, 
                        "confidence; Mean Tiles=", mean(tiles_revealed)))

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
df.comp <- d.conf.agg$TilesRevealed[d.conf.agg$Competition==0]
df.nocomp <- d.conf.agg$TilesRevealed[d.conf.agg$Competition==1]





