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

