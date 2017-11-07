rm(list=ls())

## Read in the data 
dat <- read.table('Data_for_HW_1.csv', header=TRUE, sep = ',')

## Call the states 1-4 because why use 1,2,5,6??
dat$State[dat$State == 5] <- 3
dat$State[dat$State == 6] <- 4
