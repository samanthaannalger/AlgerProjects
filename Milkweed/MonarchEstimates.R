# S. Alger, P Alexander Burnham
# Estimating Monarch population at BorderView Farm
# July 25, 2018

###### SET UP DATA SET AND VARS FOR CALCULATIONS: #######
MWm2 <- 25.5 # density of plants per 1 meter squared
plotArea <- 137088 # area of plot in meters squared
plantsSE <- 300 # number of plants sampled at each time step
eggCounts <- c(0,0,0,1,1,3,0) # vector of counts of eggs at each of 7 time steps

# calculate total number of plants at BV:
plantsInPlot <- MWm2 * plotArea

# calculate the number of estimated eggs in plot for each time step: 
totalEggsVec <- (eggCounts/plantsSE) * plantsInPlot 

# sum them up for total eggs at BV so far:
sum(totalEggsVec)
