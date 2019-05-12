# AIP Mapping
# S. Alger
# Sep. 12, 2018

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/VTApiaries/")

# read in data:

AIPReg <- read.csv("CSV_files/AIPRegistrationType.csv",                    header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)

states <- read.csv("CSV_files/USstates.csv", 
                   header=TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE)

#Get libraries
library(usmap)
library(maps)
library(ggplot2)
library(noncensus)

#download fips codes dataset
data(counties)




# merge dataframes so that states are listed with their 'fips IDs' (needed to map)
AIPReg <- merge(AIPReg, states, by= "State")

# Change column name abb -> state
names(AIPReg)[names(AIPReg) == "abb"] <- "state"


#continuing to merge dataframes so states are listed with their 'fips IDs' (needed to map)
AIPReg <- merge(AIPReg, counties, by = "state")

# Remove duplicates based on state (getting rid of excess data from fips database)
AIPReg <- AIPReg[!duplicated(AIPReg["State"]),]

#Selecting only columns I want
AIPReg <- dplyr::select(AIPReg, state, State, OverallRegistrationType, Contacted,StateID, state_fips)



#Create Map showing registration types
# reorder for plotting
AIPReg$OverallRegistrationType <- factor(AIPReg$OverallRegistrationType, levels = c("Mandatory Commercial","Mandatory Universal","Voluntary", "No Registration Program","N/A"))

# choose colors
fill.colors <- c("orange", "steelblue2", "gold", "snow4", "white")

#plot map
usmap::plot_usmap(data = AIPReg, values = "OverallRegistrationType", lines = "black") + theme(legend.position = "right", legend.text=element_text(size=15)) + guides(fill=guide_legend(title="")) + scale_fill_manual(values = fill.colors)

#create map showing whether Faith was able to contact the AIP or not
#choose colors
fill.colors <- c("grey", "orange")

#plot map
usmap::plot_usmap(data = AIPReg, values = "Contacted", lines = "black") + theme(legend.position = "right", legend.text=element_text(size=15)) + guides(fill=guide_legend(title="")) + scale_fill_manual(values = fill.colors)
