# Samantha Alger
# Negative strand analysis- selecting sample for humberto to run first PCR plate
# 2/26/18


#Preliminaries:
# Clear memory of characters
ls()
rm(list=ls())

setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")
#NegStrand <- read.csv("NegStrd_PosVirus_2015.csv", header=TRUE, stringsAsFactors=FALSE)

BombSurv <- read.csv("BombSurvNHBS.csv", header=TRUE, stringsAsFactors=FALSE)

# code to find negative bees and write out columns wanted
#NegList <- BombSurv[BombSurv$virusBINY==0,]
#NegList <- select(NegList, sample_name, target_name, Ct_mean, norm_genome_copbee)
#write.csv(NegList, file = "NegVirus_NegStrd_2015.csv")

#code to find positive bees and write out columns wanted
PosList <- BombSurv[BombSurv$virusBINY==1,]
PosList <- select(PosList, sample_name, target_name, Ct_mean, norm_genome_copbee, site, apiary_near_far,sumColonies, species)

#change the name of sample_name to ID
colnames(PosList)[1] <- "ID" 

write.csv(PosList, file = "PosVirus_NegStrd_withSites.csv")

PosVirus <- read.csv("PosVirus_NegStrd_withSites.csv", header=TRUE, stringsAsFactors=FALSE)

table(PosVirus$site, PosVirus$apiary_near_far)
table(PosVirus$target_name)

DWVOnly<- PosVirus[PosVirus$target_name == "DWV", ]

View(DWVOnly)
table(DWVOnly$site)

BQCVOnly<- PosVirus[PosVirus$target_name == "BQCV", ]

View(BQCVOnly)
table(BQCVOnly$site, BQCVOnly$species)
table(BQCVOnly$site)

# just the samples that are postive
x <- PosList$ID[!duplicated(PosList$ID)]

# read nanodrop data:
drop <- read.csv("BombSurv_RNANanodropResults.csv",
                  header=TRUE,
                  sep=",",
                  stringsAsFactors=FALSE)

# find number of samples in data set that are postive:
length(PosList$ID[!duplicated(PosList$ID)])

# cross refernce entire data set
complete <- drop[(drop$ID %in% x),]

# select columns we want from DF complete
#complete <- select(complete, ID, ng.ul, final_vol)

# merge data frames to get ng.ul in with main frame
#mergedDF <- merge(x = PosList, y = complete, by.x = "ID")

# write out .csv file:
#write.csv(mergedDF, file = "NegStrd_PosVirus_2015.csv")

# sample for USDA 1 step checking: 87, 67, 362  
#mergedDF[which(mergedDF$target_name=="DWV"),]

# calculate total volume est. left in vial
#complete$totVol <- 45 - complete$RNA_for_Dilution 

# calculate total RNA left in vial 
#complete$totRNA <- complete$totVol * complete$ng.ul

# histgoram of USDA samples ng/ul
#hist(complete$ng.ul, xlim = c(10,350))

# five number summary of the sample concentrations we brought to USDA
#summary(complete$ng.ul)

# how many are above 100 ng/ul
#z <- complete$ng.ul>150
#length(which(z==TRUE))

# which are postive for DWV
#posDWV <- PosList[PosList$target_name=="DWV",]