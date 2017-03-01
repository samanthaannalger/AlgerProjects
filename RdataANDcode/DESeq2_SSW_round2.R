# Eco Genomics - Lecture Script
# 27, February 2017
# P. Alexander Burnham


##################################################################################
# Preliminaries:
# Clear memory of characters:
ls()
rm(list=ls())

# set working directory: (files are in /222_Data)
setwd("~/EcologicalGenomics")
##################################################################################

library("DESeq2")

library("ggplot2")

countsTable <- read.delim('RawData/countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("RawData/cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)


#################### MODEL NUMBER 1: TEST EFFECT OF HEALTH CONTROLLING FOR LOCATION

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)
# In this typical model, the "health effect" represents the overall effect of health status 
# controlling for differences due to location page 26-27 manual DESeq2.pdf
# The last term in the model is what is tested.  In this case health. You need to rearrange the order 
# of the factors in the model to test for the overall effect of a diff. factor.
# This is not the same as an interaction.


dim(dds)
#[1] 13053    77

dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
#[1] 12954    77  # at > 100; we loose many fewer genes


# For practice, let's work with fewer genes so that the models can run faster...
dds <- dds[sample(nrow(dds), 1200), ]
dim(dds)

colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) #sets that "healthy is the reference

dds <- DESeq(dds) 

res <- results(dds)
res <- res[order(res$padj),]
head(res)
# log2 fold change (MAP): health S vs H 
# Wald test p-value: health S vs H 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE
# <numeric>      <numeric> <numeric>
#   TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i2_g.23139_m.23139   95.58364       2.916172 0.5141065
# TRINITY_DN45370_c1_g1_TRINITY_DN45370_c1_g1_i1_g.19244_m.19244  326.18876       1.495569 0.2910832
# TRINITY_DN45983_c5_g1_TRINITY_DN45983_c5_g1_i3_g.20837_m.20837 1292.39567       1.514289 0.2904796
# TRINITY_DN44127_c0_g1_TRINITY_DN44127_c0_g1_i1_g.16286_m.16286  250.18334       2.081641 0.4477985
# TRINITY_DN45492_c0_g1_TRINITY_DN45492_c0_g1_i5_g.19587_m.19587  276.71581       1.706664 0.3722197
# TRINITY_DN46924_c3_g1_TRINITY_DN46924_c3_g1_i1_g.23942_m.23942  253.74577       1.146243 0.2648761
# stat           pvalue
# <numeric>        <numeric>
#   TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i2_g.23139_m.23139  5.672312 0.00000001408832
# TRINITY_DN45370_c1_g1_TRINITY_DN45370_c1_g1_i1_g.19244_m.19244  5.137942 0.00000027776334
# TRINITY_DN45983_c5_g1_TRINITY_DN45983_c5_g1_i3_g.20837_m.20837  5.213065 0.00000018574588
# TRINITY_DN44127_c0_g1_TRINITY_DN44127_c0_g1_i1_g.16286_m.16286  4.648611 0.00000334177720
# TRINITY_DN45492_c0_g1_TRINITY_DN45492_c0_g1_i5_g.19587_m.19587  4.585099 0.00000453770438
# TRINITY_DN46924_c3_g1_TRINITY_DN46924_c3_g1_i1_g.23942_m.23942  4.327470 0.00001508320259
# padj
# <numeric>
#   TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i2_g.23139_m.23139 0.000006001624
# TRINITY_DN45370_c1_g1_TRINITY_DN45370_c1_g1_i1_g.19244_m.19244 0.000039442394
# TRINITY_DN45983_c5_g1_TRINITY_DN45983_c5_g1_i3_g.20837_m.20837 0.000039442394
# TRINITY_DN44127_c0_g1_TRINITY_DN44127_c0_g1_i1_g.16286_m.16286 0.000355899272
# TRINITY_DN45492_c0_g1_TRINITY_DN45492_c0_g1_i5_g.19587_m.19587 0.000386612413
# TRINITY_DN46924_c3_g1_TRINITY_DN46924_c3_g1_i1_g.23942_m.23942 0.001070907384

summary(res)

# out of 1199 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 27, 2.3% 
# LFC < 0 (down)   : 14, 1.2% 
# outliers [1]     : 31, 2.6% 
# low counts [2]   : 743, 62% 
# (mean count < 25)


#################### MODEL NUMBER 2 - INTERACTIONS

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health + location:health)

dim(dds)

dds <- dds[ rowSums(counts(dds)) > 100, ]

dds <- dds[sample(nrow(dds), 1200), ]
dim(dds)

colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S"))  #sets that "healthy is the reference

dds <- DESeq(dds, parallel=T)


resultsNames(dds)
# [1]  "Intercept"           "location_sub_vs_int" "health_S_vs_H"       "locationsub.healthS"

res <- results(dds, name="locationsub.healthS")
res <- res[order(res$padj),]
head(res)
# log2 fold change (MLE): locationsub.healthS 
# Wald test p-value: locationsub.healthS 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE
# <numeric>      <numeric> <numeric>
#   TRINITY_DN43786_c1_g1_TRINITY_DN43786_c1_g1_i1_g.15493_m.15493   123.80304      -3.271111 0.6262268
# TRINITY_DN42353_c0_g2_TRINITY_DN42353_c0_g2_i1_g.12717_m.12717 17014.16856       1.440676 0.2997342
# TRINITY_DN47127_c0_g1_TRINITY_DN47127_c0_g1_i3_g.24653_m.24653   264.88795      -1.391774 0.3537203
# TRINITY_DN41442_c3_g5_TRINITY_DN41442_c3_g5_i1_g.11211_m.11211    17.70835      -8.254507 2.0651380
# TRINITY_DN32007_c0_g2_TRINITY_DN32007_c0_g2_i1_g.4284_m.4284      10.50286      -8.011341 2.1704831
# TRINITY_DN28600_c0_g1_TRINITY_DN28600_c0_g1_i1_g.3312_m.3312      30.50985      -5.613068 1.5277455
# stat          pvalue
# <numeric>       <numeric>
#   TRINITY_DN43786_c1_g1_TRINITY_DN43786_c1_g1_i1_g.15493_m.15493 -5.223524 0.0000001755495
# TRINITY_DN42353_c0_g2_TRINITY_DN42353_c0_g2_i1_g.12717_m.12717  4.806512 0.0000015358622
# TRINITY_DN47127_c0_g1_TRINITY_DN47127_c0_g1_i3_g.24653_m.24653 -3.934674 0.0000833096630
# TRINITY_DN41442_c3_g5_TRINITY_DN41442_c3_g5_i1_g.11211_m.11211 -3.997073 0.0000641304917
# TRINITY_DN32007_c0_g2_TRINITY_DN32007_c0_g2_i1_g.4284_m.4284   -3.691041 0.0002233385194
# TRINITY_DN28600_c0_g1_TRINITY_DN28600_c0_g1_i1_g.3312_m.3312   -3.674086 0.0002387024879
# padj
# <numeric>
#   TRINITY_DN43786_c1_g1_TRINITY_DN43786_c1_g1_i1_g.15493_m.15493 0.0001353486
# TRINITY_DN42353_c0_g2_TRINITY_DN42353_c0_g2_i1_g.12717_m.12717 0.0005920749
# TRINITY_DN47127_c0_g1_TRINITY_DN47127_c0_g1_i3_g.24653_m.24653 0.0160579375
# TRINITY_DN41442_c3_g5_TRINITY_DN41442_c3_g5_i1_g.11211_m.11211 0.0160579375
# TRINITY_DN32007_c0_g2_TRINITY_DN32007_c0_g2_i1_g.4284_m.4284   0.0306732697
# TRINITY_DN28600_c0_g1_TRINITY_DN28600_c0_g1_i1_g.3312_m.3312   0.0306732697

summary(res)
# 
# out of 1199 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 2, 0.17% 
# LFC < 0 (down)   : 10, 0.83% 
# outliers [1]     : 35, 2.9% 
# low counts [2]   : 394, 33% 
# (mean count < 9)


#################### MODEL NUMBER 3 - GROUP DESIGNS can be used for contrasts of interest or interactions

colData$group <- factor(paste0(colData$location, ".", colData$health, ".", colData$score))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
dim(dds)

dds <- dds[ rowSums(counts(dds)) > 100, ]

dds <- dds[sample(nrow(dds), 1200), ]
dim(dds)

dds <- DESeq(dds, parallel=T)


resultsNames(dds)
# [1] "Intercept"    "groupint.H.0" "groupint.S.1" "groupint.S.2" "groupint.S.3" "groupint.S.4"
# [7] "groupint.S.5" "groupsub.H.0" "groupsub.S.1" "groupsub.S.2" "groupsub.S.3"

res <- results(dds, contrast=list( c("groupint.H.0","groupsub.H.0"), c("groupint.S.1","groupsub.S.1")), listValues=c(1/2, -1/2))



res <- res[order(res$padj),]
head(res)
# log2 fold change (MAP): 0.5 groupint.H.0+groupsub.H.0 vs 0.5 groupint.S.1+groupsub.S.1 
# Wald test p-value: 0.5 groupint.H.0+groupsub.H.0 vs 0.5 groupint.S.1+groupsub.S.1 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE
# <numeric>      <numeric> <numeric>
#   TRINITY_DN44744_c1_g2_TRINITY_DN44744_c1_g2_i2_g.17686_m.17686  98.83363       3.120545 0.6788041
# TRINITY_DN40954_c0_g1_TRINITY_DN40954_c0_g1_i2_g.10488_m.10488  26.44690       4.329669 0.9761004
# TRINITY_DN38692_c0_g1_TRINITY_DN38692_c0_g1_i1_g.7973_m.7973    31.70968       3.750359 0.8660999
# TRINITY_DN41666_c0_g1_TRINITY_DN41666_c0_g1_i2_g.11542_m.11542  61.78884       3.491901 0.8303922
# TRINITY_DN44993_c2_g1_TRINITY_DN44993_c2_g1_i1_g.18240_m.18240  53.01705       2.565053 0.6175604
# TRINITY_DN43116_c2_g1_TRINITY_DN43116_c2_g1_i1_g.14220_m.14220  13.16513       4.620867 1.1169674
# stat         pvalue        padj
# <numeric>      <numeric>   <numeric>
#   TRINITY_DN44744_c1_g2_TRINITY_DN44744_c1_g2_i2_g.17686_m.17686  4.597122 0.000004283664 0.004506414
# TRINITY_DN40954_c0_g1_TRINITY_DN40954_c0_g1_i2_g.10488_m.10488  4.435680 0.000009178221 0.004827744
# TRINITY_DN38692_c0_g1_TRINITY_DN38692_c0_g1_i1_g.7973_m.7973    4.330169 0.000014899488 0.005224754
# TRINITY_DN41666_c0_g1_TRINITY_DN41666_c0_g1_i2_g.11542_m.11542  4.205122 0.000026094116 0.006170206
# TRINITY_DN44993_c2_g1_TRINITY_DN44993_c2_g1_i1_g.18240_m.18240  4.153526 0.000032739074 0.006170206
# TRINITY_DN43116_c2_g1_TRINITY_DN43116_c2_g1_i1_g.14220_m.14220  4.136976 0.000035191290 0.006170206


summary(res)
# 
# out of 1200 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 67, 5.6% 
# LFC < 0 (down)   : 4, 0.33% 
# outliers [1]     : 29, 2.4% 
# low counts [2]   : 119, 9.9% 
# (mean count < 4)


#################### MODEL NUMBER 4 - TIME SERIES USING REDUCED MODELS
## based on example on bioconductor: http://www.bioconductor.org/help/workflows/rnaseqGene/#running-the-differential-expression-pipeline

ddsTS <- DESeqDataSetFromMatrix(countData = countData, colData = colData, ~ health + day + health:day)
ddsTS <- ddsTS[ rowSums(counts(ddsTS)) > 100, ]

ddsTS <- ddsTS[sample(nrow(ddsTS), 1200), ]
dim(ddsTS)

ddsTS <- DESeq(ddsTS, parallel=T, test="LRT", reduced = ~ health + day)
resTS <- results(ddsTS)
resTS$symbol <- mcols(ddsTS)$symbol

head(resTS[order(resTS$padj),],4)
summary(resTS)


#################### MODEL NUMBER 5 - EFFECT OF DISEASE SEVERITY SCORE

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ score)
dim(dds)

dds <- dds[ rowSums(counts(dds)) > 100, ]

dds <- dds[sample(nrow(dds), 1200), ]
dim(dds)

dds <- DESeq(dds, parallel=T)

res <- results(dds)
res <- res[order(res$padj),]
head(res)

# log2 fold change (MAP): score 
# Wald test p-value: score 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange      lfcSE
# <numeric>      <numeric>  <numeric>
#   TRINITY_DN46245_c3_g3_TRINITY_DN46245_c3_g3_i2_g.21719_m.21719 149.70465     -0.4474123 0.09422044
# TRINITY_DN46674_c5_g1_TRINITY_DN46674_c5_g1_i1_g.23027_m.23027  88.85428      0.6247781 0.12822890
# TRINITY_DN45492_c0_g1_TRINITY_DN45492_c0_g1_i5_g.19587_m.19587 280.76979      0.5487730 0.12095387
# TRINITY_DN40724_c3_g1_TRINITY_DN40724_c3_g1_i1_g.10214_m.10214  21.31549     -0.6056648 0.13920150
# TRINITY_DN41587_c0_g1_TRINITY_DN41587_c0_g1_i1_g.11443_m.11443 143.23701     -0.5278008 0.12146880
# TRINITY_DN44140_c0_g1_TRINITY_DN44140_c0_g1_i1_g.16300_m.16300  76.60877      0.5192346 0.12590443
# stat         pvalue
# <numeric>      <numeric>
#   TRINITY_DN46245_c3_g3_TRINITY_DN46245_c3_g3_i2_g.21719_m.21719 -4.748569 0.000002048605
# TRINITY_DN46674_c5_g1_TRINITY_DN46674_c5_g1_i1_g.23027_m.23027  4.872366 0.000001102697
# TRINITY_DN45492_c0_g1_TRINITY_DN45492_c0_g1_i5_g.19587_m.19587  4.537043 0.000005704834
# TRINITY_DN40724_c3_g1_TRINITY_DN40724_c3_g1_i1_g.10214_m.10214 -4.350993 0.000013552219
# TRINITY_DN41587_c0_g1_TRINITY_DN41587_c0_g1_i1_g.11443_m.11443 -4.345155 0.000013917726
# TRINITY_DN44140_c0_g1_TRINITY_DN44140_c0_g1_i1_g.16300_m.16300  4.124037 0.000037228827
# padj
# <numeric>
#   TRINITY_DN46245_c3_g3_TRINITY_DN46245_c3_g3_i2_g.21719_m.21719 0.0005572205
# TRINITY_DN46674_c5_g1_TRINITY_DN46674_c5_g1_i1_g.23027_m.23027 0.0005572205
# TRINITY_DN45492_c0_g1_TRINITY_DN45492_c0_g1_i5_g.19587_m.19587 0.0010344766
# TRINITY_DN40724_c3_g1_TRINITY_DN40724_c3_g1_i1_g.10214_m.10214 0.0015142486
# TRINITY_DN41587_c0_g1_TRINITY_DN41587_c0_g1_i1_g.11443_m.11443 0.0015142486
# TRINITY_DN44140_c0_g1_TRINITY_DN44140_c0_g1_i1_g.16300_m.16300 0.0033754137

summary(res)
# 
# out of 1195 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 20, 1.7% 
# LFC < 0 (down)   : 24, 2% 
# outliers [1]     : 9, 0.75% 
# low counts [2]   : 647, 54% 
# (mean count < 17)

#################################################################################
# END OF EXAMPLE CODE:




























#################### MODEL NUMBER 6 thru INFINITY - YOUR CHOICE...





################################################################
# Data quality assessment by sample clustering and visualization 

plotMA(res, main="DESeq2", ylim=c(-2,2))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN44744_c1_g2_TRINITY_DN44744_c1_g2_i2_g.17686_m.17686", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p


## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(dds, gene="TRINITY_DN46245_c3_g3_TRINITY_DN46245_c3_g3_i2_g.21719_m.21719", intgroup=(c("health","score","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
p

############## PCA plots
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("score"))
plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))
