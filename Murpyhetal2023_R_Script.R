####Murphyetal2023 Methods Paper####

setwd("/Users/kaitlynmurphy/Desktop")
getwd()

#install these packages
install.packages('nlme')
install.packages('lme4')
install.packages('lsmeans')
install.packages('ggplot2')
install_github("vqv/ggbiplot")
install.packages("ggsignif")
install.packages('vegan')
install.packages('devtools')
install.packages('tidyverse')
install.packages('dplyr')

#call on these libraries
library(nlme)
library(lme4)
library(lsmeans)
library(ggplot2)
library(ggbiplot)
library(ggsignif)
library(vegan)
library(devtools)
library(tidyverse)
library(dplyr)

#####EggSonicationExp_Project1#####

#Call on your dataset (there are multiple used in this script. This one is entitled 'Exp1_2019-2020.csv')
datum_egg1=read.csv(file.choose())
datum_egg1[datum_egg1==""] <- NA
head(datum_egg1)

######Egg/Hatchling Data######

#Egg survival

resultsglm=glmer(SURV~TREAT+EMASS+(1|CAGE),data=datum_egg1,family=binomial,na.action=na.omit)
summary(resultsglm)
plot(SURV~TREAT,data=datum_egg1)

#Convert survival data to percentage
agg <- count(datum_egg1,TREAT,SURV)
head(agg)
xtb <- xtabs(n ~ SURV + TREAT, data = agg)
xtb
xtb1 <- as.data.frame(xtb)
xtb1

#Create a percentage dataframe for egg survival figure
TREAT <- c("Control", "Control (Sonication)", "Sonication", "Swab")
SURV <- c(90.48, 94.74, 27.78, 86.96)
survival <- data.frame(TREAT, SURV)

#Figure 1A

tiff("Figure1A.tiff", width = 4, height = 4, units = 'in', res = 300)

ggplot(survival, mapping=aes(x=TREAT, y= SURV, fill= TREAT)) + 
  geom_bar(stat="identity",width=0.75,position=position_dodge()) +
  geom_signif(comparisons = list(c("Sonication", "Control")),
              map_signif_level = TRUE) +
  scale_fill_manual(values=c("Sonication" = "grey65", "Swab" = "gray80", "Control" = "grey15", "Control(Son)" = "grey34")) +
  theme_classic() +
  theme(axis.title.x = element_text(size=15,vjust=0),
        axis.text.x  = element_text(size=11,color="black"),
        axis.title.y = element_text(size=15,vjust=1),
        axis.text.y  = element_text(size=11,color="black"), 
        legend.position="none") + 
  xlab("Treatment") + ylab("Percentage egg survival")

dev.off()

#Hatchling mass
resultshatchmass=lme(MASS~TREAT+EMASS,data=datum_egg1,random=~1|CAGE,na.action=na.omit)
summary(resultshatchmass) 

#Figure 1B
datum_egg1$TREAT_2 = factor(datum_egg1$TREAT,c("Control", "Control (Sonication)", "Sonication", "Swab"))

tiff("Figure1B.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(MASS~TREAT_2, data=datum_egg1,
        col=(c("grey15", "grey34", "grey65", "gray80")),
        xlab="Treatment", ylab="Hatchling mass (g)",
        outline=FALSE)
stripchart(MASS~TREAT_2, vertical = TRUE, data = datum_egg1, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#SVL
resultssvl=lme(SVL~TREAT+EMASS,data=datum_egg1,random=~1|CAGE,na.action=na.omit)
summary(resultssvl)

tiff("Figure1C.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(SVL~TREAT_2, data=datum_egg1,
        col=(c("grey15", "grey34", "grey65", "gray80")),
        xlab="Treatment", ylab="Hatchling SVL (mm)",
        outline=FALSE)
stripchart(SVL~TREAT_2, vertical = TRUE, data = datum_egg1, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#TL
resultstl=lme(TL~TREAT+EMASS,data=datum_egg1,random=~1|CAGE,na.action=na.omit)
summary(resultstl)

#Sex
resultssex=glmer(SEX~TREAT+EMASS+(1|CAGE),data=datum_egg1,family=binomial,na.action=na.omit)
summary(resultssex)

######DNA Quantification######

results=lme(QUANTITY~TREAT,random=~1|CAGE,data=datum_egg1,na.action=na.omit)
summary(results)

datum_egg2 <- datum_egg[1:8,]

#Figure not used in manuscript
tiff("Figure1.tiff", width = 4, height = 4, units = 'in', res = 300)

plot(QUANTITY~TREAT,data=datum_egg, main="DNA Yield from Brown Anole Eggs",
     xlab="Sampling Method", ylab="DNA Quantity (ng/uL)", col=(c("darkseagreen4","darkorange2")))

dev.off()

#######Rarefaction Curve######

#Dataset with taxonomy (there are multiple used in this script. This one is entitled 'Exp1_2019_Shannons.csv')
datum=read.csv(file.choose())
head(datum)

#https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/

BCI2 <- datum[1:8, ]
raremax <- min(rowSums(BCI2))
raremax

#raremax is the minimum sample count achieved over the 26 samples. 
#We will rarefy the sample counts to this value.

col <- c("black", "grey10", "grey22", "grey34", "grey40", "grey65", "gray80", "gray90")
lty <- c("solid", "dashed", "longdash", "dotdash","solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)

#Note that I save the output from rarecurve() in object out. This object contains everything we need to draw our own version of the plot if we wish. 
#For example, we could use fewer colours and alter the line thickness1 instead to make up the required number of combinations.

tiff("rarefaction.tiff", width = 8, height = 4, units = 'in', res = 300)

out <- with(pars[1:9, ],
            rarecurve(BCI2, step = 9, col = col,
                      lty = lty, label = TRUE, abline(v = raremax)))
dev.off()

#checking for bias among treatment groups

######OTU Abund######

resultsotu=lme(OTU~TREAT,random=~1|CAGE,data=datum_egg1,na.action=na.omit)
summary(resultsotu)

#Figure 2A

#subset dataframe to only those sequenced
datum_fig2 <- datum_egg1 %>% slice(41:44, 59:62)
datum_fig2[datum_fig2==""] <- NA

#Splicing did not work, so manually went back into dataset and subsetted to only samples seqeunced
#Call on your dataset (there are multiple used in this script. This one is entitled 'Exp1_2019-2020.csv')
datum_fig2=read.csv(file.choose())
datum_egg1[datum_egg1==""] <- NA
head(datum_egg1)

#make figure
tiff("Fig_2A.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(OTU ~ TREAT,data=datum_fig2, notch=FALSE, 
        col=(c("grey80", "gray34")),
        xlab="Sampling Method", ylab="Operational taxonomic units (OTUs)")
stripchart(OTU ~ TREAT, vertical = TRUE, data = datum_fig2, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

######Shannon's Diversity######

#http://rstudio-pubs-static.s3.amazonaws.com/11473_29c2de401c55441ca5350862ebd04456.html

#Calculate Shannons diversity using the dataset subsetted for Figure 2A

##Load OTU Table in (OTUs as columns, samples as rows)
names(datum)
#for the number of species
ncol(datum)
#for the number of samples
nrow(datum)

#shannons diversity
shann <- diversity(datum)
shann

#Shannon's diversity was added to the main dataframe 

resultssh=lme(SHANNON~TREAT,random=~1|CAGE,data=datum_egg1,na.action=na.omit)
summary(resultssh)

#Figure 2B

tiff("Fig_2B.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(SHANNON ~ TREAT,data=datum_fig2, notch=FALSE, 
        col=(c("grey80", "gray34")),
        xlab="Sampling Method", ylab="Shannon's Diversity Index")
stripchart(SHANNON ~ TREAT, vertical = TRUE, data = datum_fig2, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#Faith's phylogenetic diversity

resultsfaith=lme(FAITH~TREAT,random=~1|CAGE,data=datum_egg1,na.action=na.omit)
summary(resultsfaith)

######NMDS Plot######

#https://chrischizinski.github.io/rstats/vegan-ggplot2/

#Load the metadata file
metadata=read.csv(file.choose())

#Call on the dataset (there are multiple used in this script. This one is entitled 'Exp1_2019_Shannons.csv')
datum_diversity=read.csv(file.choose())

names(datum_diversity)
#for the number of species
ncol(datum_diversity)
#for the number of samples
nrow(datum_diversity)
nrow(metadata)

#Change NA to 0
datum_diversity[is.na(datum_diversity)] <- 0

#make sure samples are in same order as metadata
head(metadata)
head(datum_diversity)

#Transpose it (Sample = rows, taxnomic=columns)
#vt_MDS= as.data.frame(t(sheet2))
#otherwise run:
Vt_MDS = datum_diversity
head (Vt_MDS)

##Actually run the nMDS, default is Bray, but you can change that (see manpage)
###########IF THERE IS A PROBLEM AT THIS STEP ABOUT NUMERIC DATA, go back to your original final and make labels in your taxonomy file numbers###
#V.nMDSAG5 <- metaMDS(Vt_MDS, k = 2, trymax = 100, trace = F)- the default (bray-curtis) is the same as this
V.nMDSAG5 <- metaMDS(Vt_MDS)
head(V.nMDSAG5)
summary(V.nMDSAG5)
stressplot(V.nMDSAG5)

##Let's put these points together with the metadata sheet so we can graph this
V.nMDSnumbersAG5=as.data.frame(V.nMDSAG5$points)
head(V.nMDSnumbersAG5)
V.nMDSnummetaAG5=cbind(V.nMDSnumbersAG5, Treatment = metadata$subject)
#V.nMDSnummetaAG5=cbind(V.nMDSnumbersAG5,metadata$subject)
head(V.nMDSnummetaAG5)

#################################plot by Method
tiff("Figure2C.tiff", width = 8, height = 4, units = 'in', res = 300)
    
ggplot(V.nMDSnummetaAG5, aes(x=MDS1, y=MDS2, color= factor(metadata$subject))) +
  coord_equal() +
  geom_polygon(data=V.nMDSnummetaAG5,aes(x=MDS1,y=MDS2,fill= factor(metadata$subject),group= factor(metadata$subject)),alpha=0.30,show.legend = NA) + # add the convex hulls
  scale_fill_manual(values=c("grey34", "black")) +
  geom_point(aes(shape = factor(metadata$subject)), size = 3.5) +
  scale_colour_manual(values=c("grey34", "black")) +
  theme_bw() +
  labs(title = "nMDS(Bray-Curtis) OTUs") 
         
dev.off()

#ANOSIM
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

m_com = as.matrix(datum_diversity)
pc=cbind(datum_diversity, Treatment = metadata$subject)
head(pc)

ano = anosim(m_com, pc$Treatment, distance = "bray", permutations = 9999)
ano

######Phyla and Genera######

#Dataset with taxonomy (there are multiple used in this script. This one is entitled 'Exp1_2019_Phyla.csv')
phyla=read.csv(file.choose())
head(phyla)

phylalm=lm(ABUND ~ TREAT, data = phyla, na.action=na.omit)
summary(phylalm)

#PieCharts

#https://www.displayr.com/how-to-make-a-pie-chart-in-r/

#subset dataframe to only those sequenced from swabbing
phyla_swab <- phyla %>% slice(81:160)
phyla_swab[phyla_swab==""] <- NA

#Figure 2D Swabbing
tiff("Swab.tiff", width = 16, height = 9, units = 'in', res = 300)

ggplot(phyla_swab, aes(x="", y=ABUND, fill=PHYLA)) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_grey(start = 0.9, end = 0) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Swabbing OTUs") +
  theme_void() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "black", size= 28),
        legend.position="right",
        legend.text = element_text(size = 7)) 

dev.off()

#subset dataframe to only those sequenced from sonication
phyla_son <- phyla %>% slice(1:80)
phyla_son[phyla_son==""] <- NA

#Figure 2D Sonication
tiff("Son.tiff", width = 16, height = 9, units = 'in', res = 300)

ggplot(phyla_son, aes(x="", y=ABUND, fill=PHYLA)) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_grey(start = 0.9, end = 0) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Sonication OTUs") +
  theme_void() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "black", size = 28),
        legend.position="right",
        legend.text = element_text(size = 7)) 

dev.off()

#####Egg_Sterilization_Project2#####

#Call on your dataset (there are multiple used in this script. This one is entitled 'Exp2_2020.csv')

datum_egg2=read.csv(file.choose())
datum_egg2[datum_egg2==""] <- NA
head(datum_egg2)

######Egg/Hatchling Data 2020######

#Egg survival
resultsglm=glmer(SURV~relevel(TREAT,ref = "20% Bleach")+EMASS+(1|CAGE),data=datum_egg2,family=binomial,na.action=na.omit)
summary(resultsglm)
plot(SURV~TREAT,data=datum_egg2)

#Convert survival data to percentage
agg1 <- count(datum_egg2,TREAT,SURV)
head(agg1)
xtb2 <- xtabs(n ~ SURV + TREAT, data = agg1)
xtb2
xtb3 <- as.data.frame(xtb2)
xtb3

#Create a percentage dataframe for egg survival figure
TREAT <- c("Control", "Control (Water)", "10% Bleach", "20% Bleach", "Ethanol")
SURV <- c(54.55, 72.73, 69.23, 75, 18.18)
survival1 <- data.frame(TREAT, SURV)

#Figure 3A

survival1$TREAT1 = factor(survival1$TREAT,c("Control", "Control (Water)", "10% Bleach", "20% Bleach", "Ethanol"))

tiff("Figure3A.tiff", width = 4, height = 4, units = 'in', res = 300)

ggplot(survival1, mapping=aes(x=TREAT1, y= SURV, fill= TREAT)) + 
  geom_bar(stat="identity",width=0.75,position=position_dodge()) +
  geom_signif(comparisons = list(c("20% Bleach", "Ethanol")),
              map_signif_level = TRUE) +
  scale_fill_manual(values=c("10% Bleach" = "grey65", "20% Bleach" = "gray80", "Ethanol" = "gray90", "Control" = "grey15", "Control(Water)" = "grey34")) +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25)) +
  theme_classic() +
  theme(axis.title.x = element_text(size=15,vjust=0),
        axis.text.x  = element_text(size=8,color="black"),
        axis.title.y = element_text(size=15,vjust=1),
        axis.text.y  = element_text(size=8,color="black"), 
        legend.position="none") + 
  xlab("Treatment") + ylab("Percentage egg survival")

dev.off()

#Hatchling mass
resultshatch=lme(MASS~relevel(TREAT,ref = "Control")+EMASS,data=datum_egg2,random=~1|CAGE,na.action=na.omit)
summary(resultshatch)
plot(MASS~TREAT,data=datum_egg2)
#non-significant interaction term

#Figure 3B
datum_egg2$TREAT_2 = factor(datum_egg2$TREAT,c("Control", "Control (Water)", "10% Bleach", "20% Bleach", "Ethanol"))

tiff("Figure3B.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(MASS~TREAT_2, data=datum_egg2,
        col=(c("grey15", "grey34", "grey65", "gray80", "grey90")),
        xlab="Treatment", ylab="Hatchling mass (g)",
        outline=FALSE)
stripchart(MASS~TREAT_2, vertical = TRUE, data = datum_egg2, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#SVL
resultssvl=lme(SVL~relevel(TREAT,ref = "Control")+EMASS,data=datum_egg2,random=~1|CAGE,na.action=na.omit)
summary(resultssvl)
#non-significant interaction term

#Figure 3C
tiff("Figure3C.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(SVL~TREAT_2, data=datum_egg2,
        col=(c("grey15", "grey34", "grey65", "gray80", "grey90")),
        xlab="Treatment", ylab="Hatchling SVL (mm)",
        outline=FALSE)
stripchart(SVL~TREAT_2, vertical = TRUE, data = datum_egg2, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#TL
resultstl=lme(TL~relevel(TREAT,ref = "Control")+EMASS,data=datum_egg2,random=~1|CAGE,na.action=na.omit)
summary(resultstl)
#non-significant interaction term

#Sex
resultssex=glmer(SEX~TREAT+EMASS+(1|CAGE),data=datum_egg2,family=binomial,na.action=na.omit)
summary(resultssex)
#significant p-value for treat

######Egg/Hatchling Data 2021-2022######

#Call on your dataset (there are multiple used in this script. This one is entitled 'Exp2_2021.csv')

datum_egg3=read.csv(file.choose())
datum_egg3[datum_egg3==""] <- NA
head(datum_egg3)

#Egg survival

resultsglm=glmer(SURV~relevel(TREAT,ref = "Control")+(1|CAGE),data=datum_egg3,family=binomial,na.action=na.omit)
summary(resultsglm)
plot(SURV~TREAT,data=datum_egg3)

#Convert survival data to percentage
agg2 <- count(datum_egg3,TREAT,SURV)
head(agg2)
xtb3 <- xtabs(n ~ SURV + TREAT, data = agg2)
xtb3
xtb4 <- as.data.frame(xtb3)
xtb4

#Create a percentage dataframe for egg survival figure
TREAT <- c("Control", "Control (Water)", "Bleach", "TSB")
SURV <- c(75, 92.31, 100, 94.12)
survival2 <- data.frame(TREAT, SURV)

#Figure 4A

survival2$TREAT1 = factor(survival2$TREAT,c("Control", "Control (Water)", "Bleach", "TSB"))

tiff("Figure4A.tiff", width = 4, height = 4, units = 'in', res = 300)

ggplot(survival2, mapping=aes(x=TREAT1, y= SURV, fill= TREAT1)) + 
  geom_bar(stat="identity",width=0.75,position=position_dodge()) +
  geom_signif(comparisons = list(c("TSB", "Control")),
              map_signif_level = TRUE) +
  scale_fill_manual(values=c("20% Bleach" = "gray80", "TSB" = "gray90", "Control" = "grey15", "Control(Water)" = "grey34")) +
  theme_classic() +
  theme(axis.title.x = element_text(size=15,vjust=0),
        axis.text.x  = element_text(size=8,color="black"),
        axis.title.y = element_text(size=15,vjust=1),
        axis.text.y  = element_text(size=8,color="black"), 
        legend.position="none") + 
  xlab("Treatment") + ylab("Percentage egg survival")

dev.off()

#Hatchling data

#Mass
resultshatch=lme(MASS~relevel(TREAT,ref = "Control")+EMASS,data=datum_egg3,random=~1|CAGE,na.action=na.omit)
summary(resultshatch)
plot(MASS~TREAT,data=datum_egg3)

#Figure 4B
datum_egg3$TREAT_2 = factor(datum_egg3$TREAT,c("Control", "Control (Water)", "Bleach", "TSB"))

tiff("Figure4B.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(MASS~TREAT_2, data=datum_egg3,
        col=(c("grey15", "grey34", "grey65", "gray80")),
        xlab="Treatment", ylab="Hatchling mass (g)",
        outline=FALSE)
stripchart(MASS~TREAT_2, vertical = TRUE, data = datum_egg3, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#SVL
resultssvl=lme(SVL~relevel(TREAT,ref = "Control")+EMASS,data=datum_egg3,random=~1|CAGE,na.action=na.omit)
summary(resultssvl)

#Figure 4C
tiff("Figure4C.tiff", width = 4, height = 4, units = 'in', res = 300)

boxplot(SVL~TREAT_2, data=datum_egg3,
        col=(c("grey15", "grey34", "grey65", "gray80")),
        xlab="Treatment", ylab="Hatchling SVL (mm)",
        outline=FALSE)
stripchart(SVL~TREAT_2, vertical = TRUE, data = datum_egg3, 
           method = "jitter", add=TRUE, pch = 20)

dev.off()

#TL
resultstl=lme(TL~relevel(TREAT,ref = "Control")+EMASS,data=datum_egg3,random=~1|CAGE,na.action=na.omit)
summary(resultstl)

#Sex
resultssex=glmer(SEX~relevel(TREAT,ref = "Control")+EMASS+(1|CAGE),data=datum_egg3,family=binomial,na.action=na.omit)
summary(resultssex)

######CFU Data 2022######

#Call on your dataset (there are multiple used in this script. This one is entitled 'Exp2_2022.csv')

datum_egg4=read.csv(file.choose())
datum_egg4[datum_egg4==""] <- NA
head(datum_egg4)

#Check to see if there was an effect of swabbing on CFU counts
#subset dataframe to only those from group A (after swabbing, before treatment)
df1 <- subset(datum_egg4,TIME=='A' )
df1[df1==""] <- NA

resultscfu=lme(COUNT~relevel(TREAT,ref = "Control"),data=df1, random=~1|CAGE,na.action=na.omit)
summary(resultscfu)

#subset dataframe to only those from group A (after swabbing)
df2 <- subset(datum_egg4,TIME=='B' )
df2[df2==""] <- NA

resultscfu1=lme(COUNT~relevel(TREAT,ref = "Control"),data=df2, random=~1|CAGE,na.action=na.omit)
summary(resultscfu1)

aggregate(df2$COUNT, list(df2$TREAT1), FUN=mean)

#Check to see if there was an effect of time between individual treatment groups
#Control (no manipulation)
df3 <- subset(datum_egg4,TREAT=='Control' )
df3 <- subset(df3, COUNT != "")

resultscfu2=lme(COUNT~TIME,data=df3, random=~1|CAGE,na.action=na.omit)
summary(resultscfu2)

aggregate(df3$COUNT, list(df3$TREAT1), FUN=mean)

#Control (water)
df4 <- subset(datum_egg4,TREAT=='Control (Water)' )
df4 <- subset(df4, COUNT != "")

resultscfu3=lme(COUNT~TIME,data=df4, random=~1|CAGE,na.action=na.omit)
summary(resultscfu3)

aggregate(df4$COUNT, list(df4$TREAT1), FUN=mean)

#Bleach
df5 <- subset(datum_egg4,TREAT=='Bleach' )
df5 <- subset(df5, COUNT != "")

resultscfu4=lme(COUNT~TIME,data=df5, random=~1|CAGE,na.action=na.omit)
summary(resultscfu4)

aggregate(df5$COUNT, list(df5$TREAT1), FUN=mean)

#TSB
df6 <- subset(datum_egg4,TREAT=='TSB' )
df6 <- subset(df6, COUNT != "")

resultscfu5=lme(COUNT~TIME,data=df6, random=~1|CAGE,na.action=na.omit)
summary(resultscfu5)

aggregate(df6$COUNT, list(df6$TREAT1), FUN=mean)

#CFUs before treatment as a covariate
resultscfu6=lme(COUNT~relevel(TREAT,ref = "Control")+TIME,data=datum_egg4,random=~1|CAGE,na.action=na.omit)
summary(resultscfu6)

#Figure 5

datum_egg4$TREAT_2 = factor(datum_egg4$TREAT,c("Control", "Control (Water)", "Bleach", "TSB"))

tiff("Figure_5.tiff", width = 8, height = 8, units = 'in', res = 300)

ggplot(datum_egg4, aes(x=TREAT_2, y=COUNT, fill=TIME)) + 
  geom_boxplot(outlier.shape = NA, notch=FALSE, width= 0.3, position = position_dodge (width = 0.4)) +
  scale_fill_grey(start = 0.9, end = 0.6) +
  labs(x = "Eggshell treatment", y ="Colony forming units (CFUs)") +
  theme_classic() +
  scale_y_continuous(limits = c(0,750)) +
  theme(axis.title.x = element_text(size=15,vjust=0),
        axis.text.x  = element_text(size=12,color="black", angle = 0, vjust = 0.5, hjust=0.4),
        axis.title.y = element_text(size=15,vjust=1),
        axis.text.y  = element_text(size=12,color="black"))

dev.off()

#Incorporating TNTC data 

#Call on your dataset (there are multiple used in this script. This one is entitled 'Supp_Fig5.csv')
datum_egg5=read.csv(file.choose())
datum_egg5[datum_egg5==""] <- NA
head(datum_egg5)

#Check for significance among treatment groups
resultscfu1=lme(CHANGE~relevel(TREAT,ref = "Control"),data=datum_egg5, random=~1|CAGE,na.action=na.omit)
summary(resultscfu1)
plot(CHANGE~TREAT, data=datum_egg5)
plot(TIME~TREAT, data=datum_egg5)

#Check to see if there was an effect of time between individual treatment groups
#Control (no manipulation)
df7 <- subset(datum_egg5,TREAT=='Control' )
df7 <- subset(df7, COUNT != "")

resultscfu2=lme(COUNT~TIME,data=df7, random=~1|CAGE,na.action=na.omit)
summary(resultscfu2)

aggregate(df7$COUNT, list(df7$TREAT1), FUN=mean)

#Control (water)
df8 <- subset(datum_egg5,TREAT=='Control (Water)' )
df8 <- subset(df8, COUNT != "")

resultscfu3=lme(COUNT~TIME,data=df8, random=~1|CAGE,na.action=na.omit)
summary(resultscfu3)

aggregate(df8$COUNT, list(df8$TREAT1), FUN=mean)

#Bleach
df9 <- subset(datum_egg5,TREAT=='Bleach' )
df9 <- subset(df5, COUNT != "")

resultscfu4=lme(COUNT~TIME,data=df9, random=~1|CAGE,na.action=na.omit)
summary(resultscfu4)

aggregate(df9$COUNT, list(df9$TREAT1), FUN=mean)

#TSB
df10 <- subset(datum_egg5,TREAT=='TSB' )
df10 <- subset(df10, COUNT != "")

resultscfu5=lme(COUNT~TIME,data=df10, random=~1|CAGE,na.action=na.omit)
summary(resultscfu5)

aggregate(df10$COUNT, list(df10$TREAT1), FUN=mean)

#Supplementary Figure 2

datum_egg5$TREAT_2 = factor(datum_egg5$TREAT,c("Control", "Control (Water)", "Bleach", "TSB"))

tiff("Supp_Figure_2.tiff", width = 8, height = 8, units = 'in', res = 300)

ggplot(datum_egg5, aes(x=TREAT_2, y=COUNT, fill=TIME)) + 
  geom_boxplot(outlier.shape = NA, notch=FALSE, width= 0.3, position = position_dodge (width = 0.4)) +
  scale_fill_grey(start = 0.9, end = 0.6) +
  labs(x = "Eggshell treatment", y ="Colony forming units (CFUs)") +
  theme_classic() +
  scale_y_continuous(limits = c(0,1038)) +
  theme(axis.title.x = element_text(size=15,vjust=0),
        axis.text.x  = element_text(size=12,color="black", angle = 0, vjust = 0.5, hjust=0.4),
        axis.title.y = element_text(size=15,vjust=1),
        axis.text.y  = element_text(size=12,color="black"))

dev.off()
