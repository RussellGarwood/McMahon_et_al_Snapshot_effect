## If packages aren't installed, install them, then load them
packages <- c("ggplot2", "Rcpp", "gridExtra", "grid","viridis", "scico")
if(length(packages[!packages %in% installed.packages()[,"Package"]]) > 0){
  install.packages(packages[!packages %in% installed.packages()[,"Package"]])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))
rm("packages")

################################################################################################################################################ 
# This version of the model does replicates by calling the C++ code in the same folder - this is good for doing replicates
################################################################################################################################################

#Clear environment  - if required
#rm(list = ls())

#This needs to be set to load the code if not running from a bash script
#setwd("~/Desktop/Programming/scripts/Simulation/McMahon_model/") #Laptop and Workstation
#outputDirectory<-"~/Desktop/"

#If running from a bash script, wd is automatically where that script is stored - which is fine
#If running from bash, so set output directory to working directory
outputDirectory<-paste(getwd(),"/results/",sep="")

################################################################################################################################################ 
# Source functions
################################################################################################################################################

#Load CPP - in the same folder as this script 
sourceCpp("McMahon_model_v9.cpp")
#All R functions are in their own folder
source("functions/graphSingleRun.R")

############################################ Single run ############################################
## Here we can graph a single run and look at the values through time plotted

#Parameters
birthChance <- 0.0020;
decayRate<-0.005;
individuals <- 100000;
repeats <- 1;
startAge <- 500; #Only used in simulation with age 
runFor <- 10000;

#Run simulation
taphonomyDFSingleRun<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,runFor)
## Now do graphing for this run
graphSingleRun(taphonomyDFSingleRun, paste("Single run - startAge ",startAge,", decayRate ",decayRate, sep=""),outputDirectory)

#Or benchmark
#benchmarks<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge, runFor, FALSE, FALSE, TRUE)
#benchmarks<-benchmarks/1000000
#diff(benchmarks)

############################################ Single run with random decay rates ############################################
## Here we can do a single run, but assign the individuals in that run their own, random, decay rate

#Parameters
birthChance <- 0.0020;
decayRate<-0.02; # This is max from decay graph
individuals <- 50000;
repeats <- 1;
startAge <- 500;
runFor <- 20000;

taphonomyDFRandomDecayRates<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge, runFor, FALSE, TRUE)
taphonomyDFRandomDecayRates$decayRateRounded<-round(taphonomyDFRandomDecayRates$decayRates, digits = 3)
write.csv(taphonomyDFRandomDecayRates, paste(outputDirectory,"taphonomyDFRandomDecayRates.csv",sep=""))

## Now do graphing for this run
ggplot(data = taphonomyDFRandomDecayRates) + geom_boxplot(mapping = aes(factor(.data[["decayRateRounded"]]), .data[["decayLevels"]])) +
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA)) +
  labs(title=paste("Taphonomy Simulator - Randomised decay rate vs preservation"), x="Individual decay rate", y="Decay level") 
ggsave(filename <- paste(outputDirectory,"Decay_rate_single_run_boxplot.pdf",sep=""))

##This code is useful for anaylsing the data - no longer used for graphing, but I have left it in
rateVector<-vector()
meanVector<-vector()
aliveVector<-vector()
for(rate in unique(taphonomyDFRandomDecayRates$decayRateRounded))
{
  rateVector<-append(rateVector,rate)
  meanVector<-append(meanVector,mean(taphonomyDFRandomDecayRates[taphonomyDFRandomDecayRates$decayRateRounded==rate,]$decayLevels))
  aliveVector<-append(aliveVector,sum(taphonomyDFRandomDecayRates$states  == 1 & taphonomyDFRandomDecayRates$decayRateRounded == rate))
}
rateMeanDataframe <- data.frame(rateVector,meanVector,aliveVector)
write.csv(rateMeanDataframe, paste(outputDirectory,"randomRateMeanDataframe.csv",sep=""))

############################################ Single run with random start ages ############################################
## Here we can do a single run, but assign the individuals in that run their own, random, decay rate

#Parameters
birthChance <- 0.00015;#This equates to about ~60k living at any time
decayRate<-0.001;#A decay rate of 100 days
individuals <- 100000;
repeats <- 1;
startAge <- 20000; #This is maximum from start age replicates
runFor <- 50000;#This is maximum from start age replicates graph

taphonomyDFRandomStartAges<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge, runFor, TRUE, FALSE)
taphonomyDFRandomStartAges$startAgeRounded<-with(taphonomyDFRandomStartAges, ifelse(startAges %in% c(0:9999),signif(taphonomyDFRandomStartAges$startAges, digits = 1), signif(taphonomyDFRandomStartAges$startAges, digits = 2)))
taphonomyDFRandomStartAges$startAgeRounded<-with(taphonomyDFRandomStartAges, ifelse(startAgeRounded %in% c(0:999),0, startAgeRounded))
write.csv(taphonomyDFRandomStartAges, paste(outputDirectory,"taphonomyDFRandomStartAges.csv",sep=""))

ageVector<-vector()
meanVector<-vector()

for(age in unique(taphonomyDFRandomStartAges$startAgeRounded))
{
  ageVector<-append(ageVector,age)
  meanVector<-append(meanVector,mean(taphonomyDFRandomStartAges[taphonomyDFRandomStartAges$startAgeRounded==age,]$decayLevels))
}
ageMeanDataframe <- data.frame(ageVector,meanVector)
write.csv(rateMeanDataframe, paste(outputDirectory,"ageMeanDataframe.csv",sep=""))

ggplot(data = taphonomyDFRandomStartAges) + geom_violin(mapping = aes(factor(.data[["startAgeRounded"]]), .data[["decayLevels"]])) +
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA)) +
  labs(title=paste("Taphonomy Simulator - Randomised life span vs preservation"), x="Individual life span", y="Decay level") 
ggsave(filename <- paste(outputDirectory,"Start_age_single_run_boxplot.pdf",sep=""))

ggplot(data = ageMeanDataframe) + geom_col(mapping = aes(ageVector,meanVector)) +
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA)) +
  labs(title=paste("Taphonomy Simulator - Randomised lifespan vs mean decay"), x="Individual life span (rounded)", y="Mean decay value")
ggsave(filename <- paste(outputDirectory,"Start_age_single_run_colplot.pdf",sep=""))

############################################ Decay rate graph using replicates ############################################
birthChance <- 0.0020;
individuals <- 5000;
repeats <- 500;
startAge <-500;

#Let's do decay rates for 2 month decay intervals between almost immediate decay (5 days) and 5 years decay time
#Max 0.02 (=0.2 per day = 5 day decay time)
#To 0.0005479452 (=1/1825 per day = 5 year decay time then divide by ten for actual rate)
decayRateV<-vector()

for(decayRate in seq(from = 0.00005479452, to = 0.02, by = (0.02-0.00005479452)/20))
{
  decayRateV<-c(decayRateV,decayRate)
}

#This will allow us to store our data
decayLevelDF<-data.frame("taphonomyValues" = 0, "decayRate" = 0);

for (decayRate in decayRateV) 
{
  print(paste("Doing decayRate", decayRate))
  if(decayRate<0.0015)taphonomyValues<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,40000)
  else taphonomyValues<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,20000)
  #To check a single run, run for 1 repeat, then as follows
  #graphSingleRun(taphonomyValues, "")
  #Otherwise for multiple repeats
  taphonomyDF<-as.data.frame(taphonomyValues)
  taphonomyDF$decayRate<-decayRate
  decayLevelDF <- (rbind (decayLevelDF, taphonomyDF))
}
decayLevelDF<-decayLevelDF[-c(1), ]
decayLevelDF$PerDay<-decayLevelDF$decayRate*10
write.csv(decayLevelDF, paste(outputDirectory,"decayLevelDF.csv",sep=""))

ggplot(data = decayLevelDF, aes(x=as.factor(PerDay), y=taphonomyValues)) + geom_boxplot() + #ylim(0.6,0.91) +
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title=paste("Taphonomy Simulator - decay rate vs preservation"), x="Decay rate (per day)", y="Average non-zero value")
ggsave(filename <- paste(outputDirectory,"Decay_rate_replicates_boxplot.pdf",sep=""))

ggplot(data = decayLevelDF, aes(x=PerDay, y=taphonomyValues)) + geom_point() + #ylim(0.6,0.91) +
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title=paste("Taphonomy Simulator - decay rate vs preservation"), x="Decay rate (per day)", y="Average non-zero value")
ggsave(filename <- paste(outputDirectory,"Decay_rate_replicates_point.pdf",sep=""))

############################################ Life span graph ############################################
birthChance <- 0.000015;#This has to be very small so as not to fill the list with long life spans
decayRate<-0.001;#A decay rate of 100 days
individuals <- 5000;
repeats <- 500;

#Lifespans of 10 days to 2000 days (5.5 years) for a decay rate of 100 days
startAgeV<-vector()
for(startAge in seq(from = 100, to = 20000, by = 500))
{
  startAgeV<-c(startAgeV,startAge)
}

decayLevelDFLifeSpan <- data.frame("taphonomyValues" = 0, "startAge" = 0);

for (startAge in startAgeV) 
{
  print(paste("Doing startAge", startAge))
  taphonomyValues<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,50000)
  #To check a single run, run for 1 repeat, then as follows
  #graphSingleRun(taphonomyValues, "")
  #Otherwise for multiple repeats
  taphonomyDF<-as.data.frame(taphonomyValues)
  taphonomyDF$startAge<-startAge
  decayLevelDFLifeSpan<-rbind(decayLevelDFLifeSpan, taphonomyDF)
}
decayLevelDFLifeSpan<-decayLevelDFLifeSpan[-c(1), ]
decayLevelDFLifeSpan$days<-decayLevelDFLifeSpan$startAge/10
write.csv(decayLevelDFLifeSpan,paste(outputDirectory,"decayLevelDFLifeSpan.csv", sep=""))

ggplot(data = decayLevelDFLifeSpan, aes(x=as.factor(days), y=taphonomyValues)) + geom_boxplot() + theme_minimal() + 
  theme(panel.border = element_rect(color="black", fill=NA), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title=paste("Taphonomy Simulator - lifespan vs preservation"), x="Lifespan (days)", y="Average non-zero value")
ggsave(filename <- paste(outputDirectory,"Life_span_boxplot.pdf",sep=""))

ggplot(data = decayLevelDFLifeSpan, aes(x=days, y=taphonomyValues)) + geom_jitter() + theme_minimal() + 
  theme(panel.border = element_rect(color="black", fill=NA), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title=paste("Taphonomy Simulator - lifespan vs preservation"), x="Lifespan (days)", y="Average non-zero value")
ggsave(filename <- paste(outputDirectory,"Life_span_points.pdf",sep=""))

############################################ Heat map ############################################
birthChance <- 0.000015;#This has to be very small so as not to fill the list with long life spans
individuals <- 10000;
repeats <- 100;

#Lifespans of 10 days to 2000 days (5.5 years), this time with log spacing
startAgeV2<-c(seq(from = 1, to = 9.5, by = 0.5) %o% 10^(2:4))
startAgeV2<-startAgeV2[startAgeV2<21000]

#Let's do decay rates for 2 month decay intervals between almost immediate decay (5 days) and 5 years decay time
#Max 0.02 (=0.2 per day = 5 day decay time)
#To 0.00005479452 (=1/1825 per day = 5 year decay time /10 to make it into dats)
#This time with log spacing
decayRateV2<-c(seq(from = 1, to = 9.5, by = 0.5) %o% 10^(-5:-2))
decayRateV2<-decayRateV2[!(decayRateV2<0.00005479452 | decayRateV2>0.02)]

decayLevelDFHeatMap <- data.frame("taphonomyValues" = 0, "decayRate" = 0, "startAge"=0);
for (decayRate in decayRateV2)
{
  for (startAge in startAgeV2) 
  {
    print(paste("Doing decayRate", decayRate, "and startAge", startAge))
    outputDirectoryLocal<-chartr(".","p",paste(outputDirectory,"DR_",decayRate,"_SA_",startAge,"/",sep=""))
    dir.create(outputDirectoryLocal)
    
    #Do single repeat to provide graphs and stats on this run
    repeats <- 1;    
    if(decayRate<0.0015 | startAge > 15000 )taphonomyValues<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,50000)
    else taphonomyValues<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,25000)
    graphSingleRun(taphonomyValues, paste("Single run - startAge ",startAge,", decayRate ",decayRate), outputDirectoryLocal)
    
    #Now do multiple repeats
    repeats <- 50;
    if(decayRate<0.0015 | startAge > 10000 )taphonomyValues<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,50000)
    else taphonomyValues<-doSimulationAge(individuals,birthChance,decayRate,repeats,startAge,25000)
    taphonomyDF<-as.data.frame(taphonomyValues)
    taphonomyDF$decayRate<-decayRate
    taphonomyDF$startAge<-startAge
    write.csv(taphonomyDF,paste(outputDirectoryLocal,"taphonomyDF.csv"))
    decayLevelDFHeatMap<-rbind(decayLevelDFHeatMap,taphonomyDF)
  }
}
decayLevelDFHeatMap <- decayLevelDFHeatMap[-c(1), ]
write.csv(decayLevelDFHeatMap,paste(outputDirectory,"decayLevelDFHeatMap.csv"))

##Now do some plots 
#Create labels
myLabels <- vector();
#Do it with integers as mod with decimals is problematic in R
for (lab in seq(from = 50, to = 600, by = 5)) 
{
  if(lab%%100==0) myLabels<-append(myLabels,as.character(lab/1000))
  else myLabels<-append(myLabels,"")
}

decayLevelDFHeatMapDays<-decayLevelDFHeatMap
decayLevelDFHeatMapDays$decayRate<-decayLevelDFHeatMapDays$decayRate*10
decayLevelDFHeatMapDays$startAge<-decayLevelDFHeatMapDays$startAge/10

ggplot(data = decayLevelDFHeatMapDays, aes(x=startAge, y=decayRate)) + geom_raster(aes(fill = taphonomyValues), hjust = 1, vjust =1)  + 
  scale_fill_distiller(palette = "YlOrRd", direction = -1) +  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA))  +
  labs(title=paste("Taphonomy Simulator - Decay rate (per day), lifespan (days) and preservation"), x="Lifespan (days)", y="Decay rate (per day)", fill = "Average non-zero value") 
ggsave(filename <- paste(outputDirectory,"Heatmap_yellow_red.pdf",sep=""))

ggplot(data = decayLevelDFHeatMapDays, aes(x=startAge, y=decayRate)) + geom_raster(aes(fill = taphonomyValues), hjust = 1, vjust =1)  + 
  scale_fill_distiller(palette = "Purples", direction = -1) +  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA))  +
  labs(title=paste("Taphonomy Simulator - Decay rate (per day), lifespan (days) and preservation"), x="Lifespan (days)", y="Decay rate (per day)", fill = "Average non-zero value")+
ggsave(filename <- paste(outputDirectory,"Heatmap_purples.pdf",sep=""))

ggplot(data = decayLevelDFHeatMapDays, aes(x=startAge, y=decayRate)) + geom_raster(aes(fill = taphonomyValues), hjust = 1, vjust =1)  + 
  scale_fill_distiller(palette = "PuBu", direction = -1) +  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA))  +
  labs(title=paste("Taphonomy Simulator - Decay rate (per day), lifespan (days) and preservation"), x="Lifespan (days)", y="Decay rate (per day)", fill = "Average non-zero value")

ggsave(filename <- paste(outputDirectory,"Heatmap_purple_blue.pdf",sep=""))

ggplot(data = decayLevelDFHeatMapDays, aes(x=startAge, y=decayRate)) + geom_raster(aes(fill = taphonomyValues), hjust = 1, vjust =1)  + 
  scale_fill_gradient2(low = "darkblue", high = "black", mid = "white", midpoint = 0.93) +  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA))  +
  labs(title=paste("Taphonomy Simulator - Decay rate (per day), lifespan (days) and preservation"), x="Lifespan (days)", y="Decay rate (per day)", fill = "Average non-zero value")

ggsave(filename <- paste(outputDirectory,"Heatmap_blue_black_white.pdf",sep=""))

ggplot(data = decayLevelDFHeatMapDays, aes(x=startAge, y=decayRate)) + geom_raster(aes(fill = taphonomyValues), hjust = 1, vjust =1)  + 
  scale_fill_gradientn(colours = rainbow(10)) +  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA))  +
  labs(title=paste("Taphonomy Simulator - Decay rate (per day), lifespan (days) and preservation"), x="Lifespan (days)", y="Decay rate (per day)", fill = "Average non-zero value")

ggsave(filename <- paste(outputDirectory,"Heatmap_rainbow.pdf",sep=""))

ggplot(data = decayLevelDFHeatMapDays, aes(x=startAge, y=decayRate)) + geom_raster(aes(fill = taphonomyValues), hjust = 1, vjust =1)  + 
  scale_fill_gradientn(colours = topo.colors(7)) +  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA))  +
  labs(title=paste("Taphonomy Simulator - Decay rate (per day), lifespan (days) and preservation"), x="Lifespan (days)", y="Decay rate (per day)", fill = "Average non-zero value")

ggsave(filename <- paste(outputDirectory,"Heatmap_topo.pdf",sep=""))