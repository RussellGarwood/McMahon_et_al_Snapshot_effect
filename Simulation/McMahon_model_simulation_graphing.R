## If packages aren't installed, install them, then load them
packages <- c("ggplot2", "Rcpp", "gridExtra", "grid","viridis", "scico")
if(length(packages[!packages %in% installed.packages()[,"Package"]]) > 0){
  install.packages(packages[!packages %in% installed.packages()[,"Package"]])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))
rm("packages")

#Set working directory if not calling from a bash script
#setwd("/SET/YOUR/PATH/HERE")

############################################ Graphing ############################################
outputDirectory<-paste(getwd(),"/graphs/",sep="")

#Load data from the CSVs if required
wd<-getwd()
paste(getwd(),"/results/",sep="") |> setwd()
decayLevelReplicates<-read.csv("decayLevelDF.csv")
randomRatesDataframe<-read.csv("taphonomyDFRandomDecayRates.csv")
lifeSpanReplicates<-read.csv("decayLevelDFLifeSpan.csv")
decayLevelDFHeatMap<-read.csv("decayLevelDFHeatMap.csv")
setwd(wd)


#Initialize plotlist
plot_list <- list()

#A function to put the Y labels to 3DP - this also serves to convert them to text, and thus factorise them 
#Note, threeDPFunction, as well as making the X axis text pretty, also forces this into making the X axis a factor, and so no need to do this separately using as.factor() 
threeDPFunction <- function(x) sprintf("%.4f", x)

#We also need to add the analytical solution to the graphs. For this, given the following:
## Tau is life span in days
## Delta is decay rate per day
#We can get the equilibrium level of decay (or mean integrity as per paper terminology) using the formula:
## equilibrium_constant_decay(tau, delta) = (2 * tau * delta + 1) / (2 * tau * delta + 2) 

############################################ Decay rate graph using replicates ############################################
#This is panel A, and experiment 1

#Simulations at evenly spaced decay rates equating to a time to total decay of 5 days through to a 5 year decay time
#We can use the same decay rates as the simulations, although we need to use the "per day" column - this is delta
#For this, sims use start age of 500 (implying a lifespan of 50 days) -  this is tau 
tau <- 50
#Calculate values for the analytical solution using above tau, and add to our dataframe
decayLevelReplicates$analyticalSolution <- ((2 * tau * decayLevelReplicates$PerDay) + 1) /  ((2 * tau * decayLevelReplicates$PerDay) + 2)
#In previous versions of this graph, I had used springf (above) to create categories for box plot
#But if we want to have an analytical line on the same graph, it makes more sense to plot the boxes on an X with the real values - hence round the data to 3DP
decayLevelReplicates$PerDay3DP <- round(decayLevelReplicates$PerDay, 3)

#We can now plot. Add to plot list for placing into panels
plot_list[[1]]<-ggplot(decayLevelReplicates) + 
  #Boxplots - Note that the width of the box plots is dictated by the continuous x axis values, and so width for the boxes needs to be set manually. Aes set here to allow different X, Y and group by for the curve and boxplots
  geom_boxplot(aes(x = PerDay3DP, y = taphonomyValues, group = factor(PerDay3DP)), outlier.alpha = 0.1, width = 0.008) +
  #Line for analytical solution - using the full perDay value and associated analytical solutio
  geom_line(aes(x = PerDay, y = analyticalSolution, colour = "Analytical solution"), inherit.aes = FALSE, alpha = 0.25, linewidth = 1) + scale_colour_manual(name = NULL, values = c("Analytical solution" = "purple")) +
  #Styling
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                          plot.margin = margin(t = 15, r = 5, b = 5, l = 15), legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
                          legend.background = element_rect(fill = alpha("white", 0.8), colour = "grey70")) + 
  labs(x="Decay rate (per day)", y="Mean integrity", tag = "A") +
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold", size = 16)) 

############################################ Single run with random decay rates ############################################
#This is panel C, and experiment 3

#No analytical solution for this
plot_list[[2]]<-ggplot(data = randomRatesDataframe, aes(x=threeDPFunction(decayRateRounded), y=decayLevels)) + geom_boxplot(outlier.alpha = 0.1) +
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.margin = margin(t = 15, r = 5, b = 5, l = 15)) +
  labs(x="Individual decay rate per day (binned)", y="Individual integrity", tag = "C") +
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold", size = 16))

############################################ Life span graph ############################################
#This is panel B, and experiment 2

#This employs decay rate of 0.001 per iteration, equivalent to 100 days to complete decay - delta is 0.01 per day
delta<- 0.01
#tau is start age in days - this is column days in our data frame
lifeSpanReplicates$analyticalSolution <- ((2 * lifeSpanReplicates$days * delta) + 1) /  ((2 * lifeSpanReplicates$days * delta) + 2)

plot_list[[3]]<-ggplot(data = lifeSpanReplicates) + geom_boxplot(aes(x=days, y=taphonomyValues, group = factor (days)), outlier.alpha = 0.1, width = 40) +
  geom_line(aes(x = days, y = analyticalSolution, colour = "Analytical solution"), inherit.aes = FALSE, alpha = 0.25, linewidth = 1) + scale_colour_manual(name = NULL, values = c("Analytical solution" = "purple")) +
  theme_minimal()  + theme(panel.border = element_rect(color="black", fill=NA), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                           plot.margin = margin(t = 15, r = 5, b = 5, l = 15), legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
                           legend.background = element_rect(fill = alpha("white", 0.8), colour = "grey70")) + 
  labs(x="Lifespan (days)", y="Mean integrity", tag = "B") +
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold", size = 16))


############################################ Heat map ############################################
#This is panel D and Experiment 4

#Now lets look at the heat map - we need to aggregate the values and calculate the mean to do the graphing
meanHeatMap <- aggregate(taphonomyValues ~ decayRate + startAge, data = decayLevelDFHeatMap,FUN = mean)
#Label the mean
names(meanHeatMap)[3] <- "meanPreservation"

#We will want to use an analtyical solution for the SI figure comparing - calculate this here
meanHeatMap$analyticalSolution <- ((2 * (meanHeatMap$startAge/10) * (meanHeatMap$decayRate*10)) + 1) /  ((2 * (meanHeatMap$startAge/10) * (meanHeatMap$decayRate*10)) + 2)
meanHeatMap$difference<-meanHeatMap$analyticalSolution-meanHeatMap$meanPreservation

# Create rectangle boundaries for the heat map - this is required as we're doing a heat map on a log scale
differences <- diff(unique(meanHeatMap$startAge))/2
differences <- append(differences, differences[length(differences)])
window <- differences[match(meanHeatMap$startAge,
                            unique(meanHeatMap$startAge))]

meanHeatMap$startAgemin <- meanHeatMap$startAge - window
meanHeatMap$startAgemax <- meanHeatMap$startAge + window

differences <- diff(unique(meanHeatMap$decayRate))/2
differences <- append(differences, differences[length(differences)])
window <- differences[match(meanHeatMap$decayRate,unique(meanHeatMap$decayRate))]

meanHeatMap$decayRatemin <- meanHeatMap$decayRate - window
meanHeatMap$decayRatemax <- meanHeatMap$decayRate + window

#plot_list[[4]]<-
  ggplot(meanHeatMap,aes(x = startAge,y = decayRate, xmin = startAgemin, xmax = startAgemax, ymin = decayRatemin, ymax = decayRatemax, fill = meanPreservation)) +
  geom_rect() +   scale_x_log10(expand = c(0, 0)) + scale_y_log10(expand = c(0, 0)) + scale_fill_viridis_c(trans = 'reverse') +
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA), legend.position = "bottom", plot.margin = margin(t = 15, r = 5, b = 5, l = 15))  +
  labs(x="Lifespan (days)", y="Decay rate (per day)", fill = "Mean integrity ", tag = "D") +
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold", size = 16))

#Now let's organise these plots using grid extra
grid.arrange(plot_list[[1]],plot_list[[3]],plot_list[[2]],plot_list[[4]], ncol = 2)

#Or save 
ggsave(filename <- paste(outputDirectory,"Taphonomy_plots_paper_figure.png",sep=""), plot = marrangeGrob(plot_list, nrow=2, ncol=2, top = NULL), width = 14, height = 8)

############################################ Stats ############################################

#Mean of replicates for slowest decay rate 
mean(decayLevelReplicates$taphonomyValues[decayLevelReplicates$decayRate==decayRateV[1]])
# 0.5133181

#Highest decay rate - note here we have floating point errors, and so need to calculate equality using a tolerance value, hence all the maths
mean(decayLevelReplicates$taphonomyValues[abs(decayLevelReplicates$decayRate-decayRateV[length(decayRateV)])<abs(decayRateV[1]-decayRateV[2])])
# 0.9519402

#Lets look at spread throughout decay rates
#Use unique here to avoid floating point errors
listSDs<-vector()
listDecays<-vector()
for (decayRate in unique(decayLevelReplicates$decayRate))
{
  stanDev<-sd(decayLevelReplicates$taphonomyValues[decayLevelReplicates$decayRate==decayRate])
  cat(stanDev,"\n")
  listSDs<-rbind(listSDs,stanDev)
  listDecays<-rbind(listDecays,decayRate)
}
max(listSDs)
formatC(max(listSDs), format="e", digits = 2)
mean(listSDs)
formatC(mean(listSDs), format="e", digits = 2)
ggplot(data = data.frame(listSDs,listDecays),aes(x=listDecays, y=listSDs)) + geom_point()

#Mean of replicates for shortest lifespan
mean(lifeSpanReplicates$taphonomyValues[lifeSpanReplicates$startAge==startAgeV[1]])
#0.5441582

#Mean of replicates for longest lifespan
mean(lifeSpanReplicates$taphonomyValues[lifeSpanReplicates$startAge==startAgeV[length(startAgeV)]])
#0.9758938

#Lets look at spread throughout start ages
#Use unique here to avoid floating point errors
listSDs<-vector()
listLifeSpan<-vector()
for (lifeSpan in unique(lifeSpanReplicates$startAge))
{
  stanDev<-sd(lifeSpanReplicates$taphonomyValues[lifeSpanReplicates$startAge==lifeSpan])
  cat(stanDev,"\n")
  listSDs<-rbind(listSDs,stanDev)
  listLifeSpan<-rbind(listLifeSpan,lifeSpan)
}
max(listSDs)
formatC(max(listSDs), format="e", digits = 2)
mean(listSDs)
formatC(mean(listSDs), format="e", digits = 2)
ggplot(data = data.frame(listSDs,listLifeSpan),aes(x=listLifeSpan, y=listSDs)) + geom_point()

#Let's look at the individual simulation with random start ages

#Mean value for minimum decay rate bin
mean(randomRatesDataframe$decayLevels[randomRatesDataframe$decayRateRounded==min(randomRatesDataframe$decayRateRounded)])
# 0.6138001

#Mean value for max decay rate bin
mean(randomRatesDataframe$decayLevels[randomRatesDataframe$decayRateRounded==max(randomRatesDataframe$decayRateRounded)])
# 0.9727146

listSDs<-vector()
for (decayRateRounded in unique(randomRatesDataframe$decayRateRounded))
{
  stanDev<-sd(randomRatesDataframe$decayLevels[randomRatesDataframe$decayRateRounded==decayRateRounded])
  cat(stanDev,"\n")
  listSDs<-rbind(listSDs,stanDev)
}
max(listSDs)
formatC(max(listSDs), format="e", digits = 2)
mean(listSDs)
formatC(mean(listSDs), format="e", digits = 2)

############################################ SI heatmap figure ############################################

#Initialize SI plotlist
SI_plot_list <- list()
#Work out the limits for our fill values using the values created above
lim <- max(abs(meanHeatMap$difference), na.rm = TRUE)

SI_plot_list[[1]]<-ggplot(meanHeatMap,aes(x = startAge,y = decayRate,xmin = startAgemin,xmax = startAgemax,ymin = decayRatemin,ymax = decayRatemax,fill = meanPreservation)) +
  geom_rect() +   scale_x_log10(expand = c(0, 0)) + scale_y_log10(expand = c(0, 0)) + scale_fill_viridis_c(trans = 'reverse') +
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA), legend.position = "bottom",  legend.title.position = "top") +
  labs(x = "Lifespan (days)",y = "Decay rate (per day)", fill = "Mean integrity in simulations")

SI_plot_list[[2]]<-ggplot(meanHeatMap,aes(x = startAge,y = decayRate,xmin = startAgemin,xmax = startAgemax,ymin = decayRatemin,ymax = decayRatemax,fill = analyticalSolution)) +
  geom_rect() +   scale_x_log10(expand = c(0, 0)) + scale_y_log10(expand = c(0, 0)) + scale_fill_viridis_c(trans = 'reverse') +
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA), legend.position = "bottom",  legend.title.position = "top") +
  labs(x = "Lifespan (days)", fill = "Analytically calculated integrity  ")

SI_plot_list[[3]]<-ggplot(meanHeatMap,aes(x = startAge,y = decayRate,xmin = startAgemin,xmax = startAgemax,ymin = decayRatemin,ymax = decayRatemax,fill = difference)) +
  geom_rect() +   scale_x_log10(expand = c(0, 0)) + scale_y_log10(expand = c(0, 0)) + 
  #scale_fill_viridis_c(option = "inferno", limits = c(-lim, lim), name = "Difference between simulation and analytical solution") + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA), legend.position = "bottom",   legend.title.position = "top") + #legend.text = element_text(angle = 45, hjust = 1)) +
  scale_fill_scico(palette = "vik", limits = c(-lim, lim),midpoint = 0, name = "Simulation mean - analytical solution") +
  labs(x = "Lifespan (days)") + guides(fill = guide_colourbar(barwidth = 15))


#Save 
#ggsave(filename <- paste(outputDirectory,"Heatmap_sim_analytical.png",sep=""), width = 7, height = 4)

#Now let's organise these plots using grid extra
grid.arrange(SI_plot_list[[1]],SI_plot_list[[2]],SI_plot_list[[3]], ncol = 3)

#Or save 
ggsave(filename <- paste(outputDirectory,"Heatmap_SI_figure.png",sep=""), plot = marrangeGrob(SI_plot_list, nrow=1, ncol=3, top = NULL), width = 14, height = 5)
