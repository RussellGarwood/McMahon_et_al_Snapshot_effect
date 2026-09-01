############################################ Graph single run function ############################################
graphSingleRun <- function(taphonomyDF, titleString, savePath = "")
{
  # Initialize plotlist
  plot_list <- list()
  
  #First the differences in the mean non zero value - this should start varying around a number near zero fairly randomly once the simulation is stationary
  #Before this time there will probably be an increase between (early) iterations
  #plot_list[[1]]<- ggplot(data = taphonomyDF) + geom_line(mapping = aes_string(y="differencesVector", x = "iterationVector")) + xlim(100,max(taphonomyDF$iterationVector)) +
  #  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA)) + 
  #  labs(title=paste("Taphonomy Simulator - Difference"), x="Time", y="Difference between iterations")
  
  #This shows the average, non zero value in the list of organisms. It wills tart near one when everything is alive, and then reduce to a ~static value when this is stationary
  plot_list[[1]]<-ggplot(data = taphonomyDF) + geom_line(mapping = aes_string(y="taphonomyVectorLocal", x = "iterationVector")) +  scale_x_continuous(expand = c(0, 0)) +
    #xlim(0,max(taphonomyDF$iterationVector)) +
    theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA)) +
    labs(title=paste("Decay level"), x="Time", y="Mean non zero taphonomy")
  
  #This shows the number of organisms alive at any time - it'll increase afte ther simulation starts, peak before things start dying, and then reduce to a static value
  plot_list[[2]]<- ggplot(data = taphonomyDF) + geom_line(mapping = aes_string(y="aliveVector", x = "iterationVector")) +  scale_x_continuous(expand = c(0, 0)) +
    #+ xlim(100,max(taphonomyDF$iterationVector)) +
    theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA)) +
    labs(title=paste("Number Alive"), x="Time", y="Number Alive")
  
  #This shows the number of decaying organisms - will start at zero when the simulation kicks off, increase, then reach a steady state
  plot_list[[3]]<-ggplot(data = taphonomyDF) + geom_line(mapping = aes_string(y="decayingVector",  x = "iterationVector")) +  scale_x_continuous(expand = c(0, 0)) +
    #+ xlim(100,max(taphonomyDF$iterationVector)) +
    theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA)) +
    labs(title=paste("Number Decaying"), x="Time", y="Number Decaying")
  
  #Save as a multipage PDF using grid extra
  if(nchar(savePath)>2)
  {
    ggsave(
      filename <- paste(savePath,"Taphonomy_plots.pdf",sep=""),
      plot = marrangeGrob(plot_list, nrow=1, ncol=3, top = titleString),
      width = 10, height = 4
    )
  }
  #Plot here
  else grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], ncol = 3, top = textGrob(titleString, gp = gpar(fontsize = 20)))
}
