## If packages aren't installed, install them, then load them
packages <- c("tidyverse", "grid","viridis", "readxl","gridExtra")
if(length(packages[!packages %in% installed.packages()[,"Package"]]) > 0){
  install.packages(packages[!packages %in% installed.packages()[,"Package"]])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))
rm("packages")

outputDirectory<-paste(getwd(),"/graphs/",sep="")

#Initialize plotlist
plot_list <- list()

#############################
## Taphonomic states
#############################

empiricalDataTaphanomicState<-read_excel("results/Arthropod taphonomic states.xlsx",range = "B4:G15",col_names = TRUE,  sheet = "Taphonomic state")
#Remove empty rows
empiricalDataTaphanomicState<-empiricalDataTaphanomicState[rowSums(!is.na(empiricalDataTaphanomicState)) > 0, ]
#Convert to a long format for graphing
empiricalDataTaphanomicStateLong <- empiricalDataTaphanomicState |>  pivot_longer(cols = starts_with("# in State"),names_to = "State",values_to = "Count")
# Reorder Taxon: put Buenellus first, then the rest in their original order
taxon_order <- c("Buenellus", setdiff(unique(empiricalDataTaphanomicStateLong$Taxon), "Buenellus"))
#Convert taxon from a character vector into a factor, then  explicitly specify the order of its levels
empiricalDataTaphanomicStateLong <- empiricalDataTaphanomicStateLong |> mutate(Taxon = factor(Taxon, levels = taxon_order))
#Identify x position for vertical line (between Buenellus and the rest)
vline_pos <- 1.5  # since Buenellus is the first factor

# Plot
plot_list[[1]]<-ggplot(empiricalDataTaphanomicStateLong, aes(x = Taxon, y = Count, fill = State)) +
  geom_col() +  scale_fill_viridis_d(option = "D", direction = -1) + labs(title = "Taphonomic State (Articulation)", y = "Count", x = "Taxon") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.8, 0.7),  legend.margin = margin(10, 10, 10, 10), legend.spacing = unit(0.5, "cm"), legend.background = element_rect(fill = "white", colour = "grey50"), panel.border = element_rect(color="black", fill=NA)) + 
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "black")

#### Normalised data
empiricalDataTaphanomicStateNormalised <- empiricalDataTaphanomicState[,-1]
empiricalDataTaphanomicStateNormalised[] <- empiricalDataTaphanomicStateNormalised/rowSums(empiricalDataTaphanomicStateNormalised, na.rm = TRUE)
empiricalDataTaphanomicStateNormalised<-add_column(empiricalDataTaphanomicStateNormalised, empiricalDataTaphanomicState[,1], .before = "# in State 5")   
empiricalDataTaphanomicStateNormalisedLong<- empiricalDataTaphanomicStateNormalised |>  pivot_longer(cols = starts_with("# in State"),names_to = "State",values_to = "Proportion")
empiricalDataTaphanomicStateNormalisedLong <- empiricalDataTaphanomicStateNormalisedLong |> mutate(Taxon = factor(Taxon, levels = taxon_order))

# Plot
plot_list[[2]]<- ggplot(empiricalDataTaphanomicStateNormalisedLong, aes(x = Taxon, y = Proportion, fill = State)) +
  geom_col() +  scale_fill_viridis_d(option = "D", direction = -1) + labs(y = "Proportion", x = "Taxon", title = " ") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", panel.border = element_rect(color="black", fill=NA)) +
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "black")


#############################
## Preservational states
#############################

empiricalDataPreservationalState<-read_excel("results/Arthropod taphonomic states.xlsx",range = "C2:F13",col_names = TRUE,  sheet = "Preservational state")
empiricalDataPreservationalState<- empiricalDataPreservationalState[rowSums(!is.na(empiricalDataPreservationalState)) > 0, ]
empiricalDataPreservationalStateLong <- empiricalDataPreservationalState |>  pivot_longer(cols = starts_with("# in State"),names_to = "State",values_to = "Count")
empiricalDataPreservationalStateLong <- empiricalDataPreservationalStateLong |> mutate(Taxon = factor(Taxon, levels = taxon_order))

# Plot
plot_list[[3]]<- ggplot(empiricalDataPreservationalStateLong, aes(x = Taxon, y = Count, fill = State)) +
  geom_col() +  scale_fill_viridis_d(option = "D", direction = -1) + labs(title = "Preservational State (Retention of labile tissues)", y = "Count", x = "Taxon") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.8, 0.7),  legend.margin = margin(10, 10, 10, 10), legend.spacing = unit(0.5, "cm"), legend.background = element_rect(fill = "white", colour = "grey50"), panel.border = element_rect(color="black", fill=NA)) + 
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "black") +
  annotate("text", x = -0.5, y = max(empiricalDataPreservationalStateLong$Count)/2, 
           label = "Preservational State (loss of labile tissues)", angle = 90, hjust = 0.5)

empiricalDataPreservationalStateNormalised <- empiricalDataPreservationalState[,-1]
empiricalDataPreservationalStateNormalised[] <- empiricalDataPreservationalStateNormalised/rowSums(empiricalDataPreservationalStateNormalised, na.rm = TRUE)
empiricalDataPreservationalStateNormalised<-add_column(empiricalDataPreservationalStateNormalised, empiricalDataPreservationalState[,1], .before = "# in State 3")   
empiricalDataPreservationalStateNormalisedLong<- empiricalDataPreservationalStateNormalised |>  pivot_longer(cols = starts_with("# in State"),names_to = "State",values_to = "Proportion")
empiricalDataPreservationalStateNormalisedLong <- empiricalDataPreservationalStateNormalisedLong |> mutate(Taxon = factor(Taxon, levels = taxon_order))

# Plot
plot_list[[4]]<-ggplot(empiricalDataPreservationalStateNormalisedLong, aes(x = Taxon, y = Proportion, fill = State)) +
  geom_col() +  scale_fill_viridis_d(option = "D", direction = -1) + labs(y = "Proportion", x = "Taxon", title = " ") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", panel.border = element_rect(color="black", fill=NA)) + 
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "black")

#############################
## Do graph
#############################

#Now let's organise these plots using grid extra
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]], ncol = 2)

# Draw line separating rows
grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(0.5, 0.5), "npc"), gp = gpar(col = "grey70", lwd = 1))

#Now save an image
paste(outputDirectory,"Empirical_columnplots.png",sep="") |> ggsave(plot = grid.grab(),  # grab the current grid including the vertical line
  width = 12, height = 12)

#############################
## Create a jitter graph
#############################

#Initialize plotlist
jitter_plot_list <- list()

#To do a jitter we need to strip this down into individuals
empiricalDataTaphanomicStateIndividuals <- empiricalDataTaphanomicStateLong |>
  #Extract state number
  mutate(StateNum = as.numeric(str_extract(State, "\\d+"))) |>
  # Expand each row by the count, removing the "Number in state" column we no longer need
  select(-State) |> uncount(Count) |>
  #Add whether it is biomineralized
  mutate(Biomineralization = if_else(Taxon == "Buenellus", "Biomineralised","Not biomineralised" ))

jitter_plot_list[[1]]<-ggplot(empiricalDataTaphanomicStateIndividuals, aes(x = Taxon, y = StateNum, colour = Biomineralization)) +
  geom_jitter(width = 0.3, height = 0.1, alpha = 0.6) + geom_vline(xintercept = vline_pos, linetype = "dotted", color = "black") +
  theme_minimal() +  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", panel.border = element_rect(color="black", fill=NA)) + scale_colour_viridis_d(option = "D", begin = 0, end = 0.6) +
  labs(y = "Taphonomic state", x = "Taxon") +  scale_y_continuous(breaks = scales::pretty_breaks(n = max(empiricalDataTaphanomicStateIndividuals$StateNum))) +
  # Add label at the top of y axis
  annotate("label", x = 0.5, y = max(empiricalDataTaphanomicStateIndividuals$StateNum) + 0.3, label = "Aritculated", hjust = 0, colour = "grey50", fill = "white") +
  annotate("label", x = 0.5, y = min(empiricalDataTaphanomicStateIndividuals$StateNum) - 0.3, label = "Disarticulated", hjust = 0, colour = "grey50") 

empiricalDataPreservationalStateIndividuals <- empiricalDataPreservationalStateLong |>
  #Extract state number
  mutate(StateNum = as.numeric(str_extract(State, "\\d+"))) |>
  # Expand each row by the count, removing the "Number in state" column we no longer need
  select(-State) |> uncount(Count) |>
  #Add whether it is biomineralized
  mutate(Biomineralization = if_else(Taxon == "Buenellus", "Biomineralised","Not biomineralised" ))

jitter_plot_list[[2]] <-ggplot(empiricalDataPreservationalStateIndividuals, aes(x = Taxon, y = StateNum, colour = Biomineralization)) +
  geom_jitter(width = 0.3, height = 0.1, alpha = 0.6) + geom_vline(xintercept = vline_pos, linetype = "dotted", color = "black") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.8, 0.8),  legend.margin = margin(10, 10, 10, 10), legend.spacing = unit(0.5, "cm"), legend.background = element_rect(fill = "white", colour = "grey50"), panel.border = element_rect(color="black", fill=NA)) + 
  scale_colour_viridis_d(option = "D", begin = 0, end = 0.6) + labs(y = "Preservational state", x = "Taxon") + scale_y_continuous(breaks = scales::pretty_breaks(n = max(empiricalDataPreservationalStateIndividuals$StateNum))) +
  annotate("label", x = 0.5, y = max(empiricalDataPreservationalStateIndividuals$StateNum) + 0.3, label = "Labile tissues remain", hjust = 0, colour = "grey50", fill = "white") +
  annotate("label", x = 0.5, y = min(empiricalDataPreservationalStateIndividuals$StateNum) - 0.3, label = "Labile tissues lost", hjust = 0, colour = "grey50") 

#Now let's organise these plots using grid extra
grid.arrange(jitter_plot_list[[1]],  jitter_plot_list[[2]], ncol = 2)

#Or save 
ggsave(filename <- paste(outputDirectory,"Empirical_jitterplots.png",sep=""),
  plot = marrangeGrob(jitter_plot_list, nrow = 1, ncol=2, top = NULL),  width = 14, height = 8)

##########################################################
## Create a jitter graph - biomineralized vs not
##########################################################

jitter_plot_list[[1]]<-ggplot(empiricalDataTaphanomicStateIndividuals, aes(x = Biomineralization, y = StateNum, colour = Biomineralization)) +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.6) + geom_vline(xintercept = vline_pos, linetype = "dotted", color = "black") +
  theme_minimal() +  theme(legend.position = "none",  axis.text.x = element_text(size = 11, color = "black"), panel.border = element_rect(color="black", fill=NA)) + scale_colour_viridis_d(option = "D", begin = 0, end = 0.6) +
  labs(y = "Taphonomic state", x = NULL) +  scale_y_continuous(breaks = scales::pretty_breaks(n = max(empiricalDataTaphanomicStateIndividuals$StateNum))) +
  # Add label at the top of y axis
  annotate("label", size = 3.5, x = 0.5, y = max(empiricalDataTaphanomicStateIndividuals$StateNum) + 0.3, label = "Articulated", hjust = 0, colour = "grey50", fill = "white") +
  annotate("label", size = 3.5, x = 0.5, y = min(empiricalDataTaphanomicStateIndividuals$StateNum) - 0.3, label = "Disarticulated", hjust = 0, colour = "grey50") 

jitter_plot_list[[2]] <-ggplot(empiricalDataPreservationalStateIndividuals, aes(x = Biomineralization, y = StateNum, colour = Biomineralization)) +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.6) + geom_vline(xintercept = vline_pos, linetype = "dotted", color = "black") +
  theme_minimal() +  theme(legend.position = "none",  axis.text.x = element_text(size = 11, color = "black"), panel.border = element_rect(color="black", fill=NA)) + scale_colour_viridis_d(option = "D", begin = 0, end = 0.6) +
  labs(y = "Preservational state", x = NULL) + scale_y_continuous(breaks = scales::pretty_breaks(n = max(empiricalDataPreservationalStateIndividuals$StateNum))) +
  # Add label at the top of y axis
  annotate("label", size = 3.5, x = 0.5, y = max(empiricalDataPreservationalStateIndividuals$StateNum) + 0.3, label = "Labile tissues remain", hjust = 0, colour = "grey50", fill = "white") +
  annotate("label", size = 3.5, x = 0.5, y = min(empiricalDataPreservationalStateIndividuals$StateNum) - 0.3, label = "Labile tissues lost", hjust = 0, colour = "grey50") 

#Now let's organise these plots using grid extra
grid.arrange(jitter_plot_list[[1]],  jitter_plot_list[[2]], ncol = 2)

#Or save 
ggsave(filename <- paste(outputDirectory,"Empirical_jitterplots_mineralization.png",sep=""),
       plot = marrangeGrob(jitter_plot_list, nrow = 1, ncol=2, top = NULL),  width = 10, height = 5)

