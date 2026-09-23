# McMahon et al. - The snapshot effect: can slower decay translate into worse average preservation of fossils? - Supplementary information

Supplementary information for the paper:

McMahon, S., Rodgers, N., Devine, L., Minter, N. and Garwood, R.J. The snapshot effect: can slower decay translate into worse average preservation of fossils?

## Contents 

This repository contains the following folders:

### Analytical

This is the analytical solution to the problem, and contains the associated python script and associated Jupyter notebook.

### Empirical

The root of this folder contains a script (Empirical.R) which loads data from the folder results, and outputs the graphs in the paper into the folder graphs. The script is fully commented, and can be run by calling from bash in a terminal pointed at this folder:

`Rscript ./Empirical.R`

### Equilibration_graphs

This folder contains a zip for each experiment, and each zip comprises the required graphs (in PDF format) to show that the simulations have reach equilibrium. 

### Simulation

This folder contains the simulation code: this is split between an R script, which can be called from bash and will run all experiments contained in the paper:

`Rscript ./McMahon_model_cpp_v9.R`

It also contains a second R script which creates all of the graphs required for the paper. The simulation is written in C++ and this code is found within the folder functions, alongside an R function for creating equilibration graphs. Outputs of the simulation are found in the folder results, and graphs output to the folder graphs (rerunning either script will overwrite the relevant outputs which are also included in the repository).

