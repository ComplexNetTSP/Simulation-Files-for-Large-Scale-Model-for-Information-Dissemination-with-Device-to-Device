SIR model with Gillespie (Tau-leap) algorithm [![Zenodo](https://zenodo.org/badge/doi/10.5281/zenodo.11759.png)](https://zenodo.org/record/11759)
===
##Software

####Usage

#####SIR Simulation
```shell
$ python simulation.py
usage: simulation.py  
  --output OUTPUT
  --duration DURATION
  [-h]
  [--tau TAU]
  [--sim-id SIM_ID]
  [--cell-id CELL_ID]
```
#####SIR Simulation with latent states
```shell
$ python simulation_latent.py
usage: simulation.py  
--output OUTPUT
--duration DURATION
[-h]
[--tau TAU]
[--sim-id SIM_ID]
[--cell-id CELL_ID]
```

#####SIR Simulation with latent states and heterogenous return probability
```shell
$ python simulation_latent_heterogeneous..py
usage: simulation.py  
--output OUTPUT
--duration DURATION
[-h]
[--tau TAU]
[--sim-id SIM_ID]
[--cell-id CELL_ID]
```

##Paper

#### References
[1] Rachit Agarwal, Vincent Gauthier, Monique Backer, Hossam Afifi, "Large Scale Model for Information Dissemination with Device to Device Communication using Call Details Records", 2014.

####Abstract

In a network of devices in close proximity such as *Device to Device* (**D2D**) communication network, we study the dissemination of public safety information at country scale. In order to provide a realistic model for the information dissemination, we extract a spatial distribution of the population of Ivory Coast from census data and determine migration pattern from the Call Detail Records (**CDR**) obtained during the *Data for Development* (**D4D**) challenge. We later apply epidemic model towards the information dissemination process based on the spatial properties of the user mobility extracted from the provided **CDR**. We then propose enhancements by adding latent states to the epidemic model in order to model more realistic user dynamics. Finally, we study dynamics of the evolution of the information spreading through the population.
