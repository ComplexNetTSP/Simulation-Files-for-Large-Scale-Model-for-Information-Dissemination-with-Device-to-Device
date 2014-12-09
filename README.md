SIR model with Gillespie (Tau-leap) algorithm [![Zenodo](https://zenodo.org/badge/doi/10.5281/zenodo.13180.svg)](http://dx.doi.org/10.5281/zenodo.13180)
===
##Paper

#### References

> *Rachit Agarwal, Vincent Gauthier, Monique Becker, Thouraya Toukabrigunes,   Hossam Afifi*, **"Large Scale Model for Information Dissemination with Device to Device Communication using Call Details Records"**, In Press, Compter Communication, 2015.

####Abstract

In a network of devices in close proximity such as *Device to Device* (**D2D**) communication network, we study the dissemination of public safety information at country scale. In order to provide a realistic model for the information dissemination, we extract a spatial distribution of the population of Ivory Coast from census data and determine migration pattern from the Call Detail Records (**CDR**) obtained during the *Data for Development* (**D4D**) challenge. We later apply epidemic model towards the information dissemination process based on the spatial properties of the user mobility extracted from the provided **CDR**. We then propose enhancements by adding latent states to the epidemic model in order to model more realistic user dynamics. Finally, we study dynamics of the evolution of the information spreading through the population.

##Simulator

####Usage

Their is three different simulation available in this package:

1. **simulation.py**: is a simple SIR simulation
2. **simulation_latent.py**: a simple SIR simulation with latent states
3. **simulation_latent_heterogeneous.py**: a simple SIR simulation with latent states and heterogeneous return probability

##### Install
```shell
$ tar zxvf archive.tar.gz
$ cd archive
$ pip install -r requirements.txt
```

##### Command line of the simulator
```shell
usage: simulation_latent_heterogeneous.py [-h]
                                          --output OUTPUT
                                          --duration DURATION
                                          [--tau TAU] [--mu MU]
                                          [--sim-id SIM_ID]
                                          [--cell-id CELL_ID]

Process SIR simulation with latent states and heterogeneous return probability.

optional arguments:
-h, --help           show this help message and exit
--output OUTPUT      output directory
--duration DURATION  simulation duration in days
--tau TAU            simulation step (fraction of day)
--mu MU              simulation mu for latent state (fraction of the
                     population)
--sim-id SIM_ID      simulation step (fraction of day)
--cell-id CELL_ID    initial cellID
```
#####Example of SIR Simulation with latent states

```shell
$ python simulation_latent.py --output ./output/latent/ --duration 7 --tau 0.1 --cell-id 0 --sim-id 1 --mu 0.3
```
