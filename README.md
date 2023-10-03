# EQRN Results
This repository contains code and reproducible results from the article

Olivier C. Pasche and Sebastian Engelke (2022). "Neural Networks for Extreme Quantile Regression with an Application to Forecasting of Flood Risk". [ArXiv:2208.07590](https://arxiv.org/abs/2208.07590) ([pdf](https://arxiv.org/pdf/2208.07590)).

## EQRN Package
A user-friendly package is under development and available in open access at [github.com/opasche/EQRN](https://github.com/opasche/EQRN). 
On the other hand, this current repository `EQRN_Results` provides the implementation for the experiments and results presented in the main paper, with a snapshot of the EQRN implementation at the time of those experiments, for convinience in reproducibility.

## Contents of this Repository
- The `R/` folder contains a snapshot of the EQRN package functions to ensure reproducibility, as well as plotting function, simulation data functions, and other helper functions.
- The `main/` folder contains the scripts necessary to reproduce every result (and more) from:
	- the iid simulation study in `main/EQRN_iid/`,
	- the sequential dependence simulation study in `main/EQRN_ts/`,
	- the analysis and forecast of river floods in Switzerland in `main/Switzerland/`, and
	- the scripts for the competitor and other models in `main/other_models/`.
- The `Results/` folder contains all the results from every script in `main/` including saved figures and tables.


By Olivier Colin PASCHE\
Research Center for Statistics, University of Geneva (CH), 2022.\
Supported by the Swiss National Science Foundation.

