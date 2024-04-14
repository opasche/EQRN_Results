# EQRN Results

This supplementary material repository contains code, data, instructions and reproducible results from the article:

> Olivier C. Pasche and Sebastian Engelke. "Neural Networks for Extreme Quantile Regression with an Application to Forecasting of Flood Risk". Published in *Annals of Applied Statistics*. 

- Published article: *To appear in 2024*.
- Preprint (2022): [ArXiv:2208.07590](https://arxiv.org/abs/2208.07590) ([pdf](https://arxiv.org/pdf/2208.07590)).


## EQRN Package

An open-source R package implementation of the proposed methods is available in open access at [github.com/opasche/EQRN](https://github.com/opasche/EQRN), which is the recommended tool to apply the methodology to new problems. 

On the other hand, this supplementary material repository, [`EQRN_Results`](https://github.com/opasche/EQRN_Results), provides the code and instructions for the experiments and results presented in the paper (and more), with a snapshot of the EQRN implementation at the time of those experiments, for convenience in reproducibility. 


## Contents of this Repository

- The `./R/` folder contains a snapshot of the EQRN package functions to help with reproducibility, as well as plotting function, simulation data functions, competitor methods wrappers or implementation, and other helper functions. 
For exact reproducibility, one can refer to `./R/SessionInfo.txt` for the exact software dependencies, including R package versions.
- The `./main/` folder contains the scripts necessary to reproduce every result (and more) from:
	- the iid simulation study in `./main/EQRN_iid/`,
	- the sequential dependence simulation study in `./main/EQRN_ts/`,
	- the analysis and forecast of river floods in Switzerland in `./main/Switzerland/`, and
	- the scripts for the competitor and other models in `./main/other_models/`.
- The `./Results/` folder contains the results output from every script in `./main/`, including saved figures and tables.
- The `./data/` folder contains the discharge data used in the paper case study application (see the data statement below for obtaining the precipitation data), and an interactive map created to explore the Swiss river and meteorological stations. It also contains intermediate computations and model saves that help significantly speed up computations when reproducing intermediate steps.


## Observational Data Availability

In the case study application, we use river discharge and precipitation data recorded in Switzerland between 1930 and 2014, in the Rhine and Aare basins. The precipitation records can be ordered online for free from [MeteoSwiss](https://gate.meteoswiss.ch/idaweb) and the discharge records from the [Swiss Federal Office for the Environment (FOEN)](https://www.bafu.admin.ch/bafu/en/home/topics/water/state/data/obtaining-monitoring-data-on-the-topic-of-water/hydrological-data-service-for-watercourses-and-lakes.html), for academic research or teaching purposes.

### River Discharge Data

For convenience in reproducing our results, the FOEN has allowed us to share our processed version of the data. Both the stations' metadata and average daily discharge records are available in `./data/data_wrangled/`. The metadata include additional information gathered by the authors. 

If you decide to use this version of the dataset, **please cite both the published article about this work (Pasche and Engelke) and the FOEN source clearly**.

### Precipitation Data

Regarding the precipitation records, we were unfortunately not granted legal permission to share them directly here. However, as stated above, they are available for free to academics and teachers, from the [MeteoSwiss portal (IDAWEB)](https://gate.meteoswiss.ch/idaweb). Filling a short registration form is required to get access to the data catalogue.

For reproducibility, the data needed are at least the daily total precipitation at stations 'BEP', 'GHS', 'THU', 'LTB', 'BRZ' and 'MER' between 1930 and 2014, in a CSV table format. It needs to be placed in `./data/Data/precip.csv`, with the stations identifiers as column names. Optionally, the station metadata used for secondary results (such as the interactive exploration map) can be placed in `./data/Data/precip_info.csv`.

### Original Sources

- Source for precipitation data: **MeteoSwiss** (<https://gate.meteoswiss.ch/idaweb>). 
- Source for river discharge data: **Swiss Federal Office for the Environment (FOEN)**  
(<https://www.bafu.admin.ch/bafu/en/home/topics/water/state/data/obtaining-monitoring-data-on-the-topic-of-water/hydrological-data-service-for-watercourses-and-lakes.html>).


____

By Olivier C. PASCHE  
Research Center for Statistics, University of Geneva (CH), 2022.  
Supported by the Swiss National Science Foundation.  

