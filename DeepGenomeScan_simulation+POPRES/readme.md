## A Deep Learning Approach for Whole Genome Scan

Qin. X, Chiang, C,Gaggiotti. EO. 2021. Detecting natural selection via deep learning.

### Datasets and Scripts for the article.

### Datasets

Simulation details are described in Capblancq et al. (2018), and the datasets can be accessed from the Dryad Digital Repository at: https://datadryad.org/review?doi=doi:10.5061/dryad.1s7v5.

The datasets used for the analyses described in this manuscript were obtained from dbGap at https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000145.v4.p2 through dbGap accession number phs000145.v4.p2.

#### Scripts used for the results reported in this paper are available through this github website.
Directory Simulations: Store the scripts for doing genome scan and making figures for simulations.

Directory POPRES: Store the scripts for doing genome scan and making figures for POPRES data. 

#### Simulations

DL_trial.R: Scripts for pre-trials selecting a better model from model list.

DLModels.R: Model lists for genome scan, these models built from different sources/libraries.

Simulations_DL_MLP.R : The scripts used for doing genome scan using DL_MLP model on 100 simulations.

Simulation_power.R: Scripts used for estimating the power of genome scan for three approaches, pcadapt, RDA, and DeepGenomeScan, as well as for making the figures.


#### POPRES

POPRES_DL_MLP_Scan: Scripts used in analysis of POPRES data, including EnvGenoemScan and SpGenomeScan.

POPRES_Manhattan: Scripts for making figures for POPRES data.

In terms of step by step examples for compiling different neural network models using DeepGenomeScan, please see [tutorials](https://xinghuq.github.io/DeepGenomeScan/index.html).

Please contact xq5@st-andrews.ac.uk if you have any problems with the paper and the package.
