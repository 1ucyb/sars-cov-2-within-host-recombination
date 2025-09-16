# sars-cov-2-within-host-recombination

_Investigating the recombination rate of SARS-CoV-2 through the ONS-CIS._

This repository contains code written for an eight-week internship project aiming to understand recombination in SARS-CoV-2 in a within-host context. Using sequences from the [Office for National Statistics Coronavirus Infection Survey](https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/methodologies/coronaviruscovid19infectionsurveyqmi), I investigated the recombination rate of SARS-CoV-2 through a modified version of the methodology from McDonald et al (2016). I then utilised the forward genetic simulation software [SLiM](https://github.com/MesserLab/SLiM) to generate similar datasets with different recombination rates for comparison.

>[!WARNING]
>The simulation aspects of this project are incomplete. Specifically, the SLiM script is not correctly parameterised, and is likely to generate a few hundred gigabytes of output if run. This work is also not published or peer reviewed.

## Running the analysis

I am unable to share data from the ONS-CIS as it is sensitive, so this analysis can only be run if you have access to that, or are able to modify the code to suit your dataset. You will also need a list of potential artefacts, which I don't know if I'm allowed to share so I'm opting not to.

1. analysis
   - `ons-cis-summary-filter.py` filters down the ONS-CIS summary document to select only useful individuals.
   - `mutation-filter.py` filters out low frequency mutations and artefacts.
   - `residuals-calc.py` calculates descriptive statistics.
   - `empirical-null.py` resamples the dataset to create an empirical null distribution.
   - `descriptive-plots.py` does what it says on the tin - creates nice plots.
2. simulation
   - `run-slim-sim.py` will run the whole simulation and analysis pathway automatically.
   - `output-summary-plot.py` needs to be called separately but will create some basic descriptive plots for your simulation outputs.

## Author

I, Lucy Back, am the author of all code in this project. I am a fourth-year MBiol Biology student at the University of Oxford. The work was carried out over an eight week internship at the Pandemic Sciences Institute, Oxford. It was supervised by Prof. Katrina Lythgoe (University of Oxford), and the idea was hers and Prof. Sarah Otto's (University of British Columbia).

## References

McDonald, M.J., Rice, D.P. & Desai, M.M. (2016). Sex speeds adaptation by altering the dynamics of molecular evolution. _Nature 531_. p233-236. [doi: 10.1038/nature17143](https://doi.org/10.1038/nature17143)
