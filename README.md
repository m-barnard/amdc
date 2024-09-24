# amdc_clustering
Code for running the method and replicating results, figures, and tables from 'Adjacency Matrix Decomposition Clustering for Human Activity Data'.

# Table of contents
* [General information](#general-information)
  
* [Replicating figures, tables and results](#reproducing-figures-tables-and-results)

    * [Simulation](#simulations-section-3)
  
    * [Daynamica analysis](#daynamica-analysis-section-4)
 
*  [Example notebooks for running code](#example-notebooks-for-running-code)
  
*  [Run analysis from R/python or shell scripts](#run-analysis-from-rpython-or-shell-scripts)
  
    * [Simulations](#simulations)
      
    * [Stability analysis](#stability-analysis)
 
    * [Computational time analysis](#computational-time-analysis)
      
* [Data Documentation](#data-see-daynamica_analysis)

## General information
All analyses were run using R 4.4.0 and Python 3.8.3. The needed python libraries can be loaded as a conda environment using nTreeClus_environ.yaml

In both src/ folders there is a file named ntrees_edit.py. The majority of this code is from https://github.com/HadiJahanshahi/nTreeClus nTreeClus.py file; however, some changes were made in order to run this method for our simulation and Daynamica data analyses. 

## Reproducing figures, tables, and results

All figures, tables, and results are reproduced in .qmd or .ipynb files. The headings in these files label which figures/tables/results are being replicated. The rendered files can be found at the [GitHub pages for this repo](https://anonymous.4open.science/w/amdc_clustering-CE3E/).

Files/results below labelled with (!) recreate tables by reading in results rather than by running the analysis (this was done for computational reasons). To actually run these analyses see the ['Run analysis from R/python or shell scripts'](#run-analysis-from-rpython-or-shell-scripts) section below. 

### Simulations (Section 3)

All files are in simulations/

* figures.qmd simulates sequences and produces Supplementary Figures 1-3
* tables.qmd produces Tables 1-2 and Supplementary Table 1 (!)

### Daynamica analysis (Section 4)

All files are in daynamica_analysis/. While the cleaned data and calculated distance matrices are already present in the repo, scripts/clean_seqs.R and scripts/get_dist_mats.R can be used to produce this data.

* run_nTreeClus_daynamica.ipynb runs nTreeClus on Daynamica day and week sequences and saves the results
* figures_day_analysis.qmd runs AMDC and hierarchical clustering on Daynamica day sequences and produces Figure 2 and Supplementary Figures 4-7
* figures_week_analysis.qmd runs AMDC and hierarchical clustering on Daynamica week sequences and produces Figure 3 and Supplementary Figures 8-9
* figures_weighted_day_analysis.qmd runs weighted AMDC on Daynamica day sequences and produces Figure 4 and Supplementary Figures 10-12
     * Note: Figure 4 can only be reproduced if you first run figures_day_analysis.qmd
* additional_analyses.qmd runs the mixed effects model and prodcues Supplementary Table 2, produces the results in Supplementary Section 1 (!), and produces the results in Section 4.4 and Supplementary Figure 13 (!)

##  Example notebooks for running code

The rendered files can be found at the [GitHub pages for this repo](https://anonymous.4open.science/w/amdc_clustering-CE3E/).

* example_run_sim.qmd and example_run_sim_nTreeClus.ipynb show how to run the simulation (from simulating the sequences to aggregating the results) for a single simulation scenario
* run_any_data.qmd shows how to run AMDC on any data source

## Run analysis from R/python or shell scripts

For running all R scripts, the working directory needs to be set to the outer folder the script is stored in (either simulations/ or daynamica_analysis/).

### Simulations

* To simulate the sequences (simulations/scripts/generate_simulated_sequences/): gen_sim_seqs.R with i = 1,2,5 or use gen_sim.sh
* To run hierarchical clustering and AMDC on the simulated sequences (simulations/scripts/run_clustering_methods): run_full_sim.R with i = 1-36 or use run_all.sh
* To run nTreeClus on the simulated sequences (simulations/scripts/run_clustering_methods): run_ntress.py with i = 0-35 or use run_all_ntrees.sh. After this is completed, run clean_ntrees.R
* To aggregate results see simulations/tables.qmd

### Stability analysis

* To run stability anaysis, run daynamica_analysis/scripts/run_all_boot.R for i = 1-20 or use run_stability.sh
* To aggregate results see daynamica_analysis/additional_analyses.qmd

### Computational time analysis

Note: since computational time depends on the computational hardware/resources used, the results in Section 4.4 cannot be perfectly replicated when running the analysis again. However, the code can be run to determine computational time on one's personal computational hardware/resources.

* To run computational time analysis for hierarchical clustering and AMDC: daynamica_analysis/scripts/run_time_msi.R
* To run computational time analysis for nTreeClus: daynamica_analysis/scripts/ntrees_time.py
* To run all computational time analysis: daynamica_analysis/scripts/assess_time.sh
* To aggregate results see daynamica_analysis/additional_analyses.qmd


## Data
> The main original, non-edited data files are found in daynamica_analysis/raw_data/processed_ucal_items_without_surveys_small.csv and daynamica_analysis/clean_data/survey_data.csv. The former is the epidsodic activity data and the later contains the response to the end of day survey question: "Indicate how much you agree with the following statement: Overall, I am concerned with having contracted Coronavirus today"
> 
> daynamica_analysis/raw_data/processed_ucal_items_without_surveys_small.csv data dictionary
>
> | Variable | Meaning |
> |---|---|
> | user_id | Unique user identifyer |
> | type | Describes episode: either ACTIVITY or TRIP |
> | subtype | Describes activity/trip category of episode: <br>HOME, CAR, PERSONAL_BUSINESS, <br>LEISURE_RECREATION, EAT_OUT,WALK, <br>SHOP, OTHER, WORK, EDUCATION, <br>BIKE, RAIL, BUS |
> | start_posix | Start date and time of the episode |
> | end_posix | End date and time of the episode |
>
> daynamica_analysis/clean_data/survey_data.csv data dictionary
>
> | Variable | Meaning |
> |---|---|
> | user_id | Unique user identifyer |
> | eod_survey_date | Date of the end of day survey |
> | contract_c19 | Response to question about COVID-19 concern<br>Options: Strongly Disagree, Disagree, Unsure,<br>Agree, Strongly Agree |
>




