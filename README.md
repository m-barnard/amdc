# amdc_clustering
Code for running the method and replicating results, figures, and tables 'Adjacency Matrix Decomposition Clustering for Human Activity Data' (https://arxiv.org/abs/2401.01949) by Martha Barnard, Yingling Fan, Julian Wolfson.

# Table of contents
[Data Documentation](#data-see-daynamica_analysis)

[Replicating figures, tables and results](#reproducing-figures-tables-and-results)

  [Simulation](#simulations-section-3)
  
  [Daynamica Analysis](#daynamica-analysis-section-4)

## Data (see daynamica_analysis/)
> The main original, non-edited data files are found in raw_data/processed_ucal_items_without_surveys_small.csv and clean_data/survey_data.csv. The former is the epidsodic activity data and the later contains the response to the end of day survey question: "Indicate how much you agree with the following statement: Overall, I am concerned with having contracted Coronavirus today"
> 
> raw_data/processed_ucal_items_without_surveys_small.csv data dictionary
>
> | Variable | Meaning |
> |---|---|
> | user_id | Unique user identifyer |
> | type | Describes episode: either ACTIVITY or TRIP |
> | subtype | Describes activity/trip category of episode: <br>HOME, CAR, PERSONAL_BUSINESS, <br>LEISURE_RECREATION, EAT_OUT,WALK, <br>SHOP, OTHER, WORK, EDUCATION, > ><br>BIKE, RAIL, BUS |
> | start_posix | Start date and time of the episode |
> | end_posix | End date and time of the episode |
>
> clean_data/survey_data.csv
>
> | Variable | Meaning |
> |---|---|
> | user_id | Unique user identifyer |
> | eod_survey_date | Date of the end of day survey |
> | contract_c19 | Response to question about COVID-19 concern<br>Options:Strongly Disagree, Disagree, Unsure,<br>Agree, Strongly Agree |
>

## Reproducing figures, tables, and results

### Simulations (Section 3)

### Daynamica Analysis (Section 4)

##  Example notebooks for running code

### Simulation

### Run AMDC on any data

## Run analysis from R/python or shell scrips

### Simulation

### Stability Analysis

### Computational time analysis


