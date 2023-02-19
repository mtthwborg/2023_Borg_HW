# Overview
This repository contains R code and datasets that together provide a reproducible example to simulate the main methodology and results for the retrospective anaysis from Borg et al. 2023, *Environmental Health Perspectives*, Current and projected heatwave-attributable occupational injuries, illnesses, and associated economic burden in Australia: a national time-series analysis. These files can be downloaded after clicking on the green button *Code* at this GitHub repository.


# System Requirements
## Hardware requirements
The code should only require a standard computer with enough RAM to support the user-defined operations. R memory usage when running this code should not exceed 4 GB of RAM.

## Software requirements
This code was tested using the following systems:
* macOS: Ventura Version 13.2

This software uses R and was tested using R version 4.2.2. The required packages can be installed using the code below within approximately 30 seconds:
`install.packages(c('readxl','data.table','devtools','lubridate','stringr','zoo','mgcv','statmod','tweedie','dlnm','mixmeta','FluMoDL','weathermetrics','pryr','MASS'))`

The versions of the attached non-base packages used during testing were:
` pryr_0.1.6 tweedie_2.3.5 statmod_1.5.0 FluMoDL_0.0.3 mvmeta_1.0.3 mixmeta_1.2.0 dlnm_2.4.7 mgcv_1.8-41 nlme_3.1-162 zoo_1.8-10 stringr_1.5.0 lubridate_1.9.` magrittr_2.0.3 data.table_1.14.6 readxl_1.4.1 weathermetrics_1.2.2 MASS_7.3-58.2`


# Files
This repository comes with the .Rproj (project) file, datasets and analysis files.

## Data
The datasets for use with this example:
  * *brambilla.max.rda* is a derived meteorological dataset sourced externally from [Brambilla et al. 2022, *Data Br.*, Hygrothermal climate analysis: An Australian dataset](https://doi.org/10.1016/j.dib.2022.108291). The code to create this dataset is included in *01_DataPrep.R*.
  * *pop.rda* is the study dataset of Australian worker population counts in Adelaide, Brisbane, Darwin, Hobart, Melbourne, Perth and Sydney, sourced from the [Australian Bureau of Statistics (ABS)](https://www.abs.gov.au/statistics/labour/employment-and-unemployment/labour-force-australia/latest-release). This includes estimates of the number of indoor and outdoor workers derived/based from [ABS Census TableBuilder Basic data](https://tablebuilder.abs.gov.au/webapi/jsf/login.xhtml) as described in the methodology by Borg et al. 2022. This dataset is stored at [https://doi.org/10.25909/21332592](https://doi.org/10.25909/21332592).
  * *public.holidays.rda* is the study list of Australian public holidays from 2004 to 2023. This dataset is sourced from [https://doi.org/10.25909/6311e7a0dcb3f67](https://adelaide.figshare.com/articles/dataset/Public_holidays_in_Australian_capital_cities_from_2004_to_2023/20732449).
  * *school.holidays.rda* is the study list of Australian school holidays from 2004 to 2023 (as a binary variable and a factor variable with levels for each school holiday period). This dataset is sourced from [https://doi.org/10.25909/6311e7b3bc760](https://adelaide.figshare.com/articles/dataset/School_holidays_in_Australian_capital_cities_from_2004_to_2023/20732173).
  * *ccia_future2.rda* is a derived meteorological dataset sourced from Climate Change in Australia and is available online [Application-ready gridded (5km) datasets - Daily time-series](https://www.climatechangeinaustralia.gov.au/en/obtain-data/download-datasets/](https://www.climatechangeinaustralia.gov.au/en/obtain-data/download-datasets/). The code to create this dataset is included in *05_Projections.R*.
  * *projpop.rda* is a study dataset of Australian projected worker population counts in each Australian capital city, derived from the online [Australian Bureau of Statistics (ABS) dataset available online]([https://www.abs.gov.au/statistics/labour/employment-and-unemployment/labour-force-australia/latest-release]). 
    * A simulated claims dataset is created in file "01_DataPrep.R"

*pop.rda*, *public.holidays.rda* and *school.holidays.rda* were also used in the main study analysis by Borg et al. 2023.
  
## Analysis
The numbered code files from *01_DataPrep.R* to *04_AnalysisStage2.R* reproduce the study analysis using the aforementioned datasets. They are designed to be run in numerical order:
  * *01_DataPrep.R* sets up the packages and datasets. It creates a simulated (fake) dataset to represent the number of claims and their associated total costs.
  * *02_AnalysisPrep.R* prepares the parameters for statistical analysis. The default parameters replicate those used for the main analysis with total costs as the outcome variable and daily maximum wet bulb globe temperature as the exposure variable. The parameters can be changed to simulate many of the supplementary analyses.
  * *03_AnalysisStage1.R* runs the Stage 1 statistical analysis, the individual models prior to multivariate meta-analysis.
  * *04_AnalysisStage2.R* runs the Stage 2 statistical analysis. This includes the multivariate meta-analysis, refitting the Stage 1 models with the best linear unbiased predictors (BLUPs), and using the refitted models to generate the main study results.
  * *05_Projections.R* projects the BLUPs derived in Stage 2 using the projected climate and population datasets.

*02_AnalysisPrep.R* includes a brief loop that can run the code for both *03_AnalysisStage1.R*, *04_AnalysisStage2.R*, and *05_Projections.R* with multiple outcome variables. Please run *01_DataPrep.R* and the preceding code in *02_AnalysisPrep.R* before running this loop.

*heat.index2.R* calls a modified function to calculate the Heat Index. It is a modified version of weathermetrics::heat.index that (1) uses a threshold of 80 instead of 79 and (2) does not automatically round results to two decimal figures. This code was only used when calculating Heat Index, which is not part of the main analysis.

# Results
*02_AnalysisPrep.R* will create the folder in the directory to store the results. Within this folder, *03_AnalysisStage1.R* will create an additional folder with the names of the outcome variable and exposure variable, as well as an additional folder called "Stage 1" inside this to replace the results that it generates. *04_AnalysisStage2.R* will create an additional folder named "Stage 2" to store its results, and *05_Projections.R* will create an additional folder to store its results, including folders assuming no climate adaptation (Adapt0) and adaptation (Adapt1).

*.gitignore* is set to exclude the folder *Results* from the working directory.

This example includes results for Melbourne and Sydney only (two models in total, instead of seven). **Because the example claims dataset was simulated and not representative of real data, this example's results are neither similar to those of the main study nor should be used to draw any formal conclusions.**

## Stage 1
The files generated are:
  * *Reduced coef and model fit.csv*. This contains, for each model, the number of days included for analysis, the reduced coefficients from the distributed lag non-linear models, the AIC and the dispersion parameters. The AIC is summed across models.
  * *Coefficients.csv*. This includes a range of descriptive statistics for the model coefficients. This was used to determine whether the included model predictors were useful predictors or not and to guide modelling choices.
  * *Tweedie shape parameters.csv*. This includes the Tweedie shape parameters from the outcome variables. It is only generated if the outcome has a Tweedie distribution. It is included in case one wishes to repeat the models without needing to re-estimate the shape parameters, which can be time-consuming with large datasets.

## Stage 2
A folder will be created to store results representing each model. Three graphs will be created here to describe the exposure-lag-response relationships prior to multivariate meta-analysis. 
The files generated are:
  * *MM tests.csv*. This contains the multivariate extension of the Cochran Q test and the I^2^ statistic for the multivariate meta-analysis.
  * *Overall e-r.png*. Generates the national overall cumulative exposure-response relationship graph.
  * *Overall lag.png*. Generates the national overall cumulative lag-response relationship graph including curves at the thresholds for heatwaves and sevre heatwaves.
  * *Fractions.csv*. This includes the proportion of the outcome attributable to heatwaves across the study period.
  * *Numbers.csv*. This includes the number of outcome values attributable to heatwaves across the study period.
  * *Numbers per year.csv*. This includes the number of outcome values attributable to heatwaves per year.
  * *City level*. This is a folder including overall cumulative exposure-response graphs for each city. Graphs wtih a "z" in front of their names include histograms representing the spread of Excess Heat Factor values. The "City vs BLUP" graph showcases the city-level curves before meta-analysis and after being refitted with their best linear unbiased predictors (BLUPs).
  
## Projections
A folder will be created to store results. This will include the file
The files generated are:
  * *Projected DMT.csv*. This contains descriptive statistics on projected meteorological variables, specifically relating to daily mean temperature (DMT).
  * *Fractions.csv*. This includes the proportion of the outcome attributable to heatwaves across the projected time period per year.
  * *Numbers.csv*. This includes the number of outcome values attributable to heatwaves across the projected time period per year.
  
*Fractions.csv* and  *Numbers.csv* will be stored in a folder named either *Adapt0* or *Adapt1*, indicating no climate adaptation and climate adaptation, respectively.
