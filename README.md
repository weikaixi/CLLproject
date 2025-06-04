# PFS prediction in patients with CLL

This repository contains the R and Python codes for implementing the fPM score and related scores.

allmodels.RData - the fPM model and the other machine learning models tested in the paper

sessioninfo.txt - the session information, workstation info and package version information 
### 🔧 Files for Model Testing

- **`demo_script.R`** – **This script runs the model using a test dataset.**
- **`demo_data.csv`** – **This is the dummy dataset used for model testing.**

prediction_result.csv - the predicted fPM score and risk group by demo_script.R

fPM_performance file - the implement of fPM score

fPM_and_otherMLmodels file - the implement of other ML models

other_prognostic_scores file - the implement of 4 other CLL prognostic scores

genetic_analysis file - the genetic analysis result

distribution_analysis file - the distribution analysis of PFS time and fPM score

permutation_test file - an example about how to produce permutaiton test in prognostic model

Contact - Weikaixin Kong at kong.weikaixin@helsinki.fi
