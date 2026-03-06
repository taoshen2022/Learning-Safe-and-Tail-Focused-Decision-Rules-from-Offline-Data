# Learning Safe and Tail-Focused Decision Rules from Offline Data
## Overview
This repository provides the R implementation for the numerical experiments in the paper ‘Learning Safe and Tail-Focused Decision Rules from Offline Data’.
## Steps to run experiments
1. Install required R dependencies
```r
install.packages("quantreg")
```

2. Load the main estimation functions
```r
# Load data generation functions
source("config.R")
# Load core estimation functions
source("Safe-IDR.R")
```

3. Run a simulation example
Define the key parameters: sample size (n), quantile level (tau), and the mean constraint (mcon). Generate synthetic data with `generate_data()` (defined in config.R), then run `qmestimate()` on the data and parameters to get the estimated decision rule.
Note: The file example.RData (R workspace file) stores results from a single run of the above simulation experiment, using parameters n=1000, tau =0.10, and mcon=3.60.
