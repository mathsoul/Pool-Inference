## Introduction
This folder contains the raw data from the M4 Competition along with code for data manipulation. 

**File Organization:**
- Small files are included directly in this folder.
- Large files include download instructions and storage guidance.
- All data manipulation procedures are provided for reproducibility.

## Data Files

### Submission Information
**File:** `Submission Info.xlsx`  
**Description:** Contains participant rankings and competition results.  
**Source:** [M4 Methods GitHub - Point Forecasts](https://github.com/Mcompetitions/M4-methods/tree/master/Point%20Forecasts)

### Training Data
**File:** `Hourly-train.csv`  
**Description:** Training dataset for hourly time series provided by competition organizers.  
**Source:** [M4 Methods GitHub - Train Dataset](https://github.com/Mcompetitions/M4-methods/tree/master/Dataset/Train)

### Test Data
**File:** `Hourly-test.csv`  
**Description:** Test dataset for hourly time series provided by competition organizers.  
**Source:** [M4 Methods GitHub - Test Dataset](https://github.com/Mcompetitions/M4-methods/tree/master/Dataset/Test)

### Series Metadata
**File:** `M4-info.csv`  
**Description:** Contains metadata about all time series including starting timestamps.  
**Source:** [M4 Methods GitHub - Main Dataset](https://github.com/Mcompetitions/M4-methods/tree/master/Dataset)

## Data Manipulation Script

### M4 Data Preparation
**File:** `M4DataManu.R`  
**Requirements:**
1. Download all submissions or the top 17 submissions (IDs: 005, 036, 039, 069, 072, 078, 104, 118, 132, 235, 237, 238, 243, 245, 250, 251, 260) from:  
   [M4 Methods GitHub - Point Forecasts](https://github.com/Mcompetitions/M4-methods/tree/master/Point%20Forecasts)
2. Unzip all of them in this folder.

**Methodology:**
- Implements data scaling using the Mean Absolute Scaled Error (MASE) methodology described in:  
  Makridakis et al. (2020). *The M4 Competition: 100,000 time series and 61 forecasting methods*.  
  International Journal of Forecasting, 36(1), 54-74.  
  DOI: [10.1016/j.ijforecast.2019.04.014](https://doi.org/10.1016/j.ijforecast.2019.04.014)