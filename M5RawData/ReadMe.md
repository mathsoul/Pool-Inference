## Introduction
This folder contains the raw data from the M5 Competition along with code for data manipulation. 

**File Organization:**
- Small files are included directly in this folder.
- Large files include download instructions and storage guidance.
- All data manipulation procedures are provided for reproducibility.

## Data Files

### Test Set
**File:** `sales_test_evaluation.csv`  
**Description:** Test set provided by the competition organizers.
**Source:** [M5 Competition Google Drive - Dataset](https://drive.google.com/drive/folders/1wxz-TAfVE7uKGCjh405eCb2Q_pG3kAm9)

### Weights of Time Series
**File:** `weights_evaluation.csv`  
**Description:** Weights of different time series provided by the competition organizers.
**Source:** [M5 Competition Google Drive - Dataset](https://drive.google.com/drive/folders/1wxz-TAfVE7uKGCjh405eCb2Q_pG3kAm9)


## Data Manipulation Script

### M5 Data Preparation
**File:** `M5DataManu.R`  
**Requirements:**
1. **Please first install RStudio and then double click the 'Pool Inference.Rproj' file to open the project. It automatically sets the right working directory.**
2. Download all submissions from:  
   [M5 Competition Google Drive - Accuracy Submissions](https://drive.google.com/drive/folders/1NZ1q8Z0gL20TED_W0Phv796MzwghOoPE)
3. Unzip all submissions directly in a subfolder called **Submissions**.
4. Download `sales_train_evaluation.csv` from [M5 Competition Google Drive - Dataset](https://drive.google.com/drive/folders/1wxz-TAfVE7uKGCjh405eCb2Q_pG3kAm9). 

**Methodology:**
- Implements data scaling using the Root Mean Squared Error (RMSE) methodology described in:  
  Makridakis et al. (2022). *M5 accuracy competition: Results, findings, and conclusions*.  
  International Journal of Forecasting, 38(4), 1346-1364.  
  DOI: [10.1016/j.ijforecast.2021.11.013](https://doi.org/10.1016/j.ijforecast.2021.11.013)