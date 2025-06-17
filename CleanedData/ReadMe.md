## Introduction  
This folder contains preprocessed datasets for the M4 and M5 competitions.  
- **`M4.Rdata`**  
  - Generated using the cleaning procedures documented in the [`M4RawData`](../M4RawData) folder.  
- **`M5.Rdata`**  
  - Generated using the cleaning procedures documented in the [`M5RawData`](../M5RawData) folder.  

## M4.Rdata
- It contains three data frames: data_test, f_data, and u_data.
   - data_test: the test data of all hourly time series with StartingDate being "1/7/15 12:00."
   - f_data: forecasts from top 17 experts for all hourly time series with StartingDate being "1/7/15 12:00."
   - u_data: forecast errors from top 17 experts for all hourly time series with StartingDate being "1/7/15 12:00."

## M5.Rdata
- It contains four types of data frames: df_true, pred_all, scale2, and df_weight.
   - df_true_Lx: the test data of level x, where $1 \leq x \leq 9$.
   - pred_Lx_all: the predictions from all 50 experts for level x, where $1 \leq x \leq 9$.
   - scale2_Lx: the scaling factor of level x, where $1 \leq x \leq 9$.
      - It is the same scaling factor as in the root mean squared scaled error (RMSSE), which used in the M5 accuracy competition.
        [(Makridakis et al. 2022)](https://doi.org/10.1016/j.ijforecast.2021.11.013)
   - df_weight: the weight for each variable on each level.
      - It is the same weight vector as in weighted RMSSE (WRMSSE), which is used in the M5 accuracy competition. [(Makridakis et al. 2022)](https://doi.org/10.1016/j.ijforecast.2021.11.013)
- It also contains the information regarding the state, store, category, and department. 
