## M4.Rdata
- It contains three data frames: data_test, f_data, and u_data.
   - data_test: the test data of all hourly time series.
   - f_data: forecasts from top 17 experts for all hourly time series.
   - u_data: forecast errors from top 17 experts for all hourly time series.

## M5.Rdata
- It contains three types of data frames: df_true, pred_all, and scale2
   - df_true_Lx: the test data of level x, where $1 \leq x \leq 9$.
   - pred_Lx_all: the predictions from all 50 experts for level x, where $1 \leq x \leq 9$.
   - scale2_Lx: the scaling factor of level x, where $1 \leq x \leq 9$.
      - It is the same scaling factor as in the root mean squared scaled error (RMSSE), which used in the M5 accuracy competition. [(Makridakis et al. 2022)](https://doi.org/10.1016/j.ijforecast.2021.11.013)
- It also contains the information regarding the state, store, category, and department. 
