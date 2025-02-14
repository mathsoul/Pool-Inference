## Reproducibility
The results for the M4 and M5 datasets can be reproduced by running the files M4.R and M5.R, respectively. The M5.R script takes approximately 10 minutes to process a single set of teams (one table), whereas M4.R requires significantly more time—several days without parallel computation and several hours with it. However, parallel execution fails unpredictably on one of our Windows machines. To circumvent this issue, we provide a demo script (M4All.R) that runs the M4 dataset using all 17 teams and 119 variables. The results confirm that P:Linear consistently outperforms other benchmarks.

The simulation results could be reproduced by running file simulation.R.

## Data Manipulation
We did not include the code for data manipulation due to the massive size of the raw data. Instead, we provide the cleaned data in the CleanedData folder and detailed information about the raw data in the RawData folder.

## File M4IdxGen.R and File M5IdxGen.R
These two files generate different selections of teams for forecast aggregation. Since solving the maximum diversity problem requires the gurobi package, which depends on a GUROBI solver license that users may not have, we have provided precomputed solutions in the TeamIndex folder.

## File func.R 
This file contains all the functions needed in this repository. We provide some documentation to explain the functions, especially the ones that avoid large matrix multiplication when the number of products is massive.


