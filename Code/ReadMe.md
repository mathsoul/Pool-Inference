## Reproducibility
The results of the M5 dataset could be reproduced directly by running ?, while the results of the M4 dataset are not. The reason is that the latter takes a long time (several days) to run without parrellel computation. Unforetunately, the parrellel computation mysterious fails on one of our computers. Thus, to avoid this issue, we only include a demo of the M4 dataset using all 17 teams and 119 variables (?). The result still shows that P:Linear dominates other benchmarks.

## Data Manipulation
We did not include the code for data manipulation due to the massive size of the raw data. Instead, we provide the cleaned data in the CleanedData folder and detailed information about the raw data in the RawData folder.

## Alternative Way of Choosing Teams
We did not include the code for an alternative method of choosing teams by solving the maximum diversity problem for two reasons. First, empirical results show that the outcomes are similar to those from selecting the top teams. Thus, we want to focus on the differences between the M4 and M5 datasets, rather than the methods of team selection. Second, we used the gurobi package to solve the maximum diversity problem, which requires a GUROBI solver license that users might not have. Additionally, installing the Gurobi package involves a non-standard procedure that could confuse users. Instead, we have included the solutions to the maximum diversity problem in the TeamIndex folder.
