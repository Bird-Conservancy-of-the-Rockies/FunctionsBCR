# FunctionsBCR
 General purpose utility functions for science and data operations at Bird Conservancy of the Rockies

## Package contents
### Analysis
```BCI``` Given a vector of MCMC posterior estimates for a parameter, provides a summary of the posterior median and credible interval<br>
```logit``` Converts a probability to the logit scale<br>
```expit``` Converts a logit probability to the probability scale<br>
```Derive_regression``` Applies linear model to series of parameters to derive their relationship with specified covariate. Typical use will be to derive trend for yearly population (e.g., occupancy or abundance) estimates<br>
```mcmcList_to_mcmcOutput``` Converts mcmcList to mcmcOutput object<br>

### Data processing & manipulation
```Impute_missing_covs_rf``` Imputes missing covariate values using a random forest informed by existing values of other covariates (maintains correlation structure of a set of covariates)<br>
```LoadFilteredWorkspace``` Builds on `base::load` to load an R workspace from .RData file while specifying objects to exclude<br>
```str_detect_any``` Detects which elements of a given vector contain any element of a second string vector<br>
```SumStats_df``` Generates summary statistics for a given set of covariates<br>
```tssr``` Generates time since sunrise value provided time at a given geographic coordinate<br>
```VIF``` Generates variance inflation factors often used to screen covariates to limit multicollinearity in regression analysis<br>
