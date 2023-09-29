# FunctionsBCR
 General purpose utility functions for science and data operations at Bird Conservancy of the Rockies

## Package contents
### Analysis
```BCI``` Given a vector of MCMC posterior estimates for a parameter, provides a summary of the posterior median and credible interval\
```logit``` Converts a probability to the logit scale\
```expit``` Converts a logit probability to the probability scale\
```RunNimbleParallel``` Wrapper function for fitting a Bayesian model in Nimble while implementing parallel processing and automated checking of convergence (Rhat) and sampling (n.effective) criteria\

### Data processing
```Impute_missing_covs_rf``` Imputes missing covariate values using a random forest informed by existing values of other covariates (maintains correlation structure of a set of covariates)\
```SumStats_df``` Generates summary statistics for a given set of covariates\
```tssr``` Generates time since sunrise value provided time at a given geographic coordinate\
```VIF``` Generates variance inflation factors often used to screen covariates to limit multicollinearity in regression analysis\
