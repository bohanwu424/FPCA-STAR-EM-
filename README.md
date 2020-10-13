# FPCA-STAR-EM-
Ensembled FPCA and STAR Model

Main File: covid_texas_fpca_fatalities

1. "face.estimator" (see file "face.estimator"): estimator implemented using fast covariance estimation
Pro: Good Fitting, Fast
Con: Log likelihood may decrease in STAR EM

2. "estimator.mle"(included in "covid_texas_fpca_fatalities") : stimator implemented using restricted MLE approach. 
Pro: Log likelihood increases in STAR EM because of MLE increase. 
Con: Slow, Bad Fitting

3. "fpca_em"(see file "fpca_em"): EM algorithm described in G.James, T.Hastie's 2001 paper. Two immediate applications of James and Hastie's idea: 
- James EM: an EM framework that takes Y and Theta (multiple parameters) as input, and return fitted Y and updated Theta. 
- James Estimator: an estimator built on James EM, takes Y and Theta and return fitted Y and updated Theta. 

EM^2
- Basic Idea:
A STAR EM framework that takes in 2 estimators: estimator_initialize, estimator_em.
In this case, estimator_initialize = face.estimator, estimator_em = james.estimator. 
Two features: 
1)Use estimator_initialize to initialize Theta. 
2)At each step S, use estimator_em to maximize the conditional expectation Theta|\hat(Y). (M step). 

The rest of EM^2 framework is the same as STAR EM. 
