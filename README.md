# FPCA-STAR-EM-
Ensembled FPCA and STAR Model

Main File: covid_texas_fpca_fatalities

1. "face.estimator" (see file "face.estimator"): estimator implemented using fast covariance estimation
Pro: Good Fitting, Fast
Con: Log likelihood may decrease in STAR EM

2. "estimator.mle"(included in "covid_texas_fpca_fatalities") : stimator implemented using restricted MLE approach. 
Pro: Log likelihood increases in STAR EM because of MLE increase. 
Con: Slow, Bad Fitting

3. "fpca_em"(see file "fpca_em"): EM algorithm described in G.James, T.Hastie's 2001 paper
*Haven't implemented in STAR_EM yet, but I have the idea of including fpca_em in STAR_EM code by feeding \theta_i at each step and get \theta_{i +1} using fpca_em. 
In other words, implemnting a double EM algorithm. 
