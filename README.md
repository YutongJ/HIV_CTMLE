# Identifying HIV sequences that escape antibody neutralization using random forests and collaborative targeted learning

## Prerequisite:
To fully replicate the research result, some R packages are required for executing the main code file:   
`dummies`
`MASS`
`plyr`
`randomForest`
`stat`

All the other dependent packages are listed at the beginning of each code file.

## Original Datasets:

The original dataset is from CATNAP database. We illustrate our method with five antibodies: VRC01, VRC26.08, PGT145, PGT121 and 10-1074.


## Code Details: 

**Main Functions:** 

*AAk_step2*\
CTMLE results from full dataset. Full data are used for both model fitting and predictions.

*cv_triplets*\
CTMLE results from cross-validation. Full data are firstly divided into training set and testing set; predicitions are obtained via cross-validated process.


*OR_AAk_tmle_ctmle*\
The estimation for one particular site of interest. This process includes two types of results: (i) full data CTMLE results and (ii) cross-validated CTMLE results. The p-values from wald test of this site are also provided.



**Nuisance Functions:** 

*fluctuate_Q*\
This function is used in TMLE to fluctuate the initial predictions from outcome regression.


*pred_fluctuate_Q*\
Given (1) estimated epsilon (from TMLE logistic regression), and (2) initial OR estimate, this function can generate updated OR estimate. It is typically useful in cross-validation for testing set. 


*waldtest*\
Wald test for all pairwise comparisons of all potential amino acids at certain site.