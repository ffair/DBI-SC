# Sparse Clustering for Customer Segmentation with High-dimensional Mixed-type Data

This folder contains data and code used in the paper [Sparse Clustering for Customer Segmentation with High-dimensional Mixed-type Data].


#### Dependencies

The code requires R >= 4.2.0. The R packages required are list in 'functions.R'.

#### Files

'co_two.cpp' is the function to calculate co-occurence distance between any two levels.
'dist_co.cpp' is the function to calculate the categorical part distance between any two samples.
'functions.R' contains the algorithm of DBI-SC method.
'examples.R' shows one simulation scenario when K=4 and the continuous variables are generated from p-generalized normal-polynomial distribution.


#### Usage

The code of DBI-SC method can be found in 'functions.R'.
The R function COC() can be called to realize DBI-SC method after sourcing the file 'functions.R'.

Explanation of parameters of function COC():

- 'x': the dataset containing continuous variables.
- 'y': the dataset containing categorical variables.
- 'K': the number of clusters.
- 'wbounds': the penalty parameter for continuous variables.
- 'vbounds': the penalty parameter for categorical variables. 
- 'group': the indicator of the relationships among dummy variables.
- 'nstart': the number of random starts for the algorithm.
- 'silent': print out progress or not.
- 'maxiter': the maximum number of iterations.


Explanation of output of function COC():

- 'vs': the weights of continuous and categorical (dummy) variables. 
- 'Cs': the clustering results.
- 'DBI.perfeature': the adjusted DBI for each continuous variable.
- 'DBI.perfeaturegroup': the adjusted DBI for each categorical variable.
- 'mDBI': the mDBI critrion of each iteration.

#### Example

In the 'example.R' file, we show one simulation scenario when K=4 and the continuous variables are generated from p-generalized normal-polynomial distribution.
The other scenarios can be reproduced by adjusting the simulation data generation function generatesimd().

#### Acknowledgements

The code for 'SAS` algorithm in this folder was adapted from code in repository of [A Simple Approach to Sparse Clustering](https://github.com/victorpu/SAS_Hill_Climb).




