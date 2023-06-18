# M2R
This is the code for our M2R project. We decided to work in RStudio, as it is easy to load data and create plots in R than in other languages that we are familiar with. Here is a brief guide outlining the steps we took in creating our models and plots:

+ Read the data into R, converting and adding noise to the times as described in the data, and filtering by station
+ Computed parameters for poisson process, and tested the fit using KS and QQ plots
+ Split data into training and test sets
+ For the self-exciting case, coded recursive form of A-Values and log-likelihood, used built-in optimisation function to compute MLE
+ Computed compensator values using optimal parameters and found p-values using time-rescaling theorem
+ Created QQ-Plots and KS Score boxplots for SE case
+ For mutually-exciting case, coded recursive form of B-Values and proceeded as in SE case
+ For self and mutually-exciting case, optimised the combined log-likelihood and proceeded as before
