## 
library(dplyr)
file_names <- c("santander_summaries/1.csv", "santander_summaries/2.csv", "santander_summaries/3.csv", "santander_summaries/4.csv", "santander_summaries/5.csv", 
                "santander_summaries/6.csv", "santander_summaries/7.csv", "santander_summaries/8.csv", "santander_summaries/9.csv", "santander_summaries/10.csv", 
                "santander_summaries/11.csv", "santander_summaries/12.csv", "santander_summaries/13.csv", "santander_summaries/14.csv", "santander_summaries/15.csv", "santander_summaries/16.csv")

#read in data
df_full = data.frame(read.table(file_names[1], header=FALSE, sep=','))  
# Create an empty list to store data frames
names(df_full) = c('src','dst','time','duration')
for (file in file_names[2:length(file_names)]) {
  df = read.table(file, header=FALSE, sep=',')
  names(df) = c('src','dst','time','duration')
  df_full <- rbind(df_full, df)
}

#convert times, sort in ascending order and filter by station
set.seed(111)
df_full$time = (df_full$time - min(df_full$time)) / 60 + runif(dim(df_full)[1])
df_full = df_full[order(df_full$time),]
df_stations = list()
for(s in sort(unique(df_full$src))){
  df_stations[[s]] = df_full[df_full$src == s,]
}

#calculate parameter for poisson process, no. of events / time period
lambda = list()
for(s in sort(unique(df_full$src))){
  df_filtered = df_stations[[s]][df_stations[[s]]$time < (60 * 24 * 7 * 12),]
  lambda[[s]] = dim(df_filtered)[1] / (60 * 24 * 7 * 12)
}

#inter-arrival times
for(s in sort(unique(df_full$src))){
  df_stations[[s]]$int_times = c(df_stations[[s]]$time[1], diff(df_stations[[s]]$time))
}

#p-values
for(s in sort(unique(df_full$src))){
  df_stations[[s]]$pp_pval = exp(-lambda[[s]] * df_stations[[s]]$int_times)
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(df_full$src))){
  f = ecdf(df_stations[[s]]$pp_pval)
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_pp = numeric()
for(s in sort(unique(df_full$src))){
  ks_scores_pp = c(ks_scores_pp, ks.test(df_stations[[s]]$pp_pval, 'punif')$statistic)
}
boxplot(ks_scores_pp, col='lightblue', pch=20)
boxplot(ks_scores_pp_test, col='lightblue', pch=20)
boxplot(ks_scores_pp_training, col='lightblue', pch=20)

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(trainingset$src))){
  f = ecdf(df_stationsh0[[s]]$pp_pval)
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')


par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(testset$src))){
  f = ecdf(df_stationsh1[[s]]$pp_pval)
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

trainingset = df_full[1:3164820,]
testset = df_full[3164821:nrow(df_full),]

df_stationsh0 = list()
for(s in sort(unique(trainingset$src))){
  df_stationsh0[[s]] = trainingset[trainingset$src == s,]
}

df_stationsh1 = list()
for(s in sort(unique(testset$src))){
  df_stationsh1[[s]] = testset[testset$src == s,]
}

T_training = max(trainingset$time) - min(trainingset$time)
lambdah0 = list()
for(s in sort(unique(trainingset$src))){
  lambdah0[[s]] = dim(df_stationsh0[[s]])[1] / T_training
}

lambdah0[[842]] = 0

for(s in sort(unique(testset$src))){
  df_stationsh1[[s]]$int_times = c(df_stationsh1[[s]]$time[1], diff(df_stationsh1[[s]]$time))
}

for(s in sort(unique(testset$src))){
  df_stationsh1[[s]]$pp_pval = exp(-lambdah0[[s]] * df_stationsh1[[s]]$int_times)
}

for(s in sort(unique(trainingset$src))){
  df_stationsh0[[s]]$int_times = c(df_stationsh0[[s]]$time[1], diff(df_stationsh0[[s]]$time))
}

for(s in sort(unique(trainingset$src))){
  df_stationsh0[[s]]$pp_pval = exp(-lambdah0[[s]] * df_stationsh0[[s]]$int_times)
}

ks_scores_pp_training = numeric()
for(s in sort(unique(trainingset$src))){
  ks_scores_pp_training = c(ks_scores_pp_training, ks.test(df_stationsh0[[s]]$pp_pval, 'punif')$statistic)
}

ks_scores_pp_test = numeric()
for(s in sort(unique(testset$src))){
  ks_scores_pp_test = c(ks_scores_pp_test, ks.test(df_stationsh1[[s]]$pp_pval, 'punif')$statistic)
}


#Self-exciting hawkes process
kernel_function <- function(t, alpha, beta) {
  return(alpha*exp(-beta*t))
}

recursive_kernel <- function(event_times, alpha, beta){
  n = length(event_times)
  A_values = rep(0, n)
  A_values[1] = 0
  for (i in 2:n) {
    A_values[i] = (kernel_function(event_times[i] - event_times[i-1], alpha, beta) / alpha) * (1 + A_values[i-1])
  }
  return(A_values)
}

negative_loglikelihood = function(pars, event_times){
  ## pars = (lambda_tilde, alpha_tilde, beta_tilde)
  lambda_tilde = pars[1]
  alpha_tilde = pars[2]
  beta_tilde = pars[3]
  ## Transform back to original scale
  lambdaop = exp(lambda_tilde)
  alphaop = exp(alpha_tilde)
  betaop = exp(beta_tilde) + alphaop
  ## Calculate the recursive sequence of A
  A_values = recursive_kernel(event_times, alphaop, betaop)
  ## First part of log-likelihood
  p1 = sum(log(lambdaop + alphaop * A_values))
  ## Second part of log-likelihood
  t_last = event_times[length(event_times)]
  p2 = lambdaop* t_last
  ## Third part of log-likelihood
  p3 = alphaop / betaop * sum(exp(-betaop * (t_last - event_times)) - 1)
  ## Return negative log-likelihood
  return(p2 - p1 - p3)
}

theta = list()
for(s in sort(unique(df_full$src))){
  print(s)
  theta[[s]] = optim(par=c(0,0,0), fn=negative_loglikelihood, 
                             event_times = df_stations[[s]]$time)$par
}

lambdaop = list()
alphaop = list()
betaop = list()
for(s in sort(unique(df_full$src))){
  lambdaop[[s]] = exp(theta[[s]][1])
  alphaop[[s]]= exp(theta[[s]][2])
  betaop[[s]] = exp(theta[[s]][3]) + exp(theta[[s]][2])
}

#compensator, using time-rescaling theorem
compavals = list()
for(s in sort(unique(df_full$src))){
  compavals[[s]] = recursive_kernel(df_stations[[s]]$time, alphaop[[s]], betaop[[s]])
}

adiff = list()
for(s in sort(unique(df_full$src))){
  adiff[[s]] = c(compavals[[s]][1], diff(compavals[[s]]))
}

compdiff = list()
for(s in sort(unique(df_full$src))){
  compdiff[[s]] = lambdaop[[s]] * df_stations[[s]]$int_times - (alphaop[[s]] / betaop[[s]]) * (adiff[[s]] - 1)
}

negcompdiff = list()
for(s in sort(unique(df_full$src))){
  negcompdiff[s] = lapply((compdiff[s]), function(x) -x)
}

hawkespvals = list()
for(s in sort(unique(df_full$src))){
  hawkespvals[[s]] = lapply((negcompdiff[s]), function(x) exp(x))
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(df_full$src))){
  f = ecdf(unlist(hawkespvals[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')


ks_scores_hawkes = numeric()
for(s in sort(unique(trainingset$src))){
  ks_scores_hawkes = c(ks_scores_hawkes, ks.test(unlist(hawkespvals[[s]]), 'punif')$statistic)
}
boxplot(ks_scores_hawkes, col='lightblue', pch=20)

lambdaoptraining = list()
alphaoptraining = list()
betaoptraining = list()
for(s in sort(unique(trainingset$src))){
  lambdaoptraining[[s]] = exp(thetatraining[[s]][1])
  alphaoptraining[[s]]= exp(thetatraining[[s]][2])
  betaoptraining[[s]] = exp(thetatraining[[s]][3]) + exp(thetatraining[[s]][2])
}

compavalstraining = list()
for(s in sort(unique(trainingset$src))){
  compavalstraining[[s]] = recursive_kernel(df_stationsh0[[s]]$time, alphaoptraining[[s]], betaoptraining[[s]])
}

adifftraining = list()
for(s in sort(unique(trainingset$src))){
  adifftraining[[s]] = c(compavalstraining[[s]][1], diff(compavalstraining[[s]]))
}

compdifftraining = list()
for(s in sort(unique(trainingset$src))){
  compdifftraining[[s]] = lambdaoptraining[[s]] * df_stationsh0[[s]]$int_times - (alphaoptraining[[s]] / betaoptraining[[s]]) * (adifftraining[[s]] - 1)
}

negcompdifftraining = list()
for(s in sort(unique(trainingset$src))){
  negcompdifftraining[s] = lapply((compdifftraining[s]), function(x) -x)
}

hawkespvalstraining = list()
for(s in sort(unique(trainingset$src))){
  hawkespvalstraining[[s]] = lapply((negcompdifftraining[s]), function(x) exp(x))
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(trainingset$src))){
  f = ecdf(unlist(hawkespvalstraining[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')


ks_scores_hawkes_training = numeric()
for(s in sort(unique(trainingset$src))){
  ks_scores_hawkes_training = c(ks_scores_hawkes_training, ks.test(unlist(hawkespvalstraining[[s]]), 'punif')$statistic)
}
boxplot(ks_scores_hawkes_training, col='lightblue', pch=20)

alphaoptest = alphaoptraining
betaoptest = betaoptraining
lambdaoptest = lambdaoptraining

alphaoptest[[842]] = 0
betaoptest[[842]] = 0
lambdaoptest[[842]] = 0

compavalstest = list()
for(s in sort(unique(testset$src))){
  compavalstest[[s]] = recursive_kernel(df_stationsh1[[s]]$time, alphaoptest[[s]], betaoptest[[s]])
}

adifftest = list()
for(s in sort(unique(testset$src))){
  adifftest[[s]] = c(compavalstest[[s]][1], diff(compavalstest[[s]]))
}

compdifftest = list()
for(s in sort(unique(testset$src))){
  compdifftest[[s]] = lambdaoptest[[s]] * df_stationsh1[[s]]$int_times - (alphaoptest[[s]] / betaoptest[[s]]) * (adifftest[[s]] - 1)
}

negcompdifftest = list()
for(s in sort(unique(testset$src))){
  negcompdifftest[s] = lapply((compdifftest[s]), function(x) -x)
}

hawkespvalstest = list()
for(s in sort(unique(testset$src))){
  hawkespvalstest[[s]] = lapply((negcompdifftest[s]), function(x) exp(x))
}

hawkespvalstest[[842]] = double(693)


par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(testset$src))){
  f = ecdf(unlist(hawkespvalstest[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_hawkes_test = numeric()
for(s in sort(unique(testset$src))){
  ks_scores_hawkes_test = c(ks_scores_hawkes_test, ks.test(unlist(hawkespvalstest[[s]]), 'punif')$statistic)
}

boxplot(ks_scores_hawkes_test, col='lightblue', pch=20)

#mutually exciting case
#sort by destinations
df_stations_a = list()
for(s in sort(unique(df_full$src))){
  df_stations_a[[s]] = df_full[df_full$dst == s,]
}
#calculate arrival times
for(s in sort(unique(df_full$src))){
  df_stations_a[[s]]$arrival_times = sort(c(df_stations_a[[s]]$time + df_stations_a[[s]]$duration/60))
}

#find intervals
intervs = list()
for(s in sort(unique(df_full$src))){
  intervs[[s]] = findInterval(df_stations[[s]]$time, df_stations_a[[s]]$arrival_times)
}

#defined kernel
mutual_kernel_function <- function(t, gamma, delta) {
  return(gamma * exp(-delta * t))
}

#recursive calculation of B-values
mutual_recursive_kernel <- function(departuretimes, eventdata, gamma, delta, intervals){
  n = length(departuretimes)
  B_values = rep(0, n)
  B_values[1] = 0
  if(intervals[1] > 0){
    B_values[1] = sum(exp(-delta * (departuretimes[1] - eventdata$arrival_times[1:intervals[1]])))
  }
  for (i in 2:n) {
    if(intervals[i-1] < intervals[i]){
      B_values[i] = exp(-delta * (departuretimes[i] - departuretimes[i-1])) * B_values[i-1] + sum(exp(-delta * (departuretimes[i] - eventdata$arrival_times[(intervals[i-1]+1):intervals[i]])))
    } else {
      B_values[i] = exp(-delta * (departuretimes[i] - departuretimes[i-1])) * B_values[i-1]
    }
  }
  return(B_values)
}

#negative log-likelihood defined similarly to the self-exciting case
negative_loglikelihood_b = function(pars, departuretimes, eventdata, intervals){
  lambdab_tilde = pars[1]
  gamma_tilde = pars[2]
  delta_tilde = pars[3]
  lambdabop = exp(lambdab_tilde)
  gammaop = exp(gamma_tilde)
  deltaop = exp(delta_tilde) + gammaop
  B_values = mutual_recursive_kernel(departuretimes, eventdata, gammaop, deltaop, intervals)
  p1b = sum(log(lambdabop + gammaop * B_values))
  t_lastb = max(c(departuretimes[[length(departuretimes)]], eventdata[[length(eventdata)]]))
  p2b = lambdabop * t_lastb
  p3b = (gammaop / deltaop) * sum(exp(-deltaop * (t_lastb - departuretimes)) - 1)
  return(p2b - p1b - p3b)
}


thetab = list()
for(s in sort(unique(df_full$src))){
  print(s)
  thetab[[s]] = optim(par = c(-3,-3,-3), fn=negative_loglikelihood_b, 
                             departuretimes = df_stations[[s]]$time, eventdata = df_stations_a[[s]],
                      intervals=intervs[[s]])$par
}

lambdabop = list()
gammaop = list()
deltaop = list()
for(s in sort(unique(df_full$src))){
  lambdabop[[s]] = exp(thetab[[s]][1])
  gammaop[[s]]= exp(thetab[[s]][2])
  deltaop[[s]] = exp(thetab[[s]][3]) + exp(thetab[[s]][2])
}

compbvals = list()
for(s in sort(unique(df_full$src))){
  print(s)
  compbvals[[s]] = mutual_recursive_kernel(df_stations[[s]]$time, df_stations_a[[s]], gammaop[[s]], deltaop[[s]], intervs[[s]])
}

bdiff = list()
for(s in sort(unique(df_full$src))){
  bdiff[[s]] = c(compbvals[[s]][1], diff(compbvals[[s]]))
}

compdiff_b = list()
for(s in sort(unique(df_full$src))){
  compdiff_b[[s]] = lambdabop[[s]] * df_stations[[s]]$int_times[2:length(df_stations[[s]]$int_times)] - (gammaop[[s]] / deltaop[[s]]) * (diff(mutual_recursive_kernel(df_stations[[s]]$time, df_stations_a[[s]], gammaop[[s]], deltaop[[s]], intervs[[s]])) -diff(intervs[[s]]))
}

hawkespvals_b = list()
for(s in sort(unique(df_full$src))){
  hawkespvals_b[s] = lapply((compdiff_b[s]), function(x) exp(-x))
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(testset$src))){
  f = ecdf(unlist(hawkespvals_b[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_hawkes_b = numeric()
for(s in sort(unique(testset$src))){
  ks_scores_hawkes_b = c(ks_scores_hawkes_b, ks.test(unlist(hawkespvals_b[[s]]), 'punif')$statistic)
}

boxplot(ks_scores_hawkes_b, col='lightblue', pch=20)

df_stations_ah0 = list()
for(s in sort(unique(trainingset$src))){
  df_stations_ah0[[s]] = trainingset[trainingset$dst == s,]
}

df_stations_ah1 = list()
for(s in sort(unique(testset$src))){
  df_stations_ah1[[s]] = testset[testset$dst == s,]
}

for(s in sort(unique(trainingset$src))){
  df_stations_ah0[[s]]$arrival_times = sort(c(df_stations_ah0[[s]]$time + df_stations_ah0[[s]]$duration/60))
}

for(s in sort(unique(trainingset$src))){
  df_stations_ah1[[s]]$arrival_times = sort(c(df_stations_ah1[[s]]$time + df_stations_ah1[[s]]$duration/60))
}

intervsh0 = list()
for(s in sort(unique(trainingset$src))){
  intervsh0[[s]] = findInterval(df_stationsh0[[s]]$time, df_stations_ah0[[s]]$arrival_times)
}

intervsh1 = list()
for(s in sort(unique(testset$src))){
  intervsh1[[s]] = findInterval(df_stationsh1[[s]]$time, df_stations_ah1[[s]]$arrival_times)
}

thetabtraining = list()
for(s in sort(unique(trainingset$src))){
  print(s)
  thetabtraining[[s]] = optim(par = c(-3,-3,-3), fn=negative_loglikelihood_b, 
                      departuretimes = df_stationsh0[[s]]$time, eventdata = df_stations_ah0[[s]],
                      intervals=intervsh0[[s]])$par
}

lambdaboptraining = list()
gammaoptraining = list()
deltaoptraining = list()
for(s in sort(unique(trainingset$src))){
  lambdaboptraining[[s]] = exp(thetabtraining[[s]][1])
  gammaoptraining[[s]]= exp(thetabtraining[[s]][2])
  deltaoptraining[[s]] = exp(thetabtraining[[s]][3]) + exp(thetabtraining[[s]][2])
}

compbvalstraining = list()
for(s in sort(unique(trainingset$src))){
  compbvalstraining[[s]] = mutual_recursive_kernel(df_stationsh0[[s]]$time, df_stations_ah0[[s]], gammaoptraining[[s]], deltaoptraining[[s]], intervsh0[[s]])
}

bdifftraining = list()
for(s in sort(unique(trainingset$src))){
  bdifftraining[[s]] = c(compbvalstraining[[s]][1], diff(compbvalstraining[[s]]))
}

compdiff_btraining = list()
for(s in sort(unique(trainingset$src))){
  compdiff_btraining[[s]] = lambdaboptraining[[s]] * df_stationsh0[[s]]$int_times[2:length(df_stationsh0[[s]]$int_times)] - (gammaoptraining[[s]] / deltaoptraining[[s]]) * (diff(mutual_recursive_kernel(df_stationsh0[[s]]$time, df_stations_ah0[[s]], gammaoptraining[[s]], deltaoptraining[[s]], intervsh0[[s]])) -diff(intervsh0[[s]]))
}

hawkespvals_btraining = list()
for(s in sort(unique(trainingset$src))){
  hawkespvals_btraining[s] = lapply((compdiff_btraining[s]), function(x) exp(-x))
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(trainingset$src))){
  f = ecdf(unlist(hawkespvals_btraining[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_hawkes_btraining = numeric()
for(s in sort(unique(trainingset$src))){
  ks_scores_hawkes_btraining = c(ks_scores_hawkes_btraining, ks.test(unlist(hawkespvals_btraining[[s]]), 'punif')$statistic)
}

boxplot(ks_scores_hawkes_btraining, col='lightblue', pch=20)

gammaoptest = gammaoptraining
deltaoptest = deltaoptraining
lambdaboptest = lambdaboptraining

gammaoptest[[842]] = 0
deltaoptest[[842]] = 0
lambdaboptest[[842]] = 0

compbvalstest = list()
for(s in sort(unique(testset$src))){
  compbvalstest[[s]] = mutual_recursive_kernel(df_stationsh1[[s]]$time, df_stations_ah1[[s]], gammaoptest[[s]], deltaoptest[[s]], intervsh1[[s]])
}

bdifftest = list()
for(s in sort(unique(testset$src))){
  bdifftest[[s]] = c(compbvalstest[[s]][1], diff(compbvalstest[[s]]))
}

compdiff_btest = list()
for(s in sort(unique(testset$src))){
  compdiff_btest[[s]] = lambdaboptest[[s]] * df_stationsh1[[s]]$int_times[2:length(df_stationsh1[[s]]$int_times)] - (gammaoptest[[s]] / deltaoptest[[s]]) * (diff(mutual_recursive_kernel(df_stationsh1[[s]]$time, df_stations_ah1[[s]], gammaoptest[[s]], deltaoptest[[s]], intervsh1[[s]])) -diff(intervsh1[[s]]))
}

hawkespvals_btest = list()
for(s in sort(unique(testset$src))){
  hawkespvals_btest[s] = lapply((compdiff_btest[s]), function(x) exp(-x))
}

hawkespvals_btest[789] = 0

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(trainingset$src))){
  f = ecdf(unlist(hawkespvals_btest[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_hawkes_btest = numeric()
for(s in sort(unique(testset$src))){
  ks_scores_hawkes_btest = c(ks_scores_hawkes_btest, ks.test(unlist(hawkespvals_btest[[s]]), 'punif')$statistic)
}

boxplot(ks_scores_hawkes_btest, col='lightblue', pch=20)

#self and mutually-exciting case
negative_loglikelihood_gen = function(pars, departuretimes, eventdata, event_times, intervals){
  lambda_tilde = pars[1]
  alpha_tilde = pars[2]
  beta_tilde = pars[3]
  gamma_tilde = pars[4]
  delta_tilde = pars[5]
  lambdaop = exp(lambda_tilde)
  alphaop = exp(alpha_tilde)
  betaop = exp(beta_tilde) + alphaop
  gammaop = exp(gamma_tilde)
  deltaop = exp(delta_tilde) + gammaop
  B_values = mutual_recursive_kernel(departuretimes, eventdata, gammaop, deltaop, intervals)
  A_values = recursive_kernel(event_times, alphaop, betaop)
  p1_gen = sum(log(lambdaop + alphaop * A_values + gammaop * B_values))
  t_last_gen = max(c(departuretimes[[length(departuretimes)]], eventdata[[length(eventdata)]], event_times[[length(event_times)]]))
  p2_gen = lambdaop * t_last_gen
  p3_gen = (alphaop / betaop) * sum(exp(-betaop * (t_last_gen - event_times)) - 1) + (gammaop / deltaop) * sum(exp(-deltaop * (t_last_gen - departuretimes)) - 1)
  return(p2_gen - p1_gen - p3_gen)
}

thetagen = list()
for(s in sort(unique(df_full$src))){
  print(s)
  thetagen[[s]] = optim(par = c(-3,-3,-3, -3, -3), fn=negative_loglikelihood_gen, 
                      departuretimes = df_stations[[s]]$time, eventdata = df_stations_a[[s]],
                      intervals=intervs[[s]], event_times = df_stations[[s]]$time)$par
}

lambdaopgen = list()
alphaopgen = list()
betaopgen = list()
gammaopgen = list()
deltaopgen = list()
for(s in sort(unique(df_full$src))){
  lambdaopgen[[s]] = exp(thetagen[[s]][1])
  alphaopgen[[s]] = exp(thetagen[[s]][2])
  betaopgen[[s]] = exp(thetagen[[s]][3]) + exp(thetagen[[s]][2])
  gammaopgen[[s]]= exp(thetagen[[s]][4])
  deltaopgen[[s]] = exp(thetagen[[s]][5]) + exp(thetagen[[s]][4])
}

compdiff_gen = list()
for(s in sort(unique(df_full$src))){
  compdiff_gen[[s]] = lambdaopgen[[s]] * df_stations[[s]]$int_times[2:length(df_stations[[s]]$int_times)] - (alphaopgen[[s]] / betaopgen[[s]]) * (diff(recursive_kernel(df_stations[[s]]$time, alphaopgen[[s]], betaopgen[[s]])) - 1) - (gammaopgen[[s]] / deltaopgen[[s]]) * (diff(mutual_recursive_kernel(df_stations[[s]]$time, df_stations_a[[s]], gammaopgen[[s]], deltaopgen[[s]], intervs[[s]])) -diff(intervs[[s]]))
}

hawkespvals_gen = list()
for(s in sort(unique(df_full$src))){
  hawkespvals_gen[s] = lapply((compdiff_gen[s]), function(x) exp(-x))
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(df_full$src))){
  f = ecdf(unlist(hawkespvals_gen[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_hawkes_gen = numeric()
for(s in sort(unique(df_full$src))){
  ks_scores_hawkes_gen = c(ks_scores_hawkes_gen, ks.test(unlist(hawkespvals_gen[[s]]), 'punif')$statistic)
}

boxplot(ks_scores_hawkes_gen, col='lightblue', pch=20)

thetagentraining = list()
for(s in sort(unique(trainingset$src))){
  print(s)
  thetagentraining[[s]] = optim(par = c(-3,-3,-3, -3, -3), fn=negative_loglikelihood_gen, 
                              departuretimes = df_stationsh0[[s]]$time, eventdata = df_stations_ah0[[s]],
                              intervals=intervsh0[[s]], event_times = df_stationsh0[[s]]$time)$par
}

lambdaopgentraining = list()
alphaopgentraining = list()
betaopgentraining = list()
gammaopgentraining = list()
deltaopgentraining = list()
for(s in sort(unique(trainingset$src))){
  lambdaopgentraining[[s]] = exp(thetagentraining[[s]][1])
  alphaopgentraining[[s]] = exp(thetagentraining[[s]][2])
  betaopgentraining[[s]] = exp(thetagentraining[[s]][3]) + exp(thetagentraining[[s]][2])
  gammaopgentraining[[s]]= exp(thetagentraining[[s]][4])
  deltaopgentraining[[s]] = exp(thetagentraining[[s]][5]) + exp(thetagentraining[[s]][4])
}

compdiff_gentraining = list()
for(s in sort(unique(trainingset$src))){
  compdiff_gentraining[[s]] = lambdaopgentraining[[s]] * df_stationsh0[[s]]$int_times[2:length(df_stationsh0[[s]]$int_times)] - (alphaopgentraining[[s]] / betaopgentraining[[s]]) * (diff(recursive_kernel(df_stationsh0[[s]]$time, alphaopgentraining[[s]], betaopgentraining[[s]])) - 1) - (gammaopgentraining[[s]] / deltaopgentraining[[s]]) * (diff(mutual_recursive_kernel(df_stationsh0[[s]]$time, df_stations_ah0[[s]], gammaopgentraining[[s]], deltaopgentraining[[s]], intervsh0[[s]])) -diff(intervsh0[[s]]))
}

hawkespvals_gentraining = list()
for(s in sort(unique(trainingset$src))){
  hawkespvals_gentraining[s] = lapply((compdiff_gentraining[s]), function(x) exp(-x))
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(trainingset$src))){
  f = ecdf(unlist(hawkespvals_gentraining[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_hawkes_gentraining = numeric()
for(s in sort(unique(trainingset$src))){
  ks_scores_hawkes_gentraining = c(ks_scores_hawkes_gentraining, ks.test(unlist(hawkespvals_gentraining[[s]]), 'punif')$statistic)
}

boxplot(ks_scores_hawkes_gentraining, col='lightblue', pch=20)

lambdaopgentest = lambdaopgentraining
alphaopgentest = alphaopgentraining
betaopgentest = betaopgentraining
gammaopgentest = gammaopgentraining
deltaopgentest = deltaopgentraining

lambdaopgentest[[842]] = 0
alphaopgentest[[842]] = 0
betaopgentest[[842]] = 0
gammaopgentest[[842]] = 0
deltaopgentest[[842]] = 0


compdiff_gentest = list()
for(s in sort(unique(testset$src))){
  compdiff_gentest[[s]] = lambdaopgentest[[s]] * df_stationsh1[[s]]$int_times[2:length(df_stationsh1[[s]]$int_times)] - (alphaopgentest[[s]] / betaopgentest[[s]]) * (diff(recursive_kernel(df_stationsh1[[s]]$time, alphaopgentest[[s]], betaopgentest[[s]])) - 1) - (gammaopgentest[[s]] / deltaopgentest[[s]]) * (diff(mutual_recursive_kernel(df_stationsh1[[s]]$time, df_stations_ah1[[s]], gammaopgentest[[s]], deltaopgentest[[s]], intervsh1[[s]])) -diff(intervsh1[[s]]))
}

hawkespvals_gentest = list()
for(s in sort(unique(testset$src))){
  hawkespvals_gentest[s] = lapply((compdiff_gentest[s]), function(x) exp(-x))
}

par(mar=c(4,4,1,1))
eval_grid = seq(0,1, length.out=250)
plot(c(0,1),c(0,1), type='n', xlab='Theoretical quantiles', ylab='Empirical quantiles')
for(s in sort(unique(testset$src))){
  f = ecdf(unlist(hawkespvals_gentest[s]))
  points(eval_grid, f(eval_grid), type='l', col='lightgray')
}
abline(0,1,col='red')

ks_scores_hawkes_gentest = numeric()
for(s in sort(unique(testset$src))){
  ks_scores_hawkes_gentest = c(ks_scores_hawkes_gentest, ks.test(unlist(hawkespvals_gentest[[s]]), 'punif')$statistic)
}

boxplot(ks_scores_hawkes_gentest, col='lightblue', pch=20)
