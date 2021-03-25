## Limmathurotsakul et al
# Data extracted from https://doi.org/10.1371/journal.pone.0012485.s004 (14/01/21)
library(dplyr)
library(rstan)
library(loo)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores())
options(shinystan.rstudio = TRUE)

path_to_data <- '/Users/SuzanneKeddie/OneDrive - London School of Hygiene and Tropical Medicine/diagnostic_test_meta_analysis/examples/Melioidosis/'
path_to_code <- '/Users/SuzanneKeddie/Documents/diagnostic-test-metaanalysis/melioidosis/'
raw_data <- read.csv(paste0(path_to_data,'melioidosis_data.csv'))


# Multinomial models set up
#################################
y <- raw_data$freq ## ordered 
M <- 5
N <- 320
C <- 2^M
# alpha_1 <- 100
# beta_1 <- 1

stan_data <- list(N=N,M=M,C=C,y=y)#,
                  # alpha_1 = alpha_1,
                  # beta_1 = beta_1)

stan_test <- paste0(path_to_code,"model0_.stan")

fit<- stan(stan_test,
           data = stan_data, 
           warmup = 500,
           iter = 4000,
           chains = 1,
           thin=1,
           #include = FALSE, pars = c('prob'),
           control=list(adapt_delta = 0.98, max_treedepth = 15)) ## may help if divergent transitions all below the diagonal

# Compute approximation to Loo-CV
log_lik_1 <- extract_log_lik(fit,merge_chains=FALSE)
r_eff <- relative_eff(exp(log_lik_1))
loo_1 <- loo(log_lik_1, r_eff = r_eff)
print(loo_1)
################################

## Individual models set up
##################################
## this is the grouped binomial data, want to make the individual level data for some models
N_rows <- sum(raw_data$freq)
N_tests <- 5

freq_data <- as.numeric(raw_data[,"freq"])
raw_data <- raw_data[,-6]

ind_data <- raw_data
ind_data <- ind_data[rep(seq_len(nrow(ind_data)), freq_data), ]

M <- 5
N <- nrow(ind_data)

stan_data <- list(N=N,M=M,t1=ind_data$a, t2=ind_data$b,t3=ind_data$c,t4=ind_data$d,t5=ind_data$e)

stan_test <- paste0(path_to_code,"model3.stan")

fit<- stan(stan_test,
           data = stan_data, 
           warmup = 500,
           iter = 4000,
           chains = 1,
           thin=1,
           #include = FALSE, pars = c('prob'),
           control=list(adapt_delta = 0.99, max_treedepth = 15)) 

fit_shiny<-as.shinystan(fit, pars=c('b1','b2','b3','b4','Sp1'))#pars arguments (=c("name1","name2")) to not include all parameters

fit_shiny<-launch_shinystan(fit_shiny)

# LOO-CV
log_lik_1 <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1))
loo_1 <- loo(log_lik_1, r_eff = r_eff)
print(loo_1)

log_lik_2 <- extract_log_lik(fit2, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_2))
loo_2 <- loo(log_lik_2, r_eff = r_eff)
print(loo_2)

comp <- loo_compare(loo_1, loo_2)
print(comp)


# Count expected frequencies
#######################################
library(rethinking)
### Want to extract grouped binomial data from this model..
## Look at an indivdiauls estimates test results for the five tests and 
## compare this against a 'profile' of test results to see which of the combinations it is
post <- extract.samples(fit3)
y_prediction <- post$y_pred
dim(y_prediction) ## 3500, 320, 5, (draws, N, M)


C <- 2^M
someData <- rep(NaN, N*C*M); 
arraymatched <- array(someData, dim=c(N,C,M))
nmatched <- array("Na_real", dim=c(N,C))
matchedprofile <- array("NA_real", dim=c(N,C)) 
freqpred <- as.data.frame(matrix("NA_real", nrow = C, ncol = 100))
profiles <- raw_data

for(i in 1:100){
  for(n in 1:N){
    for(c in 1:C){
      for(m in 1:M){
        x <- y_prediction[i,n,m]
        b <- profiles[c,m]
        arraymatched[n,c,m] <- ifelse(x == b, 1, 0)
        
      }
      xx <- as.numeric(arraymatched[n,c,])
      nmatched[n,c] <- sum(xx)
      matchedprofile[n,c] <- ifelse(nmatched[n,c] == M, 1,0)
    }
  }
  
for(c in 1:C){
  xxx <-  as.numeric(matchedprofile[,c])
  freqpred[c,i] <- sum(xxx) ## how many people matched profile 1 for example
}

print(i)

}

## now I have 3500 draws of test combination frequencies, need a summmary measure of this..
## mean?

test <- as.matrix(freqpred)
test <- apply(test, 2,as.numeric)

testmean <- apply(test, 1, mean)

profiles <- cbind(profiles, testmean)
###################################################


### box plot
testing <- freqpred
comb <- c(1:32)
test <- cbind(testing, comb)

testing2 <- freqpred
comb <- c(1:32)
test2 <- cbind(testing2, comb)

testing3 <- freqpred
comb <- c(1:32)
test3 <- cbind(testing3, comb)

library(tidyr)
data_long <- gather(test, combination, frequency, V1:V100, factor_key=TRUE)
data_long$frequency <- as.numeric(data_long$frequency)

data_long2 <- gather(test2, combination, frequency, V1:V100, factor_key=TRUE)
data_long2$frequency <- as.numeric(data_long2$frequency)

data_long3 <- gather(test3, combination, frequency, V1:V100, factor_key=TRUE)
data_long3$frequency <- as.numeric(data_long3$frequency)

obs_data <- data.frame(combination = c(1:32), frequency = c(69,6,0,0,9,0,0,1,14,3,0,5,3,0,3,6,35,15,0,5,5,6,0,7,5,18,0,25,7,11,2,60))

## combine the frequency datasets
data_long$model <- "model0"
data_long2$model <- "model1"
data_long3$model <- "model2"

names(data_long3)[2] <- "draw"
names(data_long3)[1] <- "combination"
#all_data <- rbind(data_long, data_long2)
all_data <- rbind(all_data, data_long3)
#names(all_data)[2] <- "draw"
#names(all_data)[1] <- "combination"
write.csv(all_data, paste0(path_to_data, "boxplot_data.csv"), row.names = FALSE)
# plot, maybe do a facet for each profile


g <- ggplot()+
  geom_boxplot(data = all_data, aes(x=as.factor(combination),y=frequency, fill=model), position = position_dodge(1))+
  geom_point(data =obs_data, aes(x=combination, y = frequency), shape =16, colour = "goldenrod1")+
  scale_x+
  facet_wrap(~combination)
