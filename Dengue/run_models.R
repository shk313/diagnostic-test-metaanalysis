#####  Dengue


library(cmdstanr)
library(bayesplot)
library(tidyverse)

file.path <- "C:/Users/SuzanneKeddie/OneDrive - London School of Hygiene and Tropical Medicine/diagnostic_test_meta_analysis/Reviews/"
review <- "Dengue/"
data <- read.csv(paste0(file.path,review, 'dengue_data_without_conva.csv'), check.names = FALSE)
data <- read.csv(paste0(file.path,review, 'dengue_quadas.csv'), check.names = FALSE)

##################################
# Formatting the data
names(data)
data_sub <- data[data$Include_Exclude != "Exclude" & data$Test_of_interest == "NS1 ELISA" &data$`2x2_include` != "N" & data$`Symptomatic/Asymptomatic` != "Asymptomatic"& data$`Symptomatic/Asymptomatic` != "Asymptomatic ", ]
data_sub <- data_sub[!is.na(data_sub$comparator_test_quadas),]
data_sub <- data_sub[data_sub$time_frame_group %in% c("0-4days","1-7days","1-14days"),]
data_sub <- data_sub[data_sub$time_frame_group %in% c("1-7days"),]
data_sub <- data_sub[!is.na(data_sub$comparator_test_id_all),]
#data_sub <- data
data_sub$unique_id <- 1:nrow(data_sub)
small_data <- data_sub[,c("I+/R+","I+/R-","I-/R+","I-/R-")]

small_data$test <- rowSums(small_data)
sums <- rowSums(small_data[,1:4])
sum(sums) #LQ 6519

## This data is aggregated, want it by individual
tp <- c(1,1)
fp <- c(1,0)
fn <- c(0,1)
tn <- c(0,0)

J<- length(unique(data_sub$`Article ID`)) ##study_id sometimes
l <- length(unique(data_sub$unique_id))
ref_ind <- data_sub$comparator_test_quadas ## which reference test does each study use
#ref_ind <- ifelse(ref_ind == 3,2, ifelse(ref_ind == 4, 3,1))
## Gather data by test either test 1 or 2 even with different reference tests
T1 <- data.frame(results=NA, study=NA,ref_test = NA, id = NA)#,type=NA) ## each list element contains the results from a single study for test 1
T2 <- data.frame(results=NA, study=NA,ref_test = NA, id = NA)#,type=NA)

## instead go through each row because it is no longer a separate study on each line

for(j in 1:l){
  ##for(n in 1:nrow(data)){
  for(i in 1:2){
    a <- rep(tp[i], small_data[j,1])
    b <- rep(fp[i], small_data[j,2])
    c <- rep(fn[i], small_data[j,3])
    d <- rep(tn[i], small_data[j,4])
    
    z <- c(a,b,c,d)
    study <- data_sub$unique_id[j]
    #type <- data_sub$Region_Id[j]
    #type <- data_sub$`serum/CSF`[j]
    #type <- data_sub$Symptomatic[j]
    stopifnot(length(z) == sums[j])
    t <- data.frame(results = z, study = study, ref_test = ref_ind[j], id = j)#, type = type)
    if(i == 1){
      T1 <- rbind(T1, t)
    }else if(i == 2){
      T2 <- rbind(T2,t)
    }
  }
  print(paste0(j, ",",study))
}


T1 <- T1[-1,] ## NA columns from dataframe initialising
T2 <- T2[-1,]
## Checks
sum(is.na(T1$ref_test))
stopifnot(nrow(T1) == sum(sums))
unique(T1$ref_test)
unique(T1$id)
unique(T1$study)

# ## for covariate
# pathA <- ifelse(data_sub$Region_Id == 1, 1,0)
# pathB <- ifelse(data_sub$Region_Id == 2, 1,0)
# pathC <- ifelse(data_sub$Region_Id == 3, 1,0)
# pathD <- ifelse(data_sub$Region_Id == 4, 1,0)
# pathE <- ifelse(data_sub$Region_Id == 5, 1,0)
# pathF <- ifelse(data_sub$Region_Id == 6, 1,0)
# pathG <- ifelse(data_sub$Region_Id == 7, 1,0)



path <- ifelse(data_sub$Region_Id == 1, 1, ifelse(data_sub$Region_Id == 2, 2, ifelse(data_sub$Region_Id == 3, 3, ifelse(data_sub$Region_Id == 4, 4, ifelse(data_sub$Region_Id == 5, 5,ifelse(data_sub$Region_Id == 6, 6,ifelse(data_sub$Region_Id == 7, 7,0)))))))


stan_data <- list(
  J=J,
  T1=T1$results,
  T2=T2$results,
  ind_study = T1$study, 
  ind_unique = T1$id, 
  ##study or unique id
  study= data_sub$unique_id, ##study_id for some
  #N_sums=size,
  N = nrow(T1),
  ref_test =T1$ref_test,
  RT = length(unique(T1$ref_test)),
  l = length(unique(T1$id)))#,
# regionA = pathA,
# regionB = pathB,
# regionC = pathC,
# regionD = pathD,
# regionE = pathE,
# regionF = pathF,
# regionG = pathG,
# region = path,
# R=length(unique(data_sub$Region_Id)))
#path = path,#,
#serum = pathA, ## symptomatic
#csf = pathB)#, ## not symptomatic
#pathogenAB = pathAB)#,


# normal model
stan_test <- "C:/Users/SuzanneKeddie/OneDrive - London School of Hygiene and Tropical Medicine/diagnostic_test_meta_analysis/Reviews/Dengue/hsroc_new.stan"


mod <- cmdstan_model(stan_test)
mod$print()

fit<- mod$sample(
  data = stan_data,#iter_sampling = 2000,
  chains =4,parallel_chains = 4,adapt_delta = 0.96,refresh=500)




fit$cmdstan_diagnose()

file.path.save <- 'C:/Users/SuzanneKeddie/Documents/'
fit$save_object(file =paste0(file.path.save, 'Dengue_NS1_quadas_supp.rds'))
fit <- readRDS("Dengue_PCR_all_ci_19122023.rds")
fit <- readRDS("Dengue_IgG_all_ci_15122023.rds")
# fit2 <- stan(stan_test,
#              data = stan_data, 
#              warmup = 500,
#              iter = 2600,
#              chains =4,
#              thin=2,
#          _ind_lowers <- unlist(fit$summary("Sp")[,"q5"])    #include = FALSE, pars = c('prob'),
#              control=list(adapt_delta = 0.96, max_treedepth = 12, stepsize = 0.1)) 

fit$summary(c("Se_pooled","Sp_pooled"))

##sp_ind_means <- unlist(fit$summary("Sp")[,"mean"])
sp_ind_medians <- fit$summary("Sp", ~quantile(.x, probs =c(0.5)))
sp_ind_medians <- unlist(sp_ind_medians[,"50%"])
sp_ind_lowers <- fit$summary("Sp", ~quantile(.x, probs =c(0.025)))
sp_ind_lowers <- unlist(sp_ind_lowers[,"2.5%"])
sp_ind_uppers <- fit$summary("Sp", ~quantile(.x, probs =c(0.975)))
sp_ind_uppers <- unlist(sp_ind_uppers[,"97.5%"])


se_ind_medians <- fit$summary("Se", ~quantile(.x, probs =c(0.5)))
se_ind_medians <- unlist(se_ind_medians[,"50%"])
se_ind_lowers <- fit$summary("Se", ~quantile(.x, probs =c(0.025)))
se_ind_lowers <- unlist(se_ind_lowers[,"2.5%"])
se_ind_uppers <- fit$summary("Se", ~quantile(.x, probs =c(0.975)))
se_ind_uppers <- unlist(se_ind_uppers[,"97.5%"])
## then these can be ordered like the data
specificity_est <- data.frame(sp_ind_medians=sp_ind_medians,
                              sp_ind_lowers = sp_ind_lowers,
                              sp_ind_uppers = sp_ind_uppers,
                              unique_id = data_sub$Study_no,
                              dat_test = data_sub$DAT_test)

specificity_est <- specificity_est[order(specificity_est$dat_test),]

sensitivity_est <- data.frame(se_ind_medians=se_ind_medians,
                              se_ind_lowers = se_ind_lowers,
                              se_ind_uppers = se_ind_uppers,
                              unique_id = data_sub$Study_no,
                              dat_test = data_sub$DAT_test)

sensitivity_est <- sensitivity_est[order(sensitivity_est$dat_test),]



fit$summary(c("Se_pooled","Sp_pooled"), ~quantile(.x, probs =c(0.025,0.5,0.975)))
fit$summary(c("Se_Ref1","Sp_Ref1"), ~quantile(.x, probs =c(0.025,0.5,0.975)))

fit$summary(c("Se_pred","Sp_pred"), ~quantile(.x, probs =c(0.025,0.5,0.975)))




## plots
color_scheme_set("mix-blue-pink")
p1 <- mcmc_trace(fit$draws(c("theta[1]","alpha[1]", "theta[2]", "alpha[2]")), #n_warmup = 1000,
                facet_args = list(nrow = 2, labeller = label_parsed))
p1 <-p1 + facet_text(size = 18)+theme_bw() + labs(tag="A")

p2 <- mcmc_rank_hist(fit$draws(c("theta[1]","alpha[1]", "theta[2]", "alpha[2]")),
                     facet_args = list(labeller=ggplot2::label_parsed))
p2 <- p2 + facet_text(size = 18)+theme_bw() + labs(tag = "B")


par(mfrow = c(, 2))
p
p2

grid.arrange(p1,p2,nrow=2)
ggsave(grid.arrange(p1,p2,nrow=2), filename = "convergence.png", width=7, height=6, dpi=300)
