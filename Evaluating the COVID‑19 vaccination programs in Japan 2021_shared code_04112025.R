# Code Start.........................
rm(list=ls())

# Packages:
pacman::p_load(tidyverse, magrittr, openxlsx, MASS)

################################################################################
################################################################################
# 1. Estimating protected proportion

## Condition settings
diag_bias <- 4 # Ascertainment bias of confirmed cases (1/): 1 , 2, 4, 8
vac_program <- 1 # Vaccination was implemented or not

## Population demography from statistics(2021)
age_gp <- c("a0009","a1019","a2029","a3039","a4049","a5059","a6069","a7079","a8089","a90")
pop_90over_2021 <- c(9427,10937,12642,13910,17905,17076,15260,16384,9435,2526) * 1000

pop_pref <- t(data.frame(pop_90over_2021))
rownames(pop_pref) <- NULL
colnames(pop_pref) <- age_gp

ag_labels <- colnames(pop_pref)
N <- as.vector(pop_pref)
m <- length(N)

## Immune proportion by age group
### Vaccinated individuals from VRS
#vaccinationTable_place = ".../File1.VaccinationTable.csv"
vaccinationTable_place = "./File1.VaccinationTable.csv"
vaccinationTable <- as.data.frame(read_csv(file = vaccinationTable_place, show_col_types = FALSE))

vaccinationTable %<>% dplyr::select(-status) %>% filter(date < as.Date("2021-12-01"))



vaccinationTable %>% head()
vacDataStartDate <- as.Date(vaccinationTable[1,"date"][[1]])
vacDataEndDate <- as.Date(vaccinationTable[nrow(vaccinationTable),"date"][[1]])
t_duration <- as.numeric(vacDataEndDate - vacDataStartDate + 1)

#### Figure for vaccinated individuals and proportion
vaccination_fig <- vaccinationTable %>% mutate(c0009 = cumsum(a0009),
                                               c1019 = cumsum(a1019),
                                               c2029 = cumsum(a2029),
                                               c3039 = cumsum(a3039),
                                               c4049 = cumsum(a4049),
                                               c5059 = cumsum(a5059),
                                               c6069 = cumsum(a6069),
                                               c7079 = cumsum(a7079),
                                               c8089 = cumsum(a8089),
                                               c90 = cumsum(a90)) %>% 
  dplyr::select(date,c0009,c1019,c2029,c3039,c4049,c5059,c6069,c7079,c8089,c90) %>% gather(2:(m+1), value = "vacinated", key = "age_group") 

vaccination_fig %>% ggplot() +
  geom_line(aes(x=date, y=vacinated, color=age_group), size = 0.8) +
  scale_x_date(date_labels = "%m/%d/%Y", date_breaks = "2 weeks", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x="", y="Cumulative number of vaccinated individuals", colour = "Age group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())


vaccination_fig2 <- vaccinationTable %>% mutate(c0009 = cumsum(a0009)/N[1]*100,
                                                c1019 = cumsum(a1019)/N[2]*100,
                                                c2029 = cumsum(a2029)/N[3]*100,
                                                c3039 = cumsum(a3039)/N[4]*100,
                                                c4049 = cumsum(a4049)/N[5]*100,
                                                c5059 = cumsum(a5059)/N[6]*100,
                                                c6069 = cumsum(a6069)/N[7]*100,
                                                c7079 = cumsum(a7079)/N[8]*100,
                                                c8089 = cumsum(a8089)/N[9]*100,
                                                c90 = cumsum(a90)/N[10]*100) %>% 
  dplyr::select(date,c0009,c1019,c2029,c3039,c4049,c5059,c6069,c7079,c8089,c90) 

colnames(vaccination_fig2)[2:ncol(vaccination_fig2)] <- c("0–9","10–19","20–29","30–39","40–49","50–59","60–69","70–79","80–89","≥90")
vaccination_fig2 %<>% gather(2:(m+1), value = "vacinated", key = "age_group") 

vaccination_fig2$age_group <- factor(vaccination_fig2$age_group, levels = unique(vaccination_fig2$age_group))

vaccination_fig2 %>% ggplot() +
  geom_line(aes(x=date, y=vacinated, color=age_group), size = 0.8) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "2 weeks", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  labs(x="Calendar time in 2021", y="Vaccination coverage (%)", colour = "Age group") +
  theme(text=element_text(size=12, family="sans",color="black"),
        axis.text=element_text(size=10, family="sans",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = c(0.13, 0.63))

## Vaccine efficacy profile: From Nir et al, 2022 
SeffTemporalprofile <- rep(0,300)
SeffTemporalprofile[1:34] <- seq(0,0.97,length.out=34)
SeffTemporalprofile[34+(1:300)] <- seq(0.97,0,length.out=300)
effInfProfile <- SeffTemporalprofile%*%t(rep(1,m)) 

## Protected fractions by age group
imm_ind_mx <- matrix(NA, nrow = t_duration, ncol = m)
for(t in 1:t_duration){
  
  imm_ind_mx_s <- matrix(NA, nrow = t_duration, ncol = m)
  for(s in 1:t){
    
    imm_ind_mx_s[s,] <- effInfProfile[s,] * as.matrix(vaccinationTable[t-s+1,2:(m+1)])
  }
  imm_ind_mx[t,] <- colSums(imm_ind_mx_s, na.rm = T)
}

imm_ind_df <- data.frame(date = vaccinationTable$date,
                         imm_ind_mx)
colnames(imm_ind_df)[2:(m+1)] <- ag_labels

### Figure for immune individuals 
imm_ind_fig <- imm_ind_df %>% gather(2:(m+1), value = "protected", key = "age_group") 

imm_ind_fig %>% ggplot() +
  geom_line(aes(x=date, y=protected, color=age_group), size = 0.8) +
  scale_x_date(date_labels = "%m/%d/%Y", date_breaks = "2 weeks", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x="", y="Cumulative number of protected individuals", colour = "Age group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

### Figure for immune fraction
imm_frac_df <- data.frame(date = imm_ind_df$date,
                          a0009 = imm_ind_df$a0009 / N[1],
                          a1019 = imm_ind_df$a1019 / N[2],
                          a2029 = imm_ind_df$a2029 / N[3],
                          a3039 = imm_ind_df$a3039 / N[4],
                          a4049 = imm_ind_df$a4049 / N[5],
                          a5059 = imm_ind_df$a5059 / N[6],
                          a6069 = imm_ind_df$a6069 / N[7],
                          a7079 = imm_ind_df$a7079 / N[8],
                          a8089 = imm_ind_df$a8089 / N[9],
                          a90 = imm_ind_df$a90 / N[10])

#### Figure for immune fraction 
colnames(imm_frac_df)[2:ncol(imm_frac_df)] <- c("0–9","10–19","20–29","30–39","40–49","50–59","60–69","70–79","80–89","≥90")
imm_frac_fig <- imm_frac_df %>% gather(2:(m+1), value = "protected", key = "age_group") 

imm_frac_fig$age_group <- factor(imm_frac_fig$age_group, levels = unique(imm_frac_fig$age_group))

imm_frac_fig %>% ggplot() +
  geom_line(aes(x=date, y=protected, color=age_group), size = 0.8) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "2 weeks", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  labs(x="Calendar time in 2021", y="Protected fraction due to vaccination",
       colour = "Age group") +
  theme(text=element_text(size=12, family="sans",color="black"),
        axis.text=element_text(size=10, family="sans",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = c(0.13, 0.63))

################################################################################
################################################################################
# 2. Infected cases by age group from empirical data
## Infected cases estimated from backcalculation
#hersys_case_place = ".../File2.Hersys_infection.csv"
hersys_case_place = "./File2.Hersys_infection.csv"
estim_case <- as.data.frame(read_csv(file = hersys_case_place, show_col_types = FALSE))

estim_case$Infection <- as.Date(estim_case$Infection)
estim_case[,2:(1+m)] <- round(estim_case[,2:(1+m)])

period_cases <- estim_case %>% filter(between(Infection, vacDataStartDate, vacDataEndDate))

initial_cases <- estim_case %>% filter(Infection < vacDataStartDate)

## Figure for infected cases
period_cases %>% gather(2:(m+1), value = "cases", key = "age_group") %>% 
  ggplot() +
  geom_bar(aes(x=Infection, y=cases, fill = age_group), size = 0.8, stat = "identity") +
  scale_x_date(date_labels = "%m/%d", date_breaks = "2 weeks", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,25000))+
  labs(x="Calendar time in 2021", y="Infected cases with SARS-CoV-2", fill = "Age group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

#####################
# Defined parameters

## Generation time 
### Assuming Delta variant: Nishiura, et al, 2020
gi_fit = list(shape=2.305, scale=5.452)

generation <- function(d){
  x <- 1:d
  p <- pweibull(x+1, shape = gi_fit$shape, scale = gi_fit$scale) - pweibull(x, shape = gi_fit$shape, scale = gi_fit$scale)
  p <- p/sum(p)
  return(p)
}

## Nest-generation matrix
#NGM <- read.csv('.../File3.Osaka_4th_NGM.csv')
NGM <- read.csv('./File3.Osaka_4th_NGM.csv')

## Mobility information: google mobility data
#mobility_data0 <- read.csv('.../File4.Mobility_JP.csv')
mobility_data0 <- read.csv('./File4.Mobility_JP.csv')

smooth_window <- function(dat, win_size) {
  zoo::rollapply(dat,FUN=mean,width=win_size,align="center")
}

mobility_data <- data.frame(date=as.Date(mobility_data0$date),
                            household=(1+mobility_data0$residential_percent_change_from_baseline/100),
                            work=(1+mobility_data0$workplaces_percent_change_from_baseline/100),
                            community=(1+mobility_data0$retail_and_recreation_percent_change_from_baseline/100))

mobility_data$household[4:(length(mobility_data$household)-3)] <- smooth_window(mobility_data$household,7)
mobility_data$work[4:(length(mobility_data$work)-3)] <- smooth_window(mobility_data$work,7)
mobility_data$community[4:(length(mobility_data$community)-3)] <- smooth_window(mobility_data$community,7)

mobility_data %<>% dplyr::select(date,household,work,community) %>% filter(between(date,vacDataStartDate,vacDataEndDate))

## Relative transmissibility of delta variant
#delta_0 <- read.csv('.../File5.Delta_japan.csv')
delta_0 <- read.csv('./File5.Delta_japan.csv')

delta_df <- data.frame(delta_0)
delta_df$Date <- as.Date(delta_df$Date)
delta_df$Prop <- as.numeric(delta_df$Prop)

delta_delay <- 10 
delta_df2 <- delta_df %>% filter(Date <= as.Date("2021-08-23"))
delta_df2$Date <- delta_df2$Date - delta_delay
delta_duration <- as.numeric(max(delta_df2$Date) - min(delta_df2$Date)) + 1

obs_delta <- numeric(delta_duration)
for(w in 1:nrow(delta_df2)){
  obs_delta[w*7-6] <- delta_df2$Prop[w]
}


LH <- numeric(delta_duration)
nll_func_delta <- function(param){
  phi1 <- param[1]
  phi2 <- param[2]
  phi3 <- param[3]
  
  for(t in 1:delta_duration){
    
    vac <- phi1/(1+exp(-phi3*(t-phi2)))
    expt <- obs_delta[t] 
    
    LH[t] <- (expt-vac)^2
  }
  
  LH_2 <- numeric(nrow(delta_df2))
  for(w in 1:nrow(delta_df2)){
    LH[w*7-6] -> LH_2[w]
  }
  
  LLH <- log(LH_2)
  
  return(sum(LLH))
}
param0 <- c(99,45,0.07) 
outcome_delta <- optim(param0, nll_func_delta, method = "BFGS",　control=list(maxit=10000), hessian = T)

delta_start <- as.Date("2021-05-20")

delta_lasting <- as.numeric(vacDataEndDate - delta_start) + 1
new_delta <- numeric(delta_lasting)
for(t in 1:delta_lasting){
  new_delta[t] <- outcome_delta$par[1]/(1+exp(-outcome_delta$par[3]*(t-outcome_delta$par[2])))
}

new_delta_df <- data.frame(date=seq(from = delta_start,to = vacDataEndDate, by = "day"),
                           prop=new_delta)

### Transmissibility: maximum 1.5 times higher
new_delta_df_pre <- new_delta_df
new_delta_df_pre$prop <- new_delta_df_pre$prop * 50/100

new_delta_df_v2_pre <- new_delta_df_pre %>% mutate(prop2 = 1 + prop/100 )

#### from 2/17 ~ 5/19 
new_delta_df_add <- data.frame(date = seq(vacDataStartDate,(delta_start-1),by="day"),
                               prop = 0,
                               prop2 = 1)
#### from 2/17 ~ 11/30
new_delta_df_v2 <- bind_rows(new_delta_df_add,new_delta_df_v2_pre)

## Holiday effect
#renkyu_0 <- read.csv('.../File6.Renkyu_2021.csv')
renkyu_0 <- read.csv('./File6.Renkyu_2021.csv')

holiday_df <- data.frame(renkyu_0)
holiday_df$date <- holiday_df$date

################################################################################
################################################################################
# 3. Model

## Likelihood function
b <- 1:m
previous_expect <- bind_rows(initial_cases,period_cases)

LLH_function <- function(param){
  LLH_age <- numeric(m)
  
  ja <- data.frame(matrix(NA, ncol = m, nrow = t_duration))
  colnames(ja) <- colnames(period_cases)[2:(m+1)]
  
  LLH_df_a1 <- data.frame(day = 1:t_duration,
                          expect = NA,
                          estim = NA)
  LLH_df_a2 <- LLH_df_a1
  LLH_df_a3 <- LLH_df_a1
  LLH_df_a4 <- LLH_df_a1
  LLH_df_a5 <- LLH_df_a1
  LLH_df_a6 <- LLH_df_a1
  LLH_df_a7 <- LLH_df_a1
  LLH_df_a8 <- LLH_df_a1
  LLH_df_a9 <- LLH_df_a1
  LLH_df_a10 <- LLH_df_a1
  
  LLH_df_a1$expect <- period_cases$a0009[1:t_duration]* diag_bias
  LLH_df_a2$expect <- period_cases$a1019[1:t_duration]* diag_bias
  LLH_df_a3$expect <- period_cases$a2029[1:t_duration]* diag_bias
  LLH_df_a4$expect <- period_cases$a3039[1:t_duration]* diag_bias
  LLH_df_a5$expect <- period_cases$a4049[1:t_duration]* diag_bias
  LLH_df_a6$expect <- period_cases$a5059[1:t_duration]* diag_bias
  LLH_df_a7$expect <- period_cases$a6069[1:t_duration]* diag_bias
  LLH_df_a8$expect <- period_cases$a7079[1:t_duration]* diag_bias
  LLH_df_a9$expect <- period_cases$a8089[1:t_duration]* diag_bias
  LLH_df_a10$expect <- period_cases$a90[1:t_duration]* diag_bias
  
  for(a in 1:m){
    
    for(t in 1:t_duration){ 
      
      t2 <- t + nrow(initial_cases) 
      
      s_scaling <- param[1]
      
      h_google <- mobility_data$community[t] + mobility_data$household[t] * param[2] + mobility_data$work[t] * param[3]
      
      if(new_delta_df_v2$prop2[t] == 1){
        k_delta <- 1
      } else {
        k_delta <- new_delta_df_v2$prop2[t] * param[4]
      }
      
      if(holiday_df$obon[t] == 0){
        e_holidy <- 1
      } else {
        e_holidy <- holiday_df$obon[t] * param[5]
      }
      
      all_para <- s_scaling * h_google * e_holidy * k_delta
      
      jb = colSums(previous_expect[(t2-14):(t2-1),b+1]*diag_bias*generation(14)[14:1])
      ja[t,a] = sum(as.matrix(NGM)[a,] * (1 - (imm_frac_df[t,a+1] + sum(previous_expect[1:(t2-1),a+1]*diag_bias)/N[a])) * all_para * jb) 
      
      #jb = colSums(ja[(t2-14):(t2-1),b]*generation(14)[14:1]) 
      #ja[t,a] = sum(as.matrix(NGM)[a,] * (1 - (imm_frac_df[t,a+1] + sum(ja[1:(t2-1),a])/N[a])) * all_para * jb)
    }
  }
  ja[ja<=0] <- 1e-5
  
  LLH_df_a1$estim <- ja[,1]
  LLH_df_a2$estim <- ja[,2]
  LLH_df_a3$estim <- ja[,3]
  LLH_df_a4$estim <- ja[,4]
  LLH_df_a5$estim <- ja[,5]
  LLH_df_a6$estim <- ja[,6]
  LLH_df_a7$estim <- ja[,7]
  LLH_df_a8$estim <- ja[,8]
  LLH_df_a9$estim <- ja[,9]
  LLH_df_a10$estim <- ja[,10]
  
  LLH_age[1] <- sum(-LLH_df_a1$estim + LLH_df_a1$expect * log(LLH_df_a1$estim) - lgamma(LLH_df_a1$expect + 1))
  LLH_age[2] <- sum(-LLH_df_a2$estim + LLH_df_a2$expect * log(LLH_df_a2$estim) - lgamma(LLH_df_a2$expect + 1))
  LLH_age[3] <- sum(-LLH_df_a3$estim + LLH_df_a3$expect * log(LLH_df_a3$estim) - lgamma(LLH_df_a3$expect + 1))
  LLH_age[4] <- sum(-LLH_df_a4$estim + LLH_df_a4$expect * log(LLH_df_a4$estim) - lgamma(LLH_df_a4$expect + 1))
  LLH_age[5] <- sum(-LLH_df_a5$estim + LLH_df_a5$expect * log(LLH_df_a5$estim) - lgamma(LLH_df_a5$expect + 1))
  LLH_age[6] <- sum(-LLH_df_a6$estim + LLH_df_a6$expect * log(LLH_df_a6$estim) - lgamma(LLH_df_a6$expect + 1))
  LLH_age[7] <- sum(-LLH_df_a7$estim + LLH_df_a7$expect * log(LLH_df_a7$estim) - lgamma(LLH_df_a7$expect + 1))
  LLH_age[8] <- sum(-LLH_df_a8$estim + LLH_df_a8$expect * log(LLH_df_a8$estim) - lgamma(LLH_df_a8$expect + 1))
  LLH_age[9] <- sum(-LLH_df_a9$estim + LLH_df_a9$expect * log(LLH_df_a9$estim) - lgamma(LLH_df_a9$expect + 1))
  LLH_age[10] <- sum(-LLH_df_a10$estim + LLH_df_a10$expect * log(LLH_df_a10$estim) - lgamma(LLH_df_a10$expect + 1))
  
  return(-sum(LLH_age)) 
}

## MLE
start.time <- proc.time()

#param0 =c(2.7, -0.5,  -0.08,  1.1,  1.2)
param0 =c(2.35, -0.34,  -0.2,  1.0,  1.0)

outcome2 <- optim(param0, LLH_function, method = "L-BFGS-B", lower = c(0,-5,-5,0,0), upper=c(10,5,5,5,5), control=list(maxit=10000), hessian = T)

end.time <- proc.time()
(end.time - start.time)/60 

## Bootstrap for parameter
set.seed(10)
niter <- 1000

### Variance-covariance matrix from hessian
fisher_info = solve(outcome2$hessian)

### Sampling of Rt from multivariate normal distribution
hess_sample <- mvrnorm(n=niter, mu=outcome2$par,
                       Sigma=fisher_info, tol=1e-26,
                       empirical=FALSE, EISPACK=FALSE)

### confidence interval
c(outcome2$par[1],quantile(hess_sample[,1], probs = c(0.025,0.975))) %>% round(3)
c(outcome2$par[2],quantile(hess_sample[,2], probs = c(0.025,0.975))) %>% round(3)
c(outcome2$par[3],quantile(hess_sample[,3], probs = c(0.025,0.975))) %>% round(3)
c(outcome2$par[4],quantile(hess_sample[,4], probs = c(0.025,0.975))) %>% round(3)
c(outcome2$par[5],quantile(hess_sample[,5], probs = c(0.025,0.975))) %>% round(3)


################################################################################
################################################################################
# 4. Counterfactual cases
param_est <- outcome2$par

cases_cft_optim <- data.frame(matrix(rep(0,m*t_duration),ncol=m))
colnames(cases_cft_optim) <- colnames(initial_cases)[2:(m+1)]
cases_cft_optim <- cbind(data.frame(date = imm_frac_df$date[1:t_duration]),cases_cft_optim)

## All cases 
gt_ja_pre = matrix(rep(0,m*(t_duration+nrow(initial_cases))),ncol=m)
gt_ja_pre[1:nrow(initial_cases),b] <- as.matrix(initial_cases[,b+1]) * diag_bias

for(t in 1:t_duration){ 
  t2 <- t + nrow(initial_cases) 
  s_scaling <- param_est[1]
  
  h_google <- mobility_data$community[t] + mobility_data$household[t] * param_est[2] + mobility_data$work[t] * param_est[3]
  
  if(new_delta_df_v2$prop2[t] == 1){
    k_delta <- 1
  } else {
    k_delta <- new_delta_df_v2$prop2[t] * param_est[4] 
    #k_delta <- new_delta_df_v2$prop2[t]
    #k_delta <- 1
  }
  if(holiday_df$obon[t] == 0){
    e_holidy <- 1
  } else {
    e_holidy <- holiday_df$obon[t] * param_est[5]
  }

#  x <- runif(1, min = 0.01, max = 1.99)
#  x <- rnorm(1, mean = 1, sd = 1.37)
#  if(x < 0.1){
#    x = 0.1
#  }else if(x > 1.9){
#    x = 1.9
#  }  
#  all_para <- s_scaling * h_google * e_holidy * k_delta * x
  all_para <- s_scaling * h_google * e_holidy * k_delta
  
  for(a in 1:m){
    if(vac_program == 1){
      #jb = colSums(previous_expect[(t2-14):(t2-1),b+1]*diag_bias*generation(14)[14:1]) 
      #gt_ja_pre[t2,a] = sum(as.matrix(NGM)[a,] * (1 - (imm_frac_df[t,a+1] + sum(previous_expect[1:(t2-1),a+1]*diag_bias)/N[a])) * all_para * jb)
      jb = colSums(gt_ja_pre[(t2-14):(t2-1),b]*generation(14)[14:1]) 
      gt_ja_pre[t2,a] = sum(as.matrix(NGM)[a,] * (1 - (imm_frac_df[t,a+1] + sum(gt_ja_pre[1:(t2-1),a])/N[a])) * all_para * jb)
    } else {
      jb = colSums(gt_ja_pre[(t2-14):(t2-1),b]*generation(14)[14:1]) 
      gt_ja_pre[t2,a] = sum(as.matrix(NGM)[a,] * (1 - sum(gt_ja_pre[1:(t2-1),a])/N[a]) * all_para * jb)
    }
  }
  cases_cft_optim[t,b+1] <- gt_ja_pre[t2,]
}


### infected cases considering ascertainment bias
period_cases <- estim_case %>% filter(between(Infection, vacDataStartDate, vacDataEndDate))
period_cases[,b+1] <- period_cases[,b+1] * diag_bias

#### additional packages
library(foreach); library(doParallel);
(ncore <- detectCores(logical = FALSE))

## Calculation start
start.time <- proc.time()

myCluster <- makeCluster(ncore-1, type = "PSOCK")
registerDoParallel(myCluster)

## bootstrap for each epidemic curve
boot_result_list <- list()
boot_result_list <- foreach (i = 1:niter, .packages='dplyr') %dopar% {
  
  # sample parameters
  param_est <- hess_sample[i,]
  
  cases_cft <- data.frame(matrix(rep(0,m*t_duration),ncol=m))
  colnames(cases_cft) <- colnames(initial_cases)[2:(m+1)]
  cases_cft <- cbind(data.frame(date = imm_frac_df$date),cases_cft)
  
  # All cases
  gt_ja = matrix(rep(0,m*(t_duration+nrow(initial_cases))),ncol=m)
  gt_ja[1:nrow(initial_cases),b] <- as.matrix(initial_cases[,b+1]) * diag_bias
  
  set.seed(55)
  for(t in 1:t_duration){ 
    t2 <- t + nrow(initial_cases) 
    s_scaling <- param_est[1]
    
    h_google <- mobility_data$community[t] + mobility_data$household[t] * param_est[2] + mobility_data$work[t] * param_est[3]
    
    if(new_delta_df_v2$prop2[t] == 1){
      k_delta <- 1
    } else {
      k_delta <- new_delta_df_v2$prop2[t] * param_est[4]
      k_delta <- new_delta_df_v2$prop2[t]
      #k_delta <- 1
    }
    
    if(holiday_df$obon[t] == 0){
      e_holidy <- 1
    } else {
      e_holidy <- holiday_df$obon[t] * param_est[5]
    }
    
    all_para <- s_scaling * h_google * e_holidy * k_delta
    
    for(a in 1:m){
      if(vac_program == 1){
        #jb = colSums(previous_expect[(t2-14):(t2-1),b+1]*diag_bias*generation(14)[14:1]) 
        #gt_ja[t2,a] = rpois(1,sum(as.matrix(NGM)[a,] * (1 - (imm_frac_df[t,a+1] + sum(previous_expect[1:(t2-1),a+1]*diag_bias)/N[a])) * all_para * jb)) 
        jb = colSums(gt_ja_pre[(t2-14):(t2-1),b]*generation(14)[14:1]) 
        gt_ja[t2,a] = rpois(1,sum(as.matrix(NGM)[a,] * (1 - (imm_frac_df[t,a+1] + sum(gt_ja_pre[1:(t2-1),a])/N[a])) * all_para * jb))
      } else {
        jb = colSums(gt_ja[(t2-14):(t2-1),b]*generation(14)[14:1])
        gt_ja[t2,a] = rpois(1,sum(as.matrix(NGM)[a,] * (1 - sum(gt_ja[1:(t2-1),a])/N[a]) * all_para *jb)) 
      }
    }
    cases_cft[t,b+1] <- gt_ja[t2,]
  }
  boot_result_list[[i]] <- cases_cft
}
stopCluster(myCluster)

## Making dataframe for outcome
conf_df_total <- data.frame(date = imm_frac_df$date,
                            obs = rowSums(period_cases[,2:(m+1)]),
                            estim = NA,
                            lwr95 = NA,
                            upr95 = NA,
                            med50 = NA)
conf_base_df_age <- conf_df_total

## Estimated cases by age group
estim_case_agelist <- list()
for(a in 1:m){
  conf_base_df_age$obs <- as.numeric(unlist(period_cases[,a+1]))
  conf_base_df_age$estim <- as.numeric(unlist(cases_cft_optim[,a+1]))
  boot_mx <- matrix(NA, nrow = t_duration, ncol = niter)
  for(i in 1:niter){
    boot_mx[,i] <- boot_result_list[[i]][,a+1]
  }
  for(t in 1:t_duration){
    conf_base_df_age[t,4:6] <- quantile(boot_mx[t,], probs = c(0.025,0.975,0.5))
  }
  estim_case_agelist[[a]] <- conf_base_df_age
}

## total estimated cases
conf_df_total$estim <- rowSums(cases_cft_optim[,b+1])
boot_mx_total <-  matrix(NA, nrow = t_duration, ncol = niter)
for(i in 1:niter){
  boot_mx_total[,i] <- rowSums(boot_result_list[[i]][,b+1])
}
for(t in 1:t_duration){
  conf_df_total[t,4:6] <- quantile(boot_mx_total[t,], probs = c(0.025,0.975,0.5))
}

end.time <- proc.time()
(end.time - start.time)/60 

### Figure: total
conf_df_total %>% ggplot() +
  geom_point(aes(x=date,y=obs), colour = "#DA5D17", size = 1) +
  geom_line(aes(x=date,y=estim), colour="#418E4D", size = 1) +
#  geom_ribbon(aes(ymin = lwr95, ymax = upr95, x=date), fill = "#418E4D", alpha =0.3) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "2 weeks", expand = c(0,0)) +
  scale_y_continuous(limits = c(0,max(conf_df_total[,2:6])), expand = c(0,0)) +
  labs(x="Calendar time in 2021", y="Number of infected cases with SARS-CoV-2") +
  geom_hline(yintercept = 1,linetype="dotted",col="black") +
  theme(text=element_text(size=12, family="sans",color="black"),
        axis.text=element_text(size=10, family="sans",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
## Cumulative counterfactual cases
conf_df_total_withoutV <- conf_df_total
conf_cft_df_cum <- data.frame(date = imm_frac_df$date,
                              obs = NA,
                              estim = NA,
                              lwr95 = NA,
                              upr95 = NA,
                              med50 = NA)
conf_cft_df_cum$obs <- cumsum(conf_df_total$obs)
conf_cft_df_cum$estim <- cumsum(conf_df_total$estim)
total_boot_cum <- boot_mx_total
for(i in 1:niter){
  total_boot_cum[,i] <- cumsum(boot_mx_total[,i])
}
for(t in 1:t_duration){
  conf_cft_df_cum[t,4:6] <- quantile(total_boot_cum[t,], probs = c(0.025,0.975,0.5))
}


#### by age group
conf_cft_age_cum_list <- list()
for(a in 1:m){
  conf_cft_age_cum <- conf_cft_df_cum
  conf_cft_age_cum$obs <- cumsum(estim_case_agelist[[a]]$obs)
  conf_cft_age_cum$estim <- cumsum(estim_case_agelist[[a]]$estim)
  
  boot_mx <- matrix(NA, nrow = t_duration, ncol = niter)
  for(i in 1:niter){
    boot_mx[,i] <- cumsum(boot_result_list[[i]][,a+1])
  }
  for(t in 1:t_duration){
    conf_cft_age_cum[t,4:6] <- quantile(boot_mx[t,], probs = c(0.025,0.975,0.5))
  }
  conf_cft_age_cum <- bind_cols(conf_cft_age_cum,
                                data.frame(age = ag_labels[a]))
  conf_cft_age_cum_list[[a]] <- conf_cft_age_cum
}

##### figure: total
conf_cft_df_cum %>% ggplot() +
  geom_line(aes(x=date,y=obs), colour = "#DA5D17", size = 1) +
  geom_line(aes(x=date,y=estim), colour="#418E4D", size = 1) +
#  geom_ribbon(aes(ymin = lwr95, ymax = upr95, x=date), fill = "#418E4D", alpha =0.3) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "2 weeks", expand = c(0,0)) +
  scale_y_continuous(limits = c(0,max(conf_cft_df_cum[,2:6])), expand = c(0,0)) +
  labs(x="Calendar time", y="Cumulative infected cases") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

total_cum_infected <- sum(gt_ja_pre)
print(total_cum_infected)



