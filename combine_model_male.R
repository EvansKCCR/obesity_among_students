library(tidyverse)
library(rjags)
library(rstan)
library(tidybayes)
library(cowplot)
library(patchwork)

theme_set(theme_classic())

# First indicates result of TPBF, second is result of RFM.
# NP <- Normal for BMI, obese for TPBF.
# PP <- obese for BMI, obese for TPBF
# NP <- obese for BMI, normal for TPBF.
# NP <- Normal for BMI, normal for TPBF.
# y<- c(PP, NP, PN, NN)

## Specify model parameters
male_abn <- list(y = c(66, 31, 41, 316), 
           N = 454,
           p_mu = 0.20, 
           p_kappa = 2,
           mu_se_1 = 0.955,
           mu_se_2 = 0.806, 
           mu_sp_1 = 0.889,
           mu_sp_2 = 0.889,
           
           sd_se_1 = 0.010,
           sd_se_2 = 0.019,
           sd_sp_1 = 0.015,
           sd_sp_2 = 0.015)

male_obesity <- list(y = c(25, 14, 12, 403), 
                 N = 454,
                 p_mu = 0.05, 
                 p_kappa = 2,
                 mu_se_1 = 0.818,
                 mu_se_2 = 0.818, 
                 mu_sp_1 = 0.932,
                 mu_sp_2 = 0.936,
                 
                 sd_se_1 = 0.014,
                 sd_se_2 = 0.014,
                 sd_sp_1 = 0.012,
                 sd_sp_2 = 0.011)

# To be used in plotting later
cells_male.ab <- tibble(i = 1:4,
                counts = male_abn$y,
                TPBF = c('Positive','Negative','Positive','Negative'),
                RFM= c('Positive','Positive','Negative','Negative'))

cells_male.ob <- tibble(i = 1:4,
                        counts = male_obesity$y,
                        TPBF = c('Positive','Negative','Positive','Negative'),
                        RFM= c('Positive','Positive','Negative','Negative'))

# Compile the model
model <- stan_model('combine_prevalence.stan')

# Fit the model
fit_male.ab <- sampling(model, data = male_abn, chains = 20, pars = c('SE1','SE2','SP1','SP2'), iter = 50000 , seed = 12345,include=F)

fit_male.ob <- sampling(model, data = male_obesity, chains = 20, pars = c('SE1','SE2','SP1','SP2'), iter = 50000, seed = 12345, include=F)

# Write to a little text file for easy reference later
sink('results/fit_male_ab.txt')
sink('results/fit_male_ob.txt')
print(fit_male.ab, digits=3, probs = c(0.025, 0.5, 0.975))
print(fit_male.ob, digits=3, probs = c(0.025, 0.5, 0.975))
sink()

# Make the posterior predictive check plot: Overweight/obesity
make_plot<-function(fit_male.ab,j){
  
  plotme<-fit_male.ab %>% 
    gather_draws(y_ppc[i]) %>% 
    left_join(cells_male.ab) %>% 
    filter(i==j) %>% 
    ggplot()+
    stat_histinterval(aes(.value),
                      outline_bars = T,
                      slab_color = 'gray45',
                      .width = c(0.8,0.95))+
    geom_point(aes(x = counts, y = 0.025), color = 'black', shape=25, fill = 'red')+
    labs(x='', y='')+
    theme(aspect.ratio = 1,
          plot.margin = unit(c(1,1,1,1)*0.05, "cm"),
          panel.grid.major = element_line())
  
  if(j==1){
    plotme<- plotme + 
      facet_grid(~TPBF, labeller = function(x) label_both(x,sep = ' ')) + 
      ylab('Frequency')
  }
  if(j==2){
    plotme<- plotme + facet_grid(RFM~TPBF, labeller = function(x) label_both(x,sep = ' '))
  }
  if(j==3){
    plotme<-plotme + labs(x='Cell Count', y='Frequency')
  }
  if(j==4){
    plotme<- plotme + 
      facet_grid(RFM~., labeller = function(x) label_both(x,sep = ' ')) + 
      xlab('Cell Count')
  }
  
  plotme
}

pp <- make_plot(fit_male.ab,1)
np <- make_plot(fit_male.ab, 2)
pn <- make_plot(fit_male.ab,3)
nn <- make_plot(fit_male.ab,4)

final <- (pp + np)/(pn + nn)

ggsave('plots/predictive_male_ab.png', height = 5, width = 5)

#---- Make 6x2 Figure For Paper ----

p = rstan::extract(fit_male.ab)

# Create subplots
tiff("plots/posterior_male_ab.tiff", units="in", width=8, height=5, res=300)

par(mfrow=c(2,3))

hist(p$se_1, 
     xlab = expression(delta[1]),
     main = 'Posterior TPBF Sensitivity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$sp_1, 
     xlab = expression(gamma[1]),
     main = 'Posterior TPBF Specificity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$p, 
     xlab = expression(pi),
     main = 'Posterior Prevalence',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)


hist(p$se_2, 
     xlab = expression(delta[2]),
     main = 'Posterior RFM Sensitivity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$sp_2, 
     xlab = expression(gamma[2]),
     main = 'Posterior RFM Specificity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

dev.off()



fit_male.ab %>% 
  spread_draws(p, se_2, sp_2,se_1, sp_1) %>% 
  mutate(p_survey = p*se_1 + (1-p)*(1-sp_1)) %>% 
  ggplot(aes(p))+
  geom_histogram()

## Posterior predictive plot: obesity

theme_set(theme_classic())

make_plot<-function(fit_male.ob,j){
  
  plotme<-fit_male.ob %>% 
    gather_draws(y_ppc[i]) %>% 
    left_join(cells_male.ob) %>% 
    filter(i==j) %>% 
    ggplot()+
    stat_histinterval(aes(.value),
                      outline_bars = T,
                      slab_color = 'gray45',
                      .width = c(0.8,0.95))+
    geom_point(aes(x = counts, y = 0.025), color = 'black', shape=25, fill = 'red')+
    labs(x='', y='')+
    theme(aspect.ratio = 1,
          plot.margin = unit(c(1,1,1,1)*0.05, "cm"),
          panel.grid.major = element_line())
  
  if(j==1){
    plotme<- plotme + 
      facet_grid(~TPBF, labeller = function(x) label_both(x,sep = ' ')) + 
      ylab('Frequency')
  }
  if(j==2){
    plotme<- plotme + facet_grid(RFM~TPBF, labeller = function(x) label_both(x,sep = ' '))
  }
  if(j==3){
    plotme<-plotme + labs(x='Cell Count', y='Frequency')
  }
  if(j==4){
    plotme<- plotme + 
      facet_grid(RFM~., labeller = function(x) label_both(x,sep = ' ')) + 
      xlab('Cell Count')
  }
  
  plotme
}

pp <- make_plot(fit_male.ob,1)
np <- make_plot(fit_male.ob, 2)
pn <- make_plot(fit_male.ob,3)
nn <- make_plot(fit_male.ob,4)

final <- (pp + np)/(pn + nn)

ggsave('plots/male_ob.png', height = 5, width = 5)

#---- Make 6x2 Figure For Paper ----

p = rstan::extract(fit_male.ob)

# Create subplots
tiff("plots/posterior_male_ob.tiff", units="in", width=8, height=5, res=300)

par(mfrow=c(2,3))

hist(p$se_1, 
     xlab = expression(delta[1]),
     main = 'Posterior TPBF Sensitivity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$sp_1, 
     xlab = expression(gamma[1]),
     main = 'Posterior TPBF Specificity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$p, 
     xlab = expression(pi),
     main = 'Posterior Prevalence',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)


hist(p$se_2, 
     xlab = expression(delta[2]),
     main = 'Posterior RFM Sensitivity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$sp_2, 
     xlab = expression(gamma[2]),
     main = 'Posterior RFM Specificity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

dev.off()


fit_male.ob %>% 
  spread_draws(p, se_2, sp_2,se_1, sp_1) %>% 
  mutate(p_survey = p*se_1 + (1-p)*(1-sp_1)) %>% 
  ggplot(aes(p))+
  geom_histogram()
