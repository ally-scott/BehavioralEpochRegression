cat("\f") #clear console
library(lmerTest)

library(tidyverse)
library(readxl)
library(emmeans)

library(dplyr)
library(readr)

setwd("G:/Shared drives/Richard Lab/Papers and Abstracts/Papers in Progress/Alexandra Scott/GADVPFP/Figs/fig4")

rm(list = ls())
binkernelcriteria465 <- read.csv("criteriaday465kernelbins_dff.csv")
binkernelcriteria465$rat<-gsub("_","",binkernelcriteria465$rat)
binkernelcriteria465['led']='465'

binkernelcriteria465  <-  binkernelcriteria465 %>%
  mutate(rat = parse_number(`rat`),
         event = `event`,
         binavg_b=`b_bins`,
         led=`led`) %>% 
  select(binavg_b,bin, rat, event,led)




binkernelcriteria<-binkernelcriteria465

ratinfo<- read.csv("ratinfo.csv")
names(ratinfo)[1]<- 'rat'
ratinfo  <-  ratinfo %>%
  mutate(rat = `rat`,
         sex = `sex`) %>% 
  select( rat,sex)


bins<- merge(binkernelcriteria, ratinfo, all = TRUE, by = c('rat'))

# bin_seq <- unique(bins$bin)
# summary(bins)
# ttest_fun <- function(t) {
#   browser()
#   df <- bins %>%
#     filter(bin == t)
#   test <- t.test(df$binavg_b)
#   p_value <- test$p.value
#   return(p_value)
# }
# 
# try <- lapply(bin_seq, FUN = ttest_fun(t = x))

# run linear mixed effect models
###
bin_model_DS <- lmer(binavg_b ~ (bin)*sex+ (1 | rat), data = bins%>% filter(event=="DS"))
bin_model_PE <- lmer(binavg_b ~ (bin) *sex + (1 | rat), data = bins%>% filter(event=="PE"))
bin_model_lox <- lmer(binavg_b ~ (bin) *sex + (1 | rat), data = bins%>% filter(event=="lox"))

bin_lme_DS<-summary(bin_model_DS)
bin_lme_PE<-summary(bin_model_PE)
bin_lme_lox<-summary(bin_model_lox)

capture.output(bin_lme_DS, file = "bin_lme_DS_dff.doc")
capture.output(bin_lme_PE, file = "bin_lme_PE_dff.doc")
capture.output(bin_lme_lox, file = "bin_lme_lox_dff.doc")

bin_lmeanova_DS<- anova(bin_model_DS)
bin_lmeanova_PE<- anova(bin_model_PE)
bin_lmeanova_lox<- anova(bin_model_lox)

capture.output(bin_lmeanova_DS, file = "bin_lmeanova_DS_dff.doc")
capture.output(bin_lmeanova_PE, file = "bin_lmeanova_PE_dff.doc")
capture.output(bin_lmeanova_lox, file = "bin_lmeanova_lox_dff.doc")




#capture.output(bin_lme, file = "binkernel_lme_results.doc")
#capture.output(bin_lmeanova, file = "binkernel_lmeanova_results.doc")




# Pairwise comparisons between time points at each group levels
# Paired t-test is used because we have repeated measures by time
## need  to write a for loop for each bin? 
stat.test <- bins %>%
  group_by(bins) %>%
  t_test(binavg_b ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test



stat.testDS <- t.test(binavg_b~bin,data=bins%>% filter(event=="DS"), paired=FALSE ) 
stat.testPE <- t.test(binavg_b~sex,data=bins%>% filter(event=="PE"), paired=FALSE) 
stat.testlox <- t.test(binavg_b~sex,data=bins%>% filter(event=="lox"), paired=FALSE ) 






library(ggplot2)


ggplot(bins, aes(x=bin, y=binavg_b,color=led)) + 
  stat_summary(fun.data = mean_se,geom = "ribbon")+      
  stat_summary(fun = "mean", geom = "line")+
  
  facet_grid(vars(sex), vars(event)) + 
  
  #stat_summary(fun.data = mean_se,geom = "ribbon")+      
  #stat_summary(fun = "mean", geom = "line")+
  
  
  #(aes(group=rat))+
  #stat_summary(fun.data = mean_se,geom = "area")+
 
  labs(title="time bins criteria day kernels",x ="event type", y = "b coeff from kernel ")  
ggsave("ribbon_bin criteria day kernels_dff.pdf")



ggplot(bins, aes(x=bin, y=binavg_b,color=led)) + 
  stat_summary(fun.data = mean_se,geom = "ribbon")+      
  stat_summary(fun = "mean", geom = "line")+
  
  facet_grid(vars(event)) + 
  
  #stat_summary(fun.data = mean_se,geom = "ribbon")+      
  #stat_summary(fun = "mean", geom = "line")+
  
  
  #(aes(group=rat))+
  #stat_summary(fun.data = mean_se,geom = "area")+
  
  labs(title="time bins criteria day kernels",x ="event type", y = "b coeff from kernel ")  
ggsave("ribbon_bin criteria day kernels_nosex_dff.pdf")
