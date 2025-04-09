cat("\f")  ## Clear the console; comment this if you do not use RStudio
##' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##'
##'             Discrete Choice: Binary and Multinomial Models
##'
##' NOTE: 
##' 1. No need to specify the path; this script will run as long as you put the 
##'    data in the same folder as the script.
##' 2. This script automatically installs and imports packages.
##' 3. Text is commented with ##', alternative codes with #.
##'
##' INPUT:
##' 1. Online Data
##'
##' OUTPUT:
##' 1. 
##'
##' AUTHOR: Conghan Zheng
##' LAST UPDATED: 10 Oct 2024 (https://conghanzheng.github.io)
##'
##' TO-DO:
##' - Nested Logit: add regressors that affect the upper nest choice
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

rm(list = ls())  ## clear the environments

## Auto install and import all packages listed in the vector:
if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c('rstudioapi',  # if you use RStudio
              # 'this.path',  # if you don't use RStudio
              'dplyr','tidyr','stringr','haven','data.table','ggplot2','car',
              'stargazer','margins','nnet','marginaleffects','mlogit','effects')

pacman::p_load(packages, character.only = TRUE)
invisible(lapply(packages, require, character.only = TRUE, quietly = TRUE))

##' List packages called in this script (package 'NCmisc' required)
##' - This is used to check if some package dependencies have been forgotten.
# if (!require('NCmisc')) install.packages('NCmisc')
# library('NCmisc')
# list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)

## Change the work directory to (a subdirectory under) that of the current script:
## - if you do not use RStudio as your IDE:
# scriptpath <- this.dir() %>% setwd()
## - if you are using RStudio:
scriptpath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
scriptpath %>% setwd()

options(dplyr.summarise.inform = FALSE)  ## mute the summarise() info

## Start timer
time_start <- Sys.time()

cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] is running. \n ==== \n')  ## Showing this in console helps when you are working on a project of many scripts, RStudio required

## Load Data ----

data1_url <- "https://raw.githubusercontent.com/conghanzheng/conghanzheng.github.io/master/assets/R/DiscreteChoice.dta"

download.file(data1_url, destfile = "DiscreteChoice.dta")

data1 <- haven::read_dta("DiscreteChoice.dta") %>%
  data.table::setDT()

## BINARY OUTCOME MODELS ----

### LPM vs. Logit vs. Probit ----

## LPM (Linear Probability Model)
lpm <- lm(lfp ~ age + I(age^2) + married + educ + black + nchild + citiz, 
          data = data1)

## Logit
logit_mod <- glm(lfp ~ age + I(age^2) + married + educ + black + nchild + citiz, 
                 data = data1, 
                 family = binomial(link = "logit"))

## Odds ratios for logit
exp(coef(logit_mod))

## Probit
probit_mod <- glm(lfp ~ age + I(age^2) + married + educ + black + nchild + citiz, 
                  data = data1, family = binomial(link = "probit"))

## Compare
stargazer(lpm, logit_mod, probit_mod, 
          type = "text",
          keep.stat = c("n","rsq",'ll'),
          column.labels = c("LPM","Logit","Probit"),
          model.names = F)

##' For robust errors:
##' use packages 'sandwich', 'lmtest'
# coeftest(logit_mod, vcov = vcovHC(logit_mod, type = "HC1"))

### Marginal Effects ----

## Marginal effects at means (MEM)

data1 <- data1 %>%
  mutate(across(c(age, married, educ, black, nchild, citiz), as.numeric))

mean_vals <- data1 %>% 
  summarize(
    age = mean(age, na.rm = TRUE),
    age2 = mean(age^2, na.rm = TRUE),
    married = mean(married, na.rm = TRUE),
    educ = mean(educ, na.rm = TRUE),
    black = mean(black, na.rm = TRUE),
    nchild = mean(nchild, na.rm = TRUE),
    citiz = mean(citiz, na.rm = TRUE)
  )

mem_logit <- margins(logit_mod, at = as.list(mean_vals))
summary(mem_logit)

mem_probit <- margins(probit_mod, at = as.list(mean_vals))
summary(mem_probit)

## Marginal effects at a representative value (MER)

## At age = 20, age2 = 400, etc.
rep_vals <- data.frame(age = 20, age2 = 400, married = 1, educ = 73, black = 1, nchild = 2, citiz = 1)
mer_logit <- margins(logit_mod, at = rep_vals)
summary(mer_logit)

## Average marginal effects (AME)
ame_logit <- margins(logit_mod) 
summary(ame_logit)

ame_probit <- margins(probit_mod)
summary(ame_probit)

## Example: Margins by Education

## Factorization by education level
educ_lvl <- c("Less than high-school","High-school diploma",
              "Some college","College degree",
              "Master/Professional/PhD")

data1 <- data1 %>%
  mutate(educ_1 = case_when(
    educ == 73 ~ 1,
    educ == 81 ~ 2,
    educ >= 91 & educ <= 111 ~ 3,
    educ > 111 ~ 4,
    TRUE ~ 0
  ),
  educ_1 = factor(educ_1, levels = 0:4,
                  labels = educ_lvl))

logit_edu <- glm(lfp ~ age + I(age^2) + married + educ_1 + black + nchild + citiz,
                 data = data1, family = binomial("logit"))
mem_logit_edu <- margins(logit_edu, at = mean_vals)
summary(mem_logit_edu)

## Evaluate margins at specific factor levels
margins_edu_levels <- avg_slopes(logit_edu, variables = "educ_1", at = mean_vals) %>%
  ## Extract the first factor in the contrast (the "treatment" category)
  mutate(treatment = str_extract(contrast, "(?<=mean\\().+?(?=\\))")) %>%
  ## Convert to factor with desired order
  mutate(treatment = factor(treatment, levels = educ_lvl)) %>%
  ## Arrange by the treatment factor
  arrange(treatment)

margins_edu_levels

## Plotting predicted probabilities by education
preds_edu <- data.frame(educ_1 = levels(data1$educ_1))

pred_link <- predict(logit_edu, 
                     newdata = data.frame(age = mean_vals$age,
                                          age2 = mean_vals$age2,
                                          married = mean_vals$married,
                                          educ_1 = preds_edu$educ_1,
                                          black = mean_vals$black,
                                          nchild = mean_vals$nchild,
                                          citiz = mean_vals$citiz),
                     type = "link", se.fit = TRUE)

preds_edu <- preds_edu %>%
  mutate(pred = plogis(pred_link$fit),
         lower = plogis(pred_link$fit - 1.96 * pred_link$se.fit),
         upper = plogis(pred_link$fit + 1.96 * pred_link$se.fit))

## Plot with CI
ggplot(preds_edu, aes(x = sort(educ_1), y = pred)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  ylab("Prob.(lfp)") + xlab("Education Level") +
  ggtitle("Predicted Probability of lfp by Education") +
  theme_minimal()

### Goodness of fit ---

## Classification
cutoff <- 0.5
phat_logit <- predict(logit_mod, type = "response")
pred_logit <- ifelse(phat_logit > cutoff,1,0)
tab_logit <- table(Predicted = pred_logit, Actual = data1$lfp)

## Sensitivity & specificity 
sensitivity <- tab_logit[2,2]/(tab_logit[2,2] + tab_logit[1,2])
specificity <- tab_logit[1,1]/(tab_logit[1,1] + tab_logit[2,1])
c(sensitivity = sensitivity, specificity = specificity)

## Compare predictions (logit vs probit vs OLS)
phat_probit <- predict(probit_mod, type = "response")
phat_ols <- predict(lpm, type = "response")

## Plot predicted probabilities by age
age_sort <- data1 %>% group_by(age) %>%
  summarize(m_lfp = mean(lfp,na.rm = TRUE),
            m_logit = mean(phat_logit[age == .$age],na.rm = TRUE),
            m_probit = mean(phat_probit[age == .$age],na.rm = TRUE),
            m_ols = mean(phat_ols[age == .$age],na.rm = TRUE)) %>%
  ungroup()

ggplot(age_sort, aes(x = age)) +
  theme_minimal() +
  geom_line(aes(y = m_lfp, color = "Actual")) +
  geom_line(aes(y = m_ols, color = "OLS"), linetype = "longdash") +
  geom_line(aes(y = m_logit, color = "Logit"), linetype = "dashed") +
  geom_line(aes(y = m_probit, color = "Probit"), linetype = "dotdash") +
  labs(x = 'Age', y = "Prob.(lfp = 1)", title = "Actual vs Predicted Choice Prob. Across Models by Age", colour = NULL) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"))

### Specification Tests ----

data1_test <- data1 %>%
  mutate(
    age_married = age * married,
    age2_married = I(age^2) * married
  )

##' Wald test 
##' Testing if the coefficients of the interaction terms are jointly zero
logit_wald <- glm(lfp ~ age + I(age^2) + married + educ_1 + black + nchild + citiz + 
                    age_married + age2_married, 
                  data = data1_test, 
                  family = binomial(link = "logit"))

wald_test <- linearHypothesis(logit_wald, c("age_married = 0", "age2_married = 0"))
print(wald_test)

## Likelihood-ratio test

## Full model (with interactions)
logit_full <- glm(lfp ~ age + I(age^2) + married + educ_1 + black + nchild + citiz + 
                    age_married + age2_married, 
                  data = data1_test, 
                  family = binomial(link = "logit"))

## Restricted model (without interactions)
logit_restricted <- glm(lfp ~ age + I(age^2) + married + educ_1 + black + nchild + citiz,
                        data = data1_test, 
                        family = binomial(link = "logit"))

lr_test <- lrtest(logit_restricted, logit_full)
print(lr_test)

## Lagrange multiplier test

xbhat <- predict(logit_restricted, type = "link")

## Add squared linear predictor
data1_test <- data1_test %>%
  mutate(xbhat2 = xbhat^2)

## Fit model with squared linear predictor
logit_lm <- glm(lfp ~ age + I(age^2) + married + educ_1 + black + nchild + citiz + xbhat2,
                data = data1_test, 
                family = binomial(link = "logit"))

## Test significance of squared term
lm_test <- linearHypothesis(logit_lm, "xbhat2 = 0")
print(lm_test)

## MULTINOMIAL MODELS ----

data1_mnl <- data1 %>%
  mutate(
    sector = case_when(
      lfp == 0 ~ 0,
      classwkr %in% c(13, 14) ~ 1,
      classwkr %in% c(22, 23) ~ 2,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(sector)) %>%
  mutate(sector = factor(sector, 
                         levels = 0:2,
                         labels = c("Not participating", "Self-employed", "Private sector employee")
  ))

### Multinomial Logit ----

mnl_model <- multinom(sector ~ age + age2 + married + educ_1 + black + nchild + citiz,
                      data = data1_mnl, 
                      trace = FALSE)
mnl_coef <- coef(mnl_model)

## Get relative risk ratios
rrr <- exp(coef(mnl_model))

## Predicted probabilities
mnl_probs <- predict(mnl_model, type = "probs")
colnames(mnl_probs) <- paste0("pmlogit", 1:3)

## Add predictions to data
data1_mnl <- cbind(data1_mnl, mnl_probs)

n <- nrow(data1_mnl)

## Calculate marginal effects

calc_me_numeric <- function(var) {
  delta <- sd(data1_mnl[[var]]) / 1000
  temp_data <- data1_mnl
  temp_data[[var]] <- temp_data[[var]] + delta
  pred_delta <- predict(mnl_model, newdata = temp_data, type = "probs")
  me <- (pred_delta - mnl_probs) / delta
  me_means <- colMeans(me)
  me_ses <- apply(me, 2, sd) / sqrt(n)
  return(list(me = me_means, se = me_ses))
}

calc_me_binary <- function(var) {
  temp_data_1 <- temp_data_0 <- data1_mnl
  temp_data_1[[var]] <- 1
  temp_data_0[[var]] <- 0
  pred_1 <- predict(mnl_model, newdata = temp_data_1, type = "probs")
  pred_0 <- predict(mnl_model, newdata = temp_data_0, type = "probs")
  me <- pred_1 - pred_0
  me_means <- colMeans(me)
  me_ses <- apply(me, 2, sd) / sqrt(n)
  return(list(me = me_means, se = me_ses))
}

calc_me_factor <- function(var) {
  temp_data <- data1_mnl
  levels <- levels(data1_mnl[[var]])
  me_list <- list()
  
  for (lev in levels[-1]) {  # Compare to reference level
    temp_data_1 <- temp_data_0 <- temp_data
    temp_data_1[[var]] <- factor(lev, levels = levels)
    temp_data_0[[var]] <- factor(levels[1], levels = levels)
    pred_1 <- predict(mnl_model, newdata = temp_data_1, type = "probs")
    pred_0 <- predict(mnl_model, newdata = temp_data_0, type = "probs")
    me <- pred_1 - pred_0
    me_means <- colMeans(me)
    me_ses <- apply(me, 2, sd) / sqrt(n)
    me_list[[lev]] <- list(me = me_means, se = me_ses)
  }
  return(me_list)
}

numeric_vars <- c("age",'age2',"nchild")
binary_vars <- c("age",'age2',"nchild","married", "black", "citiz")
factor_vars <- c("educ_1")

mnl_me <- list()
for (var in numeric_vars) mnl_me[[var]] <- calc_me_numeric(var)
for (var in binary_vars) mnl_me[[var]] <- calc_me_binary(var)
for (var in factor_vars) mnl_me[[var]] <- calc_me_factor(var)

mnl_me <- mnl_me %>% as.data.frame()

### Conditional Logit ----

data("Fishing", package = "mlogit")
data2_cl <- dfidx(Fishing, varying = 2:9, shape = "wide", choice = "mode") %>%
  rename(p = price,
         c = catch)

## Basic model
cl1 <- mlogit(mode ~ p + c | 0, data = data2_cl)

## Model with income
cl2 <- mlogit(mode ~ p + c | income, data = data2_cl)

## Calculate ME
calc_cl_me <- function(model, covariate) {
  me <- effects(model, covariate = covariate)
  
  if(covariate == "income") {
    var_terms <- paste0("income:", c("boat", "charter", "pier"))
    var_me <- sum(diag(vcov(model)[var_terms, var_terms]))
  } else {
    var_me <- vcov(model)[covariate, covariate]
  }
  
  se <- sqrt(var_me)
  
  data.frame(
    Variable = covariate,
    ME = me[1],
    SE = se,
    t = me[1]/se, 
    p = 2*pnorm(-abs(me[1]/se))
  )
}

cl1_me_p <- calc_cl_me(cl1, "p")
cl1_me_c <- calc_cl_me(cl1, "c")
cl1_me_list <- list(cl1_me_p, cl1_me_c)
cl1_me_table <- do.call(rbind, cl1_me_list)

cl2_me_p <- calc_cl_me(cl2, "p")
cl2_me_c <- calc_cl_me(cl2, "c")
cl2_me_inc <- calc_cl_me(cl2, "income")
cl2_me_list <- list(cl2_me_p, cl2_me_c, cl2_me_inc)
cl2_me_table <- do.call(rbind, cl2_me_list)

### Nested Logit ----

data2_nl <- data2_cl

nl1 <- mlogit(mode ~ p + c | 0, 
              data2_nl,
              nests = list(coast = c('beach','pier'), 
                          water = c('charter','boat')))

nl2 <- mlogit(mode ~ p + c | income, 
              data2_nl,
              nests = list(coast = c('beach','pier'), 
                           water = c('charter','boat')))

## Compare CL and NL

stargazer(cl1, nl1, cl2, nl2,
          type = "text",
          keep.stat = c("n",'ll'),
          column.labels = c("CL1","NL1","CL2","NL2"),
          model.names = F)

## Closing ----

#file.remove("DiscreteChoice.dta")

## Output
# save(list = c(),
#      file = file.path(scriptpath,'DiscreteChoice_results.RData'))  ## or you can save the whole image

## Ends timer
time_end <- Sys.time()
cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] takes',format(time_end - time_start),'\n ==== \n')  ## RStudio required