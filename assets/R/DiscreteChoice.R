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
##' 1. DiscreteChoice.dta
##'
##' OUTPUT:
##' 1. DiscreteChoice_results.RData
##'
##' AUTHOR: Conghan Zheng
##' LAST UPDATED: 10 Oct 2024 (https://conghanzheng.github.io)
##'
##' TO-DO:
##' - At-means marginal effects for 2.3
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

rm(list = ls())  ## clear the environments

## Auto install and import all packages listed in the vector:
if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c('rstudioapi',  # if you use RStudio
              # 'this.path',  # if you don't use RStudio
              'dplyr','tidyr','haven','data.table','ggplot2',
              'stargazer','margins','nnet','marginaleffects','mlogit')

pacman::p_load(packages, character.only = TRUE)
invisible(lapply(packages, require, character.only = TRUE, quietly = TRUE))

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

## Load data

data_url <- "https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/DiscreteChoice.dta"

download.file(data2_url, destfile = "DiscreteChoice.dta")

data <- haven::read_dta("DiscreteChoice.dta") %>%
  data.table::setDT()

## Exercise 1: Binary Choice Models ----

### Exercise 1.1: Logit, Probit, and Odds Ratio for Logit -----

## Generate variables
data <- data %>%
  mutate(age2 = age^2,
         same_ind = as.integer(ind_prior == ind_after))

## Logit
blogit <- glm(same_ind ~ exper + unempl_dur + age + age2 + female + ln_wage + veneto_resid,
                   data = data, family = binomial(link = 'logit'))

## Probit
bprobit <- glm(same_ind ~ exper + unempl_dur + age + age2 + female + ln_wage + veneto_resid,
                    data = data, family = binomial(link = 'probit'))

## Compare models
stargazer(blogit, bprobit,
          type = 'text',
          title = "Logit vs Probit Models",
          dep.var.caption = "Dependent variable: same_ind",
          dep.var.labels = c(""),
          model.numbers = FALSE,
          column.labels = c("Logit", "Probit"),
          covariate.labels = c("Experience", "Unemployment Duration", "Age", "Age squared", "Female", "Log(Wage)", "Veneto Resident"),
          omit.stat = c("aic"),
          single.row = TRUE)

## Report odds ratios for logit model
logit_or <- exp(coef(blogit))
logit_or_confint <- exp(confint(blogit))

odds_ratios <- cbind(Odds_Ratio = logit_or, logit_or_confint)
print(odds_ratios)

### Exercise 1.2: Marginal Effects ----

mean_values <- data.frame(t(colMeans(data[, c("exper", "unempl_dur", "age", "age2", "female", "ln_wage", "veneto_resid")], na.rm = TRUE)))

## Marginal effects at means for the logit model

logit_margins_atmeans <- margins(blogit, at = mean_values)

summary(logit_margins_atmeans)

## Marginal effects at means for the probit model

probit_margins_atmeans <- margins(bprobit, at = mean_values)

summary(probit_margins_atmeans)

## Exercise 2: Multinomial and Conditional Logit Models ----

### Exercise 2.1: Multinomial Logit ----

## Generate choice_ind
data <- data %>%
  mutate(choice_ind = NA_integer_,
         choice_ind = case_when(ind_after %in% c(1,2,3) ~ 1,
                                ind_after %in% c(4,5,6) ~ 2,
                                ind_after %in% c(7,8) ~ 3),
         choice_ind = factor(choice_ind, levels = c(1,2,3),
                             labels = c("Manufacturing", "Services", "Public sector")))

## Multinomial Logit
MNL <- nnet::multinom(choice_ind ~ exper + unempl_dur + age + age2 + female + ln_wage + veneto_resid,
                          data = data, base = 1)
summary(MNL)

## Extract coefficients and compute relative risk ratios
MNL_coefs <- coef(MNL)
MNL_rrr <- exp(MNL_coefs)
print(MNL_rrr)

## Exercise 2.2: At-Means Marginal Effects for MNL ----

MNL_me <- marginaleffects::slopes(MNL, newdata = datagrid()) %>%
  filter(term %in% c("exper", "unempl_dur", "female"))

## Exercise 2.3: Conditional Logit ----

alternatives_finer <- c("manuf_light", "manuf_heavy", "manuf_elecon",
                        "serv_sales", "serv_fin", "serv_other",
                        "public_adm", "public_hlth")

data <- haven::read_dta(file.path(scriptpath, "DiscreteChoice.dta")) %>% 
  data.table::setDT() %>%
  mutate(
    age2 = age^2,
    choice_ind = as.integer(ind_after)
  )

## Map choice_ind to alternative labels
alt_labels <- setNames(alternatives_finer, 1:8)
data <- data %>%
  mutate(
    alternative_chosen = alt_labels[as.character(choice_ind)]
  )

## Create alternative-specific choice dummies (d_)
for (alt in alternatives_finer) {
  data[[paste0("d_", alt)]] <- as.integer(data$choice_ind == match(alt, alternatives_finer))
}

## Lists of vars
d_vars <- grep("^d_", names(data), value = TRUE)
e_vars <- grep("^e_", names(data), value = TRUE)
indiv_vars <- c("unempl_dur", "age", "age2", "female", "ln_wage", "veneto_resid")

## The conditional logit long form data
data_long <- data %>%
  pivot_longer( ## Reshape longer
    cols = all_of(c(d_vars, e_vars)),
    names_to = c(".value", "alternative"),
    names_pattern = "(d|e)_(.*)"
  ) %>% 
  mutate(
    choice = d,
    alternative_num = match(alternative, alternatives_finer)
  )

data_clogit <- mlogit::mlogit.data(
  data_long,
  choice = "choice",
  shape = "long",
  alt.var = "alternative",
  id.var = "id"
)

## Model formula: choice ~ alternative-specific variables | individual-specific variables
clogit <- mlogit::mlogit(
  choice ~ e | unempl_dur + age + age2 + female + ln_wage + veneto_resid,
  data = data_clogit
)

summary(clogit)

## Compute Marginal Effects

beta_e <- coef(clogit)['e']
data_clogit_df <- model.frame(clogit)
index_data <- index(clogit$model)
data_clogit_df$chid <- index_data$chid
data_clogit_df$alternative <- index_data$alt

if (!'e' %in% names(data_clogit_df)) {
  stop("Variable 'e' is not found in the data frame.")
}

data_clogit_df$e <- as.numeric(data_clogit_df$e)
data_clogit_df$e[is.na(data_clogit_df$e)] <- 0

data_clogit_df <- data_clogit_df %>%
  mutate(V_ij = beta_e * e)

data_clogit_df <- data_clogit_df %>%
  mutate(
    exp_V_ij = exp(V_ij)
  ) %>%
  group_by(chid) %>%
  mutate(
    denominator = sum(exp_V_ij),
    P_ij = exp_V_ij / denominator
  ) %>%
  ungroup()

P_ij_df <- data_clogit_df %>%
  select(chid, alternative, P_ij)

me_data <- data_clogit_df %>%
  rename(e_alternative = alternative, P_ej = P_ij) %>%
  select(chid, e_alternative, P_ej)

me_data_full <- me_data %>%
  inner_join(P_ij_df, by = "chid") %>%
  rename(alternative = alternative, P_ik = P_ij)

me_data_full <- me_data_full %>%
  mutate(
    ME = ifelse(
      e_alternative == alternative,
      beta_e * P_ik * (1 - P_ik),
      -beta_e * P_ik * P_ej
    )
  )

## Average marginal effects
avg_me <- me_data_full %>%
  group_by(e_alternative, alternative) %>%
  summarise(
    Average_Marginal_Effect = mean(ME, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = alternative,
    values_from = Average_Marginal_Effect
  )

## Exercise 2.4: Nested Logit ----

nests <- list(
  Manufacturing = c("manuf_light", "manuf_heavy", "manuf_elecon"),
  Services = c("serv_sales", "serv_fin", "serv_other"),
  Public = c("public_adm", "public_hlth")
)

nested_clogit <- mlogit::mlogit(
  choice ~ e | unempl_dur + age + age2 + female + ln_wage + veneto_resid,
  data = data_clogit,
  nests = nests,
  reflevel = "manuf_light"
)

summary(nested_clogit)

## Closing ----

file.remove("DiscreteChoice.dta")

## Output
# save(list = c('blogit','bprobit','MNL','clogit','nlogit'),
#      file = file.path(scriptpath,'DiscreteChoice_results.RData'))  ## or you can save the whole image

## Ends timer
time_end <- Sys.time()
cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] takes',format(time_end - time_start),'\n ==== \n')  ## RStudio required