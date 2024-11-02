# clear the console; comment this if you do not use Rstudio
cat("\f") 
##' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##'
##'                      Panel Data
##'
##' NOTE: 
##' 1. No need to specify the path, this script will run as long as you put the 
##'    data in the same folder as the script.
##' 2. This script automatically installs and imports packages
##' 3. Text is commented with ##, alternative codes with #.
##'
##' INPUT:
##' 1. online data
##'
##' OUTPUT:
##' 1. 
##'
##' AUTHOR: Conghan Zheng (https://conghanzheng.github.io)
##' LAST UPDATED: 05 Oct 2024
##'
##' TO-DO:
##' - Translate Stata code to R code.
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

rm(list = ls()) ## clear the environments

## Auto install and import all packages listed in the vector:
if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c('rstudioapi', ## if you use Rstudio
              # 'this.path', ## if you don't use Rstudio
              'dplyr','haven','zoo','plm','stargazer','data.table','ggplot2',
              'patchwork', 'lfe', 'AER')

pacman::p_load(packages, character.only = TRUE)
invisible(lapply(packages, require, character.only = TRUE, quietly = TRUE))

## Set working directory to the script's location
scriptpath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
scriptpath %>% setwd()

options(dplyr.summarise.inform = FALSE) ## mute the summarise() info

## Start timer
time_start <- Sys.time()

cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] is running. \n ==== \n') ## Showing this in console helps when you are working on a project of many scripts, Rstudio required

data_url <- "https://raw.githubusercontent.com/conghanzheng/conghanzheng.github.io/master/assets/R/PanelData.dta"

download.file(data_url, destfile = "PanelData.dta")

##' PART I: MANIPULATING PANEL DATA ----

## Load data
data <- haven::read_dta("PanelData_3.dta") %>% data.table::setDT()

## Reshape: wide -> long
# Assuming the data is in wide format with variables like n1, n2, ..., y1, y2, ..., for years
# We will reshape it to long format
data_long <- data %>%
  tidyr::pivot_longer(
    cols = matches("^(n|w|k|y)[0-9]+$"),
    names_to = c(".value", "year"),
    names_pattern = "([a-z]+)([0-9]+)"
  ) %>%
  mutate(year = as.integer(year))

## Sort
data_long <- data_long %>%
  arrange(firm, year)

## Overview
summary(data_long)

## Drop missing (in all variables)
data_long <- data_long %>%
  filter(complete.cases(n, w, k, y))

## Exclude single observations (firms that only appear once in the panel)
data_long <- data_long %>%
  group_by(firm) %>%
  mutate(count = n()) %>%
  ungroup()

## Check
table(data_long$count)

## Drop single observations (if any)
data_long <- data_long %>%
  filter(count > 1) %>%
  select(-count)

## Exclude duplicates
data_long <- data_long %>%
  distinct(firm, year, .keep_all = TRUE)

## Generate order variable
data_long <- data_long %>%
  group_by(firm, year) %>%
  mutate(order = row_number()) %>%
  ungroup()

##' PART II: STATIC MODELS ----

## II.0 Panel Setting
# We will use the plm package for panel data models
# Create a pdata.frame
data_long <- pdata.frame(data_long, index = c("firm", "year"))

## II.1 FE

# Fixed Effects model
FE_model <- plm(n ~ k + w + y, data = data_long, model = "within", effect = "individual")
# Clustered standard errors
FE_vcov <- vcovHC(FE_model, method = "arellano", type = "HC0", cluster = "group")
FE_se <- sqrt(diag(FE_vcov))

# Fixed Effects model with time dummies
FE_twoway_model <- plm(n ~ k + w + y + factor(year), data = data_long, model = "within", effect = "individual")
FE_twoway_vcov <- vcovHC(FE_twoway_model, method = "arellano", type = "HC0", cluster = "group")
FE_twoway_se <- sqrt(diag(FE_twoway_vcov))

## II.2 Least Squares Dummy Variables (LSDV) Estimator

# OLS with firm dummies
OLS_model <- lm(n ~ k + w + y + factor(firm), data = data_long)
OLS_se <- sqrt(diag(vcovHC(OLS_model, type = "HC0")))

## Comparison
stargazer(OLS_model, FE_model, FE_twoway_model,
          type = "text",
          se = list(OLS_se, FE_se, FE_twoway_se),
          dep.var.labels = "n",
          covariate.labels = c("k", "w", "y"),
          keep = c("k", "w", "y"),
          model.names = FALSE,
          column.labels = c("OLS", "FE", "FE with Year Dummies"))

## II.3 Large number of fixed effects

# One fixed effect: firm
areg_model <- felm(n ~ k + w + y | firm, data = data_long)

# Extract the estimated fixed effects
fe_firm <- getfe(areg_model)

# Two fixed effects: firm and year
reghdfe_model <- felm(n ~ k + w + y | firm + year, data = data_long)

# Save fixed effects
fe_firm_year <- getfe(reghdfe_model)

## II.4 First-Differenced Least Squares (FDLS)

# First differences
data_long <- data_long %>%
  group_by(firm) %>%
  arrange(year) %>%
  mutate(
    D_n = n - lag(n),
    D_k = k - lag(k),
    D_w = w - lag(w),
    D_y = y - lag(y)
  ) %>%
  ungroup()

# FD model with intercept
FD_model <- lm(D_n ~ D_k + D_w + D_y, data = data_long)
FD_se <- sqrt(diag(vcovHC(FD_model, type = "HC0")))

# FD model without intercept
FD_nocons_model <- lm(D_n ~ D_k + D_w + D_y + 0, data = data_long)
FD_nocons_se <- sqrt(diag(vcovHC(FD_nocons_model, type = "HC0")))

## Compare with previous approaches
stargazer(FE_model, FD_model, FD_nocons_model,
          type = "text",
          se = list(FE_se, FD_se, FD_nocons_se),
          dep.var.labels = "Change in n",
          covariate.labels = c("k", "w", "y"),
          keep = c("k", "w", "y"),
          model.names = FALSE,
          column.labels = c("FE", "FD", "FD No Const"))

## II.5 Random Effects Model (FGLS estimator)

RE_model <- plm(n ~ k + w + y, data = data_long, model = "random")
RE_se <- sqrt(diag(vcovHC(RE_model, method = "arellano", type = "HC0", cluster = "group")))

## Compare with FE estimator
stargazer(FE_model, RE_model,
          type = "text",
          se = list(FE_se, RE_se),
          dep.var.labels = "n",
          covariate.labels = c("k", "w", "y"),
          keep = c("k", "w", "y"),
          model.names = FALSE,
          column.labels = c("FE", "RE"))

## II.6 FE or RE

# Hausman Test
hausman_test <- phtest(FE_model, RE_model)
print(hausman_test)

## II.7 Panel IV

# Fixed Effects IV
FE_IV_model <- plm(y ~ k + n | k + w, data = data_long, model = "within")
FE_IV_se <- sqrt(diag(vcovHC(FE_IV_model, method = "arellano", type = "HC0", cluster = "group")))

# First Difference IV
FD_IV_model <- plm(y ~ k + n | k + w, data = data_long, model = "fd")
FD_IV_se <- sqrt(diag(vcovHC(FD_IV_model, method = "arellano", type = "HC0", cluster = "group")))

# Random Effects IV
RE_IV_model <- plm(y ~ k + n | k + w, data = data_long, model = "random")
RE_IV_se <- sqrt(diag(vcovHC(RE_IV_model, method = "arellano", type = "HC0", cluster = "group")))

# Compare IV models
stargazer(FE_IV_model, FD_IV_model, RE_IV_model,
          type = "text",
          se = list(FE_IV_se, FD_IV_se, RE_IV_se),
          dep.var.labels = "y",
          covariate.labels = c("k", "n"),
          keep = c("k", "n"),
          model.names = FALSE,
          column.labels = c("FE IV", "FD IV", "RE IV"))

##' PART III: DYNAMIC MODELS ----

## III.1 Anderson and Hsiao (1981, 1982)

# Generate lags
data_long <- data_long %>%
  group_by(firm) %>%
  arrange(year) %>%
  mutate(
    nL1 = lag(n, 1),
    nL2 = lag(n, 2),
    kL1 = lag(k, 1),
    kL2 = lag(k, 2),
    wL1 = lag(w, 1),
    wL2 = lag(w, 2),
    yL1 = lag(y, 1),
    yL2 = lag(y, 2)
  ) %>%
  ungroup()

# First differences
data_long <- data_long %>%
  mutate(
    D_n = n - nL1,
    D_nL1 = nL1 - nL2,
    D_nL2 = nL2 - lag(nL2, 1),
    D_w = w - wL1,
    D_wL1 = wL1 - wL2,
    D_k = k - kL1,
    D_kL1 = kL1 - kL2,
    D_kL2 = kL2 - lag(kL2, 1),
    D_y = y - yL1,
    D_yL1 = yL1 - yL2,
    D_yL2 = yL2 - lag(yL2, 1)
  )

# Generate year dummies
data_long$year_factor <- as.factor(data_long$year)

# Anderson and Hsiao IV Regression
AH81_model <- ivreg(D_n ~ D_nL1 + D_nL2 + D_w + D_wL1 + D_k + D_kL1 + D_kL2 + D_y + D_yL1 + D_yL2 + year_factor | nL2 + D_nL2 + D_w + D_wL1 + D_k + D_kL1 + D_kL2 + D_y + D_yL1 + D_yL2 + year_factor, data = data_long)
summary(AH81_model)

## III.2 Arellano and Bond (1991): Difference GMM

# Note: This may take a while to run

AB_model <- pgmm(n ~ lag(n, 1) + lag(n, 2) + w + lag(w, 1) + lag(k, 0:2) + lag(y, 0:2) | lag(n, 2:99),
                 data = data_long, 
                 index = c("firm", "year"),
                 effect = "individual", 
                 model = "twosteps", 
                 transformation = "d",
                 collapse = TRUE)

summary(AB_model)

## III.3 Blundell and Bond (1998): System GMM

AB_Sys_model <- pgmm(n ~ lag(n, 1) + lag(n, 2) + w + lag(w, 1) + lag(k, 0:2) + lag(y, 0:2) | lag(n, 2:99),
                     data = data_long, 
                     index = c("firm", "year"),
                     effect = "individual", 
                     model = "twosteps", 
                     transformation = "ld",
                     collapse = TRUE)

summary(AB_Sys_model)

## Compare models
stargazer(AB_model, AB_Sys_model,
          type = "text",
          dep.var.labels = "n",
          covariate.labels = c("L1.n", "L2.n", "w", "L1.w", "k", "L1.k", "L2.k", "y", "L1.y", "L2.y"),
          model.names = FALSE,
          column.labels = c("Difference GMM", "System GMM"))

## Closing ----

## Ends timer
time_end <- Sys.time()
cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] takes',format(time_end - time_start),'\n ==== \n') ## Rstudio required