cat("\f") ## clear the console; comment this if you do not use Rstudio
##' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##'
##'                  Panel Data: Static and Dynamics Models
##'
##' NOTE: 
##' 1. No need to specify the path, this script will run as long as you put the 
##'    data in the same folder as the script.
##' 2. This script automatically installs and imports packages
##' 3. Text is commented with ##, alternative codes with #.
##'
##' INPUT:
##' 1. Online Data
##'
##' OUTPUT:
##'
##' AUTHOR: Conghan Zheng (https://conghanzheng.github.io)
##' LAST UPDATED: 05 Oct 2024
##'
##' TO-DO:
##' -
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

rm(list = ls()) ## clear the environments

## Auto install and import all packages listed in the vector:
if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c('rstudioapi', ## if you use Rstudio
  # 'this.path', ## if you don't use Rstudio
  'dplyr','plm','stargazer','haven','data.table','ggplot2',
  'tidyr','lmtest','fixest')

pacman::p_load(packages, character.only = TRUE)
invisible(lapply(packages, require, character.only = TRUE, quietly = TRUE))

##' List packages called in this script (package 'NCmisc' required)
##' - This is used to check if some package dependencies have been forgotten.
# if (!require('NCmisc')) install.packages('NCmisc')
# library('NCmisc')
# list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)

## Change the work directory to (a subdirectory under) that of the current script:
## - if you do not use Rstudio as your IDE:
# scriptpath <- this.dir() %>% setwd()
## - if you are using Rstudio:
scriptpath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
scriptpath %>% setwd()

options(dplyr.summarise.inform = FALSE) ## mute the summarise() info

## Start timer
time_start <- Sys.time()

cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] is running. \n ==== \n') ## Showing this in console helps when you are working on a project of many scripts, Rstudio required

##' TA ----

data_url <- "https://raw.githubusercontent.com/conghanzheng/conghanzheng.github.io/master/assets/R/PanelData.dta"

download.file(data_url, destfile = "PanelData.dta")

data_raw <- haven::read_dta("PanelData.dta") %>%
  data.table::setDT() ## things get quicker with data.table

## PART I: MANIPULATING PANEL DATA ----

## The data has variables n, w, k, y with repeated firm-year observations.

data_long <- data_raw %>%
  tidyr::pivot_longer(
    cols = -firm, 
    names_to = c(".value", "year"),
    names_pattern = "([nwky])([0-9]{4})"
  ) %>%
  mutate(year = as.integer(year)) %>%
  data.table::setDT() %>%
  arrange(firm, year) %>% ## order by firm and year
  filter(if_all(c(n, w, k, y), ~ !is.na(.))) %>% ## Drop missing (in all variables n, w, k, y)
  group_by(firm) %>%
  mutate(count_byfirm = n()) %>%
  ungroup() %>%
  filter(count_byfirm > 1) %>% ## Exclude single observations (firms that only appear once)
  select(-count_byfirm) %>%
  distinct(firm, year, .keep_all = TRUE) ## Exclude duplicates: keep first appearances per firm-year

## Overview
summary(data_long)

## Panel setting is done via plm by specifying index=c("firm","year").

##' Patterns in panel data
##' check T for each firm
plm::pdim(data_long, index = c("firm","year"))

## Within/Between Summary

## Function: within and between variation decomposition
xtsum <- function(data, id, time, vars) {
  results <- lapply(vars, function(var) {
    ## Filter out NA values for that variable
    df <- data %>% filter(!is.na(.data[[var]]))
    
    ## Overall mean and standard deviation
    overall_mean <- mean(df[[var]], na.rm = TRUE)
    overall_sd   <- sd(df[[var]], na.rm = TRUE)
    
    ## Between variation: mean by group (id)
    firm_means <- df %>%
      group_by(.data[[id]]) %>%
      summarize(mean_i = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
    
    between_sd <- sd(firm_means$mean_i, na.rm = TRUE)
    
    ##' Within variation:
    ##' Compute deviations from by-id means
    df_within <- df %>%
      left_join(firm_means, by = id) %>%
      mutate(diff_within = .data[[var]] - mean_i)
    
    within_var <- mean(df_within$diff_within^2, na.rm = TRUE)
    within_sd  <- sqrt(within_var)
    
    data.frame(
      variable = var,
      overall_mean = overall_mean,
      overall_sd   = overall_sd,
      between_sd   = between_sd,
      within_sd    = within_sd
    )
  })
  
  do.call(rbind, results)
}

## Within/Between Variation decomposition
variation_desomposition <- xtsum(data_long, id = "firm", time = "year", vars = c("n", "w", "k", "y"))
print(variation_desomposition)

## PART II: STATIC MODELS ----

## II.1 Fixed Effects (FE)
fe_model <- plm(n ~ k + w + y,
                data = data_long,
                index = c("firm","year"),
                model = "within") # FE model
summary(fe_model)

## Clustered standard errors by firm can be done using vcovHC:
fe_model_clse <- coeftest(fe_model, vcov = vcovHC(fe_model, method = "arellano", type = "HC1", cluster = "group"))

## Including year dummies (two-way FE):
fe_twoway_model <- plm(n ~ k + w + y + factor(year),
                       data = data_long,
                       index = c("firm","year"),
                       model = "within")
fe_twoway_clse <- coeftest(fe_twoway_model, vcov = vcovHC(fe_twoway_model, cluster = "group"))

## II.2 LSDV Estimator (include firm dummies manually)

## lm with factor(firm) gives the LSDV
ols_model <- lm(n ~ k + w + y + factor(firm), data = data_long)
summary(ols_model)

## Compare FE and OLS using stargazer:
stargazer(fe_model, ols_model, type = "text",
          column.labels = c("FE","LSDV"),
          model.names = F,
          keep = c('k','w','y'),
          keep.stat = c("n","rsq"),
          title = "Comparison of FE and OLS with firm dummies")

## II.3 Large number of fixed effects

## lm with factor(firm) or fixest
fe_fixest <- fixest::feols(n ~ k + w + y | firm, data = data_long)
summary(fe_fixest)

## Two-way FEs: firm and year
fe_fixest_2way <- fixest::feols(n ~ k + w + y | firm + year, data = data_long)
summary(fe_fixest_2way)

## II.4 First Differences (FD)

## Create first differences
data_long <- data_long %>%
  group_by(firm) %>%
  arrange(year) %>%
  mutate(Dn = n - lag(n),
         Dk = k - lag(k),
         Dw = w - lag(w),
         Dy = y - lag(y)) %>%
  ungroup()

fd_model <- lm(Dn ~ Dk + Dw + Dy, data = data_long)
summary(fd_model)

## Compare FE and FD
stargazer(fe_model, fd_model, keep.stat = c("n","rsq"), type = "text",
          model.names = F,
          column.labels = c("FE","FD"))

## II.5 Random Effects Model
re_model <- plm(n ~ k + w + y,
                data = data_long,
                index = c("firm","year"),
                model = "random")
summary(re_model)

## Compare FE and RE
stargazer(fe_model, re_model, type = "text",
          keep.stat = c("n","rsq"),
          column.labels = c("FE","RE"),
          model.names = F)

## II.6 FE or RE: Hausman test
hausman_test <- phtest(fe_model, re_model)
hausman_test

##' II.7 Panel IV
##' plm with Instrument variable approach: 'pvcm' or 'panelvar' or use iv within plm

##' If w is instrument for n, specify properly
##' Actually for IV: if n is endogenous and w is instrument:
##' formula: y ~ k + w + n | k + w + some_instrument_that_excludes n

##' No exact instruments:
##' y ~ k + n | k + w means w is excluded instrument for n
iv_fe <- plm(y ~ k + n | k + w,
             data = data_long,
             index = c("firm","year"),
             model = "within")
summary(iv_fe)

## FD IV
iv_fd <- plm(y ~ k + n | k + w,
             data = data_long,
             index = c("firm","year"),
             model = "fd")
summary(iv_fd)

## RE IV
iv_re <- plm(y ~ k + n | k + w,
             data = data_long,
             index = c("firm","year"),
             model = "random")
summary(iv_re)

stargazer(iv_fe, iv_fd, iv_re, type = "text", 
          keep.stat = c("n","rsq"),
          model.names = F,
          column.labels = c("FE-IV",'FD-IV','RE-IV'))

## PART III: DYNAMIC MODELS ----

## Lags
data_long <- data_long %>%
  group_by(firm) %>%
  arrange(year) %>%
  mutate(Ln = lag(n,1),
         Lk = lag(k,1),
         Lw = lag(w,1),
         Ly = lag(y,1),
         L2n = lag(n,2)) %>%
  ungroup()

##' Anderson and Hsiao (1981,1982):
##' 2SLS in differences: 
##' Δn_t  =  Δn_{t-1} instrumented by n_{t-2}
##' pgmm with transformation = "d"
AH_model <- pgmm(n ~ lag(n,1) + k + w + y | lag(n,2),
                 data = data_long,
                 index = c("firm","year"),
                 effect = "individual", model = "onestep", transformation = "d")
summary(AH_model)

##' Arellano and Bond (1991) - Difference GMM
##' pgmm with transformation = "d"
##' Suppose we want to use lagged values of n starting from 2 periods back as GMM-type instruments.

AB_1step <- pgmm(
  n ~ lag(n, 1) + k + w + y | lag(n, 2:99),
  data = data_long,
  index = c("firm","year"),
  effect = "individual",
  model = "onestep", ## the default
  transformation = "d",
  collapse = TRUE
)
summary(AB_1step)

AB_2step <- pgmm(
  n ~ lag(n, 1) + k + w + y | lag(n, 2:99),
  data = data_long,
  index = c("firm","year"),
  effect = "individual",
  model = "twosteps", 
  transformation = "d",
  collapse = TRUE
)
summary(AB_2step)

##' Blundell and Bond (1998) - System GMM
##' transformation = "ld" for system GMM in pgmm
BB_1step <- pgmm(
  n ~ lag(n, 1) + k + w + y | lag(n, 2:99),
  data = data_long,
  index = c("firm","year"),
  effect = "individual",
  model = "onestep", ## the default
  transformation = "ld",  # System GMM
  collapse = TRUE
)
summary(BB_1step)

BB_2step <- pgmm(
  n ~ lag(n, 1) + k + w + y | lag(n, 2:99),
  data = data_long,
  index = c("firm","year"),
  effect = "individual",
  model = "twosteps",
  transformation = "ld",  # System GMM
  collapse = TRUE
)
summary(BB_2step)

## Compare dynamic estimators:
stargazer(AH_model, AB_1step, AB_2step, BB_1step, BB_2step, type = "text",
          column.labels = c("A-H","A-Bond 1S","A-Bond 2S","A-Bover 1S","A-Bover 2S"),
          model.names = F)

## Closing ----

## Output
# save(list = c('ols_model','fe_model','re_model','fd_model','AH_model','AB_model','BB_model'), file = file.path(scriptpath,'PanelData_results.RData')) ## or you can save the whole image

## Ends timer
time_end <- Sys.time()

## End of script
cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] takes',format(time_end - time_start),'\n ==== \n') ## Rtudio required