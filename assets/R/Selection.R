cat("\f")  ## Clear the console; comment this if you do not use RStudio
##' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##'
##'                 Censoring, Truncation, and Selection
##'
##' NOTE: 
##' 1. No need to specify the path; this script will run as long as you put the 
##'    data in the same folder as the script.
##' 2. This script automatically installs and imports packages.
##' 3. Text is commented with ##, alternative codes with #.
##'
##' INPUT:
##' 1. Selection.dta
##'
##' OUTPUT:
##' 1. Selection_results.RData
##'
##' AUTHOR: Conghan Zheng (https://conghanzheng.github.io)
##' LAST UPDATED: 16 Oct 2024
##'
##' TO-DO:
##' - Find a more efficient way of Selection FIML model (chunk the data / parallel computing)
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

rm(list = ls())  ## clear the environment

## Auto install and import all packages listed in the vector:
if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c('rstudioapi',  ## if you use RStudio
              # 'this.path',  ## if you don't use RStudio
              'dplyr','tidyr','haven','data.table','ggplot2',
              'stargazer','texreg','AER','sampleSelection','margins','censReg',
              'survival','future','caret')

pacman::p_load(packages, character.only = TRUE)
invisible(lapply(packages, require, character.only = TRUE, quietly = TRUE))

## Change the working directory to (a subdirectory under) that of the current script:
## - if you do not use RStudio as your IDE:
# scriptpath <- this.dir() %>% setwd()
## - if you are using RStudio:
scriptpath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(scriptpath)

options(dplyr.summarise.inform = FALSE)  ## mute the summarise() info

## Start timer
time_start <- Sys.time()

cat('\n ==== \n Script [', rstudioapi::getSourceEditorContext()$path, '] is running. \n ==== \n')

## Load data

data_url <- "https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/Selection.dta"

download.file(data_url, destfile = "Selection.dta")

data_raw <- haven::read_dta("Selection.dta") %>%
  data.table::setDT()

file.remove("Selection.dta")

## Label definitions
data <- data %>%
  mutate(
    educ = factor(educ,
                  levels = c(1, 2, 3, 4),
                  labels = c("Less than high-school graduated",
                             "Some years in college",
                             "College graduated",
                             "Higher education"))
  )

## Replace educ_sp = 0 if married == 0
data <- data %>%
  mutate(
    educ_sp = ifelse(married == 0, 0, educ_sp),
    educ_sp = factor(educ_sp,
                     levels = c(0, 1, 2, 3, 4),
                     labels = c("Non-married",
                                "SP: Less than high-school graduated",
                                "SP: Some years in college",
                                "SP: College graduated",
                                "SP: Higher education"))
  )

## Exercise 1: Tobit Model ----

## Log wages and dummy for positive wages
data <- data %>%
  mutate(
    lny = log(earn_h_r + 1e-12),
    dy = as.integer(earn_h_r > 0)
  )

## Censoring point
gamma <- min(data$lny, na.rm = TRUE)
cat("Censoring point of log(y) / gamma = ", gamma, "\n")
data %>% filter(lny < gamma + 0.01) %>% count(lny)
data %>% filter(earn_h_r < 0.02) %>% count(earn_h_r)

## Define factors
data <- data %>%
  mutate(
    cohort = factor(cohort),
    region = factor(region),
    educ = factor(educ),
    educ_sp = factor(educ_sp)
  )

## Tobit model for real hourly wages (left-censored at 0)

tobit_earn_h_r <- survreg(Surv(earn_h_r, earn_h_r > 0, type = "left") ~ married + educ * cohort + region * time + educ:cohort + educ:time + cohort:time + region:time + educ * time * cohort,
                          data = data,
                          dist = "gaussian",
                          weights = weight)

## Tobit model for log real hourly wages (left-censored at gamma)

tobit_lny <- survreg(Surv(lny, lny > gamma, type = "left") ~ married + educ * cohort + region * time + educ:cohort + educ:time + cohort:time + region:time + educ * time * cohort,
                     data = data,
                     dist = "gaussian",
                     weights = weight)

lltobit <- logLik(tobit_lny)

## Prediction
xb <- predict(tobit_lny, type = "linear")
sigma <- tobit_lny$scale

threshold <- (gamma - xb) / sigma
yhat_tobit <- exp(xb + 0.5 * sigma^2) * (1 - pnorm((gamma - xb - sigma^2) / sigma))
ytrunchat_tobit <- ifelse(data$dy == 1, yhat_tobit / (1 - pnorm(threshold)), NA)

## Compare models
screenreg(
  list(leveltobit = tobit_earn_h_r, logtobit = tobit_lny),
  stars = c(0.01, 0.05, 0.1),
  custom.model.names = c("Level Tobit", "Log Tobit"),
  digits = 3  # Increase the number of decimal places
)

### Exercise 1.2: Marginal effects on the left-censored mean of log wages ----

data_prl <- data

tobit_lny_AER <- tobit(lny ~ married + educ * cohort + region * time + educ:cohort + educ:time + cohort:time + region:time + educ * time * cohort, 
                       left = gamma, 
                       right = Inf, 
                       data = data_prl,
                       weights = data$weight)

mean_values <- data_prl[, .(
  married = mean(married, na.rm = TRUE),
  time = mean(time, na.rm = TRUE),
  region = names(sort(table(region), decreasing = TRUE)[1]),  # Most common region
  cohort = names(sort(table(cohort), decreasing = TRUE)[1]),  # Most common cohort
  educ = names(sort(table(educ), decreasing = TRUE)[1])  # Most common education level
)]

xb <- predict(tobit_lny_AER, type = "linear")
sigma <- tobit_lny_AER$scale
marginal_effect_time <- coef(tobit_lny_AER)["time"] * (1 - pnorm(-xb / sigma))

summary(marginal_effect_time)

### Exercise 1.3: Two-part model of log real hourly wages ----

## Part 1: Probit model
twopart_probit <- glm(dy ~ married + educ * cohort + region * time + educ:cohort + educ:time + cohort:time + region:time + educ * time * cohort,
                      data = data,
                      family = binomial(link = 'probit'),
                      weights = weight)
llprobit <- logLik(twopart_probit)

## Part 2: OLS regression for positive wages
twopart_ols <- lm(lny ~ married + educ * cohort + region * time + educ:cohort + educ:time + cohort:time + region:time + educ * time * cohort,
                  data = data[data$dy == 1, ],
                  weights = data$weight[data$dy == 1])

## Compute log-likelihood for OLS regression
residuals <- residuals(twopart_ols)
sigma2 <- sum(data$weight[data$dy == 1] * residuals^2) / sum(data$weight[data$dy == 1])
lllognormal <- sum(data$weight[data$dy == 1] * dnorm(residuals, mean = 0, sd = sqrt(sigma2), log = TRUE))

## Log-likelihood of the two-part model
lltwopart <- as.numeric(llprobit) + lllognormal

cat("The log-likelihood of the Tobit model is", lltobit, "\n")
cat("The log-likelihood of the two-part model is", lltwopart, "\n")

## Compare coefficients
screenreg(list(logtobit = tobit_lny, twopart_ols = twopart_ols), 
          stars = c(0.01, 0.05, 0.1),
          custom.model.names = c("Log Tobit", "Two-Part OLS"), 
          digits = 3)

## Exercise 2: Heckman Correction Method for Selection ----

### Exercise 2.1: Heckman two-step estimation without weights ----

heckman_liml <- selection(selection = dy ~ educ_sp + benefits + married + educ * cohort + region * time + educ:cohort + educ:time + cohort:time + region:time + educ * time * cohort,
                          outcome = lny ~ married + educ * cohort + region * time + educ:cohort + educ:time + cohort:time + region:time + educ * time * cohort,
                          data = data,
                          method = "2step")

summary(heckman_liml)

## Extract parameters from the Heckman model
sigma <- heckman_liml$sigma
rho <- heckman_liml$rho
sig2sq <- sigma^2
sig12sq <- (rho * sigma)^2

## Generate model matrices
x1b1 <- model.matrix(heckman_liml$probit) %*% coef(heckman_liml$probit)
x2b2 <- model.matrix(heckman_liml$lm) %*% coef(heckman_liml$lm)

x1b1_filtered <- x1b1[match(rownames(model.matrix(heckman_liml$lm)), rownames(model.matrix(heckman_liml$probit)))]

probpos <- pnorm(x1b1_filtered)
yhat_heckman_liml <- exp(x2b2 + 0.5 * sig2sq) * (1 - pnorm(-x1b1_filtered - sig12sq))
ytrunchat_heckman_liml <- yhat_heckman_liml / probpos

## Closing ----

## Output
# save(list = c('tobit_lny', 'twopart_probit', 'twopart_ols', 'heckman_liml', 'heckman_fiml', 'truncols'),
#      file = file.path(scriptpath,'Selection_results.RData'))  ## or you can save the whole image

## Ends timer
time_end <- Sys.time()
cat('\n ==== \n Script [', rstudioapi::getSourceEditorContext()$path, '] takes', format(time_end - time_start), '\n ==== \n')