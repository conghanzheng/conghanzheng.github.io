cat("\f")  ## Clear the console; comment this if you do not use RStudio
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
##                               DURATION ANALYSIS
##
## NOTE: 
## 1. No need to specify the path; this script will run as long as you put the 
##    data in the same folder as the script.
## 2. This script automatically installs and imports packages.
## 3. Text is commented with ##, alternative codes with #.
##
## INPUT:
## 1. Online Data
##
##' AUTHOR: Conghan Zheng (https://conghanzheng.github.io)
## LAST UPDATED: 02 Nov 2024
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

rm(list = ls())  ## Clear the environment

## Auto install and import all packages listed in the vector:
if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c('rstudioapi',  ## if you use RStudio
              # 'this.path',  ## if you don't use RStudio
              'dplyr','tidyr','haven','data.table','ggplot2','coxme','multiwayvcov',
              'survival','survminer','stargazer','texreg','reshape2','sandwich','lmtest')

pacman::p_load(packages, character.only = TRUE)
invisible(lapply(packages, require, character.only = TRUE, quietly = TRUE))

## Change the working directory to that of the current script:
## - if you do not use RStudio as your IDE:
# scriptpath <- this.path::this.dir() %>% setwd()
## - if you are using RStudio:
scriptpath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(scriptpath)

options(dplyr.summarise.inform = FALSE)  ## Mute the summarise() info

## Start timer
time_start <- Sys.time()

cat('\n ==== \n Script is running. \n ==== \n')

data1_url <- "https://raw.githubusercontent.com/conghanzheng/conghanzheng.github.io/master/assets/R/Duration_1.dta"
download.file(data1_url, destfile = "Duration_1.dta")

data2_url <- "https://raw.githubusercontent.com/conghanzheng/conghanzheng.github.io/master/assets/R/Duration_2.dta"
download.file(data2_url, destfile = "Duration_2.dta")

## Empirical Survival Function ----

data2 <- haven::read_dta("Duration_2.dta") %>% data.table::setDT()

## Round hours to one decimal point
data2 <- data2 %>%
  mutate(rounded_hours = round(hours, 1))

## Generate failure (the driver stops) variable
data2 <- data2 %>%
  mutate(exit = 1)

## Generate shift identifier
data2 <- data2 %>%
  mutate(shid = id * 100000 + date)

## Create survival object
surv_obj <- with(data2, Surv(time = rounded_hours, event = exit))

## Describe data
summary(surv_obj)

## Histogram of rounded_hours
ggplot(data2, aes(x = rounded_hours)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Rounded Hours", x = "Rounded Hours", y = "Frequency")

## Kaplan-Meier estimate
km_fit <- survfit(surv_obj ~ 1, data = data2)

## Plot Kaplan-Meier Survival Estimate
ggsurvplot(km_fit, data = data2, 
           xlab = "Time (Rounded Hours)", ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Estimate")

## Impact of High-Wage Shift on Working Hours ----

## Get the median wage
medwage <- median(data2$wage, na.rm = TRUE)

## Dummy for receiving a wage above the median
data2 <- data2 %>%
  mutate(high_wage = ifelse(!is.na(wage) & wage > medwage, 1, 0))

## The Cox PH model with one regressor
cox_model <- coxph(surv_obj ~ high_wage, data = data2)
summary(cox_model)

## Extract baseline hazard
basehaz_cox <- basehaz(cox_model, centered = FALSE)

## Extract beta coefficient
beta <- coef(cox_model)["high_wage"]

## Compute hazards for high_wage = 0 and high_wage = 1
hazard_df <- basehaz_cox
hazard_df$hazard0 <- hazard_df$hazard * exp(0)
hazard_df$hazard1 <- hazard_df$hazard * exp(beta)

## Reshape data for plotting
hazard_plot_df <- hazard_df %>%
  select(time, hazard0, hazard1) %>%
  melt(id.vars = "time", variable.name = "high_wage", value.name = "hazard")

## Plot hazard estimates
ggplot(hazard_plot_df, aes(x = time, y = hazard, color = high_wage)) +
  geom_line() +
  labs(title = "Hazard Estimates from Cox Model at Different Covariate Values",
       x = "Time", y = "Hazard", color = "High Wage") +
  scale_color_manual(labels = c("No", "Yes"), values = c("blue", "red")) +
  xlim(2, 12) +
  ylim(0, 1.5) +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  theme_minimal()

## Cox Regressions and Test of PH Assumption ----

## Model 1
cox_model1 <- coxph(surv_obj ~ high_wage + rainfall + snowfall, data = data2)
summary(cox_model1)

## Test the PH assumption
test_ph1 <- cox.zph(cox_model1)
print(test_ph1)
ggcoxzph(test_ph1)

## Model 2 with cluster(id)
cox_model2 <- coxph(surv_obj ~ high_wage + rainfall + snowfall + cluster(id), data = data2)
summary(cox_model2)

## Test the PH assumption
test_ph2 <- cox.zph(cox_model2)
print(test_ph2)
# Note: cox.zph may not work with cluster(), use caution interpreting results

## Compare coefficient estimates
screenreg(list(Model1 = cox_model1, Model2 = cox_model2),
          stars = c(0.01, 0.05, 0.1),
          custom.model.names = c("Model 1", "Model 2"),
          digits = 3)

## Plot to compare empirical survival function with the predicted one

## Kaplan-Meier estimate by high_wage
km_fit_hw <- survfit(surv_obj ~ high_wage, data = data2)

## Cox model predicted survival curves
new_data <- data.frame(high_wage = c(0, 1),
                       rainfall = mean(data2$rainfall, na.rm = TRUE),
                       snowfall = mean(data2$snowfall, na.rm = TRUE))

cox_fit_hw <- survfit(cox_model1, newdata = new_data)

## Plot comparison
ggsurvplot_list <- list(
  "Kaplan-Meier (No High Wage)" = survfit(surv_obj ~ 1, data = data2[data2$high_wage == 0, ]),
  "Kaplan-Meier (High Wage)" = survfit(surv_obj ~ 1, data = data2[data2$high_wage == 1, ]),
  "Cox Model (No High Wage)" = survfit(cox_model1, newdata = new_data[1, ]),
  "Cox Model (High Wage)" = survfit(cox_model1, newdata = new_data[2, ])
)

ggsurvplot_combine <- ggsurvplot(ggsurvplot_list,
                                 data = data2,
                                 combine = TRUE,
                                 legend.title = "Model",
                                 legend.labs = c("KM - No High Wage", "KM - High Wage", "Cox - No High Wage", "Cox - High Wage"),
                                 ggtheme = theme_minimal(),
                                 title = "Comparison of Empirical and Predicted Survival Functions by High Wage")

## Unobserved Heterogeneity of Drivers ----

## Cox model with frailty (shared frailty by id) using coxph
cox_model_re <- coxph(Surv(rounded_hours, exit) ~ high_wage + rainfall + snowfall + frailty(id), data = data2)
summary(cox_model_re)

## Extract frailty terms
frailty_terms <- cox_model_re$frail

## Ensure names are correctly assigned to frailty terms
# The names of frailty_terms should correspond to unique IDs
# If names are missing, assign them using unique IDs
if (is.null(names(frailty_terms))) {
  unique_ids <- as.character(unique(data2$id))
  names(frailty_terms) <- unique_ids
}

## Create the frailty data frame
frailty_df <- data.frame(id = as.numeric(names(frailty_terms)), lnnu = frailty_terms)

## Predicted frailties: low means that the chance of stopping is lower

## Identify drivers with highest and lowest frailty
max_frailty <- frailty_df[which.max(frailty_df$lnnu), ]
min_frailty <- frailty_df[which.min(frailty_df$lnnu), ]

cat("Driver with highest frailty (laziest):", max_frailty$id, "\n")
cat("Driver with lowest frailty (most diligent):", min_frailty$id, "\n")

## Compute baseline survival function
base_surv <- survfit(cox_model_re)
time <- base_surv$time
S0 <- base_surv$surv

## Compute survival functions for drivers
S_low <- S0 ^ exp(min_frailty$lnnu)
S_high <- S0 ^ exp(max_frailty$lnnu)

## Create data frame for plotting
surv_plot_df <- data.frame(
  time = time,
  S0 = S0,
  S_low = S_low,
  S_high = S_high
)

## Melt data for plotting
surv_melt <- melt(surv_plot_df, id.vars = "time", variable.name = "Survival", value.name = "Probability")

## Limit time to less than 200
surv_melt_filtered <- surv_melt %>% filter(time < 200)

## Plot survival functions
ggplot(surv_melt_filtered, aes(x = time, y = Probability, color = Survival)) +
  geom_line() +
  labs(title = "Survival Function for Drivers with Highest and Lowest Frailty",
       x = "Time", y = "Survival Probability") +
  theme_minimal()

## Empirical Hazard of Stopping by Hour and Income ----

data_ps4_1 <- haven::read_dta("Duration_1.dta") %>% data.table::setDT()

## Discretize total hours at the end of the trip using intervals
data_ps4_1 <- data_ps4_1 %>%
  mutate(
    int_hour = case_when(
      tot_hrs < 3 ~ 0,
      tot_hrs >= 3 & tot_hrs < 6 ~ 1,
      tot_hrs >= 6 & tot_hrs < 7 ~ 2,
      tot_hrs >= 7 & tot_hrs < 8 ~ 3,
      tot_hrs >= 8 & tot_hrs < 9 ~ 4,
      tot_hrs >= 9 & tot_hrs < 10 ~ 5,
      tot_hrs >= 10 & tot_hrs < 11 ~ 6,
      tot_hrs >= 11 & tot_hrs < 12 ~ 7,
      tot_hrs >= 12 ~ 8
    ),
    int_hour = factor(int_hour, levels = 0:8, labels = c("<3", "3-5", "6", "7", "8", "9", "10", "11", ">=12"))
  )

## Discretize total income at the end of the trip using intervals
data_ps4_1 <- data_ps4_1 %>%
  mutate(
    int_income = case_when(
      tot_inc < 25 ~ 0,
      tot_inc >= 25 & tot_inc < 50 ~ 1,
      tot_inc >= 50 & tot_inc < 75 ~ 2,
      tot_inc >= 75 & tot_inc < 100 ~ 3,
      tot_inc >= 100 & tot_inc < 125 ~ 4,
      tot_inc >= 125 & tot_inc < 150 ~ 5,
      tot_inc >= 150 & tot_inc < 175 ~ 6,
      tot_inc >= 175 & tot_inc < 200 ~ 7,
      tot_inc >= 200 & tot_inc < 225 ~ 8,
      tot_inc >= 225 ~ 9
    ),
    int_income = factor(int_income, levels = 0:9, labels = c("<25", "25-49", "50-74", "75-99", "100-124", "125-149", "150-174", "175-199", "200-224", ">=225"))
  )

## Hazards of stopping by hour
probit_model_hr <- glm(stop ~ int_hour, family = binomial(link = "probit"), data = data_ps4_1)

## Clustered standard errors by shid
cluster_se_hr <- coeftest(probit_model_hr, vcov = vcovCL(probit_model_hr, cluster = ~ shid))

## Predict probabilities
data_ps4_1$prhazard_hr <- predict(probit_model_hr, type = "response")

## Compute mean predicted probabilities by int_hour
mean_prhazard_hr <- data_ps4_1 %>%
  group_by(int_hour) %>%
  summarise(mean_prhazard_hr = mean(prhazard_hr, na.rm = TRUE))

print(mean_prhazard_hr)

## Hazards of stopping by income
probit_model_inc <- glm(stop ~ int_income, family = binomial(link = "probit"), data = data_ps4_1)

## Clustered standard errors by shid
cluster_se_inc <- coeftest(probit_model_inc, vcov = vcovCL(probit_model_inc, cluster = ~ shid))

## Predict probabilities
data_ps4_1$prhazard_inc <- predict(probit_model_inc, type = "response")

## Compute mean predicted probabilities by int_income
mean_prhazard_inc <- data_ps4_1 %>%
  group_by(int_income) %>%
  summarise(mean_prhazard_inc = mean(prhazard_inc, na.rm = TRUE))

print(mean_prhazard_inc)

## Logit Discrete Duration Model ----

## Load data
data_ps4_1 <- haven::read_dta("Duration_1.dta") %>% data.table::setDT()

## Day of the week
data_ps4_1 <- data_ps4_1 %>%
  mutate(
    weekday = weekdays(as.Date(date)),
    weekday = factor(weekday, levels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))
  )

## Discretize total hours for each trip
data_ps4_1 <- data_ps4_1 %>%
  mutate(
    int_hour_num = floor(tot_hrs),  ## Integer part of tot_hrs
    int_hour_num = ifelse(int_hour_num >= 12, 12, int_hour_num)  ## Cap at 12
  )

## Create labels without special characters
label_int_hour <- paste0("h", 0:12)
data_ps4_1$int_hour <- factor(data_ps4_1$int_hour_num, levels = 0:12, labels = label_int_hour)

## Collapse data
collapsed_data <- data_ps4_1 %>%
  group_by(shid, int_hour_num) %>%
  summarise(
    stop = sum(stop, na.rm = TRUE),
    id = mean(id, na.rm = TRUE),
    tot_hrs = max(tot_hrs, na.rm = TRUE),
    tot_inc = max(tot_inc, na.rm = TRUE),
    weekday = first(weekday),
    raining = max(raining, na.rm = TRUE),
    snowing = max(snowing, na.rm = TRUE),
    nonmanh = max(nonmanh, na.rm = TRUE),
    hour = min(hour, na.rm = TRUE),
    trip = n(),
    int_hour = first(int_hour)  # Include int_hour here
  ) %>%
  ungroup()

## Durations are censored at time = 12
collapsed_data <- collapsed_data %>%
  mutate(
    stop = ifelse(tot_hrs >= 12, 0, stop)
  )

## Change the scale of income
collapsed_data <- collapsed_data %>%
  mutate(
    tot_inc_100 = tot_inc / 100
  )

## Define duration structure of the data
collapsed_data <- collapsed_data %>%
  mutate(
    event = stop,  ## stop is the event indicator
    time = tot_hrs  ## tot_hrs is the time variable
  )

## Factor covariates
collapsed_data$hour <- factor(collapsed_data$hour)
collapsed_data$weekday <- factor(collapsed_data$weekday)
collapsed_data$int_hour <- factor(collapsed_data$int_hour)

## Prepare the formula
formula <- event ~ -1 + int_hour + tot_inc_100 + raining + snowing + hour + weekday + nonmanh

## Fit logistic regression model with cluster-robust standard errors
logit_model <- glm(formula, family = binomial(link = "logit"), data = collapsed_data)

## Compute clustered standard errors using vcovCL
cluster_se <- vcovCL(logit_model, cluster = collapsed_data$shid)

## Use coeftest to apply the clustered standard errors
logit_model_cl <- coeftest(logit_model, vcov = cluster_se)

## Extract coefficients for int_hour dummies
coefficients <- coef(logit_model)

## Get the coefficients for int_hour dummies
int_hour_coeffs <- coefficients[grep("^int_hour", names(coefficients))]

## Create data frame of int_hour_num and their variable names
coeff_df <- data.frame(
  int_hour_num = 0:12,
  coeff_name = paste0('int_hour', label_int_hour),
  stringsAsFactors = FALSE
)

## Map coefficients to int_hour_num
coeff_df$igamma <- int_hour_coeffs[coeff_df$coeff_name]

## Merge igamma into collapsed_data
collapsed_data <- collapsed_data %>%
  left_join(coeff_df[, c('int_hour_num', 'igamma')], by = 'int_hour_num')

## Compute hazard
collapsed_data$haz_free <- 1 / (1 + exp(-collapsed_data$igamma))

## Plot

## Remove observations with int_hour_num > 12 (if applicable)
plot_data <- collapsed_data %>%
  filter(as.numeric(as.character(int_hour_num)) <= 12)

## Aggregate haz_free by int_hour_num
plot_data_summary <- plot_data %>%
  group_by(int_hour_num) %>%
  summarise(haz_free = mean(haz_free, na.rm = TRUE))

## Convert int_hour_num to numeric for plotting
plot_data_summary$int_hour_num <- as.numeric(as.character(plot_data_summary$int_hour_num))

## Plot the hazard function using step function
ggplot(plot_data_summary, aes(x = int_hour_num, y = haz_free)) +
  geom_step(direction = "hv") +
  scale_x_continuous(breaks = seq(0, 12, 1), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x = "Time (Hours)", y = "Hazard", title = "Discrete Hazard Function") +
  theme_minimal()

## Closing ----

#file.remove("Duration_1.dta")
#file.remove("Duration_2.dta")

## Ends timer
time_end <- Sys.time()
cat('\n ==== \n Script takes', format(time_end - time_start), '\n ==== \n')