# cat("\f") ## clear the console; comment this if you do not use Rstudio
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
##' 1. PS1_1.dta
##' 2. PS1_2.dta
##'
##' OUTPUT:
##' 1. PS1_results.RData
##'
##' AUTHOR: Conghan Zheng (https://conghanzheng.github.io)
##' LAST UPDATED: 05 Oct 2024
##'
##' TO-DO:
##' - Find a more efficient package for the dynamic panel models.
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

rm(list = ls()) ## clear the environments

## Auto install and import all packages listed in the vector:
if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c('rstudioapi', ## if you use Rstudio
  # 'this.path', ## if you don't use Rstudio
  'dplyr','haven','zoo','plm','stargazer','data.table','ggplot2',
  'patchwork')

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

##' Exercise 1: Static Panel ----
##' Reference: Borjas (2003)

data1_raw <- haven::read_dta(file.path(scriptpath,"PS1_1.dta")) %>%
  data.table::setDT() ## things get quicker with data.table

## 1.1 - 1.2 FE vs RE vs OLS

data1 <- data1_raw %>%
  mutate(year_month = zoo::as.yearmon(paste0(year, month), "%Y %m"))

data1_h <- data1 %>% filter(h_skill == 1)

## Random effects (h_immigr == immigr for the high-skilled subsample)
RE <- plm(ln_wage  ~ h_immigr + age + age2 + female + black + asian, 
          data = data1_h,
          index = c("id","year_month"), 
          model = "random") 
summary(RE)

## Fixed effects
FE <- plm(ln_wage  ~ h_immigr + age + age2 + female + black + asian, 
          data = data1_h,
          index = c("id","year_month"), 
          model = "within") 
summary(FE)

## OLS
OLS <- lm(ln_wage ~ h_immigr + age + age2 + female + black + asian,  
          data = data1_h)
summary(OLS)

## Comparison

##' Note: 
##' 1. You can produce latex output by specifying "for (fmt in c('text','latex'))" in the for loop.
##' 2. "text" means display the table in the console

for (fmt in c('text')) {
  
  if (fmt == "latex") {
    tmppath <- file.path(scriptpath,"Ex1.tex")
  } else {
    tmppath <- NULL
  }
  
  stargazer(RE, FE, OLS,
            type = fmt,
            out = tmppath,
            title = paste0("Regression Results on Wage, Sample: Skilled Workers"),
            dep.var.caption = "\\textit{Dependent Variable}: Log(Wage)",
            dep.var.labels = c(""),
            order = paste0("^", c("h_immigr",'age','age2','female','black','asian'
            ) , "$"),
            covariate.labels = c("Migrant",'Age','Age-sq','Female','Black','Asian'
            ),
            star.cutoffs = c(0.05, 0.01, 0.001),
            #notes = c(""),
            #notes.append = FALSE,
            #notes.align = 'l',
            #notes.label = '',
            model.numbers = FALSE,
            keep.stat = c("n","rsq"),
            #column.sep.width = "0.5pt",
            align = FALSE,
            column.separate = c(1,1,1),
            column.labels = c("RE",'FE','OLS'),
            omit.table.layout = "m" ## Omit model-type labels
  )
}

## 1.3 Hausman Test

hausman_test <- phtest(FE, RE)

## 1.4 Lower Frequency Data by-group jobtype*locaiton*time with Multi-way FEs

data1_annual <- data1 %>%
  filter(hours_worked > 0) %>%
  mutate(jobtype = 1,
         jobtype = ifelse(routine == 1 & cognitive == 0, 2, jobtype),
         jobtype = ifelse(routine == 0 & cognitive == 1, 3, jobtype),
         jobtype = ifelse(routine == 1 & cognitive == 1, 4, jobtype),
         hours_worked = hours_worked * 4.348125, ## hours per week to hours per month
         wage = exp(ln_wage),
         ) %>%
  group_by(jobtype, statefip, year) %>%
  summarise(
    wage = weighted.mean(wage, weight),
    hours_worked = weighted.mean(hours_worked, weight),
    ph = sum(h_immigr, na.rm = TRUE) / n(),
    share_hskill = weighted.mean(h_skill, weight),
    share_female = weighted.mean(female, weight),
    share_black = weighted.mean(black, weight),
    share_asian = weighted.mean(asian, weight)
  ) %>%
  mutate(ln_wage = log(wage),
         job_state = interaction(jobtype, statefip),
         job_year = interaction(jobtype, year),
         state_year = interaction(statefip, year))

FE2_wage <- plm(ln_wage ~ ph + share_hskill + share_female + share_black + share_asian 
                + job_state + job_year + state_year,
                data = data1_annual,
                index = c("jobtype",'statefip',"year"),
                model = 'within')

FE2_hour <- plm(hours_worked ~ ph + share_hskill + share_female + share_black + share_asian 
                + job_state + job_year + state_year,
                data = data1_annual,
                index = c("jobtype",'statefip',"year"),
                model = 'within')

FE2_wage_nc <- plm(ln_wage ~ ph + job_state + job_year + state_year,,
                   data = data1_annual,
                   index = c("jobtype",'statefip',"year"),
                   model = 'within')

FE2_hour_nc <- plm(hours_worked ~ ph + job_state + job_year + state_year,,
                   data = data1_annual,
                   index = c("jobtype",'statefip',"year"),
                   model = 'within')

stargazer(FE2_wage, FE2_wage_nc, FE2_hour, FE2_hour_nc,
          type = 'text',
          out = tmppath,
          title = paste0("Regression Results on Yearly Wage"),
          dep.var.caption = "\\textit{Dependent Variable}:",
          dep.var.labels = c("Log(Wage)",'Hours Worked'),
          order = paste0("^", c('ph',"ph","share_hskill","share_female","share_black","share_asian"
          ) , "$"),
          covariate.labels = c("Share: H-skill migrant",'Share: H-skill','Share: Female','Share: Black','Share: Asian'
          ),
          star.cutoffs = c(0.05, 0.01, 0.001),
          #notes = c(""),
          #notes.append = FALSE,
          #notes.align = 'l',
          #notes.label = '',
          model.numbers = FALSE,
          keep.stat = c("n","rsq"),
          #column.sep.width = "0.5pt",
          align = FALSE,
          column.separate = c(1,1,1,1),
          column.labels = c("(1)","(2)","(3)","(4)"),
          omit = c("job_state", "job_year", "state_year"),
          add.lines = list(  # Manually add lines for the FEs
            c("FE(Job)", "Yes", "Yes", "Yes", "Yes"),
            c("FE(State)", "Yes", "Yes", "Yes", "Yes"),
            c("FE(Year)", "Yes", "Yes", "Yes", "Yes"),
            c("FE(Job X State)", "Yes", "Yes", "Yes", "Yes"),
            c("FE(Job X Year)", "Yes", "Yes", "Yes", "Yes"),
            c("FE(State X Year)", "Yes", "Yes", "Yes", "Yes")
          )
)

##' Exercise 2: Dynamic Panel ----
##' Reference: Guner et al. (2018)

data2_raw <- haven::read_dta(file.path(scriptpath,"PS1_2.dta")) %>%
  data.table::setDT()

data2 <- data2_raw %>%
  mutate(across(c("agecat",'married'), as.factor)) %>%
  bind_cols(model.matrix(~ agecat - 1, data = .)) %>%
  bind_cols(model.matrix(~ agecat:married - 1, data = .)) %>%
  select(id, year, wgt, everything())

colnames(data2) <- gsub(":", "_", colnames(data2))

agecat_vars <- grep("^agecat[0-9]+(.*[^0])$", names(data2), value = TRUE) ## regex used here
children_vars <- grep("^children", names(data2), value = TRUE)
birth_vars <- grep("^birthdum", names(data2), value = TRUE)

## 2.1 FE

WG21 <- plm(as.formula(paste("healthy ~", paste(c(agecat_vars, children_vars, birth_vars, "college", "taxincome"), collapse = " + "))),
            data = data2,
            index = c("id", "year"),
            weights = wgt,
            model = "within")

## 2.2 Arellano-Bond and Arellano-Bover

##' Note:
##' 1. Writing the formula manually without using `paste` makes it faster.
##' 2. Some regressors (child dummies, birth dummies) are omitted in favor of speed.
##' 3. The following code takes 8 mins on Mac OS with the M1 chip.

ABondL1 <- pgmm(healthy ~ lag(healthy, 1) + agecat27 + agecat32 + agecat37 + agecat42 + agecat47 + agecat52 + agecat57 + agecat62 + agecat22_married1 + agecat27_married1 + agecat32_married1 + agecat37_married1 + agecat42_married1 + agecat47_married1 + agecat52_married1 + agecat57_married1 + agecat62_married1 + college + taxincome| lag(healthy, 1),
                data = data2, 
                index = c('id','year'),
                effect = "individual", 
                model = "onestep", 
                transformation = "d", ## A-Bond
                collapse = TRUE, # nodiffsargan, collapse instruments
)

ABondL4 <- pgmm(healthy ~ lag(healthy, 1) + agecat27 + agecat32 + agecat37 + agecat42 + agecat47 + agecat52 + agecat57 + agecat62 + agecat22_married1 + agecat27_married1 + agecat32_married1 + agecat37_married1 + agecat42_married1 + agecat47_married1 + agecat52_married1 + agecat57_married1 + agecat62_married1 + college + taxincome| lag(healthy, 1:4),
                data = data2, 
                index = c('id','year'),
                effect = "individual", 
                model = "onestep", 
                transformation = "d", ## A-Bond
                collapse = TRUE, # nodiffsargan, collapse instruments
)

ABoverL1 <- pgmm(healthy ~ lag(healthy, 1) + agecat27 + agecat32 + agecat37 + agecat42 + agecat47 + agecat52 + agecat57 + agecat62 + agecat22_married1 + agecat27_married1 + agecat32_married1 + agecat37_married1 + agecat42_married1 + agecat47_married1 + agecat52_married1 + agecat57_married1 + agecat62_married1 + college + taxincome| lag(healthy, 1),
                data = data2, 
                index = c('id','year'),
                effect = "individual", 
                model = "onestep",
                transformation = "ld", ## A-Bover
                collapse = TRUE, # nodiffsargan, collapse instruments
)

## Table
stargazer(ABondL1, ABondL4, ABoverL1,
          type = "text",
          dep.var.caption = "Health Status",
          dep.var.labels = c(""),
          covariate.labels = c("L.health", 
                               "Age 25-29", "Age 30-34", "Age 35-39",
                               "Age 40-44", "Age 45-49", "Age 50-54", 
                               "Age 55-59", "Age 60-64",
                               "Age 20-24 X Married", "Age 25-29 X Married",
                               "Age 30-34 X Married", "Age 35-39 X Married", 
                               "Age 40-44 X Married", "Age 45-49 X Married", 
                               "Age 50-54 X Married", "Age 55-59 X Married", 
                               "Age 60-64 X Married", 
                               "College", "Taxable Income"),
          keep = c("L.healthy", "agecat*", "agecat*_married*"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          keep.stat = c("n","rsq"),
          digits = 3,
          model.numbers = FALSE,
          column.labels = c("ABond 1L", "ABond 4L", "ABoverL1 1L")
          )

## 2.3 Marriage Health Gap Plot

## Function: to extract coefficients and standard errors from plm and pgmm models
extract_coef <- function(model, pattern) {
  # Get summary of the model
  model_summary <- summary(model)
  
  # Extract coefficients
  coef_table <- model_summary$coefficients
  
  # Convert to data frame and filter by pattern (agecat and married)
  coef_df <- as.data.frame(coef_table, stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("term") %>%
    filter(grepl(pattern, term)) %>%
    select(term, estimate = Estimate, std.error = `Std. Error`)
  
  return(coef_df)
}

## Extract coef and se
coef_wg <- extract_coef(WG21, "agecat.*married") %>%
  mutate(age = as.numeric(gsub("agecat(\\d{2}).*", "\\1", term))) 

coef_abond <- extract_coef(ABondL1, "agecat.*married") %>%
  mutate(age = as.numeric(gsub("agecat(\\d{2}).*", "\\1", term)))

coef_ABoverL1 <- extract_coef(ABoverL1, "agecat.*married") %>%
  mutate(age = as.numeric(gsub("agecat(\\d{2}).*", "\\1", term))) %>%
  filter(age > 22)

## Set common y-axis limits based on the overall range of the data
y_min <- min(c(coef_wg$estimate - coef_wg$std.error,
               coef_abond$estimate - coef_abond$std.error,
               coef_ABoverL1$estimate - coef_ABoverL1$std.error))

y_max <- max(c(coef_wg$estimate + coef_wg$std.error,
               coef_abond$estimate + coef_abond$std.error,
               coef_ABoverL1$estimate + coef_ABoverL1$std.error))

y_limits <- c(y_min, y_max)  # Set dynamic y-axis limits

## Plot for WG21
plot_wg <- ggplot(coef_wg, aes(x = age, y = estimate)) +
  geom_line(size = 1.5) +  # Main line for the estimate
  geom_line(aes(y = estimate - std.error), linetype = "dashed") +  # Lower CI
  geom_line(aes(y = estimate + std.error), linetype = "dashed") +  # Upper CI
  geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = coef_wg$age) + # Set x-axis breaks to Age
  scale_y_continuous(limits = y_limits) +  # Set shared y-axis limits
  labs(x = "Age", y = "Probability gap married vs. single, β(a)", title = "A. Fixed Effects") +
  theme_minimal()

## Plot for ABondL1
plot_abond <- ggplot(coef_abond, aes(x = age, y = estimate)) +
  geom_line(size = 1.5) +  # Main line for the estimate
  geom_line(aes(y = estimate - std.error), linetype = "dashed") +  # Lower CI
  geom_line(aes(y = estimate + std.error), linetype = "dashed") +  # Upper CI
  geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = coef_abond$age) + # Set x-axis breaks to Age
  scale_y_continuous(limits = y_limits) +  # Set shared y-axis limits
  labs(x = "Age", y = "Probability gap married vs. single, β(a)", title = "B. Difference GMM") +
  theme_minimal()

## Plot for ABoverL1
plot_ABoverL1 <- ggplot(coef_ABoverL1, aes(x = age, y = estimate)) +
  geom_line(size = 1.5) +  # Main line for the estimate
  geom_line(aes(y = estimate - std.error), linetype = "dashed") +  # Lower CI
  geom_line(aes(y = estimate + std.error), linetype = "dashed") +  # Upper CI
  geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = coef_ABoverL1$age) + # Set x-axis breaks to Age
  scale_y_continuous(limits = y_limits) +  # Set shared y-axis limits
  labs(x = "Age", y = "Probability gap married vs. single, β(a)", title = "C. System GMM") +
  theme_minimal()

plot_mhgap <- plot_wg + plot_abond + plot_ABoverL1 + plot_layout(nrow = 1)
plot_mhgap

##' Note: You can also save the combined plot in a plot object, and then use the 
##' "ggsave" command to save it to the local disk.

## Closing ----

## Output
save(list = c('RE','FE','OLS','FE2_wage','FE2_wage_nc','FE2_hour','FE2_hour_nc','ABondL1','ABondL4','ABoverL1','plot_mhgap'), file = file.path(scriptpath,'PS1_results.RData')) ## or you can save the whole image

## Ends timer
time_end <- Sys.time()
cat('\n ==== \n Script [',rstudioapi::getSourceEditorContext()$path,'] takes',format(time_end - time_start),'\n ==== \n') ## Rtudio required