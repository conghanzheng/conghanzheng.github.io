rm(list = ls())
# cat("\f") ## clear the console; comment this line if you are not working with Rstudio
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
##              Automated Data Collection from eurostat open data
##
##
## DESCRIPTION: 
## Download selected business statistics from eurostat through API.
##
## INPUT:  
## 1. "manually input country names (iso code, in [country_chr]), time range ([data_start] and [data_end]), and names of variables to be requested ([dsid]) in the script"
##
## OUTPUT: 
## 1. data and TOC (csv, rds, and RData files by lines at the end of the script)
##
## LAST UPDATED: 26 APR 2023
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PRELIMINARIES ----

time_start <- Sys.time()

if (!require('pacman')) install.packages('pacman')
library(pacman)

packages <- c(# 'rstudioapi', ## if you use Rstudio
              'this.path', ## if you don't use Rstudio
              'tidyr','dplyr','zoo','countrycode','data.table','eurostat','stringr')

pacman::p_load(packages, character.only = TRUE)
invisible(lapply(packages, require, character.only = TRUE, quietly = TRUE))

## Change the work directory to (a subdirectory under) that of the current script:
## - if you do not use Rstudio as your IDE:
scriptpath <- this.dir()
## - if you are using Rstudio:
# scriptpath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path))

scriptpath %>% setwd()

inputpath  <- file.path(scriptpath, 'input')
outputpath <- file.path(scriptpath, 'output/eurostatB')
dir.create(outputpath,showWarnings = FALSE)

## list packages called in this script (package 'NCmisc' required)
if (!require('NCmisc')) install.packages('NCmisc')
library('NCmisc')
list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)

cat('\n ==== \n Script [',this.path::this.path(),'] is running. \n ==== \n')

## Set country and time ----

## time range of the current dataset: 
##              Euro area - 19 countries (Quarterly data, from 2015/2018); 
##              European Union - 27 countries (Monthly data, from 2020)

## start date
data_start <- format(as.Date('2000-01-01'),'%Y-%m')
## end date
data_end <- format(Sys.Date(), '%Y-%m')

## country code

country_chr <- c('AD','BA','CA') ## iso code

ALLcountrycode <- data.frame(iso2c = country_chr, 
                             name = countrycode(country_chr, origin = 'iso2c', destination = 'country.name'))

## Eurostat Business Statistics ----
## https://ec.europa.eu/eurostat/web/main/data/database
##
## - Production in industry: t_sts_ind_prod
## - Labor input in industry: t_sts_ind_labo
## - Construction, building, and civil engineering: t_sts_cons
## - Producer prices in industry: t_sts_ind_pric
## - Import prices in industry: t_sts_ind_impi

## eurostat json and unicode web services: 
## https://ec.europa.eu/eurostat/web/json-and-unicode-web-services

## Variables to be Requested

dsid <- c('teiis080','teiis090','teiis100','teiis110','teiis120','teiis130','teiis140', ## production in industry
          'teiis400', ## labor input in industry
          'teiis500','teiis550','teiis510', ## construction, building, and civil engineering
          'teiis010','teiis020','teiis030','teiis040','teiis050','teiis060','teiis070', ## producer prices in industry
          'teiis011','teiis012','teiis013','teiis014','teiis015','teiis016' ## import prices in industry
          )

## Prepare for Table of Contents

toc_raw <- eurostat::get_eurostat_toc() %>%
  filter(code %in% dsid,
         !grepl("Q",`data start`)) ## some of the requested variables are not available at the monthly basis

count <- 1

## Retrieve Data

for (i in toc_raw$code) {
  
  tmpdata <- eurostat::get_eurostat(id = i) %>%
    mutate(indicator = toc_raw$code[count],
           iso2c = countrycode(geo, origin = "iso2c", destination = "iso2c"),
           year_month = time,
           value = values) %>%
    group_by(iso2c,year_month) %>%
    arrange(unit, .by_group = T) %>%
    filter(complete.cases(iso2c), iso2c %in% ALLcountrycode$iso2c, complete.cases(values), row_number() == 1) %>%
    select(c('iso2c','year_month','indicator','value'))
  
  if (count == 1) {
    data_raw <- tmpdata
  } else {
    data_raw <- rbind(data_raw, tmpdata)
  }
  
  count <- count + 1
}

## Reshape

data_raw$year_month <- data_raw$year_month %>%
  str_replace('M','-') %>%
  as.yearmon() %>%
  format('%Y-%m')

data_wide <- data_raw %>%
  setDT %>% 
  dcast(iso2c + year_month ~ indicator, value.var = 'value')

## Make TOC ----

eurostat_toc <- toc_raw[!duplicated(toc_raw$code),c('code','title')] %>%
  as.data.frame()

colnames(eurostat_toc) <- c('indicator','description')

num = 2 ## how many columns are there in the data table before all the indicators;

eurostat_toc <- eurostat_toc[eurostat_toc$indicator %in% colnames(data_wide)[(num + 1):ncol(data_wide)],]

## Get the time range for each indicator

startchr <- character(ncol(data_wide) - num)
endchr <- character(ncol(data_wide) - num)

for (i in (num + 1):ncol(data_wide)) {
  tmpdata <- data.frame(data_wide[,1:num],data_wide[,..i]) %>%
    na.omit(cols = ncol(tmpdata))
  
  startchr[i - num] <- min(tmpdata$year_month)
  endchr[i - num] <- max(tmpdata$year_month)
}

eurostat_toc <- eurostat_toc %>%
  mutate(data_start = startchr,
         data_end = endchr,
         frequency = rep('Monthly',nrow(eurostat_toc)),
         source = rep('eurostatB',nrow(eurostat_toc)),
         link = rep('https://ec.europa.eu/eurostat/web/main/data/database',nrow(eurostat_toc)),
         subject = rep('Short-term business statistics',nrow(eurostat_toc)),
         type = rep(NA,nrow(eurostat_toc)),
         min = rep(NA,nrow(eurostat_toc)),
         max = rep(NA,nrow(eurostat_toc))) %>%
  setDT() %>%
  setcolorder(c('indicator','type','min','max','description','frequency','data_start','data_end','source','link','subject'))

## Variable type/min/max

for (i in 1:nrow(eurostat_toc)) {
  eurostat_toc$type[i] <- class(data_wide[[which(colnames(data_wide) == eurostat_toc$indicator[i])]])
  eurostat_toc$min[i] <- min(data_wide[[which(colnames(data_wide) == eurostat_toc$indicator[i])]], na.rm = T)
  eurostat_toc$max[i] <- max(data_wide[[which(colnames(data_wide) == eurostat_toc$indicator[i])]], na.rm = T)
}

num2 <- ncol(eurostat_toc)

## Count the number of observations of each variable

eurostat_toc$NumObs <- rep(0,nrow(eurostat_toc))

for (i in eurostat_toc$indicator) {
  tmpdata <- data_wide[, c('iso2c','year_month',i), with = F]
  tmpdata <- tmpdata[complete.cases(tmpdata),]
  
  tmpdata_o <- tmpdata[tmpdata$iso2c %in% origin_list,]
  tmpdata_d <- tmpdata[tmpdata$iso2c %in% destination_list,]
  
  eurostat_toc[which(eurostat_toc$indicator == i),'NumObs'] <- nrow(tmpdata)
  eurostat_toc[which(eurostat_toc$indicator == i),'NumObs_o'] <- nrow(tmpdata_o)
  eurostat_toc[which(eurostat_toc$indicator == i),'NumObs_d'] <- nrow(tmpdata_d)
}

## Add country coverage

country_coverage <- data.frame(indicator = names(data_wide)[(num + 1):ncol(data_wide)])

for (i in 1:nrow(country_coverage)) {
  tmpcol1 <- data_wide[, c(which(colnames(data_wide) == 'iso2c'),which(colnames(data_wide) == country_coverage$indicator[i])), with = F]
  tmpcol2 <- unique(tmpcol1[complete.cases(tmpcol1[,2]),1])
  
  for (j in tmpcol2$iso2c) {
    if (!j %in% colnames(country_coverage)) {
      country_coverage[,ncol(country_coverage) + 1] <- 0
      colnames(country_coverage)[ncol(country_coverage)] <- j
      country_coverage[which(country_coverage$indicator == country_coverage$indicator[i]),ncol(country_coverage)] <- 1
    } else {
      country_coverage[which(country_coverage$indicator == country_coverage$indicator[i]),which(colnames(country_coverage) == j)] <- 1
    }
  }
}

eurostat_toc <- merge(eurostat_toc,country_coverage,by = 'indicator')

num3 <- num2 + 3

eurostat_toc_tmp <- eurostat_toc

eurostat_toc_tmp_o <- eurostat_toc_tmp[,c(1:num3,which(colnames(eurostat_toc_tmp) %in% origin_list)), with = F]
eurostat_toc_tmp_o[,'NumOriginsCovered'] <- rowSums(eurostat_toc_tmp_o[,(num3 + 1):ncol(eurostat_toc_tmp_o)])

eurostat_toc_tmp_d <- eurostat_toc_tmp[,c(1:num3,which(colnames(eurostat_toc_tmp) %in% destination_list)), with = F]
eurostat_toc_tmp_d[,'NumDestinationsCovered'] <- rowSums(eurostat_toc_tmp_d[,(num3 + 1):ncol(eurostat_toc_tmp_d)])

eurostat_toc_tmp <- merge(eurostat_toc_tmp,eurostat_toc_tmp_o,by = colnames(eurostat_toc_tmp)[colnames(eurostat_toc_tmp) %in% colnames(eurostat_toc_tmp_o)])
eurostat_toc_tmp <- merge(eurostat_toc_tmp,eurostat_toc_tmp_d,by = colnames(eurostat_toc_tmp)[colnames(eurostat_toc_tmp) %in% colnames(eurostat_toc_tmp_d)])

ind_delete <- eurostat_toc_tmp[eurostat_toc_tmp$NumOriginsCovered == 0 & eurostat_toc_tmp$NumDestinationsCovered == 0,'indicator']

if (length(ind_delete) > 0) {
  data_wide <- data_wide[,colnames(data_wide)[!colnames(data_wide) %in% ind_delete],with = F]
  eurostat_toc <- eurostat_toc[!eurostat_toc$indicator %in% ind_delete,]
}

## OUTPUT ----

write.csv(data_wide,file.path(outputpath,'eurostatB.csv'))

#saveRDS(data_wide,file.path(outputpath,'eurostatB.rds'))
#saveRDS(eurostat_toc, file.path(outputpath,'eurostatB_toc.rds'))
#save(list = c('data_wide','eurostat_toc'), file = file.path(outputpath,'eurostatB.RData'))

time_end <- Sys.time()
cat('\n ==== \n Script [',this.path::this.path(),'] takes ',format(time_end - time_start),'\n ==== \n eurostat short-term business data is updated to ',max(data_wide$year_month, na.rm = T),'\n')