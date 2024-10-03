** Problem Set 1: Panel Data
** Microeconometrics with Joan Llull, IDEA, FALL 2024
** TA: Conghan Zheng
**
** - Inputs:  PS1_1.dta
**            PS1_2.dta
** - Output: PS1_2.csv

cls, clear all
capture set more off

cd "/Users/zheng/Documents/02 IDEA_PhD/Teaching/TA_Microeconometrics_Fall_IDEA/2024/Part I/PS/PS1/code and data"

** Packages to be installed:	
** 	1. estout: Tools for making regression tables
* cap ssc install estout
** 	2. xtabond2: Arellano-Bond and Arellano-Bover estimators 
* cap ssc install xtabond2
** 	3. coefplot: Plot regression coefficients and other results
* cap ssc install coefplot 

** Exercise 1 ------------------------------------------------------------------
		
use "PS1_1.dta", clear

** generate a date variable combining year and month
gen time = ym(year, month)
format time %tm
order id year month time weight
	
** declare dataset as panel
xtset id time
	
** Exercise 1.1 ----
	
eststo clear
	
** RE
xtreg ln_wage h_immigr age age2 female black asian if h_skill==1, re
eststo RE

** FE
xtreg ln_wage h_immigr age age2 female black asian if h_skill==1, fe
eststo FE
	
** Compare
esttab RE FE, se ar2 replace nogaps compress mtitle

** Exercise 1.2 ----
	
** OLS
reg ln_wage h_immigr age age2 female black asian if h_skill==1
eststo OLS
	
** Compare
esttab RE FE OLS, se star stats(N r2 r2_a) nogaps compress mtitle
	
** Exercise 1.3 ----

** test FE versus RE
hausman FE RE

hausman FE RE, constant

hausman FE RE, constant sigmamore
	
** Exercise 1.4 ----

** Define job types 
gen job_type = 1 if routine==0 & cognitive==1 
replace job_type = 2 if routine==1 & cognitive==1 
replace job_type = 3 if routine==0 & cognitive==0 
replace job_type = 4 if routine==1 & cognitive==0 
	
** Labelling
label define labeljt 1 "NRC" 2 "RC" 3 "NRM" 4 "RM", replace
label values job_type labeljt

** Wages
gen wage = exp(ln_wage)
		
** Drop cases in which hours worked vary (check it by tabulating)
** After this we are left with individuals having hours worked above 0

drop if hours_worked==-8
	
** Hours_worked now are 'hours per week'. We need 'hours per month'
gen hours_worked_m = hours_worked*4.35
label var hours_worked_m "Hours worked per month"
	
** Dummy: Native
gen native = (immigr==0)
	
/* - Generating immigrant and native dummies before aggregating the data.

   The dummies are generated for the first entry of each (id,year) by-group, 
   they will work as indicators for that by-group when we later aggregate the 
   data.
*/
foreach var in immigr native h_immigr h_skill {
		bysort id year: gen `var'_d1=(`var'==1 & _n==1)
		
		** In case someone change occupation or location during the year:
		bysort id year: replace `var'_d1=1 ///
		if `var'==1 & _n>1 & (statefip!=statefip[_n-1] | job_type!=job_type[_n-1])
	}
	
/* - Aggregate data by job type, state and year.

   After collapse the data, we obtain sums of dummies, which are no longer
   binary. They are the number of according individuals in each (job,state,time)
   by-group, to be used for calculating the share later.
*/
collapse (sum) immigr_d1 native_d1 h_immigr_d1 h_skill_d1 ///
		 (mean) wage hours_worked_m [pw=weight], ///	
		 by(job_type statefip year)
	
order job_type statefip year
			
** Log wages
	gen ln_wage = ln(wage)	
	
/* Share of high-skilled defined by Equation (3). Here we omit suffix d1 because
   according to Stata's matching method, it will catch the most closed one in 
   existing variables. 
*/
gen phm = h_immig/(immigr+native)
gen ph = h_skill/(immigr+native)
	
** Job type, state and year interactions
egen jobxstate = group(job_type statefip) 
egen jobxyear = group(job_type year)
egen statexyear = group(statefip year)
	
/* Estimate the model without Alaska, Montana, and Wyoming, where there are no 
   high-skilled immigrants. 
*/

set matsize 1000
	
** 1st way of coding estimation
areg ln_wage phm ph i.job_type i.year i.job_type#i.statefip ///
	 i.job_type#i.year i.statefip#i.year ///
	 if !inlist(statefip, 2, 30, 56), absorb(statefip)
eststo m141

areg ln_wage phm i.job_type i.year i.job_type#i.statefip ///
	 i.job_type#i.year i.statefip#i.year ///
	 if !inlist(statefip, 2, 30, 56), absorb(statefip)
eststo m144

** 2nd way of coding
areg ln_wage phm ph i.job_type i.jobxyear i.jobxstate if !inlist(statefip, 2, 30, 56), absorb(statexyear)
estimates store m142

/* 3rd way of coding (we don't include job and state fixed effects because they 
are captured by the individual FE from id variable jobxstate) */
xtset jobxstate year

xtreg ln_wage phm ph i.year i.job_type#i.year i.statefip#i.year ///
      if !inlist(statefip, 2, 30, 56), fe
eststo m143
	
** Effect on hours worked
areg hours_worked_m phm ph i.job_type i.jobxyear i.jobxstate ///
     if !inlist(statefip, 2, 30, 56), absorb(statexyear)
eststo m145

areg hours_worked_m phm i.job_type i.jobxyear i.jobxstate ///
     if !inlist(statefip, 2, 30, 56), absorb(statexyear)
eststo m146

label var hours_worked_m "Hours"
label var ph "% H Workers"
label var phm "% H Migrants"

** Compare
esttab m141 m144 m145 m146, label se star stats(N r2 r2_a) nogaps compress keep(phm ph) depvar

** Exercise 2 ------------------------------------------------------------------
	
use "PS1_2.dta", clear

xtset id year
	
** Generate dummies for age categories
xi, noomit prefix("i_") gen i.agecat
	
** Generate dummies for interactions between married and age categories
xi, noomit prefix("i_") gen i.agecat*married
	
order id year wgt
	
** Exercise 2.1 ----
	
eststo clear
	
xtreg healthy i_agecat_27-i_agecat_62 i_ageXmarr_22-i_ageXmarr_62 college taxincome children* birthdum* [pw=wgt], fe
eststo WG21 

** Exercise 2.2 ----

/* - 'nodiffsargan' option in 'xtabond2' command
	
   All of the tests that xtabond2 would report are weak when the instrument 
   count is high.  Difference-in-Sargan/Hansen tests are are computationally 
   intensive since they involve re-estimating the model for each test; the 
   nodiffsargan option is available to prevent them.
   
   - the choice between one step, two-step, non-robust, and robust options
   
   In using xtabond2 command, for the choice between one step, two-step, 
   non-robust, and robust options, it doesn't matter if you understand what's 
   happening and could interpret it, because none of them is good enough. 
   
   Two-step robust is asymptotically more efficient than one-step robust, 
   especially for system GMM, but the reported std. err. tend to be severely 
   downward biased for two-step.
   
   In this exercise, you can use the default option and not specify any extra. 
   
   When you are doing your own research, for one-step robust estimation and for 
   all two-step estimations, the Hansen J statistics is reported by the command,
   it's the minimized two-step GMM criterion. Since Hansen J is robust but could
   be weakened by adding more instruments, you can also complement it with the 
   Sargan statistics, which is not robust to heteroskedasticity but doesn't 
   have the instruments problem.
*/

** On error: command xtabond2 is unrecognized
* cap ssc install xtabond2, replace
	
** Arellano-Bond
xtabond2 healthy l.healthy i_agecat_27-i_agecat_62 i_ageXmarr_22-i_ageXmarr_62 ///
		 college taxincome children* birthdum* [pw=wgt], nodiffsargan ///
		 gmmstyle(l.healthy) ///
		 ivstyle(i_agecat_27-i_agecat_62 i_ageXmarr_22-i_ageXmarr_62 college taxincome children* birthdum*) ///
		 noleveleq 
eststo ABond01	

** Four lags
xtabond2 healthy l.healthy i_agecat_27-i_agecat_62 i_ageXmarr_22-i_ageXmarr_62 ///
		 college taxincome children* birthdum* [pw=wgt], nodiffsargan ///
		 gmmstyle(l.healthy, lag(1 4)) ///
		 ivstyle(i_agecat_27-i_agecat_62 i_ageXmarr_22-i_ageXmarr_62 college taxincome children* birthdum*) ///
		 noleveleq 
eststo ABond02
	
** Arellano-Bover 
xtabond2 healthy l.healthy i_agecat_27-i_agecat_62 i_ageXmarr_22-i_ageXmarr_62 ///
		 college taxincome children* birthdum* [pw=wgt], nodiffsargan ///
		 gmmstyle(l.healthy) ///
		 ivstyle(i_agecat_27-i_agecat_62 i_ageXmarr_22-i_ageXmarr_62 college taxincome children* birthdum*)
eststo ABover01
	
** Compare
esttab ABond01 ABond02 ABover01, se star stats(N r2 r2_a) nogaps compress ///
	   keep(L.healthy i_agecat_* i_ageXmarr_*) mtitle

** Exercise 2.3
	
** Change labels of age-married interactions for the graphs
foreach ii in 22 27 32 37 42 47 52 57 62 {
		lab var i_ageXmarr_`ii' "`ii'"
}
	
** Make three plots and combine them
coefplot (WG21, label(WG)), ///
		 recast(line) lwidth(*2) ciopts(recast(rline) lpattern(dash)) ///
		 keep(i_ageXmarr_*) yline(0) vertical ///
		 xtitle("Age") ytitle("Probability gap married vs. single, {&beta}(a)") ///
		 title("{it: A. Fixed Effects}", color(black)) ///
		 name(plotwg, replace)
		
coefplot (ABond01, label(Differenced GMM)), ///
		 recast(line) lwidth(*2) ciopts(recast(rline) lpattern(dash)) ///
		 keep(i_ageXmarr_*) yline(0) vertical ///
		 xtitle("Age") ytitle("Probability gap married vs. single, {&beta}(a)") ///
		 title("{it: B. Difference GMM}", color(black)) ///
		 name(plotabond, replace)
		
coefplot (ABover01, label(System GMM)), ///
		 recast(line) lwidth(*2) ciopts(recast(rline) lpattern(dash)) ///
		 keep(i_ageXmarr_*) yline(0) vertical ///
		 xtitle("Age") ytitle("Probability gap married vs. single, {&beta}(a)") ///
		 title("{it: C. System GMM}", color(black)) ///
		 name(plotabover, replace)
		
graph combine plotwg plotabond plotabover, ycommon graphregion(color(white)) ///
	  rows(1) name(plotcombined, replace) 
graph export "PS1.png", as(png) replace 
