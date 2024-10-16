** =============================================================================
** TA session 4 - CHAPTER 4: DURATION MODELS
** Microeconometrics, IDEA, FALL 2024
** TA: Conghan Zheng

/* Contents

I Continuous Duration
	
	I.0 Duration Data
	
	I.1 Nonparametric Approach 
	
	I.2 Semi-parametric Approach
		I.2.1 Cox's Proportional Hazard Model
		I.2.2 Model diagnostics
		I.2.3 Stratified analysis (optional)
		
	I.3 Parametric Approach

	I.4 Unobserved heterogeneity
	
	
II Discrete Duration  

	II.1 Logit

	II.2 Complementary Log Log 

	II.3 Unobserved Heterogeneity
	
*/

/* Dataset
	
	TA4.dta: college dropouts data 
*/	

cls, clear all
	
cd "/Users/zheng/Documents/02 IDEA_PhD/Teaching/TA_Microeconometrics_Fall_IDEA/2024/Part I/TA/TA4"
	
** Packages to be installed	

** estout - Making regression tables
* ssc install estout 
** dthaz - Estimates the hazard and survival probabilities of the population 
* ssc install dthaz
** pgmhaz8 - Estimate discrete time (grouped data) proportional hazard models
* ssc install pgmhaz8 


** PART I: CONTINUOUS DURATION -------------------------------------------------

** I.0 Duration Data ----

use "TA4.dta", clear
describe 
	
/* stset: to tell Stata we are using duration data

   We set the data as duration data based on variable 'span'. 
	
   Other options of 'stset':
	
	origin() defines analysis time, when a subject becomes at risk (here 
	date0). Without specifying 'origin()', observations before 1 January 1960 
	are ignored by Stata. 
	   
	time0() specifies the beginning of the period spanned by each record. 
	
	scale() helps to adapt to duration in weeks, months, years, etc. Specify 7 
    for weeks, 30 for months, 365.25 for years.
	
	exit() specifies the time when individuals exit the study. 'id()' helps in
	cases when we have more than one observation per individual (e.g., repeated 
	unemployment is possible in practice).
	
	enter() specifies the time when new individuals enter the study.
*/
stset duration, failure(event=1) id(id)
list id duration event _t0 _t _d _st if id<=20, noobs sepby(id)

** Describe duration data
stdescribe

/* Report whether values of variables within subject vary over time and reports 
   their pattern of missing values. A zero value in the varying column means 
   that within subjects, the variable doesn't vary. 
*/
stvary

** I.1 Nonparametric Approach ----

** Generate married dummy (1 if married, 0 otherwise) using time of marriage 
gen married = (marriage != 99)

** a) Kaplan-Meier estimator 
	
** Report the Kaplan-Meier survivor function
sts list
sts list, by(female) compare

** Graph the Kaplan-Meier survivor function
sts graph, survival ci risktable legend(pos(8) col(1) ring(0))
* graph export "ta4_np1.png", as(png) replace
		
** Graph by gender
sts graph, by(female) survival ci risktable legend(pos(8) col(1) ring(0))
* graph export "ta4_np2.png", as(png) replace
	
** b) Hazard function estimation
sts graph, hazard kernel(gaussian) width(4 5)
* graph export "ta4_np3.png", as(png) replace
	
** I.2 Semi-parametric Approach ----

** I.2.1 Cox's Proportional Hazard Model

global x female grade part_time lag stm married
	
/* Estimate a PH model to study time until college dropout. 

   There is no need to specify the outcome variable, because Stata knows it 
   after using 'stset'. 
   
   'nohr' option reports coefficients, not hazard ratios (proportional change in
   hazard when the regressor increases by one unit) 
*/
stcox $x, nohr
	
/*  Part-time students are 3.35 times more likely to drop out than 
	full-time students. 
	   
	Another example: starting college one month later multiplies the hazard rate
	by 1.11 (decreases chances of college dropout)
*/
** Part-time students are 3.35 times more likely to drop out than other students
display exp(_b[part_time])
/* Starting college one month later multiplies the hazard rate by 1.11 
   (decreases chances of college dropout).
*/
display exp(_b[stm])

stcox $x , vce(robust)
	
** Prediction
		
** Estimate a PH model with only one regressor
stcox part_time, vce(robust)
		
** Cumulative hazard 

** Baseline H(t, x = 0)
predict H0, basechazard
line H0 _t, c(J) sort 
* graph export "ta4_base1.png", as(png) replace
		
** Visualize how being part-time student influences the hazard
gen H1 = H0*exp(_b[part_time])
label variable H0 H0 
line H1 H0 _t, c(J J) sort
* graph export "ta4_base2.png", as(png) replace
		
** Survivor function

predict S0, basesurv	
line S0 _t, c(J) sort	
* graph export "ta4_base3.png", as(png) replace

gen S1 = S0^exp(_b[part_time])
label variable S0 S0
line S1 S0 _t, c(J J) sort
* graph export "ta4_base4.png", as(png) replace
		
** Hazard function h(t, x)

stcurve, hazard at(part_time=0)
* graph export "ta4_base5.png", as(png) replace
stcurve, hazard at1(part_time=0) at2(part_time=1)
* graph export "ta4_base6.png", as(png) replace

** Using a different kernel smoothing function:
stcurve, hazard at1(part_time=0) at2(part_time=1) kernel(gaussian) width(4)
* graph export "ta4_base7.png", as(png) replace

/* Difference between commands 'sts graph' and 'stvurve': 

Command sts graph, hazard by() plots the the non-parametric estimates (the product-limit), and the by() option means the estimations are done on two subsamples.

Command stcurve is used after command stcox, streg, etc.  Options at1(), at2() means the semi-parametrically estimated function is plotted for different values of the covariate.
*/

/* I.2.2 Model diagnostics: 

    Test the proportional hazard assumption
	  
	Link for Manual: 
	https://www.stata.com/manuals13/ststcoxph-assumptiontests.pdf
*/

** College dropout data
use "TA4.dta", clear
stset duration, failure(event) id(id)
gen married = (marriage != 99)
global x female grade part_time lag stm married

stcox $x , vce(robust)
		
** 1. Statistical test: based on scaled Schoenfeld residuals 
estat phtest, detail
		
** 2. Graphical test 1: Log–log plot
stphplot, by(part_time)
*graph export "ta4_log_log_prt.png", as(png) replace
stphplot, by(female)
*graph export "ta4_log_log_fml.png", as(png) replace
		
/* 3. Graphical test 2：

   The fitted survivor function from Cox regression and the (nonregression) 
   Kaplan-Meier estimate of the survivor function should be similar if PH holds.
*/
stcoxkm, by(part_time) // separate legend(cols(1)) 
*graph export "ta4_stcoxkm_prt.png", as(png) replace
stcoxkm, by(female) // separate legend(cols(1)) 
*graph export "ta4_stcoxkm_fml.png", as(png) replace

/* I.2.3 Stratified analysis (optional) ----

   Assumption: The baseline hazard h{0i}(t) can be group-specific, while the 
   coefficients from the exp(x*beta) are the same for all individuals. 

   Below we estimate the PH model without the regressor female, with the 
   regressor female, only for men and only for women, and compare the stratified 
   estimates.
*/
	
qui {
	stcox grade part_time lag stm married, hr vce(robust)
	estimates store original

	stcox grade part_time lag stm married female, hr vce(robust)
	estimates store by_gender

	stcox grade part_time lag stm married if female==0, hr vce(robust)
	estimates store males

	stcox grade part_time lag stm married if female==1, hr vce(robust)
	estimates store females
}

stcox grade part_time lag stm married, hr strata(female) vce(robust)
estimates store strata

estimates table original by_gender males females strata, star stats(N)
		
** Cumulative hazard
cap drop H0 H1
predict H, basech
gen H0 = H if female==0
gen H1 = H if female==1
line H0 H1 _t, c(J J) sort legend(label(1 "Men") label(2 "Women"))
* graph export "ta4_strata.png", as(png) replace

** I.3 Parametric Approach ----

/* The 'streg' commands allows for a parametric hazard function: exponential, 
   Weibull, Gompertz, ...
*/
	
use "TA4.dta", clear
stset duration, failure(event) id(id)
gen married = (marriage != 99)
global x female grade part_time lag stm married
	
** Model estimates

qui stcox $x, vce(robust)
eststo Cox
stcurve, hazard title("Cox PH") saving(hazard_cox, replace)

qui streg $x, dist(exponential) vce(robust) time
eststo Exponential
stcurve, hazard title("Exponential") saving(hazard_exponential, replace)

qui streg $x, dist(weibull) vce(robust) time
eststo Weibull
stcurve, hazard title("Weibull") saving(hazard_weibull, replace)

qui streg $x, dist(loglogistic) vce(robust)
eststo Loglogit
stcurve, hazard title("Loglogistic") saving(hazard_loglogit, replace)

qui streg $x, dist(lognormal) vce(robust)
eststo Lognormal
stcurve, hazard title("Lognormal") saving(hazard_lognormal, replace)

esttab Cox Exponential Weibull Loglogit Lognormal, se keep($x _cons) compress nogap mtitle stat(ll aic bic N)

** Fitted hazards (evaluated at the mean of the regressors)

sts graph, hazard title("Smoothed Kaplan-Meier") saving(hazard_nonparametric, replace)

gr combine hazard_nonparametric.gph hazard_cox.gph hazard_exponential.gph hazard_weibull.gph hazard_loglogit.gph hazard_lognormal.gph
*graph export "ta4_parametric_hazard.png", as(png) replace

** - Postestimation -

use "TA4.dta", clear
stset duration, failure(event) id(id)
gen married = (marriage != 99)
global x female grade part_time lag stm married
		
** Predicting time of failure 

streg $x, dist(weibull) vce(robust)
		
predict t_mean, time mean // Mean survival time
predict t_median, time median // Median survival time 
predict lnt_mean, lntime mean // Mean ln(survival time)
predict lnt_median, lntime median // Median ln(survival time)

list _t $x t_mean t_median lnt_mean lnt_median in L, abbrev(10)

** Predicting hazard 

qui streg $x, dist(weibull) vce(robust)
predict h, hazard
list id _t0 _t $x _d h if inlist(id,10,21)
predict s, surv	// Predict S(t|t0)
predict cs, csurv // Cumulative survival: S(t|earliest t0 for subject)
list id _t0 _t $x _d s cs if inlist(id,10,21)
		
/* More on 'stcurve': plot the hazard and survivor function for different groups
   according to being-part time student and high-school grades.
*/
qui streg $x, dist(weibull) vce(robust)
stcurve, hazard at1(part_time=0 grade=1) at2(part_time=0 grade=2) ///
                at3(part_time=0 grade=3) at4(part_time=0 grade=4) ///
				at5(part_time=0 grade=5)
stcurve, hazard at1(part_time=1 grade=1) at2(part_time=1 grade=2) ///
                at3(part_time=1 grade=3) at4(part_time=1 grade=4) ///
				at5(part_time=1 grade=5)
stcurve, survival at(part_time=0 grade=5)

** I.4 Unobserved heterogeneity ----

use "TA4.dta", clear
stset duration, failure(event) id(id)
gen married = (marriage != 99)
global x female part_time lag stm married

** Cox PH with Gamma-distributed RE

stcox $x, nohr shared(id)
predict nu, effects // nu = log(alpha), the RE
sum nu // very small in magnitude and not varying much

** Weibull model with frailty of inverse-Gaussian form
streg, dist(weibull) frailty(invgau) vce(robust) nolog nohr

** Visualize the effect of the frailty on the fitted hazard rate

stcurve, hazard alpha1 title("Conditional on frailty")
*graph export "ta4_frailty1.png", as(png) replace
stcurve, hazard unconditional title("Unconditional on frailty")
*graph export "ta4_frailty0.png", as(png) replace

** The model with regressors can't converge, so Weibull model may not be a good choice.

** PART II: DISCRETE DURATION -------------------------------------------------
	
use "TA4.dta", clear
gen married = (marriage!=99)
global x female grade part_time stm married

/* Expand each id-duration record to #duration records, time-indexed by _period 
   which ranges from 1 to $duration.

   Now each row refers to each person-period.
   
   The unit of duration variable is months.
*/
cap ssc install dthaz
prsnperd id duration

** Group periods (t: month --> a: quarter)
gen quarter = ceil(_period/3)
	
** Binary dependent variable: whether exit happens in each period or not
gen y = (duration==_period & event==1)
	
** Logit
qui logit y i.quarter $x, nocons vce(cluster id)
eststo logit
	
/* Complementary Log-Log:
   Dummie for each time period are included as additional regressors, equivalent
   to the Cox PH model (see Cameron and Trivedi, 2005, Chapter 17.10). */
qui cloglog y i.quarter $x, nocons vce(robust)
eststo cloglog	
**The predicted probability is the hazard 
predict hazard, pr			
	
** Alternative: Generalized linear model (cloglog family)
qui glm y i.quarter $x, family(binomial) link(cloglog) nocons vce(robust)
eststo glm

** Comparison 
esttab logit cloglog, se mtitle compress nogap keep($x) stat(ll aic bic N)
	
** Unobserved individual heterogeneity (optional)
	
/* LR test results sugggests that there is significant unobserved heterogeneity. 

   Gamma variance: from the gamma mixture distribution of the unobserved 
   individual heterogeneity (frailty).
*/
cap ssc install pgmhaz8
pgmhaz8 quarter $x, id(id) seq(_period) dead(y) nocons