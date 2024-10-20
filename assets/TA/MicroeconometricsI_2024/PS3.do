** PROBLEM SET 3: CENSORING, TRUNCATION, AND SELECTION
** Microeconometrics, IDEA, FALL 2024
** TA: Conghan Zheng

cls, clear all
capture set more off

cd "/Users/zheng/Documents/02 IDEA_PhD/Teaching/TA_Microeconometrics_Fall_IDEA/2024/Part I/PS/PS3/code and data"
	
** Packages to be installed 
* estout // Tools for making regression tables
	
use PS3.dta, clear

label define Education 1"Less than high-school graduated" ///
					   2"Some years in college" ///
					   3"College graduated" ///
					   4"Higher education"
label value educ Education

replace educ_sp = 0 if married == 0
label define Education_sp 0"Non-married" ///
					   	  1"SP: Less than high-school graduated" ///
					      2"SP: Some years in college" ///
					      3"SP: College graduated" ///
					      4"SP: Higher education"
label value educ_sp Education_sp

** Exercise 1.1: Tobit model ----
	
**Log wages
gen lny= ln(earn_h_r)

** Dummy for positive wages
gen dy = (earn_h_r > 0)
	
** New threshold
summarize lny
scalar gamma = r(min)		
disp "gamma = " gamma
scalar gamma01 = gamma - 0.0000001
replace lny = gamma01 if lny ==.
tab earn_h_r if earn_h_r < 0.02
tab lny if lny < gamma + 0.01 

** Regressors
global x married educ i.cohort i.reg time ///
       i.educ#i.cohort i.educ#c.time i.cohort#c.time i.reg#c.time ///
	   i.educ##c.time##i.cohort
	
** Tobit model for real hourly wages (left-censored at 0) 
tobit earn_h_r $x [pw=weight], ll(0)		
eststo leveltobit

** Tobit model for log real hourly wages (left-censored at gamma01)
tobit lny $x [pw=weight], ll(gamma01) vce(robust)
eststo logtobit

scalar lltobit = e(ll)

** Prediction
predict xb, xb
matrix btobit = e(b)
scalar sigma = sqrt(btobit[1,colsof(btobit)]) // the last element of the estimates matrix
gen threshold = (gamma01-xb)/sigma
gen yhat_tobit = exp(xb+0.5*sigma^2)*(1-normal((gamma01-xb-sigma^2)/sigma))
gen ytrunchat_tobit = yhat_tobit/(1 - normal(threshold)) if dy == 1

esttab leveltobit logtobit, star scalars(ll_0 ll) aic noomitted wide se /// 
       title("Exercise 1.1"\label{table1})	nonumbers
	   
** Exercise 1.2: Marginal effects on the left-censored mean of log wages ----
	
qui tobit lny $x [pw=weight], ll(gamma01)
eststo m1: margins, dydx(time) atmeans predict(ystar(gamma01, .))
esttab m1 // this stores the tobit (the active estimation in this case)
	
qui tobit lny $x [pw=weight], ll(gamma01)
eststo m2: margins, dydx(time) atmeans predict(ystar(gamma01, .)) post
esttab m2 // with the 'post' option, now the active estimation is the marginal effects

** Exercise 1.3: Two-part model of log real hourly wages ----
	
** Part 1: Probit for observing wages larger than 0
probit dy $x [pw=weight]
eststo twopart_probit
scalar llprobit = e(ll)
	
** Part 2: Regression of log wages
reg lny $x [pw=weight] if dy==1
eststo twopart_ols
scalar lllognormal = e(ll)
	
cap drop res_lny
predict res_lny, residuals
	
** Log-likelihood of the two-part model 
scalar lltwopart = llprobit + lllognormal
display "The log-likelihood of the Tobit model is " lltobit
display "The log-likelihood of the two-part model is lltwopart = " lltwopart

** Compare coefficients
esttab logtobit twopart_ols, se compress mtitle
	   
** Exercise 2.1: Heckman correction method for selection ----

/* 
In Stata, when using the Heckman two-step estimation (twostep option), you cannot combine it with certain options like weights. 

This is because these options involve adjustments to the standard errors or estimation procedure that are not compatible with the two-step process in Heckman’s model.

So we do the Heckman twostep estimation without weights, and then we use the FIML estimation with weights to compare the results.
*/
heckman lny $x, select(dy = i.educ_sp benefit $x) twostep
eststo heckman_liml
	
** Predicttion
predict probpos, psel
predict x1b1, xbsel
predict x2b2, xb
scalar sig2sq = e(sigma)^2
scalar sig12sq = (e(rho)*e(sigma))^2
gen yhat_heckman_liml = exp(x2b2 + 0.5*(sig2sq))*(1-normal(-x1b1-sig12sq))
gen ytrunchat_heckman_liml = yhat_heckman_liml/probpos

** Exercise 2.2: FIML estimation of the selection model ----

** FIML with sampling weights
heckman lny $x [pw=weight], select(dy = i.educ_sp benefit $x)
eststo heckman_fiml

** Predicttion
cap drop probpos x1b1 x2b2
predict probpos, psel
predict x1b1, xbsel
predict x2b2, xb
scalar sig2sq = e(sigma)^2
scalar sig12sq = (e(rho)*e(sigma))^2
gen yhat_heckman_fiml = exp(x2b2 + 0.5*(sig2sq))*(1-normal(-x1b1-sig12sq))
gen ytrunchat_heckman_fiml = yhat_heckman_fiml/probpos

** FIML without sampling weights
heckman lny $x, select(dy = i.educ_sp benefit $x)
eststo heckman_fiml_nowgt

** Predicttion
cap drop probpos x1b1 x2b2
predict probpos, psel
predict x1b1, xbsel
predict x2b2, xb
scalar sig2sq = e(sigma)^2
scalar sig12sq = (e(rho)*e(sigma))^2
gen yhat_heckman_fiml_nowgt = exp(x2b2 + 0.5*(sig2sq))*(1-normal(-x1b1-sig12sq))
gen ytrunchat_heckman_fiml_nowgt = yhat_heckman_fiml_nowgt/probpos
	
** Comparison of estimates
esttab heckman_liml heckman_fiml_nowgt, se wide mtitle
	
** Comparison of prediced wages
display "yhat_heckman_LIML = " yhat_heckman_liml ", yhat_heckman_FIML_woWeights = " = yhat_heckman_fiml_nowgt

** Exercise 2.3: Truncated OLS, Tobit and two-step selection models ----
	
** Truncated OLS (part two of the two-part model)
reg lny $x [pw=weight] if dy==1
eststo truncols

** Prediction
predict mu, xb
gen ytrunchat_ols = exp(mu+0.5*e(rmse)^2)

** Compare prediced wages
display "yhat_truncated_OLS = " ytrunchat_ols ", yhat_Tobit = " ytrunchat_tobit ", yhat_heckman_FIML = " ytrunchat_heckman_fiml  "
	
** Compare the estimates
esttab truncols logtobit heckman_fiml, scalars(ll) aic wide noomitted replace se ///
	   title("Exercise 2"\label{table2}) nonumbers ///
	   mtitles("Truncated OLS" "Tobit" "Selection FIML") 	
	
** Summarize observations and predictions for truncated wages
summarize earn_h_r ytrunchat_ols ytrunchat_tobit ytrunchat_heckman_fiml if dy==1

** Plot truncated means 
preserve
	collapse (mean) earn_h_r ytrunchat_ols ytrunchat_tobit ytrunchat_heckman_fimlif dy==1, by(year)
	
	twoway (line earn_h_r year) (line ytrunchat_ols year) (line ytrunchat_tobit year) (line ytrunchat_heckman_fiml year), ///
		legend(label(1 "Observed") label(2 "Truncated OLS") label(3 "Tobit") label(4 "Selection FIML")) ///
		ytitle(" ") xtitle("Year") ylabel(,angle(horizontal)) ///
		subtitle("Truncated mean: E(y|y>0)", position(11) justification(left)) ///
		title("Average Real Hourly Wage, 2001-2015")
	graph export "PS3_comparison.png", as(png) replace
restore