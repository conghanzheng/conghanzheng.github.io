** TA Session - TREATMENT EFFECTS: RD and DiD
** Microeconometrics, IDEA, FALL 2024
** TA: Conghan Zheng

/* Contents
	
I) Regression Discontinuity

	- Sharp design 
	- Fuzzy design 
	
II) Difference in Differences 	
	
	- Repeated cross-sections
	- Panel data 
*/	

/* Data

TA2_1.dta: An artificial cross section of 20 observations.

TA2_1.dta: An artificial panel of 12 individuals and two periods.
*/	

** Packages to be installed ----

** estout - Making regression tables
* ssc install estout 
	
** RD
* ssc install rd 
	
** DID
* ssc install diff

cls, clear all
cap set more off

cd "..."
	
** I REGRESSION DISCONTINUITY ----

** I.1 Sharp design

** Simulate some data (N = 1000)
clear all
set obs 1000
gen s = 10 + 5*invnormal(uniform())
	
** Generate X
global s_star = 10
gen x = s - $s_star
	
** Generate D (sharp design)
gen d = 1 if s > $s_star
replace d = 0 if s <= $s_star
	
** Define y. The true value of the treatment effect will be E[Y1] - E[Y0] = 400
gen y1 = 600 + 6.5*x - 2*x^2 + 0.001*x^3 + 300*invnorm(uniform()) 
gen y0 = 200 + 6.5*x - 0.20*x^2 + 0.01*x^3 + 300*invnorm(uniform()) 
	
** Generate the observable outcome 
gen y = y0 + d*(y1 - y0)
	
** Visualize outcome discontinuity 
twoway histogram y if s > $s_star, barw(60) color(green%50) || ///
	hist y if s < $s_star, barw(60) color(red%30) ///
	legend(order(1 "Right side" 2 "Left side") pos(11) col(1) ring(0)) ///
	xtitle() ytitle(Frequency) ylabel()
	
** Polynomials for the control function
gen dx = d*x
gen dx2 = d*(x^2)
gen dx3 = d*x^3
gen x2 = x^2
gen x3 = x^3
	
** Estimation without high order polynomials 
reg y d x dx 
cap drop y_hat_1
predict y_hat_1, xb
graph twoway (scatter y s if s>=$s_star , clstyle(p1)) ///
	(scatter y s if s<=$s_star , clstyle(p1)) ///
	(line y_hat_1 s if s>=$s_star , msymbol(o)) ///
	(line y_hat_1 s if s<=$s_star , msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Sharp–RDD – Parametric linear regression") ///
	legend( label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right Prediction") label(4 "Left Prediction"))
		
graph export "ta6_sharprd1.png", as(png) replace
	
** Improve the fit by including higher order polynomials 
reg y d x x2 x3 dx dx2 dx3
cap drop y_hat_2
predict y_hat_2, xb
graph twoway (scatter y s if s>=$s_star , clstyle(p1)) ///
	(scatter y s if s<=$s_star , clstyle(p1)) ///
	(scatter y_hat_2 s if s>=$s_star , msymbol(o)) ///
	(scatter y_hat_2 s if s<=$s_star , msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Sharp–RDD – Parametric Polynomial Regression") ///
	legend(label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right Prediction") label(4 "Left Prediction"))

graph export "ta6_sharprd2.png", as(png) replace
	
/* Local polynomial functions f0 and f1: We can select different kernels, 
   bandwiths and degrees. 
*/
global bandw 5 
cap drop f0 f1
lpoly y s if s <= $s_star, gen(f0) at(s) kernel(triangle) bwidth($bandw) degree(3) nograph
lpoly y s if s > $s_star, gen(f1) at(s) kernel(triangle) bwidth($bandw) degree(3) nograph

graph twoway (scatter y s if s >= $s_star , clstyle(p1)) ///
	(scatter y s if s <= $s_star , clstyle(p1)) ///
	(scatter f0 s if s < $s_star, msize(medsmall) msymbol(o)) ///
	(scatter f1 s if s >= $s_star, msize(medsmall) msymbol(o)), ///
	xline($s_star, lpattern(dash)) ///
	title("Sharp RDD – Local polynomial regression (LPR)") ///
	legend(label(1 "Right actual data") label(2 "Left actual data") ///
	label(3 "Right LPR prediction") label(4 "Left LPR prediction")) ///
	note(Bandwidth = $bandw)

graph export "ta6_sharprd3.png", as(png) replace

** Calculate the treatment effect 
gen z=$s_star
cap drop f0 f1 
qui lpoly y s if s<=$s_star, gen(f0) at(z) k(tri) bw($bandw) deg(3) nogr
qui lpoly y s if s>$s_star, gen(f1) at(z) k(tri) bw($bandw) deg(3) nogr
scalar rdef=f1[1]-f0[1]
display rdef

** Obtain standard errors for the local polynomial estimates 
cap program drop rdd_s
program rdd_s, rclass
	version 16 
	args deg ker band cut
	cap drop f0 f1 z
	gen z = `cut'
	qui lpoly y s if s < `cut', gen(f0) at(z) k(`ker') bw(`band') deg(`deg') nogr
	qui lpoly y s if s >= `cut', gen(f1) at(z) k(`ker') bw(`band') deg(`deg') nogr
	return scalar rdef = f1[1]-f0[1]
end
	
/* 3-degree polynomial with triangular kernel, bandwith of 5, and discontinuity 
   cut equal to 10 
*/ 
rdd_s 3 tri 5 10
return list 
bootstrap r(rdef), reps(50): rdd_s 3 tri 5 10
	
/* Regression discontinuity estimates using the 'rd' command 

   'rd' estimates local linear or kernel regression models on both sides of the
   cutoff.  Estimates are sensitive to the choice of bandwidth. 
*/ 
rd y s, z0($s_star) bwidth($bandw)
	
** Various bandwiths
rd y s, z0($s_star) bwidth($bandw) mbw	

** I.2 Fuzzy design
	
** Simulate some data 
clear all
set obs 1000
	
gen s = -1 + 2*runiform()
gen Z = (s >= 0)
gen v = rnormal(0,1)
gen d = (-0.5 + Z + s + v >= 0)
gen y1 = 2 + s + s^2 + 3*s^3 + invnorm(uniform())
gen y0 = 1 + s + s^2 + 3*s^3 + invnorm(uniform())
gen y = y0 + d*(y1 - y0)

gen s2 = s^2
gen s3 = s^3

** Under the fuzzy design, the treatment is endogeneous
ivregress 2sls y s s2 s3 (d=Z), first 
	
** Nonparametric estimation using a 3rd-degree polynomial
global s_star = 0
global bandw = 5 
cap drop f0 f1 
lpoly y s if s<$s_star, gen(f0) at(s) k(tri) bw($bandw) deg(3) nogr
lpoly y s if s>=$s_star, gen(f1) at(s) k(tri) bw($bandw) deg(3) nogr

** discountinuity in y
graph twoway (scatter y s if s>=$s_star , clstyle(p1)) ///
	(scatter y s if s<=$s_star , clstyle(p1)) ///
	(scatter f0 s if s<$s_star, msize(medsmall) msymbol(o)) ///
	(scatter f1 s if s>=$s_star, msize(medsmall) msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Fuzzy-RDD - Outcome Non-parametric Local Linear Regression", size(medlarge)) ///
	legend(label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right LLR Prediction") label(4 "Left LLR Prediction")) ///
	note(Bandwidth = $bandw)
graph export "ta6_fuzzyrd1.png", as(png) replace
	
capture drop g0 g1
lpoly d s if s<$s_star, gen(g0) at(s) k(tri) bw($bandw) deg(3) nogr
lpoly d s if s>=$s_star, gen(g1) at(s) k(tri) bw($bandw) deg(3) nogr

** discontinuity in P(D=1)
graph twoway ///
	(scatter d s if s>=$s_star & d==1, clstyle(p1)) ///
	(scatter d s if s<=$s_star & d==0, clstyle(p1)) ///
	(scatter g0 s if s<$s_star, msize(medsmall) msymbol(o)) ///
	(scatter g1 s if s>=$s_star, msize(medsmall) msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Fuzzy-RDD - Probability Non-parametric Local Linear Regression") ///
	legend(label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right LLR Prediction") label(4 "Left LLR Prediction")) ///
	note(Bandwidth = $bandw)
graph export "ta6_fuzzyrd2.png", as(png) replace

** Obtain standard errors for the local polynomial estimates by bootstrapping
cap program drop rdd_f
program rdd_f, rclass
	version 16
	args deg ker band cut
	* Outcome discontinuity
	cap drop z f0 f1 g0 g1
	gen z = `cut'
	qui lpoly y s if s<`cut', gen(f0) at(z) k(`ker') bw(`band') deg(`deg') nogr
	qui lpoly y s if s>=`cut', gen(f1) at(z) k(`ker') bw(`band') deg(`deg') nogr
	scalar disc_y = f1[1]-f0[1]
	* Probability discontinuity
	cap drop g0 g1
	qui lpoly d s if s<`cut', gen(g0) at(z) k(`ker') bw(`band') deg(`deg') nogr
	qui lpoly d s if s>=`cut', gen(g1) at(z) k(`ker') bw(`band') deg(`deg') nogr
	scalar disc_d = g1[1]-g0[1]
	return scalar rddef = disc_y/disc_d
end

cap drop z 
cap drop f0 f1 g0 g1
bootstrap r(rddef), reps(10): rdd_f 3 tri 5 0
	
** Or we can use the 'rd' command 
rd y d s, z0(0) bwidth(5)
	
** II) DIFFERENCE IN DIFFERENCES ----
	
** II.1 Repeated cross-sections

** We assume that the sample composition does not vary over time
use "TA2_1.dta", clear
	
** Method 1: DID estimate for ATE is the coefficient estimate of DT
reg Y D T DT 
	
** Method 2: step by step: 

** Average over treated units in T = 0
qui summarize Y if D==1 & T==0
scalar treat0=r(mean)
	
** Average over treated units in T = 1
qui summarize Y if D==1 & T==1 
scalar treat1=r(mean)
	
** Average over untreated units in T = 0 
qui summarize Y if D==0 & T==0
scalar untreat0=r(mean)
	
** Average over untreated units in T = 1
qui summarize Y if D==0 & T==1
scalar untreat1=r(mean)	
	
** Diff-in-diff coefficient
scalar didbeta = treat1 - treat0 - (untreat1 - untreat0)
di "The DID coefficient is " didbeta
	
** Method 3: the 'diff' command
diff Y, treated(D) period(T)

** II.2 Panel data

** By construction, the sample composition is the same across time
use "TA2_2.dta", clear

xtset id year 
	
** Lags
gen y_1 = L.y 
gen d_1 = L.d 
gen x1_1 = L.x1 
gen x2_1 = L.x2 

** Differences
gen delta_y = D.y
gen delta_x1 = D.x1 
gen delta_x2 = D.x2 
gen delta_d = D.d 

** Method 1: estimate for ATE is the coefficient estimate for d
reg delta_y d if d_1==0

** Method 2: step by step
qui sum delta_y if d == 1 & d_1==0 
scalar mean_t = r(mean)
qui sum delta_y if d == 0 & d_1==0
scalar mean_c = r(mean)
scalar did = mean_t - mean_c 
di "The DID coefficient is " did 
	
** If fixed effects are considered, the first-differenced model can be used  
reg delta_y delta_d, noconst
	
** WG estimator: the coefficient estimate for d*dyear_2000 is the DID estimate
xi, noomit prefix("d") gen i.year
xtreg y d##dyear_2000, fe		

** DID with covariates 
reg delta_y d delta_x1 delta_x2