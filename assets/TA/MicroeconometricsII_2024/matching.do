** TREATMENT EFFECTS: Matching
** Microeconometrics, IDEA FALL 2024
** Conghan Zheng

/* Contents

1 Covariates matching	
2 Propensity score matching

*/	

/* Datasets

matching.dta: Ham and LaLonde (1996), link: https://www.jstor.org/stable/2171928.
		   Source: National Supported Work (NSW) labor market experiment.
		   Treatment: on-the-job training lasting between 9 months and one year (1976-1977).
		   
		   y = real earnings in 1978; d = on-the-job training.
*/	

** Package Dependencies ----

** 1. estout - Making regression tables
* ssc install estout 
	
** 2. psmatch2 - Perform propensity score matching 
* ssc install psmatch2
	
** 3. pscore - Estimation of average treatment effects by Becker and Ichino
* help st0026_2
	
** 4. Module to perform coarsened exact matching 
* ssc install cem

cls, clear all
cap set more off

cd "..."
	
** MATCHING ----

** matching.dta: y = real earnings in 1978; d = on-the-job training
use "matching.dta", clear

gen re74sq = re74^2
gen u74hisp = unem74*hisp
gen educsq = educ^2

global x_ps age agesq educ marr nodegree black hisp re74 re75  

** 1 Covariates matching ----
	
** Nearest-neighbor (NN) matching 
teffects nnmatch (y $x_ps) (d)
	
** Specify the minimum number of matches per observation (default is 1)
teffects nnmatch (y $x_ps) (d), nneighbor(3)

** NN matching with exact matching on some binary variables 
teffects nnmatch (y $x_ps) (d), nneighbor(3) ematch(hisp black)	

** Bias adjustment for large samples
teffects nnmatch (y $x_ps) (d), biasadj(re74 re75)	

** 2 Propensity score matching ----

/* Stata commands: 
   
   - estimating propensity score: pscore
   
   - matching propensity score: psmatch2
   
   - estimating treatment effect: teffects psmatch; att*
*/
	
** NN matching, propensity scores are predicted from a logit of D on X
teffects psmatch (y) (d $x_ps)
	
** Assess the overalap graphically 
qui teffects psmatch (y) (d $x_ps), generate(near_obs)
** Nearest neighbors 
sum near_*
teffects overlap
graph export "ta1_psoverlap.png", as(png) replace

/* Command 'pscore': estimates the propensity score and tests the Balancing 
   Hypothesis.

   Link: https://www.sscnet.ucla.edu/soc/faculty/mason/readings/becker_ichino_pscore_sj_2002.pdf
   
   - Option 'pscore' is a compulsory option and specifies the variable name for 
   the estimated propensity score.
   
   - Option 'comsup' restricts the analysis of the balancing property to all 
   treated plus those controls in the region of common support. A dummy variable
   named comsup is added to the dataset to identify the observations in the 
   common support.
   
   - Option 'blockid' specifies the variable name for the block number of the 
   estimated propensity score.
   
   - What is a block? 
   Blocks defined over intervals of propensity score. Stratification or interval
   matching is based on the idea of dividing the range of variation of the 
   propensity score in intervals such that within each interval, the treated and
   control units have, on the average, the same propensity score.
   
   The overall treatment effect is a weighted average of block-specific 
   treatment effects, where the weight for each block is given by the 
   corresponding fraction of treated units.
   
   - Option 'level(#)' specifies the significance level of the tests of the 
   balancing property.
   
   - Option 'logit' specifies that a logit model to estimate the propensity 
   score be used instead of the default probit model
*/
pscore d $x_ps, pscore(myscore) comsup blockid(myblock) level(0.005) logit

/** 'psmatch2': matching the propensity score

   Link: http://repec.org/bocode/p/psmatch2.html
*/

** k-Nearest neighbors matching: neighbor(k)
psmatch2 d, out(y) pscore(myscore) neighbor(3) common 
	
** Test of the balancing property 
pstest myscore
pstest $x_ps
	
** pstest results for treated and untreated 
summarize myscore [aweight=_weight] if d==0
summarize myscore [aweight=_weight] if d==1 

** Graph: quality of our matching
label define tstatus 0 "Control sample" 1 "Treated sample"
label values d tstatus
label variable d "Treatment Status"
	
** Before matching 
qui graph twoway (kdensity myscore if d==1, msize(small)) ///
	(kdensity myscore if d==0, msize(small) ///
	lpattern(shortdash_dot)), ///
	subtitle(, bfcolor(none)) ///
	xtitle("propensity–score (Before)", size(medlarge)) xscale(r(0.1 0.7)) /// 
	ytitle("Density", size(medlarge)) yscale(r(0 9)) ///
	legend(pos(12) ring(0) col(1)) ///
	legend( label(1 "Treated") label(2 "Untreated")) saving(BEFORE, replace)

** After matching 
qui graph twoway (kdensity myscore [aweight=_weight] if d==1,  msize(small)) ///
	(kdensity myscore [aweight=_weight] if d==0, msize(small) ///
	lpattern(shortdash_dot)), ///
	subtitle(, bfcolor(none)) ///
	xtitle(" propensity–score (After)", size(medlarge)) xscale(r(0.1 0.7)) ///
	ytitle("Density", size(medlarge)) yscale(r(0 9)) ///
	legend(pos(12) ring(0) col(1)) ///
	legend( label(1 "Treated") label(2 "Untreated")) saving(AFTER, replace)
		
** Combine graphs
graph combine BEFORE.gph AFTER.gph
graph export "ta1_match.png", as(png) replace

/* Commands 'att*': estimate the ATT based on the propensity score

   - 'attnd' and 'attnw': NN matching, standard errors are obtained analytically
   or by bootstrapping
   
   Difference: 
   'attnw' gives equal weight to the groups of forward and backward 
   matches; 
   'attnd' randomly draws either the forward or backward matches.
   
   
   - 'attr': radius matching
   
   - 'attk': kernel matching, default: Gaussian kernel, standard errors are 
   obtained by bootstrapping.
   
   - 'atts': stratification method

*/

** ATT estimate based on NN matching
attnd y d $x_ps, comsup

** ATT estimate based on radius matching
attr y d $x_ps, comsup logit
	
attr y d $x_ps, comsup logit radius(0.001)

attr y d $x_ps, comsup logit radius(0.0001)

** ATT estimate based on statification matching 
atts y d, pscore(myscore) blockid(myblock) comsup

** ATT estimate based on Kernel matching, bootstrap repetitions: 50
attk y d $x_ps, comsup boot reps(50) dots logit