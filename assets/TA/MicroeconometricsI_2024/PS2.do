** Problem Set 2: Discrete Choice
** Microeconometrics, IDEA, FALL 2024
** TA: Conghan Zheng
**
** - Inputs: PS2.dta

cls, clear all

cd "..."

use "PS2.dta", clear

** Exercise 1.1 ----------------------------------------------------------------
	
gen age2 = age^2
	
gen same_ind = (ind_prior == ind_after)
	
eststo clear
	
** Logit
logit same_ind exper unempl_dur age age2 female ln_wage veneto_resid
eststo Logit
	
** Probit
probit same_ind exper unempl_dur age age2 female ln_wage veneto_resid
eststo Probit
	
** Compare logit and probit 
esttab Logit Probit, se compress nogaps stat(ll N) mtitle

** Report odds ratios	
logit same_ind exper unempl_dur age age2 female ln_wage veneto_resid, or

** Exercise 1.2 ----------------------------------------------------------------
	
** At-means marginal effects for the logit model
qui logit same_ind exper unempl_dur age age2 female ln_wage veneto_resid 
margins, dydx(*) atmeans

** At-means marginal effects for the probit model
qui probit same_ind exper unempl_dur age age2 female ln_wage veneto_resid 
margins, dydx(*) atmeans
	
** Exercise 2.1 ----------------------------------------------------------------
	
** Generate an categorical outcome variable for the choice of industry after displacement
gen choice_ind = 1 if inlist(ind_after, 1, 2, 3)
replace choice_ind = 2 if inlist(ind_after, 4, 5, 6)
replace choice_ind = 3 if inlist(ind_after, 7, 8)
label define lbl_choice_ind 1 "Manufacturing" 2 "Services" 3 "Public sector"
label value choice_ind lbl_choice_ind

** Multinomial Logit
mlogit choice_ind exper unempl_dur age age2 female ln_wage veneto_resid, baseoutcome(1) rrr 

** Exercise 2.2 ----------------------------------------------------------------
	
** At means marginal effects
qui mlogit choice_ind exper unempl_dur age age2 female ln_wage veneto_resid, ///
		baseoutcome(1) rrr 
		
margins, dydx(*) predict(outcome(1)) atmeans // manufacturing
margins, dydx(*) predict(outcome(2)) atmeans // services
margins, dydx(*) predict(outcome(3)) atmeans // public sector

** Exercise 2.3 ----------------------------------------------------------------

** Reload dataset and do different cleansing
use "PS2.dta", clear

gen age2 = age^2
	
** Outcome variable
gen choice_ind = ind_after 
	
** Choice dummies 
gen d_manuf_light = (choice_ind == 1)
gen d_manuf_heavy = (choice_ind == 2)
gen d_manuf_elecon = (choice_ind == 3)
gen d_serv_sales = (choice_ind == 4)
gen d_serv_fin = (choice_ind == 5)
gen d_serv_other = (choice_ind == 6)
gen d_public_adm = (choice_ind == 7)
gen d_public_hlth = (choice_ind == 8)

** Reshape dataset to accommodate alternative-varying regressors
reshape long d_ e_, i(id) j(industry manuf_light manuf_heavy manuf_elecon ///
		serv_sales serv_fin serv_other public_adm public_hlth) string
		
order id industry year choice_ind d_ e_

** Conditional logit
asclogit d_ e_, case(id) alternatives(industry) ///
		 casevars(unempl_dur age age2 female ln_wage veneto_resid)
eststo CLogit8
	
** At-means marginal effects

estat mfx, varlist(e_)
	
gen type = 1 if inlist(industry, "manuf_elecon", "manuf_heavy", "manuf_light")
replace type = 2 if inlist(industry, "serv_sales", "serv_fin", "serv_other")
replace type = 3 if inlist(industry, "public_adm", "public_hlth")
	
label define lbl_choice_ind 1 "Manufacturing" 2 "Services" 3 "Public sector"
label value type lbl_choice_ind
	
** Predict choice probabilities
predict p
bysort type: summ p

** Exercise 2.4 ----------------------------------------------------------------
	
drop type
	
** Tree structure 
nlogitgen type = industry(manufacturing: manuf_light | manuf_heavy | manuf_elecon, ///
		  services: serv_sales | serv_fin | serv_other, public: public_adm | public_hlth)

nlogittree industry type, choice(d_)
	
** Nested logit 
nlogit d_ e_ || type:, || industry: unempl_dur age age2 female ln_wage veneto_resid, ///
	   noconst case(id) 
eststo NestedLogit
	
** Compare with conditional logit
esttab CLogit8 NestedLogit, se stat(ll N) mtitle nogaps compress