***********************************************
*Title: Analysis for Cell Phone RCT 
*Date Last Modfied: October, 2018
*Author(s): Arun Chandrasekhar, Devika Lakhote
***********************************************

clear 
set more off
version 12

*setting date variables
loc curr_month = lower(substr(subinstr("`c(current_date)'"," ","",.),-7,.))

*defining and setting paths 
global path_main "~/Dropbox/Gossip Centrality/Final_Data"
global path_data "$path_main/data"
global path_tables "$path_main/tables"


cd "$path_data"
set scheme s2mono

use karnataka_cell_rct.dta


*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* FIGURE 1 */
*-------------------------------------------------------------------------------------

preserve
	cd "$path_tables"
	drop if flierVillage == 1 // checkme - flier village
	keep if inj>0 // drop villages where the expt didn't take off
	* CDFs
	sort CallsReceived
	cumul CallsReceived if random==1, gen(cuml_random_noc)
	label var cuml_random "Random"
	sort CallsReceived
	cumul CallsReceived if gossip==1, gen(cuml_gossip_noc)
	label var cuml_gossip "Gossip"
	sort CallsReceived
	cumul CallsReceived if elder==1, gen(cuml_elder_noc)
	label var cuml_elder "Elder"
	sort CallsReceived

	stack cuml_random_noc CallsReceived cuml_gossip_noc CallsReceived cuml_elder_noc CallsReceived, into(c temp) wide clear
	label var temp "Calls Received" 
	line cuml_random_noc cuml_gossip_noc cuml_elder_noc temp, lpattern(dot solid dash) ///
	legend(order(1 "Random" 2 "Gossip" 3 "Elder")) ///
	lwidth(thin thin thin) ///
	ytitle("F(x)") ///
	graphregion(color(white))
	graph export CDF_RF.pdf, replace
restore



*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* FIGURE 2 */
*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
preserve
	cd "$path_tables"
	set scheme s2mono
	
	loc gossip_count_control gossip_miss gossip_ct both_miss
	loc elder_count_control  elder_miss  elder_ct

	qreg CallsRec gossip elder num_seeds `gossip_count_control' `elder_count_control' if flierVillage == 0, level(90) 
	grqreg gossip elder, ci format(%9.0f) 
	graph export Quantile_Reg_RF.pdf, replace


	qreg CallsRec atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control' if flierVillage == 0, quantile(.5) level(90)
	grqreg atleast_gossip_1 atleast_elder_1, ci format(%9.0f) level(90) 
	graph export Quantile_Reg_OLS_v2.pdf, replace 
restore


*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* TABLE 1*/
*-------------------------------------------------------------------------------------
* Panel A
*-------------------------------------------------------------------------------------


loc gossip_count_control gossip_miss gossip_ct both_miss
loc elder_count_control  elder_miss  elder_ct


reg CallsRec gossip elder num_seeds `gossip_count_control' `elder_count_control'  if flierVillage ==0, ro //level(95)
//su CallsRec if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test gossip == elder
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_2_1_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop(RF) replace ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

*-------------------------------------------------------------------------------------

reg CallsRec atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control' if flierVillage ==0, ro level(90)
//su CallsRec if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0 
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 					
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_2_1_v2", keep(atleast_gossip_1 atleast_elder_1) label tex(frag) nocons noaster nor2 nonotes cttop(OLS) append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')

*-------------------------------------------------------------------------------------

reg atleast_gossip_1 gossip elder num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_gossip_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 					
di "Control Mean is `control_mean'"
test gossip == elder
loc testpval = round(`r(p)',0.01) 
outreg2 using "$path_tables/Table_RCT_A_2_1_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop("IV 1: First Stage") append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')


*-------------------------------------------------------------------------------------

reg atleast_elder_1 gossip elder num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_elder_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 					
di "Control Mean is `control_mean'" 
test gossip == elder
loc testpval = round(`r(p)',0.01) 
outreg2 using "$path_tables/Table_RCT_A_2_1_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop("IV 2: First Stage") append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

*-------------------------------------------------------------------------------------

ivregress 2sls CallsRec num_seeds `gossip_count_control' `elder_count_control' (atleast_gossip_1 atleast_elder_1 = gossip elder) if flierVillage ==0, ro
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 				
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_2_1_v2", keep(atleast_gossip_1 atleast_elder_1) label tex(frag) nocons noaster nor2 nonotes cttop("IV: Second Stage") append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')

*-------------------------------------------------------------------------------------
* Panel B
*-------------------------------------------------------------------------------------
reg NormalizedCalls_Ex gossip elder num_seeds `gossip_count_control' `elder_count_control'  if flierVillage ==0, ro //level(90)
//su NormalizedCalls_Ex if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test gossip == elder
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_2_1_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop(RF) replace ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

*-------------------------------------------------------------------------------------

reg NormalizedCalls_Ex atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control' if flierVillage ==0, ro //level(90)
//su NormalizedCalls_Ex if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0 
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test atleast_gossip_1 == atleast_elder_1
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_2_1_v2", keep(atleast_gossip_1 atleast_elder_1) label nocons tex(frag) noaster nor2 nonotes cttop(OLS) append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')

*-------------------------------------------------------------------------------------

reg atleast_gossip_1 gossip elder num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_gossip_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test gossip == elder
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_2_1_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop("IV 1: First Stage") append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

*-------------------------------------------------------------------------------------

reg atleast_elder_1 gossip elder num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_elder_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test gossip == elder
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_2_1_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop("IV 2: First Stage") append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

*-------------------------------------------------------------------------------------

ivregress 2sls NormalizedCalls_Ex num_seeds `gossip_count_control' `elder_count_control' (atleast_gossip_1 atleast_elder_1 = gossip elder) if flierVillage ==0, ro 
//su NormalizedCalls_Ex if e(sample) & random == 1 
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test atleast_gossip_1 == atleast_elder_1
loc testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_2_1_v2", keep(atleast_gossip_1 atleast_elder_1) label nocons tex(frag) noaster nor2 nonotes cttop("IV: Second Stage") append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')


*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* TABLE D.1*/
*-------------------------------------------------------------------------------------
* Panel A
*-------------------------------------------------------------------------------------


reg CallsRec gossip gossip5 elder elder5 num_seeds `gossip_count_control' `elder_count_control'  if flierVillage ==0, ro
//su CallsRec if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_A_2_2", keep(gossip gossip5 elder elder5) label tex(frag) nocons noaster nor2 nonotes cttop(RF) replace ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

reg CallsRec atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control' if flierVillage ==0, ro
//su CallsRec if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0 
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_A_2_2", keep(atleast_gossip_1 atleast_elder_1) label tex(frag) nocons noaster nor2 nonotes cttop(OLS) append ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

reg atleast_gossip_1 gossip elder gossip5 elder5 num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_gossip_1 if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_A_2_2", keep(gossip gossip5 elder elder5) label tex(frag) nocons noaster nor2 nonotes cttop("IV 1: First Stage") append ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

reg atleast_elder_1 gossip elder gossip5 elder5 num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_elder_1 if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_A_2_2", keep(gossip gossip5 elder elder5) label tex(frag) nocons noaster nor2 nonotes cttop("IV 2: First Stage") append ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

ivregress 2sls CallsRec num_seeds `gossip_count_control' `elder_count_control' (atleast_gossip_1 atleast_elder_1 = gossip elder gossip5 elder5) if flierVillage ==0, first ro 
//su CallsRec if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0 
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_A_2_2", keep(atleast_gossip_1 atleast_elder_1) label tex(frag) nocons noaster nor2 nonotes cttop("IV: Second Stage") append ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------
* Panel B
*-------------------------------------------------------------------------------------


reg NormalizedCalls_Ex gossip gossip5 elder elder5 num_seeds `gossip_count_control' `elder_count_control'  if flierVillage ==0, ro
//su NormalizedCalls_Ex if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_B_2_2", keep(gossip gossip5 elder elder5) label tex(frag) nocons noaster nor2 nonotes cttop(RF) replace ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

reg NormalizedCalls_Ex atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control' if flierVillage ==0, ro
//su NormalizedCalls_Ex if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0 
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_B_2_2", keep(atleast_gossip_1 atleast_elder_1) label tex(frag) nocons noaster nor2 nonotes cttop(OLS) append ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

reg atleast_gossip_1 gossip elder gossip5 elder5 num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_gossip_1 if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_B_2_2", keep(gossip elder gossip5 elder5) label tex(frag) nocons noaster nor2 nonotes cttop("IV 1: First Stage") append ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

reg atleast_elder_1 gossip elder gossip5 elder5 num_seeds `gossip_count_control' `elder_count_control' if flierVillage== 0, ro 
//su atleast_elder_1 if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_B_2_2", keep(gossip elder gossip5 elder5) label tex(frag) nocons noaster nor2 nonotes cttop("IV 2: First Stage") append ///
addstat("Control Group Mean", `control_mean')

*-------------------------------------------------------------------------------------

ivregress 2sls NormalizedCalls_Ex num_seeds `gossip_count_control' `elder_count_control' (atleast_gossip_1 atleast_elder_1 = gossip gossip5 elder elder5) if flierVillage ==0, ro 
//su NormalizedCalls_Ex if e(sample) & gossip == 0 & gossip5 == 0 & elder == 0 & elder5 == 0 
loc cons = _b[_cons]
di "Constant is `cons'"
foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}
loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
outreg2 using "$path_tables/Table_RCT_B_2_2", keep(atleast_gossip_1 atleast_elder_1) label tex(frag) nocons noaster nor2 nonotes cttop("IV: Second Stage") append ///
addstat("Control Group Mean", `control_mean')



*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* TABLE E.1*/
*-------------------------------------------------------------------------------------
* Panel A
*-------------------------------------------------------------------------------------


reg CallsRec gossip elder num_seeds `gossip_count_control' `elder_count_control', ro
//su CallsRec if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test gossip == elder 
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_flyer_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop(RF) replace ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

reg CallsRec atleast_gossip_1 atleast_elder_1  num_seeds `gossip_count_control' `elder_count_control', ro
//su CallsRec if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1  == 0 
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1 
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_flyer_v2", keep(atleast_gossip_1 atleast_elder_1 ) label tex(frag) nocons noaster nor2 nonotes cttop(OLS) append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')

reg atleast_gossip_1 gossip elder num_seeds `gossip_count_control' `elder_count_control', ro 
//su atleast_gossip_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test gossip == elder 
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_flyer_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop("IV 1: First Stage") append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

reg atleast_elder_1 gossip elder num_seeds `gossip_count_control' `elder_count_control', ro 
//su atleast_elder_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test gossip == elder 
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_flyer_v2", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes cttop("IV 2: First Stage") append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

ivregress 2sls CallsRec num_seeds `gossip_count_control' `elder_count_control' (atleast_gossip_1 atleast_elder_1 = gossip elder), ro 
//su CallsRec if e(sample) & random == 1 
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1 
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_A_flyer_v2", keep(atleast_gossip_1 atleast_elder_1 ) label tex(frag) nocons noaster nor2 nonotes cttop("IV: Second Stage") append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')


*-------------------------------------------------------------------------------------
* Panel B
*-------------------------------------------------------------------------------------


reg NormalizedCalls_Ex gossip elder num_seeds `gossip_count_control' `elder_count_control', ro
//su NormalizedCalls_Ex if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test gossip == elder
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_flyer_v2", keep(gossip elder) label tex(frag) nocons noaster nor2  nonotes cttop(RF)  replace ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

reg NormalizedCalls_Ex atleast_gossip_1 atleast_elder_1  num_seeds `gossip_count_control' `elder_count_control', ro
//su NormalizedCalls_Ex if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1  == 0
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_flyer_v2", keep(atleast_gossip_1 atleast_elder_1 ) label tex(frag) nocons noaster nor2 nonotes cttop(OLS) append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')

reg atleast_gossip_1 gossip elder num_seeds `gossip_count_control' `elder_count_control', ro 
//su atleast_gossip_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test gossip == elder
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_flyer_v2", keep(gossip elder) label tex(frag) nocons noaster nor2  nonotes cttop("IV 1: First Stage")  append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')


reg atleast_elder_1 gossip elder num_seeds `gossip_count_control' `elder_count_control', ro 
//su atleast_elder_1 if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test gossip == elder
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_flyer_v2", keep(gossip elder) label tex(frag) noaster nor2 nonotes cttop("IV 2: First Stage") append ///
addstat("Control Group Mean", `control_mean', "Gossip Treatment=Elder Treatment (pval.)", `testpval')

ivregress 2sls NormalizedCalls_Ex num_seeds `gossip_count_control' `elder_count_control' (atleast_gossip_1 atleast_elder_1 = gossip elder), ro 
//su NormalizedCalls_Ex if e(sample) & random == 1
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1
local testpval = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_RCT_B_flyer_v2", keep(atleast_gossip_1 atleast_elder_1 ) label tex(frag) nocons noaster nor2 nonotes cttop("IV: Second Stage") append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval')



*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* TABLE 9 */
*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
replace flyer = flierVillage 
replace num_seeds = NumberofSeeds


reg CallsRec atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control' if flyer == 0 & villageid != . & villageid != 62,  ro
//su CallsRec if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01) 
outreg2 using "$path_tables/Table_SeedDC_OLS_v2", keep(atleast_gossip_1 atleast_elder_1) nocons noaster tex(frag)  nonotes nor2 label replace ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval1')

reg CallsRec atleast_gossip_1 atleast_elder_1 dcAtleast1 num_seeds `gossip_count_control' `elder_count_control' if flyer == 0 & villageid != . & villageid != 62,  ro
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01)
test atleast_gossip_1 == dcAtleast1 
loc testpval2 = round(`r(p)',0.01)
test atleast_elder_1 == dcAtleast1
loc testpval3 = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_SeedDC_OLS_v2", keep(atleast_gossip_1 atleast_elder_1 dcAtleast1) nocons noaster tex(frag) nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean', ///
"At least 1 Gossip=At least 1 Elder (pval.)", `testpval1', ///
"At least 1 Gossip=At least 1 High \textit{DC} Seed (pval.)", `testpval2', ///
"At least 1 Elder=At least 1 High \textit{DC} Seed (pval.)", `testpval3')

 

reg CallsRec dcAtleast1 num_seeds `gossip_count_control' `elder_count_control' if flyer == 0 & villageid != . & villageid != 62,  ro 
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
outreg2 using "$path_tables/Table_SeedDC_OLS_v2", keep(dcAtleast1) noaster tex(frag) nocons nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean')

 

reg NormalizedCalls_Ex atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control' if flyer == 0 & villageid != . & villageid != 62,  ro
//su NormalizedCalls_Ex if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_SeedDC_OLS_v2", keep(atleast_gossip_1 atleast_elder_1) nocons noaster tex(frag)  nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval1')

reg NormalizedCalls_Ex atleast_gossip_1 atleast_elder_1 dcAtleast1 num_seeds `gossip_count_control' `elder_count_control' if flyer == 0 & villageid != . & villageid != 62,  ro
//su NormalizedCalls_Ex if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0 & dcAtleast1 == 0
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01)
test atleast_gossip_1 == dcAtleast1 
loc testpval2 = round(`r(p)',0.01)
test atleast_elder_1 == dcAtleast1
loc testpval3 = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_SeedDC_OLS_v2", keep(atleast_gossip_1 atleast_elder_1 dcAtleast1) nocons noaster tex(frag) nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean', ///
"At least 1 Gossip=At least 1 Elder (pval.)", `testpval1', ///
"At least 1 Gossip=At least 1 High \textit{DC} Seed (pval.)", `testpval2', ///
"At least 1 Elder=At least 1 High \textit{DC} Seed (pval.)", `testpval3')
 

reg NormalizedCalls_Ex dcAtleast1 num_seeds `gossip_count_control' `elder_count_control'  if flyer == 0 & villageid != . & villageid != 62,  ro 
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
outreg2 using "$path_tables/Table_SeedDC_OLS_v2", keep(dcAtleast1) noaster tex(frag) nocons nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean')


*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* TABLE E.2 */
*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------

reg CallsRec atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control'  if villageid != . & villageid != 62,  ro
//su CallsRec if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01) 
outreg2 using "$path_tables/Table_SeedDC_OLS_flyer_v2", keep(atleast_gossip_1 atleast_elder_1) nocons noaster tex(frag)  nonotes nor2 label replace ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval1')

reg CallsRec atleast_gossip_1 atleast_elder_1 dcAtleast1 num_seeds `gossip_count_control' `elder_count_control'  if villageid != . & villageid != 62,  ro
//su CallsRec if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0 & dcAtleast1 == 0
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01) 
test atleast_gossip_1 == dcAtleast1 
loc testpval2 = round(`r(p)',0.01)
test atleast_elder_1 == dcAtleast1
loc testpval3 = round(`r(p)',0.01)
outreg2 using "$path_tables/Table_SeedDC_OLS_flyer_v2", keep(atleast_gossip_1 atleast_elder_1 dcAtleast1) nocons noaster tex(frag) nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean', ///
"At least 1 Gossip=At least 1 Elder (pval.)", `testpval1', ///
"At least 1 Gossip=At least 1 High \textit{DC} Seed (pval.)", `testpval2', ///
"At least 1 Elder=At least 1 High \textit{DC} Seed (pval.)", `testpval3')
 

reg CallsRec dcAtleast1 num_seeds `gossip_count_control' `elder_count_control' if villageid != . & villageid != 62,  ro 
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
outreg2 using "$path_tables/Table_SeedDC_OLS_flyer_v2", keep(dcAtleast1) noaster tex(frag) nocons nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean')

 

reg NormalizedCalls_Ex atleast_gossip_1 atleast_elder_1 num_seeds `gossip_count_control' `elder_count_control'  if villageid != . & villageid != 62,  ro
//su NormalizedCalls_Ex if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"
	
	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
di "Control Mean is `control_mean'"
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01) 
outreg2 using "$path_tables/Table_SeedDC_OLS_flyer_v2", keep(atleast_gossip_1 atleast_elder_1) nocons noaster tex(frag)  nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean', "At least 1 Gossip=At least 1 Elder (pval.)", `testpval1')

reg NormalizedCalls_Ex atleast_gossip_1 atleast_elder_1 dcAtleast1 num_seeds `gossip_count_control' `elder_count_control'  if villageid != . & villageid != 62,  ro
//su NormalizedCalls_Ex if e(sample) & atleast_gossip_1 == 0 & atleast_elder_1 == 0 & dcAtleast == 0
test atleast_gossip_1 == atleast_elder_1 
loc testpval1 = round(`r(p)',0.01)
test atleast_gossip_1 == dcAtleast1 
loc testpval2 = round(`r(p)',0.01)
test atleast_elder_1 == dcAtleast1
loc testpval3 = round(`r(p)',0.01) 
outreg2 using "$path_tables/Table_SeedDC_OLS_flyer_v2", keep(atleast_gossip_1 atleast_elder_1 dcAtleast1) nocons  noaster nor2 tex(frag) nonotes label append ///
addstat("Control Group Mean", `control_mean', ///
"At least 1 Gossip=At least 1 Elder (pval.)", `testpval1', ///
"At least 1 Gossip=At least 1 High \textit{DC} Seed (pval.)", `testpval2', ///
"At least 1 Elder=At least 1 High \textit{DC} Seed (pval.)", `testpval3')


reg NormalizedCalls_Ex dcAtleast1 num_seeds `gossip_count_control' `elder_count_control' if villageid != . & villageid != 62,  ro 
loc cons = _b[_cons]
di "Constant is `cons'"

foreach var in num_seeds `gossip_count_control' `elder_count_control' {
	loc coef_`var' = _b[`var']
	di "Coefficent of `var' is `coef_`var''"

	qui su `var' if flierVillage ==0 & villageid != . & villageid != 62
	loc mean_`var' `r(mean)'
	di "Mean of `var' is `mean_`var''"

	loc est_`var' = `coef_`var'' * `mean_`var''
	di "est of `var' is `est_`var''"
}

loc control_mean = `cons' + `est_num_seeds' + `est_gossip_miss' + `est_elder_miss' + `est_both_miss' + `est_gossip_ct' + `est_elder_ct' 
					
di "Control Mean is `control_mean'" 
outreg2 using "$path_tables/Table_SeedDC_OLS_flyer_v2", keep(dcAtleast1) noaster tex(frag) nocons nonotes nor2 label append ///
addstat("Control Group Mean", `control_mean')





*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------
/* TABLE G.1 */
*-------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------

use "$path_data/all_hh_data", clear 

su scst_dummy if gossip==0 & elder==0 
loc constant = r(mean)
reg scst_dummy gossip elder, cluster(villageid)
test gossip = elder
loc pval1 = r(p)
outreg2 using "$path_tables/Table_GossipvsElder", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes ctitle("SCST") replace addstat("Random Household Mean", `constant', "Gossip = Elder p-val", `pval1')


su labor_dummy if gossip==0 & elder==0 
loc constant = r(mean)
reg labor_dummy gossip elder, cluster(villageid)
test gossip = elder
loc pval1 = r(p)
outreg2 using "$path_tables/Table_GossipvsElder", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes ctitle("Laborer") append addstat("Random Household Mean", `constant', "Gossip = Elder p-val", `pval1')


su land_owner if gossip==0 & elder==0 
loc constant = r(mean)
reg land_owner gossip elder, cluster(villageid)
test gossip = elder
loc pval1 = r(p)
outreg2 using "$path_tables/Table_GossipvsElder", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes ctitle("Land Owner") append addstat("Random Household Mean", `constant', "Gossip = Elder p-val", `pval1')


su have_electricity if gossip==0 & elder==0 
loc constant = r(mean)
reg have_electricity gossip elder, cluster(villageid)
test gossip = elder
loc pval1 = r(p)
outreg2 using "$path_tables/Table_GossipvsElder", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes ctitle("Electrified") append addstat("Random Household Mean", `constant', "Gossip = Elder p-val", `pval1')


su private_electricity if gossip==0 & elder==0 
loc constant = r(mean)
reg private_electricity gossip elder, cluster(villageid)
test gossip = elder
loc pval1 = r(p)
outreg2 using "$path_tables/Table_GossipvsElder", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes ctitle("Private Electricity") append addstat("Random Household Mean", `constant', "Gossip = Elder p-val", `pval1')


su owned_house if gossip==0 & elder==0 
loc constant = r(mean)
reg owned_house gossip elder if num_rooms <= 9, cluster(villageid)
test gossip = elder
loc pval1 = r(p)
outreg2 using "$path_tables/Table_GossipvsElder", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes ctitle("Own House") append addstat("Random Household Mean", `constant', "Gossip = Elder p-val", `pval1')


su num_rooms if gossip==0 & elder==0 
loc constant = r(mean)
reg num_rooms gossip elder if num_rooms <=9, cluster(villageid)
test gossip = elder
loc pval1 = r(p)
outreg2 using "$path_tables/Table_GossipvsElder", keep(gossip elder) label tex(frag) nocons noaster nor2 nonotes ctitle("No. of Rooms") append addstat("Random Household Mean", `constant', "Gossip = Elder p-val", `pval1')







