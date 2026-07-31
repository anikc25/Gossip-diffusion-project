***********************************************
*Title: Analysis for Cell Phone RCT 
*Date Last Modfied: October, 2018
*Author(s): Arun Chandrasekhar, Devika Lakhote, Francisco Munoz
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



use karnataka_gossip_processed.dta


*********************************************************************************
// REGRESSION TABLES 
*********************************************************************************
cd "$path_tables"

************
*EVENT
************
replace leader = 0 if leader==. // fix missing leadership data as 0s

gen Y = Nnom_event
label var Y "Event"
gen Yround = round(Y)

*1. Table 6, Panel A - Poisson
xi: poisson Yround DC7_norm, cluster(village)
outreg2 using tablePoiA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)   label tex(frag)  excel bdec(3) replace   nonotes nocons noaster nor2
xi: poisson Yround degree_norm , cluster(village)
outreg2  using tablePoiA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons  noaster nor2
xi: poisson Yround eig_norm, cluster(village)
outreg2  using tablePoiA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons  noaster nor2
xi: poisson Yround leader, cluster(village)
outreg2 using tablePoiA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons  noaster nor2
xi: poisson Yround gps_sd  , cluster(village)
outreg2 using tablePoiA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons  noaster nor2

*2. Table C.1, Panel A - OLS
xi:reg Y DC7_norm  , cluster(village)
outreg2  using tableOLSA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) replace   nonotes nocons   noaster nor2
xi:reg Y degree_norm  , cluster(village)
outreg2  using tableOLSA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y eig_norm  , cluster(village)
outreg2  using tableOLSA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y leader  , cluster(village)
outreg2  using tableOLSA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)   label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y gps_sd  , cluster(village)
outreg2 using tableOLSA_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2

*3. Table C.2, Panel A - OLS + Lasso
xi:reg Y DC7_norm degree_norm, cluster(village)
outreg2 using tableOLSC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) replace   nonotes nocons   noaster nor2 
xi:reg Y DC7_norm eig_norm , cluster(village)
outreg2 using tableOLSC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y DC7_norm  leader  , cluster(village)
outreg2 using tableOLSC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y DC7_norm  gps_sd  , cluster(village)
outreg2 using tableOLSC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y DC7_norm  degree_norm eig_norm , cluster(village)
outreg2 using tableOLSC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
** post-lasso (Belloni, Chernozhukov and Hansen)  
/*
. di 2*1.1*sqrt(6466)*invnormal(1-0.05/(2*5))
455.67753
*/
lassoShooting Y DC7_norm eig_norm degree_norm leader gps_sd, lasiter(100) lambda(455.67753)
/*
          Selected|          LASSO     Post-LASSO
------------------+--------------------------------
          DC7_norm|      .01049831       .2930833
*/
xi: reg Y DC7_norm, cluster(village)
outreg2 using tableOLSC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append addtext(Post-LASSO, \checkmark  ) nocons nonotes  noaster nor2



*4.1 Table 7, Panel A - Poisson + Lasso
xi: poisson Yround DC7_norm degree_norm, cluster(village)
outreg2 using tablePoiC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) replace   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm eig_norm   , cluster(village)
outreg2 using tablePoiC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson  Yround DC7_norm  leader   , cluster(village)
outreg2 using tablePoiC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)   label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm  gps_sd  , cluster(village)
outreg2 using tablePoiC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm degree_norm eig_norm  , cluster(village)
outreg2 using tablePoiC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
  * using OLS Lasso for first stage selection (above)
xi: poisson Yround DC7_norm , cluster(village)
outreg2 using tablePoiC_diamE_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append addtext(Post-LASSO,  \checkmark ) nonotes nocons   noaster nor2

*4.2 Panel C - Only Poisson (No Lasso - for main body) 
xi: poisson Yround DC7_norm degree_norm, cluster(village)
outreg2 using tablePoiC_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) replace   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm eig_norm   , cluster(village)
outreg2 using tablePoiC_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson  Yround DC7_norm  leader   , cluster(village)
outreg2 using tablePoiC_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)   label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm  gps_sd  , cluster(village)
outreg2 using tablePoiC_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm degree_norm eig_norm  , cluster(village)
outreg2 using tablePoiC_diamE, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Event)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2


***************
*LOAN
***************
drop Y Yround
gen Y = Nnom_loan
label var Y "Loan"
gen Yround = round(Y)

*1. Table 6, Panel B - Poisson 
xi: poisson Yround DC7_norm, cluster(village)
outreg2 using tablePoiA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) replace  nonotes nocons   noaster nor2
xi: poisson Yround degree_norm , cluster(village)
outreg2  using tablePoiA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes nocons   noaster nor2
xi: poisson Yround eig_norm, cluster(village)
outreg2  using tablePoiA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes nocons   noaster nor2
xi: poisson Yround leader, cluster(village)
outreg2 using tablePoiA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes nocons   noaster nor2
xi: poisson Yround gps_sd  , cluster(village)
outreg2 using tablePoiA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes nocons   noaster nor2

*2. Table C.1, Panel B - OLS
xi:reg Y DC7_norm  , cluster(village)
outreg2  using tableOLSA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Loan)  label tex(frag)  excel bdec(3) replace  nonotes nocons  noaster nor2
xi:reg Y degree_norm  , cluster(village)
outreg2  using tableOLSA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append nonotes nocons  noaster nor2 
xi:reg Y eig_norm  , cluster(village)
outreg2  using tableOLSA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes nocons   noaster nor2
xi:reg Y leader  , cluster(village)
outreg2  using tableOLSA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes nocons  noaster nor2 
xi:reg Y gps_sd  , cluster(village)
outreg2 using tableOLSA_diamL, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append nonotes nocons   noaster nor2


*3. Table C.2, Panel B - OLS + Lasso 
xi:reg Y DC7_norm degree_norm, cluster(village)
outreg2 using tableOLSC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) replace   nonotes nocons   noaster nor2
xi:reg Y DC7_norm eig_norm , cluster(village)
outreg2 using tableOLSC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y DC7_norm  leader  , cluster(village)
outreg2 using tableOLSC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y DC7_norm  gps_sd  , cluster(village)
outreg2 using tableOLSC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi:reg Y DC7_norm degree_norm eig_norm  , cluster(village)
outreg2 using tableOLSC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes  nocons   noaster nor2
** post-lasso (Belloni, Chernozhukov and Hansen)
/*
. di 2*1.1*sqrt(6466)*invnormal(1-0.05/(2*5))
455.67753
*/
lassoShooting Y DC7_norm eig_norm degree_norm leader gps_sd, lasiter(100) lambda(455.67753)
/*
          Selected|          LASSO     Post-LASSO
------------------+--------------------------------
          DC7_norm|      .05871772      .32468549
       degree_norm|      .07529575      .09224629
*/
xi: reg Y DC7_norm degree_norm, cluster(village)
outreg2 using tableOLSC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) append addtext(Post-LASSO,  \checkmark ) nocons   noaster nor2


*4 Table 7, Panel B - Poisson + Lasso 
xi: poisson Yround DC7_norm degree_norm, cluster(village)
outreg2 using tablePoiC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd)  ctitle(Loan)  label tex(frag)  excel bdec(3) replace  nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm eig_norm   , cluster(village)
outreg2 using tablePoiC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Loan)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson  Yround DC7_norm  leader   , cluster(village)
outreg2 using tablePoiC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Loan)   label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm  gps_sd  , cluster(village)
outreg2 using tablePoiC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Loan)  label tex(frag)  excel bdec(3) append   nonotes nocons   noaster nor2
xi: poisson Yround DC7_norm degree_norm eig_norm , cluster(village)
outreg2 using tablePoiC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Loan)  label tex(frag)  excel bdec(3) append  nonotes nocons   noaster nor2
  * using OLS Lasso for first stage selection (above)
xi: poisson Yround DC7_norm degree_norm , cluster(village)
outreg2 using tablePoiC_diamL_lasso, keep(DC7_norm eig_norm degree_norm leader gps_sd) ctitle(Loan)  label tex(frag)  excel bdec(3) append addtext(Post-LASSO,  \checkmark  ) nonotes nocons   noaster nor2


 

*********************************************************************************
//  FIGURES 4 & 5
*********************************************************************************



************
*EVENT
************


cumul DC7_norm if leader_Dnom_e == 1, gen(cumul_leader_nom_e)
sort DC7_norm

cumul DC7_norm if Not_leader_Dnom_e == 1, gen(cumul_notleader_nom_e)
sort DC7_norm

cumul DC7_norm if leader_Not_Dnom_e == 1, gen(cumul_leader_notnom_e)
sort DC7_norm

cumul DC7_norm if Not_leader_Not_Dnom_e == 1, gen(cumul_notleader_notnom_e)
sort DC7_norm

label var cumul_leader_nom_e "Nominated, Leader"
label var cumul_notleader_nom_e "Nominated, Not Leader"
label var cumul_leader_notnom_e "Not Nominated, Leader"
label var cumul_notleader_notnom_e "Not Nominated, Not Leader"

preserve 
	stack cumul_leader_nom_e DC7_norm ///
	cumul_notleader_nom_e DC7_norm ///
	cumul_leader_notnom_e DC7_norm ///
	cumul_notleader_notnom_e DC7_norm , into(c temp) wide clear
	
	label var temp "Diffusion Centrality" 
	twoway line cumul_leader_nom_e cumul_notleader_nom_e cumul_leader_notnom_e cumul_notleader_notnom_e temp, lpattern(solid solid dash dash) ///
	lwidth(thick thin thick thin) ///
	legend(order(1 "Nominated, Leader" 2 "Nominated, Not Leader" 3 "Not Nominated, Leader" 4 "Not Nominated, Not Leader")) legend(size(*.6)) ///
	ytitle("F(x)") ///
graphregion(color(white))
	graph export CDFs_NomLeader_Event.pdf, replace
restore


************
*LOAN
************

cumul DC7_norm if leader_Dnom_l == 1, gen(cumul_leader_nom_l)
sort DC7_norm

cumul DC7_norm if Not_leader_Dnom_l == 1, gen(cumul_notleader_nom_l)
sort DC7_norm

cumul DC7_norm if leader_Not_Dnom_l == 1, gen(cumul_leader_notnom_l)
sort DC7_norm

cumul DC7_norm if Not_leader_Not_Dnom_l == 1, gen(cumul_notleader_notnom_l)
sort DC7_norm

label var cumul_leader_nom_l "Nominated, leader"
label var cumul_notleader_nom_l "Nominated, not leader"
label var cumul_leader_notnom_l "Not nominated, leader"
label var cumul_notleader_notnom_l "Not nominated, not leader"

preserve 
	stack cumul_leader_nom_l DC7_norm ///
	cumul_notleader_nom_l DC7_norm ///
	cumul_leader_notnom_l DC7_norm ///
	cumul_notleader_notnom_l DC7_norm , into(c temp) wide clear
	
	label var temp "Diffusion Centrality" 
	line cumul_leader_nom_l cumul_notleader_nom_l cumul_leader_notnom_l cumul_notleader_notnom_l temp, lpattern(solid solid dash dash) ///
	lwidth(thick thin thick thin) ///
	legend(order(1 "Nominated, Leader" 2 "Nominated, not leader" 3 "Not nominated, leader" 4 "Not nominated, not leader")) legend(size(*.6)) ///
	ytitle("F(x)") ///
graphregion(color(white))
	graph export CDFs_NomLeader_Loan.pdf, replace
restore







*********************************************************************************
//  TABLE 8
*********************************************************************************

use "$path_data/nomination_data.dta", clear

su nominated
loc loc_mean = r(mean)

xi3: reg nominated perc_H_ji, cluster(village)
outreg2 using table_network_gossip, keep(perc_H_ji) ctitle(Nominated)   label tex(frag)  excel bdec(3) replace  addstat("Dep. var mean", `loc_mean') nonotes nocons noaster nor2
xi3: reghdfe nominated perc_H_ji, absorb(unique_id) cluster(village)
outreg2 using table_network_gossip, keep(perc_H_ji) ctitle(Nominated)   label tex(frag)  excel bdec(3) append addtext(Respondent FE, \checkmark  ) addstat("Dep. var mean", `loc_mean')   nonotes nocons noaster nor2
xi3: reg nominated perc_H_ji z*z*z*z*z, cluster(village)
outreg2 using table_network_gossip, keep(perc_H_ji) ctitle(Nominated)   label tex(frag)  excel bdec(3) append addtext(Flexible Controls for DC, \checkmark  ) addstat("Dep. var mean", `loc_mean')   nonotes nocons noaster nor2

xi3: reghdfe nominated perc_H_ji z*z*z*z*z, absorb(unique_id) cluster(village)
outreg2 using table_network_gossip, keep(perc_H_ji) ctitle(Nominated)   label tex(frag)  excel bdec(3) append addtext(Respondent FE, \checkmark, Flexible Controls for DC, \checkmark  ) addstat("Dep. var mean", `loc_mean')   nonotes nocons noaster nor2
xi3: reghdfe nominated perc_H_ji, absorb(unique_id_j) cluster(village)
outreg2 using table_network_gossip, keep(perc_H_ji) ctitle(Nominated)   label tex(frag)  excel bdec(3) append addtext(Rankee FE, \checkmark  )  addstat("Dep. var mean", `loc_mean')  nonotes nocons noaster nor2
xi3: reghdfe nominated perc_H_ji, absorb(unique_id unique_id_j) cluster(village)
outreg2 using table_network_gossip, keep(perc_H_ji) ctitle(Nominated)   label tex(frag)  excel bdec(3) append addtext(Respondent FE, \checkmark, Rankee FE, \checkmark  ) addstat("Dep. var mean", `loc_mean')    nonotes nocons noaster nor2





*********************************************************************************
//  TABLE F.1
*********************************************************************************

gen p98H = H_ji
gen p99H = H_ji
gen DC_i = H_ji
collapse (mean) village (sum) nominated DC_i   (p99) p99H (p98) p98H , by(unique_id)

label var p99 "99th percentile of $ NG_{ji} $"
label var p98 "98th percentile of $ NG_{ji} $"

xi3: logit nominated p99 DC_i if DC_i~=0, cluster(village)
outreg2 using table_dontknow, keep(p99) ctitle(Nominated Anyone)   label tex(frag)  excel bdec(3) replace  addtext(99th Percentile, \checkmark  )  nonotes nocons noaster nor2
mfx
xi3: clogit nominated p99 DC_i if DC_i~=0, cluster(village) group(village)
outreg2 using table_dontknow, keep(p99) ctitle(Nominated Anyone)   label tex(frag)  excel bdec(3) append addtext(99th Percentile, \checkmark, Village FE, \checkmark   )   nonotes nocons noaster nor2
xi3: logit nominated p98 DC_i if DC_i~=0, cluster(village)
outreg2 using table_dontknow, keep(p98) ctitle(Nominated Anyone)   label tex(frag)  excel bdec(3) append addtext(98th Percentile, \checkmark )  nonotes nocons noaster nor2
mfx
xi3: clogit nominated p98 DC_i if DC_i~=0, cluster(village) group(village)
outreg2 using table_dontknow, keep(p98) ctitle(Nominated Anyone)   label tex(frag)  excel bdec(3) append addtext(98th Percentile, \checkmark, Village FE, \checkmark   )  nonotes nocons noaster nor2





*********************************************************************************
//  TABLE F.2
*********************************************************************************

use "$path_data/nominating_features.dta", clear

areg didnominate_loan DC_std Nnom_loan Nnom_event leader scst have_el private_el own_house num_rooms land labor business gps female if num_rooms <=8, absorb(village) cluster(village)
outreg2 using table_drive_nomination, keep(DC_std Nnom_loan Nnom_event leader scst have_el private_el own_house num_rooms land labor business gps female) ctitle(Nominates Someone (Loan))   label tex(frag)  excel bdec(3) replace   nonotes nocons noaster nor2

areg didnominate_event DC_std Nnom_loan Nnom_event leader scst have_el private_el own_house num_rooms land labor business gps female if num_rooms <=8, absorb(village) cluster(village)
outreg2 using table_drive_nomination, keep(DC_std Nnom_loan Nnom_event leader scst have_el private_el own_house num_rooms land labor business gps female) ctitle(Nominates Someone (Event))   label tex(frag)  excel bdec(3) append   nonotes nocons noaster nor2




