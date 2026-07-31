

clear all
set more off
set scheme s1color
cap log close
pause off

global yourdropbox "C:/Users/asankar/Dropbox (MIT)"
global output "${yourdropbox}/HIM DATASETS/Analysis/ArunEstherPaper/tables"

use "${yourdropbox}/HIM DATASETS/Analysis/ArunEstherPaper/data/haryana_rct_Table2_seedslevel.dta"



levelsof commtreat_arm, loc(lvls_commt)
cap file close sumstats
file open sumstats using "${output}/Table2.tex", write replace
file write sumstats"\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \small \begin{tabular}{@{\extracolsep{2.7pt}}{l}*{4}{>{\centering\arraybackslash}m{2.0cm}}@{}}  \toprule" _n      
*file write sumstats"\addlinespace[2mm]" _n
file write sumstats "& Random Seed & Gossip Seed & Trusted Seed  & Trusted Gossip Seed \\" _n
file write sumstats " & (1) & (2) &  (3) & (4)  \\" _n
file write sumstats "\midrule" _n 
file write sumstats "\multicolumn{2}{c}{\emph{Nominations Statistics (per village)}}  & & & \\" _n

preserve
loc vil_outcomes "nnoms_byvillage nomstop6"
use "${yourdropbox}/HIM DATASETS/Analysis/ArunEstherPaper/data/haryana_rct_Table2_villagelevel.dta", clear 
 foreach y in `vil_outcomes' {
	local varlab: variable label `y' 
 	foreach j in `lvls_commt' {
 		if `j'!=1 {
	 		sum `y' if commtreat_arm==`j' 
	 		loc mean_`j' = string(r(mean), "%9.3f")
	 		loc sd_`j' =  string(r(sd), "%9.3f") 
 		}

 	}
 		file write sumstats " `varlab' & . & `mean_3' & `mean_4'  & `mean_5' \\ " _n     
 		file write sumstats "  & . & (`sd_3') & (`sd_4')  & (`sd_5') \\ " _n     
 
 } 


restore 

file write sumstats "\multicolumn{2}{c}{\emph{Seed Characteristics}}  & & & \\" _n
loc outcomes "refused seedage female education ownland  asset_index hindu muslim caste_stsc caste_obc  member_panchayat num_or_chauk  interact_1  parti_1  f7_aware_immu_camp_1  f1_anm_aware_1 f4_aware_asha_1 "

 foreach y in `outcomes' {
 	 eststo clear 
 	 // (1)
	 local varlab: variable label `y' 

 	foreach j in `lvls_commt' {
 		sum `y' if commtreat_arm==`j' 
 		loc mean_`j' = string(r(mean), "%9.3f")
 		loc sd_`j' =  string(r(sd), "%9.3f") 

 	}
 		file write sumstats " `varlab' & `mean_1' & `mean_3' & `mean_4'  & `mean_5' \\ " _n     
 		file write sumstats "  & (`sd_1') & (`sd_3') & (`sd_4')  & (`sd_5') \\ " _n     
 
 } 

 	foreach j in `lvls_commt' {
 	 	count if commtreat_arm==`j'  
 		loc obs_`j' = string(`r(N)', "%9.0f")
 	}

 	file write sumstats "\addlinespace[1mm]" _n           
 	file write sumstats	"\midrule" _n 
 	file write sumstats " Observations & `obs_1' & `obs_3' & `obs_4'  & `obs_5' \\ " _n     



file write sumstats "\addlinespace[1mm]" _n
file write sumstats "\bottomrule" _n              
file write sumstats "\end{tabular}" _n 
file close sumstats
