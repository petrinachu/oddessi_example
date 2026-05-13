/*
	ICTMC 2026 POSTER: ODDESSI EXAMPLE
	AUTHOR: PETRINA CHU
	STATA VERSION 18
	Version 0.0 23jan2026
*/
cap log close
log using "ictmc2026_oddessi_log", text replace

** USING EXAMPLE OF ANY RE-REFFERAL TO CRISIS SERVICES OUTCOME
use "C:\Users\k2150120\OneDrive - King's College London\Documents\ODDESSI\out\rr_analysis2.dta", clear
logit any_reref i.arm ib5.catchment list_size imd, or
*get marginal predicted probabilities of any re-referral in each arm
margins, at(arm=(0 1)) post
	*if odds ratio is the target estimand:
	nlcom (_b[2._at]/(1-_b[2._at]))/(_b[1._at]/(1-_b[1._at]))
		/*
			PLEASE NOTE THAT NLCOM P-VALUES FOR RATIOS WILL NOT MAKE SENSE
			NULL HYPOTHESIS IS THAT TRANSFORMED COEFFICIENT = 0
			(see page 9: https://www.stata.com/manuals13/rnlcom.pdf)
			If a p-value is desired, you need to use 'testnl'
		*/
	*if risk difference is the target estimand:
	nlcom (_b[2._at]) - (_b[1._at])
	*or risk ratio:
	nlcom (_b[2._at])/(_b[1._at])
		*in case we want to compare against glm for r
		//glm any_reref i.arm ib5.catchment list_size imd, fam(bin) link(log) eform
		//WARNING: CONVERGENCE NOT ACHIEVED
	
** USING EXAMPLE OF TOTAL DAYS IN PSYCH INPATIENT OUTCOME
	// ZERO-INFLATED PORTION TECHNICALLY GIVES ODDS OF ZERO DAYS
	// BUT TEAM WAS INTERESTED IN ODDS OF ANY DAYS
use "C:\Users\k2150120\OneDrive - King's College London\Documents\ODDESSI\out\rr_analysis2.dta", clear
	*gen exposure variable
	//(date of index crisis [RR_01] to last record review [RR_05])
	gen timein = RR_05 - RR_01
zinb indays_tot i.arm ib5.catchment list_size imd, exposure(timein) irr inflate(i.arm ib5.catchment list_size imd)
*get marginal predicted probabilities of ZERO psych inpatient days in each arm
margins, at(arm=(0 1)) predict(pr) post
	*if you want mOR of ZERO psych inpatient days:
	nlcom (_b[2._at]/(1-_b[2._at]))/(_b[1._at]/(1-_b[1._at]))
	*then get marginal predicted probabilities of ANY psych inpatient days
	nlcom (1-_b[1._at])
	nlcom (1-_b[2._at])
		//can present these probabilities if desired, or extend as needed
		*mOR for ANY psych inpatient days:
		nlcom ((1-_b[2._at])/_b[2._at])/((1-_b[1._at])/_b[1._at])
		*risk difference for ANY psych inpatient days:
		nlcom (1-_b[2._at]) - (1-_b[1._at])
	
log close