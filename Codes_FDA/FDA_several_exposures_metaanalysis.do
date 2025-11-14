

local title "title"
use "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/boot_data_shannon_mos_xJ01ZZ234.dta",clear

foreach type in invsimpson richness shannon{
	foreach anti_type in J01CE J01CF J01CA J01A J01FF J01MA{
	if "`title'"=="title"{
		if "`anti_type'"=="J01CE" local name "Penicillin V"
		if "`anti_type'"=="J01CF" local name "Flucloxacillin"
		if "`anti_type'"=="J01CA" local name "Extended-spectrum penicillin"
		if "`anti_type'"=="J01A" local name "Tetracyline"
		if "`anti_type'"=="J01FF" local name "Clindamycin"
		if "`anti_type'"=="J01MA" local name "Fluoroquinolones"
	}
		
		foreach var in mos scapis simpler{
			clear
			save "/home/ulfha881/Desktop/Ulf/Antibiotics/Results//`type'_`var'_meta.dta",emptyok replace
			forv i=1/250{
				cap use "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/boot_data_`type'_`var'_x`anti_type'`i'.dta",clear
				gen month=mod(_n,98)
				replace month=98 if month==0
				append using "/home/ulfha881/Desktop/Ulf/Antibiotics/Results//`type'_`var'_meta.dta"
				sleep 1000
				save "/home/ulfha881/Desktop/Ulf/Antibiotics/Results//`type'_`var'_meta.dta",replace
			}
				egen stderr=sd(V1),by(month)
				keep month stderr
				duplicates drop
				gen cohort="`var'"
				sleep 1000
				save "/home/ulfha881/Desktop/Ulf/Antibiotics/Results//`type'_`var'_meta.dta",replace
		}

		cap postclose myfile
		postfile myfile double month mean se using "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/meta_analysis_`type'_`anti_type'.dta",replace

		foreach var in scapis simpler mos{
			use "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/main_data_`var'_`type'_x`anti_type'.dta",clear
			gen cohort="`var'"
			gen month=_n
			save "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/main_data_`var'_`type'2.dta",replace
		}

		use "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/main_data_scapis_`type'2.dta",clear
		append using "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/main_data_simpler_`type'2.dta"
		append using "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/main_data_mos_`type'2.dta"
		save "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/main_data_`type'.dta",replace

		use "/home/ulfha881/Desktop/Ulf/Antibiotics/Results//`type'_scapis_meta.dta"
		append using "/home/ulfha881/Desktop/Ulf/Antibiotics/Results//`type'_mos_meta.dta"
		append using "/home/ulfha881/Desktop/Ulf/Antibiotics/Results//`type'_simpler_meta.dta"
		merge 1:1 cohort month using "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/main_data_`type'.dta",nogenerate keep(match)

		sum month
		forv i=1/`=r(max)'{
			preserve
			keep if month==`i'
			meta set y_vals stderr
			meta summarize,common(invvariance)
			post myfile (`i') (r(theta)) (r(se))
			restore
		}
		postclose myfile

		use "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/meta_analysis_`type'_`anti_type'.dta",clear
		gen hci=mean+1.96*se
		gen lci=mean-1.96*se
		
		gen year=month/12

		if "`type'"=="shannon"{
			twoway (line mean year,sort lcolor(black)) (line hci year,lpattern(longdash) lcolor(black) sort) (line lci year,lpattern(longdash) lcolor(black) sort), ///
			legend(off) yline(0,lcolor(gs12)) ytitle("Shannon-Weaver Index") xtitle("Time since last antibiotic (years)") saving("`type'`anti_type'",replace) xline(1,lcolor(gs12) lpattern(longdash)) xline(4,lcolor(gs12) lpattern(longdash)) xline(8,lcolor(gs12) lpattern(longdash)) xlabel(1/8) title("`name'")
			*graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI_`type'_`name'.tif",width(2500) replace
		}
		if "`type'"=="richness"{
			twoway (line mean year,sort lcolor(black)) (line hci year,lpattern(longdash) lcolor(black) sort) (line lci year,lpattern(longdash) lcolor(black) sort), ///
			legend(off) yline(0,lcolor(gs12)) ytitle("Species richness") xtitle("Time since last antibiotic (years)")  saving("`type'`anti_type'",replace) xline(1,lcolor(gs12) lpattern(longdash)) xline(4,lcolor(gs12) lpattern(longdash)) xline(8,lcolor(gs12) lpattern(longdash)) xlabel(1/8) title("`name'")
			*graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI_`type'_`name'.tif",width(2500) replace
		}
		if "`type'"=="invsimpson"{
			twoway (line mean year,sort lcolor(black)) (line hci year,lpattern(longdash) lcolor(black) sort) (line lci year,lpattern(longdash) lcolor(black) sort), ///
			legend(off) yline(0,lcolor(gs12)) ytitle("Inverse Simpson Index") xtitle("Time since last antibiotic (years)")  saving("`type'`anti_type'",replace) xline(1,lcolor(gs12) lpattern(longdash)) xline(4,lcolor(gs12) lpattern(longdash)) xline(8,lcolor(gs12) lpattern(longdash)) xlabel(1/8) title("`name'")
			*graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI_`type'_`name'.tif",width(2500) replace
		}
	}
}
graph combine "shannonJ01CE" "shannonJ01CF" "shannonJ01CA" "shannonJ01A" "shannonJ01FF" "shannonJ01MA",ycommon xcommon
graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI/CI_shannon`title'.tif",width(2500) replace
graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI/CI_shannon`title'.pdf",replace

graph combine "richnessJ01CE" "richnessJ01CF" "richnessJ01CA" "richnessJ01A" "richnessJ01FF" "richnessJ01MA",ycommon xcommon
graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI/CI_richness`title'.tif",width(2500) replace
graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI/CI_richness`title'.pdf",replace
graph combine "invsimpsonJ01CE" "invsimpsonJ01CF" "invsimpsonJ01CA" "invsimpsonJ01A" "invsimpsonJ01FF" "invsimpsonJ01MA",ycommon xcommon
graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI/CI_invsimpson`title'.tif",width(2500) replace
graph export "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/CI/CI/CI_invsimpson`title'.pdf",replace
