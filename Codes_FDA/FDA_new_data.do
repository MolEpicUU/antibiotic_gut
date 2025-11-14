*cd "/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data"
*import delimited "/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/scapis_dataset_for_fda_ulf__full.tsv",clear

cd "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/"
foreach cohort in scapis mos simpler{
	foreach type in All J01CE J01CF J01CA J01CR J01FA J01D J01A J01FF J01MA J01E J01XE01{
		import delimited "/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data//`cohort'_dataset_for_fda_ulf__full.tsv",clear
		
		if "`cohort'"=="mos"{ 
			rename lopnrmos subject
			rename extraction_plate site_plate
			rename smoking smokestatus
			gen placebirth="Scandinavia"
			rename polypharmacy_cat polypharmacy12m_cat
		}	
		replace polypharmacy12m_cat=subinstr(polypharmacy12m_cat,"1-4","1_4",.)
		replace polypharmacy12m_cat=subinstr(polypharmacy12m_cat,"5 or more","5",.)
		
		if "`cohort'"=="simpler"{
			rename simpkey subject
			rename aliquoting_plate site_plate
			replace sex="1" if sex=="M"
			replace sex="2" if sex=="F"
			destring sex,replace
			foreach var in ppi statins metformin betablock ssri{
				replace `var'="0" if `var'=="yes"
				replace `var'="1" if `var'=="no"
				destring `var',replace
			}
			replace site_plate="frst_shipment" if site_plate=="1st shipment"
			replace smokestatus="former" if smokestatus=="ex"
		}
		
		gen dat_tmp=date(edatum,"YMD")
		bysort subject (dat_tmp): gen n=_n
		drop dat_tmp
		reshape wide class atc edatum,i(subject) j(n)
		isid subject

		local multiplier=1 //6 och vi får den richness som liknar PEK-et.
		local rang=30*`multiplier'
		local loopnr=98/`multiplier'
		
		//doublecheck so this is correct in MOS
		if "`cohort'"=="scapis"{ 
			replace sex="1" if sex=="male"
			replace sex="2" if sex=="female"
			destring sex,replace
			
			foreach var in ppi statins metformin betablock ssri{
				replace `var'="0" if `var'=="yes"
				replace `var'="1" if `var'=="no"
				destring `var',replace
			}
		}
		
		*confirm variable antipsycho
		
		replace site_plate=subinstr(site_plate,"ö","o",.)
		replace education=subinstr(education,"education_university","education_University",.)
		replace education=subinstr(education,"education_Upper Secondary","education_upper_secondary",.)
		replace site_plate=subinstr(site_plate,"Box","box",.)
		unab edat: edatum*
		local i=0
		foreach var of local edat{
			local i=`i'+1
			gen ed`i'=date(`var',"YMD")
			format ed`i' %td
		}
		
		gen tmp_date=date(visit1,"YMD")
		local i=0
		foreach var of local edat{
			local i=`i'+1
			gen diff`i'=(tmp_date-ed`i')
			gen diff_all`i'=.
			gen diff_`type'`i'=.
		}
		local n=`i'

		/*
		forv i=1/`loopnr'{
			gen all_anti`i'=0
			forv j=1/`n'{
				replace all_anti`i'=1 if inrange(diff`i',`rang'*(`i'-1),`rang'*(`i'))==1
			}
		}
		*/
		*assert all_anti99==0
		*if _rc==0 drop all_anti99
		
		forv i=1/`loopnr'{
			gen anti`i'=0
			forv j=1/`n'{
				gen ATC=(regexm(atc`j',"^`type'"))
				replace ATC=0 if "`type'"=="J01D" & !((strpos(atc`j',"J01DB")) | (strpos(atc`j',"J01DC")) | (strpos(atc`j',"J01DD")))
				replace ATC=1 if "`type'"=="All" & antibiotic==1
				replace anti`i'=1 if inrange(diff`j',`rang'*(`i'-1),`rang'*(`i'))==1 & ATC==1
				drop ATC
			}
		}
		noi count
		local N=r(N)
		*drop if diff_`type'==. & diff_all!=.
		noi di "`type' `cohort'"
		*noi tab ATC
		count
		assert `=r(N)'==`N' if "`type'"=="All"
		drop tmp_date
		*drop all_anti* 
		*replace richness=richness-50 if ATC==1 //just a check so the results make sense
		*count if ATC==1
		*drop if antibiotic==0
		
		*bsample _N
		isid subject
		rename antipsycho psych
		replace psych="1" if psych=="yes"
		replace psych="0" if psych=="no"
		destring psych,replace

		save "FDA_`cohort'_`type'.dta",replace
	}
}

cd "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/"

*J01CE J01CF J01CA J01CR J01FA J01D J01A J01FF J01MA J01E J01XE01
foreach cohort in scapis mos simpler{
	use "FDA_`cohort'_All.dta",clear
	drop eda* ed? ed?? diff* atc* class*
	rename anti* All*
	foreach type in J01CE J01CF J01CA J01CR J01A J01FF J01D J01XE01 J01E J01MA J01FA{
		*"J01CE","J01CF","J01CA","J01A","J01FF"
		noi di "`type' `cohort'"
		merge 1:1 subject using "FDA_`cohort'_`type'.dta",nogenerate
		rename anti* `type'*
		drop eda* ed? ed?? diff* atc* class*
	}
	
	forv i=1/98{
		egen J01ZZ`i'=rowmax(J01D`i' J01XE01`i' J01E`i' J01FA`i')
	}
	drop *biotic All* J01D* J01XE01* J01E* J01FA*
	save "FDA_`cohort'_combined.dta",replace
}
use "FDA_scapis_combined.dta",clear
