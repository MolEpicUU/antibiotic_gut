/*
We retrieved information on all dispensed prescriptions with the Anatomical
Therapeutic Chemical code J01 (antibacterials for systemic use) and classified them as betalactamase-sensitive penicillins (J01CE), beta-lactamase-resistant penicillins (J01CF), extendedspectrum penicillins (J01CA), penicillins combined with beta-lactamase inhibitors (J01CR),
cephalosporins (J01DB, J01DC, J01DD), macrolides (J01FA), tetracyclines (J01A),
lincosamides (J01FF), fluoroquinolones (J01MA), sulfonamides and trimethoprim (J01E), and
nitrofurantoin (J01XE01). 
*/

 /*simpler*/
set seed 476
cd "/proj/sens2019512/boot/"
foreach cohort in simpler scapis mos{
forv j=1/250{
		use "/home/ulfha881/Desktop/Ulf/Antibiotics/Results/FDA_`cohort'_combined.dta",clear
		bsample _N
		save "FDA_`cohort'`j'.dta",replace
	}
}
