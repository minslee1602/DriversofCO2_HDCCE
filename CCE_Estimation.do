log using "CCE_Estimations.log"

clear all
set more off

* Import Excel file
import excel "Cabon_Emm_disaggregate.xlsx", firstrow clear

* Check data
describe
summarize

* Panel Setup
encode CountryName, gen(country_id)

* Install packages
capture which xtcsd
if _rc ssc install xtcsd

capture which xtcse2
if _rc net install xtdcce2, from("https://janditzen.github.io/xtdcce2/")

* Convert date variable to a simple year integer
gen year_int = year(Year)

* Declare panel
xtset country_id year_int

* Create groups for heterogeneity check
destring G7 BRICS CarbonTax, replace force

gen group = .
replace group = 1 if G7 == 1
replace group = 2 if BRICS == 1
replace group = 3 if CarbonTax == 1 & G7 == 0 & BRICS == 0
replace group = 4 if G7 == 0 & BRICS == 0 & CarbonTax == 0

label define grouplbl 1 "G7" 2 "BRICS" 3 "Carbon tax (non-G7)" 4 "Other"
label values group grouplbl

* Demeaning GDP
sum gdp_pc
gen gdp_pc_dm = gdp_pc - r(mean)
gen gdpsq_pc_dm = gdp_pc_dm^2

* Declare the panel using the integer year
xtset country_id year_int

* Rename long variable names
rename fossil_fuel_production_pc 	fossil_prod
rename fossil_fuel_consumption_pc 	fossil_cons
rename renewables_consumption_pc 	ren_cons
rename renewables_production_pc 	ren_prod

* CCEMG Estimator
xtdcce2 carbon_emissions_pc ///
    gdp_pc_dm gdpsq_pc_dm fd_index ///
    fossil_prod fossil_cons ///
    ren_prod ren_cons, ///
    crosssectional(_all) cr_lags(0) mgmissing

* Test CD on residuals
predict resid_ccemg, residuals
xtcd2 resid_ccemg
xtcse2 resid_ccemg


* CCE Pooled Estimator
xtdcce2 carbon_emissions_pc ///
    gdp_pc_dm gdpsq_pc_dm fd_index ///
    fossil_prod fossil_cons ren_prod ren_cons, ///
    crosssectional(_all) cr_lags(0) ///
    pooled(gdp_pc_dm gdpsq_pc_dm fd_index ///
        fossil_prod fossil_cons ren_prod ren_cons) ///
    mgmissing
estimates store ccep_full
	
* Test CD 
predict resid_ccep, residuals
xtcd2 resid_ccep
xtcse2 resid_ccep


log close










