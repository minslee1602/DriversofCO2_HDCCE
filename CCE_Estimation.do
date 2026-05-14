*log using "CCE_Estimations.log", replace
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
capture which estout
if _rc ssc install estout

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


* Declare the panel using the integer year
xtset country_id year_int

* Rename long variable names
rename fossil_fuel_production_pc 	fossil_prod
rename fossil_fuel_consumption_pc 	fossil_cons
rename renewables_consumption_pc 	ren_cons
rename renewables_production_pc 	ren_prod

* Fixed effects panel
xtreg carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    fossil_prod fossil_cons ///
    ren_prod ren_cons, ///
    fe vce(cluster country_id)
estimates store fe_full

* Test CD on FE residuals
predict resid_fe, e
xtcd2 resid_fe
xtcse2 resid_fe


* CCEMG Estimator
xtdcce2 carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    fossil_prod fossil_cons ///
    ren_prod ren_cons, ///
    crosssectional(_all) cr_lags(0) mgmissing

* Test CD on residuals
predict resid_ccemg, residuals
xtcd2 resid_ccemg
xtcse2 resid_ccemg

* CCE Pooled Estimator
xtdcce2 carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    fossil_prod fossil_cons ren_prod ren_cons, ///
    crosssectional(_all) cr_lags(0) ///
    pooled(gdp_pc gdpsq_pc fd_index ///
        fossil_prod fossil_cons ren_prod ren_cons) ///
    mgmissing
estimates store ccep_full
	
* Test CD on CCEP residuals
predict resid_ccep, residuals
xtcd2 resid_ccep
xtcse2 resid_ccep

* ---------------------- HETEROGENEITY ANALYSIS ------------------------------
* Compute cross-sectional averages from the FULL 43-country sample.
* We import these averages into the CCEP equation and run OLS filtered by each country subgroup.
foreach var of varlist carbon_emissions_pc gdp_pc gdpsq_pc fd_index ///
                       fossil_prod fossil_cons ren_prod ren_cons {
    bysort year_int: egen csa_`var' = mean(`var')
}

* Group-wise OLS augmented with full-sample CS averages.

* Group 1: G7
xtreg carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    fossil_prod fossil_cons ren_prod ren_cons ///
    csa_carbon_emissions_pc csa_gdp_pc csa_gdpsq_pc csa_fd_index ///
    csa_fossil_prod csa_fossil_cons csa_ren_prod csa_ren_cons ///
    if group == 1, fe vce(cluster country_id)
estimates store ols_g7

* Group 2: BRICS
xtreg carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    fossil_prod fossil_cons ren_prod ren_cons ///
    csa_carbon_emissions_pc csa_gdp_pc csa_gdpsq_pc csa_fd_index ///
    csa_fossil_prod csa_fossil_cons csa_ren_prod csa_ren_cons ///
    if group == 2, fe vce(cluster country_id)
estimates store ols_brics

* Group 3: Carbon tax (non-G7/BRICS)
xtreg carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    fossil_prod fossil_cons ren_prod ren_cons ///
    csa_carbon_emissions_pc csa_gdp_pc csa_gdpsq_pc csa_fd_index ///
    csa_fossil_prod csa_fossil_cons csa_ren_prod csa_ren_cons ///
    if group == 3, fe vce(cluster country_id)
estimates store ols_ctax

* Group 4: Other
xtreg carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    fossil_prod fossil_cons ren_prod ren_cons ///
    csa_carbon_emissions_pc csa_gdp_pc csa_gdpsq_pc csa_fd_index ///
    csa_fossil_prod csa_fossil_cons csa_ren_prod csa_ren_cons ///
    if group == 4, fe vce(cluster country_id)
estimates store ols_other

* Comparison table 
esttab ccep_full ols_g7 ols_brics ols_ctax ols_other, ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(gdp_pc gdpsq_pc fd_index fossil_prod fossil_cons ///
         ren_prod ren_cons) ///
    mtitles("Full CCEP" "G7" "BRICS" "Carbon Tax" "Other") ///
    title("Group-wise OLS Estimates Augmented with Cross-Sectional Averages") ///
    stats(N N_g r2_w, labels("Observations" "Countries" "Within R²")) ///


* CD tests on residuals from each group regression
display _newline "CD TESTS ON GROUP RESIDUALS"
foreach grp in g7 brics ctax other {
    estimates restore ols_`grp'
    predict resid_`grp', e
    display _newline "--- `grp' ---"
    xtcd2 resid_`grp'
}

*log close
