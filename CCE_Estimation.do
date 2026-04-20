log using "CCE_Estimation.log"

clear all
set more off

* Import Excel file
import excel "Cabon_Emm_disaggregate.xlsx", firstrow clear

* Check data
describe
summarize

* Panel Setup
encode CountryName, gen(country_id)

* Declare panel
xtset country_id Year

* Install packages
capture which xtcsd
if _rc ssc install xtcsd

capture which xtcse2
if _rc net install xtdcce2, from("https://janditzen.github.io/xtdcce2/")

* Convert date variable to a simple year integer
gen year_int = year(Year)

* Create groups for heterogeneity check
destring G7 BRICS CarbonTax, replace force

gen group = .
replace group = 1 if G7 == 1                                    // G7: N = 7
replace group = 2 if BRICS == 1                                 // BRICS: N = 5
replace group = 3 if CarbonTax == 1 & G7 == 0 & BRICS == 0      // CarbonTax non-G7/BRICS: N = 15
replace group = 4 if G7 == 0 & BRICS == 0 & CarbonTax == 0      // Other: N = 16

label define grouplbl 1 "G7" 2 "BRICS" 3 "Carbon tax (non-G7)" 4 "Other"
label values group grouplbl

* Verify country counts
preserve
bysort country_id: keep if _n == 1
tab group
restore


* Rename variables with long names
rename coal_consumption_pc         coal_cons_pc
rename oil_consumption_pc          oil_cons_pc
rename gas_consumption_pc          gas_cons_pc
rename renewables_consumption_pc   ren_cons_pc
rename coal_production_pc          coal_prod_pc
rename oil_production_pc           oil_prod_pc
rename gas_production_pc           gas_prod_pc
rename solar_production_pc         solar_prod_pc
rename wind_production_pc          wind_prod_pc
rename hydro_production_pc         hydro_prod_pc
rename other_renewables_production_pc othren_prod_pc

* Declare the panel using the integer year
xtset country_id year_int

* Pesaran CCE Test with cross-sectional averages on all. This CCE will break down.
// xtdcce2 carbon_emissions_pc ///
//     gdp_pc gdpsq_pc fd_index ///
//     coal_cons_pc oil_cons_pc gas_cons_pc ren_cons_pc ///
//     coal_prod_pc oil_prod_pc gas_prod_pc ///
//     solar_prod_pc wind_prod_pc hydro_prod_pc othren_prod_pc, ///
//     crosssectional( ///
//         carbon_emissions_pc ///
//         gdp_pc gdpsq_pc ///
//         fd_index ///
//         coal_cons_pc oil_cons_pc gas_cons_pc ren_cons_pc ///
//         coal_prod_pc oil_prod_pc gas_prod_pc ///
//         solar_prod_pc wind_prod_pc hydro_prod_pc othren_prod_pc ///
//     ) cr_lags(0) mgmissing
	
* Pesaran CCE Test cross-sectional averages only on carbon_emissions_pc gdp_pc fd_index
xtdcce2 carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    coal_cons_pc oil_cons_pc gas_cons_pc ren_cons_pc ///
    coal_prod_pc oil_prod_pc gas_prod_pc ///
    solar_prod_pc wind_prod_pc hydro_prod_pc othren_prod_pc, ///
    crosssectional(carbon_emissions_pc gdp_pc fd_index) cr_lags(0) mgmissing

estat ebistructure

* Test if CD still exists
predict resid_cce, residuals
xtcd2 resid_cce
xtcse2 resid_cce

* Carbon tax non-G7/BRICS (N = 15)
xtdcce2 carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    coal_cons_pc oil_cons_pc gas_cons_pc ren_cons_pc ///
    coal_prod_pc oil_prod_pc gas_prod_pc ///
    solar_prod_pc wind_prod_pc hydro_prod_pc othren_prod_pc ///
    if group == 3, ///
    crosssectional(carbon_emissions_pc gdp_pc fd_index) cr_lags(0)
estimates store cce_tax

* Non Carbon Tax, Non G7, Non BRICS(N = 16)
xtdcce2 carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    coal_cons_pc oil_cons_pc gas_cons_pc ren_cons_pc ///
    coal_prod_pc oil_prod_pc gas_prod_pc ///
    solar_prod_pc wind_prod_pc hydro_prod_pc othren_prod_pc ///
    if group == 4, ///
    crosssectional(carbon_emissions_pc gdp_pc fd_index) cr_lags(0)
estimates store cce_other


* G7 (N = 7)
xtdcce2 carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    coal_cons_pc oil_cons_pc gas_cons_pc ren_cons_pc ///
    coal_prod_pc oil_prod_pc gas_prod_pc ///
    solar_prod_pc wind_prod_pc hydro_prod_pc othren_prod_pc ///
    if group == 1, ///
    crosssectional(carbon_emissions_pc gdp_pc fd_index) cr_lags(0) 
estimates store cce_g7


* BRICS (N = 5)
xtdcce2 carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    coal_cons_pc oil_cons_pc gas_cons_pc ren_cons_pc ///
    coal_prod_pc oil_prod_pc gas_prod_pc ///
    solar_prod_pc wind_prod_pc hydro_prod_pc othren_prod_pc ///
    if group == 2, ///
    crosssectional(carbon_emissions_pc gdp_pc fd_index) cr_lags(0) 
estimates store cce_brics


* Table for comparison 
estimates table cce_g7 cce_brics cce_tax cce_other, ///
    b(%9.4f) se(%9.4f) ///
    title("CCEMG by country group")
	
log close