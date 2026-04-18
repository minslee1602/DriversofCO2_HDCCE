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

* Declare the panel using the integer year
xtset country_id year_int

* Estimate fixed effects model
xtreg carbon_emissions_pc ///
    gdp_pc gdpsq_pc fd_index ///
    coal_consumption_pc oil_consumption_pc gas_consumption_pc ///
    renewables_consumption_pc ///
    coal_production_pc oil_production_pc gas_production_pc ///
    solar_production_pc wind_production_pc hydro_production_pc ///
    other_renewables_production_pc, fe

* Re-run diagnostics
predict epsilon2, e
xtcsd, pesaran abs
xtcse2 epsilon2