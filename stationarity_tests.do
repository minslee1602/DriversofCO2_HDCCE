*==============================================================================

* Panel root tests to verify nonstationarity
*
* Output: xtcips_results_all.dta, xtcips_results.xlsx (sheets Stata_summary,
* Stata_all_tests, and Compare_to_R if the R workbook is found).
*==============================================================================

clear all
set more off
capture log close
log using "xtcips_results.log", replace text

*--- 0. Settings  --------------------------------------------------
local datafile    "Final_Cabon_Emm_disaggregate_zscored.xlsx"
local idvar       "CountryCode"
local yearvar     "Year"
local taxvar      "CarbonTax"
local outcome     "carbon_emissions_pc"
local exclude     "carbon_emissions_pc_z"
local nyears      30     
local dropcountry "CRI"
local maxlags_list "1 2"
local bglags      2      // only used by xtcips for its own BG diagnostics

* Groups and deterministic specifications 
local g1name "All"
local g1cond "1==1"
local g2name "Carbon tax"
local g2cond "`taxvar'==1"
local g3name "No carbon tax"
local g3cond "`taxvar'==0"

local d1name "Intercept"
local d1opt  ""
local d2name "Intercept + trend"
local d2opt  "trend"

*--- 1. Install xtcips if needed -----------------------------------------------
capture which xtcips
if _rc ssc install xtcips

*--- 2. Load and prepare the data ----------------------------------------------
import excel using "`datafile'", firstrow clear


foreach v in `idvar' `yearvar' `taxvar' `outcome' {
    capture confirm variable `v'
    if _rc {
        display as error "Variable `v' not found. Edit the settings block. Variables present:"
        ds
        exit 111
    }
}

foreach v in `yearvar' `taxvar' `outcome' {
    capture confirm string variable `v'
    if !_rc {
        display as text "Converting `v' from text to numeric"
        destring `v', replace force
    }
}


quietly levelsof `yearvar', local(yrs)
local nyr : word count `yrs'
local firstyear : word `=`nyr' - `nyears' + 1' of `yrs'
local lastyear  : word `nyr' of `yrs'
display as text "Sample window: `firstyear'-`lastyear'"

keep if inrange(`yearvar', `firstyear', `lastyear')
drop if `idvar' == "`dropcountry'"
quietly count if !inlist(`taxvar', 0, 1)
if r(N) > 0 {
    display as error "`taxvar' has values other than 0/1 (or missing) after the sample restriction:"
    tab `taxvar', missing
    exit 198
}
egen id = group(`idvar')
* Use a plain 1..T time index for xtset. Year was imported with a date format
* (%td), and a daily-date time variable would make D. (first differences) look
* for the previous DAY instead of the previous year.
egen tidx = group(`yearvar')
xtset id tidx
if "`r(balanced)'" != "strongly balanced" {
    display as error "Panel is not strongly balanced. Diagnostics:"
    xtdescribe
    bysort id: gen __n = _N
    display as error "Countries without exactly `nyears' observations in the window:"
    tab `idvar' if __n != `nyears'
    exit 459
}
quietly summarize id
display as text "Countries in panel: " r(max)

ds *_z
local zvars `r(varlist)'
foreach v in `zvars' {
    capture confirm string variable `v'
    if !_rc {
        display as text "Converting `v' from text to numeric"
        destring `v', replace force
    }
}
local drivers : list zvars - exclude
local testvars `outcome' `drivers'
local ndrv : word count `drivers'
assert `ndrv' == 22

*--- 3. Run xtcips for every combination and collect the results ----------------
tempname H
tempfile res
postfile `H' str15 group str40 variable str16 transformation str18 deterministics ///
    int maxlags int n_used int n_total double cips cv01 cv05 cv10 using `res', replace

foreach v of local testvars {
    foreach tr in "Level" "First difference" {

        capture drop __x __nmis __sd __ok
        if "`tr'" == "Level" {
            quietly gen double __x = `v'
            local expr "`v'"
        }
        else {
            quietly gen double __x = D.`v'
            local expr "D.`v'"
        }

        * Country-level eligibility: no missing values and a non-constant series
        quietly {
            bysort id: egen __nmis = total(missing(`v'))
            bysort id: egen __sd   = sd(__x)
            gen byte __ok = (__nmis == 0 & __sd > 1e-8 & __sd < .)
            sort id tidx
        }

        forvalues g = 1/3 {
            local gname "`g`g'name'"
            local gcond "`g`g'cond'"

            quietly {
                capture drop __t1 __t2
                egen __t1 = tag(id) if `gcond'
                egen __t2 = tag(id) if `gcond' & __ok == 1
                count if __t1 == 1
                local n_total = r(N)
                count if __t2 == 1
                local n_used = r(N)
            }

            forvalues d = 1/2 {
                local dname "`d`d'name'"
                local dopt  "`d`d'opt'"

                foreach ml of local maxlags_list {

                    local cips = .
                    local c1   = .
                    local c5   = .
                    local c10  = .

                    if `n_used' >= 5 {
                        capture quietly xtcips `expr' if `gcond' & __ok == 1, ///
                            maxlags(`ml') bglags(`bglags') `dopt'
                        if _rc == 0 {
                            local cips = r(cips)
                            matrix CV  = r(cv)
                            matrix CVv = vec(CV)
                            * Assumed order of the critical values: 1%, 5%, 10%.
                            * Check once with -matrix list r(cv)- after a manual
                            * xtcips run and edit here if the layout differs.
                            if rowsof(CVv) == 3 {
                                local c1  = CVv[1,1]
                                local c5  = CVv[2,1]
                                local c10 = CVv[3,1]
                            }
                            else {
                                display as error "Unexpected r(cv) layout for `v' / `gname'"
                            }
                        }
                    }

                    post `H' ("`gname'") ("`v'") ("`tr'") ("`dname'") (`ml') ///
                        (`n_used') (`n_total') (`cips') (`c1') (`c5') (`c10')
                }
            }
        }
    }
    display as text "done: `v'"
}
capture drop __x __nmis __sd __ok __t1 __t2
postclose `H'

*--- 4. Decisions and classification -------------------------------------------
use `res', clear
gen byte countries_dropped = n_total - n_used
gen byte reject_01 = cips < cv01 if !missing(cips, cv01)
gen byte reject_05 = cips < cv05 if !missing(cips, cv05)
gen byte reject_10 = cips < cv10 if !missing(cips, cv10)
label define yn 0 "No" 1 "Yes"
label values reject_01 reject_05 reject_10 yn
save "xtcips_results_all.dta", replace

* One row per group x variable x deterministics x maxlags (levels vs differences)
preserve
keep if transformation == "Level"
keep group variable deterministics maxlags n_used n_total cips cv05 reject_05
rename (cips cv05 reject_05) (cips_level cv05_level rej05_level)
tempfile lvl
save `lvl'
restore

keep if transformation == "First difference"
keep group variable deterministics maxlags cips cv05 reject_05
rename (cips cv05 reject_05) (cips_diff cv05_diff rej05_diff)
merge 1:1 group variable deterministics maxlags using `lvl', nogenerate

gen str60 class_stata = "Not available"
replace class_stata = "Stationary in levels (I(0)-like)" if rej05_level == 1
replace class_stata = "Unit root in levels, stationary in differences (I(1)-like)" ///
    if rej05_level == 0 & rej05_diff == 1
replace class_stata = "Unit root in levels and differences (inconclusive)" ///
    if rej05_level == 0 & rej05_diff == 0
order group variable deterministics maxlags n_used n_total cips_level cv05_level ///
    cips_diff cv05_diff class_stata
sort group deterministics maxlags variable

display as text _newline "Classification counts (5% level):"
tab class_stata group if maxlags == 1 & deterministics == "Intercept"
tab class_stata group if maxlags == 1 & deterministics == "Intercept + trend"

save "xtcips_summary.dta", replace
export excel using "xtcips_results.xlsx", sheet("Stata_summary") firstrow(variables) replace
use "xtcips_results_all.dta", clear
export excel using "xtcips_results.xlsx", sheet("Stata_all_tests") firstrow(variables) sheetreplace


log close
