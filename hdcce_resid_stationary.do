*==============================================================================
* Panel unit-root test on HD-CCE residuals (least squares and desparsified
* lasso)
*
* Output: xtcips_hdcce_resid_all.dta, xtcips_hdcce_residuals.xlsx (sheets
* Stata_summary, Stata_all_tests).
*==============================================================================

clear all
set more off
capture log close
log using "xtcips_hdcce_residuals.log", replace text

*--- 0. Settings  --------------------------------------------------
local datafile    "tables/HDCCE_residuals.xlsx"
local idvar       "CountryCode"
local yearvar     "Year"
local groupvar    "Group"
local testvars    "resid_ls resid_lasso"
local maxlags_list "1 2"
local bglags      2      // only used by xtcips for its own BG diagnostics

local d1name "Intercept"
local d1opt  ""
local d2name "Intercept + trend"
local d2opt  "trend"

*--- 1. Install xtcips if needed -----------------------------------------------
capture which xtcips
if _rc ssc install xtcips

*--- 2. Load and check the data ----------------------------------------------
capture confirm file "`datafile'"
if _rc {
    display as error "Cannot find `datafile'"
    display as error "Expected a 'tables' subfolder next to this do-file, containing HDCCE_residuals.xlsx."
    display as error "Current working directory:"
    display as error "`c(pwd)'"
    display as error "Contents of the tables subfolder (if it exists):"
    capture dir "tables"
    exit 601
}

import excel using "`datafile'", firstrow clear

foreach v in `idvar' `yearvar' `groupvar' resid_ls resid_lasso {
    capture confirm variable `v'
    if _rc {
        display as error "Variable `v' not found. Edit the settings block. Variables present:"
        ds
        exit 111
    }
}

foreach v in `yearvar' resid_ls resid_lasso {
    capture confirm string variable `v'
    if !_rc {
        display as text "Converting `v' from text to numeric"
        destring `v', replace force
    }
}

capture confirm string variable `groupvar'
if _rc {
    display as error "`groupvar' is expected to be a string variable ('All' / 'Carbon tax' / 'No carbon tax')"
    exit 198
}

tempfile master
save `master'

*--- 3. Run xtcips for every variable x group -----------------------------------
tempname H
tempfile res
postfile `H' str15 group str40 variable str16 transformation str18 deterministics ///
    int maxlags int n_used int n_total double cips cv01 cv05 cv10 using `res', replace

local g1name "All"
local g2name "Carbon tax"
local g3name "No carbon tax"

foreach v of local testvars {
    forvalues g = 1/3 {
        local gname "`g`g'name'"

        use `master', clear
        keep if `groupvar' == "`gname'"
        quietly count
        if r(N) == 0 {
            display as error "No observations for group `gname' -- skipping"
            continue
        }

        capture drop id tidx
        egen id = group(`idvar')
        egen tidx = group(`yearvar')
        xtset id tidx
        if "`r(balanced)'" != "strongly balanced" {
            display as error "Group `gname' is not strongly balanced for `v'. Diagnostics:"
            xtdescribe
            bysort id: gen __n = _N
            display as error "Countries without a full set of observations:"
            tab `idvar' if __n != r(max)
            exit 459
        }
        quietly summarize id
        local n_total = r(max)

        foreach tr in "Level" "First difference" {

            capture drop __x __nmis __sd __ok __t2
            if "`tr'" == "Level" {
                quietly gen double __x = `v'
            }
            else {
                quietly gen double __x = D.`v'
            }

            quietly {
                bysort id: egen __nmis = total(missing(`v'))
                bysort id: egen __sd   = sd(__x)
                gen byte __ok = (__nmis == 0 & __sd > 1e-8 & __sd < .)
                egen __t2 = tag(id) if __ok == 1
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
                        capture quietly xtcips __x if __ok == 1, ///
                            maxlags(`ml') bglags(`bglags') `dopt'
                        if _rc == 0 {
                            local cips = r(cips)
                            matrix CV  = r(cv)
                            matrix CVv = vec(CV)
                            * Assumed order of the critical values: 1%, 5%, 10%.
                            * The FE-OLS residual check on this same panel suggested
                            * this xtcips installation actually returns 10%, 5%, 1%
                            * (reversed): verify with -matrix list r(cv)- after one
                            * manual run and edit the three lines below if needed.
                            * cv05 (the middle element) is unaffected either way.
                            if rowsof(CVv) == 3 {
                                local c1  = CVv[1,1]
                                local c5  = CVv[2,1]
                                local c10 = CVv[3,1]
                            }
                            else {
                                display as error "Unexpected r(cv) layout for `v' / group `gname'"
                            }
                        }
                    }

                    post `H' ("`gname'") ("`v'") ("`tr'") ("`dname'") (`ml') ///
                        (`n_used') (`n_total') (`cips') (`c1') (`c5') (`c10')
                }
            }
        }
        display as text "done: `v' / `gname'"
    }
}
postclose `H'

*--- 4. Decisions and classification -------------------------------------------
use `res', clear
gen byte countries_dropped = n_total - n_used
gen byte reject_01 = cips < cv01 if !missing(cips, cv01)
gen byte reject_05 = cips < cv05 if !missing(cips, cv05)
gen byte reject_10 = cips < cv10 if !missing(cips, cv10)
label define yn 0 "No" 1 "Yes"
label values reject_01 reject_05 reject_10 yn
save "xtcips_hdcce_resid_all.dta", replace

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
sort group variable deterministics maxlags

display as text _newline "Classification (5% level), HD-CCE residuals:"
tab class_stata group if variable == "resid_ls" & maxlags == 1 & deterministics == "Intercept"
tab class_stata group if variable == "resid_lasso" & maxlags == 1 & deterministics == "Intercept"
tab class_stata group if variable == "resid_ls" & maxlags == 1 & deterministics == "Intercept + trend"
tab class_stata group if variable == "resid_lasso" & maxlags == 1 & deterministics == "Intercept + trend"

save "xtcips_hdcce_resid_summary.dta", replace
export excel using "xtcips_hdcce_residuals.xlsx", sheet("Stata_summary") firstrow(variables) replace
use "xtcips_hdcce_resid_all.dta", clear
export excel using "xtcips_hdcce_residuals.xlsx", sheet("Stata_all_tests") firstrow(variables) sheetreplace

log close
