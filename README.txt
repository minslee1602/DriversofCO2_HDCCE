Drivers of Carbon Emissions in OECD and BRICS Countries - A High-Dimensional CCE Panel Approach
==================================================================================================

Replication materials for the HD-CCE panel analysis of per capita carbon emissions
across 42 OECD and BRICS countries, 1990-2019.

Final_Cabon_Emm_disaggregate_zscored.xlsx
    Source data: panel of 42 countries, 1990-2019, with the outcome and 22
    drivers (raw and z-scored)

HDCCE_LASSOLS_CODE.R
    Estimates the HD-CCE least-squares and desparsified-lasso models, with
    bandwidth sensitivity, model fit, and residual diagnostics

HDCCE_LASSOLS_RESULTS.xlsx
    Output of the R script above

stationarity_tests.do
    Runs the CIPS panel unit-root test on the outcome and all 22 drivers

xtcips_results.xlsx, xtcips_results_all.dta, xtcips_summary.dta, xtcips_results.log
    Output of the do-file above

Monte Carlo simulations are provided under releases: "Monte Carlo Replication Code" --> "hdcce_Simulations.zip"
