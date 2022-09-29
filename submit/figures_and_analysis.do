********************************************************************************
* Title: Time Series Project @ UT Austin MA in Econ Opt III			           *
* Author: Lutfi Sun                                                            *
* Date: Thursday April 25, 2021                                                *
********************************************************************************

capture log close
clear all
drop _all
macro drop _all 

cd "/Users/lutfisun/Desktop/turkey_taylor_project"

log using tsproject.log, replace

*-------------------------------------------------------------------------------
* Overnight Interest Rates
*-------------------------------------------------------------------------------

import delimited "data/cb_oni.csv"

* data on overnight interest rates is missing for some most months and 
* not available quarterly (could not find) so, i will take averages of 
* observations within the same quarter

gen edate = date(date, "DM20Y")
format edate %d 
gen quarter = quarter(edate)
gen year = year(edate)
gen yq = yq(year, quarter)

egen oni_b = mean(oni_borrowing), by(yq)
egen oni_l = mean(oni_lending ), by(yq)

label variable edate "date"
label variable oni_b "borrowing"
label variable oni_l "lending"

graph twoway line oni_b oni_l edate, name(cb_oni)  ///
	title("CBRT Overnight Interest Rates") xline(15600 17000 18500 20700) ///
	text(58 15602  "implicit targeting" 58 17002 "inflation targeting", placement(east) size(.25cm) color(maroon)) ///
	text(58 18503  "wide corridor system" 58 20703 "not wide corridor", placement(east) size(.25cm) color(maroon))

graph export "figures/cb_oni.pdf", name(cb_oni) replace

drop oni_borrowing oni_lending edate quarter year
duplicates drop yq, force

save "data/cb_oni.dta", replace

*-------------------------------------------------------------------------------
* Inflation Target
*-------------------------------------------------------------------------------

clear
import delimited "data/cb_inflation.csv"

graph twoway line target realization year, name(cb_inflation) xsc(r(2002 2021)noextend)  ///
	title("Year to Year % Changes in CPI") xline(2002 2006 2010 2016) ///
	text(38 2002.1  "implicit targeting" 38 2006.06 "inflation targeting", placement(east) size(.25cm) color(maroon)) ///
	text(38 2010.1  "wide corridor system" 38 2016.1 "not wide corridor", placement(east) size(.25cm) color(maroon))

graph export "figures/cb_inflation.pdf", name(cb_inflation) replace

gen quarter = 1
gen yq = yq(year, quarter)

drop year quarter
duplicates drop yq, force

save "data/cb_inflation.dta", replace

*-------------------------------------------------------------------------------
* Merge All
*-------------------------------------------------------------------------------

clear
import delimited "data/evds_clean.csv", parselocale(en_US)

gen yq = quarterly(date, "YQ")
format yq %tq

merge m:m yq using "data/cb_oni.dta", generate(merg1)
merge m:m yq using "data/cb_inflation.dta", generate(merg2)

drop merg*
duplicates drop yq, force
sort yq
tset yq

ssc install fillmissing, replace
fillmissing oni_b oni_l target, with(previous)

save "data/merged_clean.dta", replace

*-------------------------------------------------------------------------------
* Figures
*-------------------------------------------------------------------------------

* output

gen ln_gdp = log(gdp) 
sort yq
gen d_ln_gdp = D.ln_gdp

gen gdp_usd = gdp / usd_buy
gen ln_gdp_usd = log(gdp_usd )

label variable yq "Quarter"
label variable ln_gdp "log gdp try"
label variable ln_gdp_usd "log gdp usd"

graph twoway line ln_gdp ln_gdp_usd yq, name(gdp) title("Natural Log of Gross Domestic Product")
graph export "figures/gdp.pdf", name(gdp) replace

gen quarter = quarter(dofq(yq))
gen period2 = _n * _n
reg ln_gdp_usd yq i.quarter
predict ln_gdp_usd_hat

predict output_gap, resid
gen output_gap_10 = 10*output_gap  
label variable output_gap_10 "detrended gdp"
label variable oecd_cli_dtd "detrended cli"

graph twoway (line output_gap_10 oecd_cli_dtd yq), name(output_gap) title("Output Gap Using GDP and Composite Leading Indicator")
graph export "figures/output_gap.pdf", name(output_gap) replace

* inflation

gen ln_cpi95 = log(cpi95_all)
gen d_ln_cpi95 = D.ln_cpi95

graph twoway line d_ln_cpi95 yq, name(d_ln_cpi95) title("Changes in Log Cost of Living 95")
graph export "figures/d_ln_cpi95.pdf", name(d_ln_cpi95) replace

gen ln_cpi03 = log(cpi03_all)
gen d_ln_cpi03 = D.ln_cpi03

graph twoway line d_ln_cpi03 yq, name(d_ln_cpi03) title("Changes in Log Cost of Living 03")
graph export "figures/d_ln_cpi03.pdf", name(d_ln_cpi03) replace

* forex

graph twoway line usd_buy yq, name(usd_fx_rates) title("USD to TRY Exchange Rates")
graph export "figures/usd_fx_rates.pdf", name(usd_fx_rates) replace

*-------------------------------------------------------------------------------
* Short Run SVAR Analysis (USD)
*-------------------------------------------------------------------------------

gen OG_01 = output_gap_10
gen dP_01 = D.ln_cpi95
gen oi_01 = oni_b

matrix A1 = (1,0,0 \ .,1,0 \ .,.,1)
matrix B1 = (.,0,0 \ 0,.,0 \ 0,0,.)

svar OG_01 dP_01 oi_01, lags(1/8) aeq(A1) beq(B1)
irf create svar1, set(var_diff1, replace)
irf cgraph (svar1 OG_01 OG_01 sirf) (svar1 OG_01 dP_01 sirf) (svar1 OG_01 oi_01 sirf) ///
           (svar1 dP_01 OG_01 sirf) (svar1 dP_01 dP_01 sirf) (svar1 dP_01 oi_01 sirf) ///
		   (svar1 oi_01 OG_01 sirf) (svar1 oi_01 dP_01 sirf) (svar1 oi_01 oi_01 sirf), ///
		   title ("Irfs of VAR Model: dtd_gdp_usd d_ln_cpi95 oni", size(small))
		   
graph export "figures/svar1.pdf", replace

*
gen Y_02 = ln_gdp_usd
gen P_03 = ln_cpi03
gen i_02 = oni_b

svar Y_02 P_03 i_02, lags(1/8) aeq(A1) beq(B1)
irf create svar2, set(var_ln2, replace)
irf cgraph (svar2 Y_02 Y_02 sirf) (svar2 Y_02 P_03 sirf) (svar2 Y_02 i_02 sirf) ///
           (svar2 P_03 Y_02 sirf) (svar2 P_03 P_03 sirf) (svar2 P_03 i_02 sirf) ///
		   (svar2 i_02 Y_02 sirf) (svar2 i_02 P_03 sirf) (svar2 i_02 i_02 sirf), ///
		   title ("Irfs of VAR Model: ln_gdp ln_cpi95 oni", size(small))

graph export "figures/svar2.pdf", replace

save "data/merged_clean_post.dta", replace


*-------------------------------------------------------------------------------
* Short Run SVAR Analysis (TRY and PPI)
*-------------------------------------------------------------------------------

gen ln_gdp_try = log(gdp)
reg ln_gdp_try yq i.quarter
predict og_try, resid

gen ln_ppi03 = log(ppi_03)

gen OG_02 = og_try
gen dP_02 = D.ln_ppi03
gen ex_02 = D.usd_buy
gen oi_02 = oni_l

matrix A2 = (1,0,0,0 \ .,1,0,0 \ .,.,1,0 \ .,.,.,1)
matrix B2 = (.,0,0,0 \ 0,.,0,0 \ 0,0,.,0 \ 0,0,0,.)

svar OG_02 dP_02 ex_02 oi_02, lags(1/8) aeq(A2) beq(B2)
irf create svar4, set(var_diff2, replace)
irf cgraph (svar4 OG_02 OG_02 sirf) (svar4 OG_02 dP_02 sirf) (svar4 OG_02 ex_02 sirf) (svar4 OG_02 oi_02 sirf) ///
           (svar4 dP_02 OG_02 sirf) (svar4 dP_02 dP_02 sirf) (svar4 dP_02 ex_02 sirf) (svar4 dP_02 oi_02 sirf) ///
		   (svar4 ex_02 OG_02 sirf) (svar4 ex_02 dP_02 sirf) (svar4 ex_02 ex_02 sirf) (svar4 ex_02 oi_02 sirf) ///
		   (svar4 oi_02 OG_02 sirf) (svar4 oi_02 dP_02 sirf) (svar4 oi_02 ex_02 sirf) (svar4 oi_02 oi_02 sirf), ///
		   title ("Irfs of VAR Model: dtd_gdptry dln_ppi03 d_usdx oni_lending", size(small))
		   
graph export "figures/svar4_try.pdf", replace

*
gen ln_usd = log(usd_buy)

gen OG_03 = ln_gdp_try
gen dP_03 = ln_ppi03
gen ex_03 = ln_usd
gen oi_03 = oni_l

svar OG_03 dP_03 ex_03 oi_03, lags(1/8) aeq(A2) beq(B2)
irf create svar4, set(var_diff2, replace)
irf cgraph (svar4 OG_03 OG_03 sirf) (svar4 OG_03 dP_03 sirf) (svar4 OG_03 ex_03 sirf) (svar4 OG_03 oi_03 sirf) ///
           (svar4 dP_03 OG_03 sirf) (svar4 dP_03 dP_03 sirf) (svar4 dP_03 ex_03 sirf) (svar4 dP_03 oi_03 sirf) ///
		   (svar4 ex_03 OG_03 sirf) (svar4 ex_03 dP_03 sirf) (svar4 ex_03 ex_03 sirf) (svar4 ex_03 oi_03 sirf) ///
		   (svar4 oi_03 OG_03 sirf) (svar4 oi_03 dP_03 sirf) (svar4 oi_03 ex_03 sirf) (svar4 oi_03 oi_03 sirf), ///
		   title ("Irfs of VAR Model: lngdptry lnppi03 lnusdx oni_lending", size(small))
		   

graph export "figures/svar4_try2.pdf", replace

save "data/merged_clean_post.dta", replace

/*
*-------------------------------------------------------------------------------
* Long Run SVAR Analysis 1
*-------------------------------------------------------------------------------

*Assume 1) in the long run, demand shocks (inflation shocks) and policy shocks have no effect on the output. 
*And the policy shocks have no effect on the inflation level (inflation rate is the real variable).

matrix C = (.,0,0 \ .,.,0 \ .,.,.) //no 1's on the diagonal

svar d_Y_01 d_P_01 i_01, lags(1/12) lreq(C) 
irf create lr, set(lrirf1) step(40) replace
*view IRF (not cumulative )
irf graph sirf, yline(0,lcolor(black)) xlabel(0(4)40) byopts(yrescale) set(lrirf1)
graph export "figures/lvar_irf.pdf", replace

*create cumulative irf graphs for responses of output variable 
use lrirf1.irf, clear
sort irfname impulse response step
gen csirf = sirf
by irfname impulse: replace csirf = sum(sirf) if response== "d_Y_01"
order irfname impulse response step sirf csirf
save lrirf3.irf, replace

irf set lrirf3.irf
irf graph csirf, yline(0,lcolor(black)) noci xlabel(0(4)40) byopts(yrescale)

*create cumulative irf graphs for responses of all variables
by irfname impulse: replace csirf = sum(sirf) if response != "d_Y_01"
save lrirf3.irf, replace

irf set lrirf3.irf
irf graph csirf, yline(0,lcolor(black)) noci xlabel(0(4)40) byopts(yrescale)
graph export "figures/lvar_cumulative.pdf", replace
*/
*-------------------------------------------------------------------------------
* Long Run SVAR Analysis 2
*-------------------------------------------------------------------------------

matrix C = (.,0,0,0 \ .,.,0,0 \ .,.,.,0 \ .,.,.,.) //no 1's on the diagonal

svar OG_02 dP_03 ex_03 oi_03, lags(1/12) lreq(C) 
irf create lr, set(lrirf1) step(40) replace
*view IRF (not cumulative )
irf graph sirf, yline(0,lcolor(black)) xlabel(0(4)40) byopts(yrescale) set(lrirf1)
graph export "figures/lvar_irf3.pdf", replace

*create cumulative irf graphs for responses of output variable 
use lrirf1.irf, clear
sort irfname impulse response step
gen csirf = sirf
by irfname impulse: replace csirf = sum(sirf) if response== "OG_02"
order irfname impulse response step sirf csirf
save lrirf3.irf, replace

irf set lrirf3.irf
irf graph csirf, yline(0,lcolor(black)) noci xlabel(0(4)40) byopts(yrescale)

*create cumulative irf graphs for responses of all variables
by irfname impulse: replace csirf = sum(sirf) if response != "OG_02"
save lrirf3.irf, replace

irf set lrirf3.irf
irf graph csirf, yline(0,lcolor(black)) noci xlabel(0(4)40) byopts(yrescale)
graph export "figures/lvar_cumulative3.pdf", replace

*-------------------------------------------------------------------------------
use "data/merged_clean_post.dta", clear


/* a. Test for a unit root in r1 and r3 series using several augmented Dickey Fuller tests. 
gen t = _n
tset t
graph twoway line r1 r3 t

dfuller r1
dfuller r1, trend
dfuller r1, drift
*/


log close
