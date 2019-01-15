*! Hello, I'm med4way.ado
*! v2.3.0 - XXxxxXXXX

/* 
Previous versions:
v2.2.3 - 25nov2018
v2.2.2 - 14jun2018
v2.2.1 - 19sep2017
v2.2.0 - 31jul2017
v2.1.1 - 28jul2017
v2.1.0 - 26jul2017
v2.0.0 - 27mar2017
*/

capture program drop med4way 
program define med4way, eclass
	version 11.0

	if replay() {
		if ("`e(cmd)'" != "med4way") error 301
		
		if "`e(vce)'" == "bootstrap" {
			local vcetype "boostrap"
		}
		else {
			local vcetype "delta method"
		}
		display _n(1) "-> 4-way decomposition: `vcetype'" _n(1)
		Display `0'
		exit
	}

	syntax varlist(min=2 numeric) [if] [in], /*
	*/	/*
	*/ a0(real) /*
	*/ a1(real) /*
	*/ m(real) /*
	*/ yreg(string) /*
	*/ mreg(string) /*
	*/	/*
	*/ [ c(string) /*
	*/ CASEControl /*
	*/ FULLoutput /*
	*/ NODELTAMethod /*
	*/ ROBUST /* undocumented - superseeded by yregoptions()
	*/ NOLEGEND /* undocumented
	*/ NOWARNing /* undocumented
	*/ level(cilevel) /*
	*/ BOOTstrap /*
	*/ reps(integer 1000) /*
	*/ seed(passthru) /*
	*/ SAving(passthru) /*
	*/ YREGOPTions(passthru) /*
	*/ MREGOPTions(passthru) /*
	*/ ESTSTOre /* undocumented
	*/ BCA  * ]

	//[if] [in] marksample
	marksample touse
	
	// Step 1===================================================================
	// Checks and parsing of program options
	local wrnngtxt 0 
		// initialize wrnngtxt local. The idea behind wrnngtxt is
		// that, in case of Error in Step 1, no Warnng msgs are issued. 
		// Plus, it is possible to customize the order of the warnngs

	//parse dioptions
	_get_diopts diopts, `options'

	//parse yreg
	gettoken yreg_1 yreg_2 : yreg, parse(",")
	gettoken comma opt : yreg_2,   parse(",")
// 	if length("`comma'") > 1 { unnecessary
//  		local comma = trim("`comma'")
// 	}
	
	local yreg = trim(`"`yreg_1'"')
    local dist = trim(`"`opt'"')

	
	local l = length("`yreg'")
	if substr("linear", 1, max(3,`l')) == "`yreg'" {
		local yreg "linear"
	}
	if substr("logistic", 1, max(4,`l')) == "`yreg'" {
		local yreg "logistic"
	}
	if substr("logbinomial", 1, max(4,`l')) == "`yreg'" {
		local yreg "logbinomial"
	}
	if substr("poisson", 1, max(3,`l')) == "`yreg'" {
		local yreg "poisson"
	}
	if substr("negbinomial", 1, max(4,`l')) == "`yreg'" {
		local yreg "negbinomial"
	}
	
	//parse mreg
	local l = length("`mreg'")
	if substr("linear", 1, max(3,`l')) == "`mreg'" {
		local mreg "linear"
	}
	if substr("logistic", 1, max(4,`l')) == "`mreg'" {
		local mreg "logistic"
	}
	
	//validate yreg and mreg 
	local yregtypes linear logistic logbinomial poisson negbinomial cox aft
	local nyreg : list posof "`yreg'" in yregtypes
	if !`nyreg' {
		display as error "Error: yreg must be chosen from: `yregtypes'."
		error 198
	} 
	else {
		local yreg : word `nyreg' of `yregtypes'
	}
	
	local mregtypes linear logistic
	local nmreg : list posof "`mreg'" in mregtypes
	if !`nmreg' {
		display as error "Error: mreg must be chosen from: `mregtypes'."
		error 198		
	}
	else {
		local mreg : word `nmreg' of `mregtypes'
	}
	
	//parse dist
	if ("`yreg'"=="aft") {
		local l = length("`dist'")
		if `l' == 0 local dist "exponential"
		
		if substr("exponential", 1, max(1,`l')) == "`dist'" {
			local dist "exponential"
		}
		if substr("weibull", 1, max(1,`l')) == "`dist'" {
			local dist "weibull"
		}	
	}
	else {
		local dist ""
	}
	
	//validate dist
	if ("`yreg'"=="aft") {
		local disttypes exponential weibull
			local ndist : list posof "`dist'" in disttypes
			if !`ndist' {
				display as error "Error: distribution for aft models must be chosen from: `disttypes'."
				error 198
			} 	
	}
	
	
	//bootsrap
	if ("`bootstrap'" !="") {
		local bootstrap true
	}
	else {
		local bootstrap false
	}
	
	//suppress delta method
	if ("`nodeltamethod'" !="") {
		local deltamethod false
	}
	else {
		local deltamethod true
	} 

	//full output or reduced output
	if ("`fulloutput'" != "") {
		local output full
	}
	else {
		local output reduced
	}

	//casecontrol or not
	if ("`casecontrol'" !="") {
		local casecontrol true
	}
	else {
		local casecontrol false
	}
	
	//legend or not
	if ("`nolegend'" !="") {
		local legend false
	}
	else {
		local legend true
	}
	if ("`deltamethod'"=="false") & ("`bootstrap'"=="false") { 
		// if deltamethod not requested and bootstrap not requested, don't print the legend (there is no legend to print!)
		local legend false 
	}

	//warnings or not
	if ("`nowarning'" !="") {
		local warning false
	}
	else {
		local warning true
	}
	
	//store estimates
	if ("`eststore'" != "") {
		local eststore true
	}
	else {
		local eststore false
	}
	
	//survival outcome
	if ("`yreg'"=="cox") | ("`yreg'"=="aft") {
		local survoutcome true
	}
	else {
		local survoutcome false
	}
	
	
	//check stset when appropriate
	if ("`survoutcome'"=="true") {
		st_is 2 analysis
	}
	
	//tokenize main variables	
	if ("`survoutcome'"=="true") {
		gettoken avar varlist	: varlist
		gettoken mvar cvars		: varlist
	}
	else if ("`survoutcome'"=="false") {	
		gettoken yvar varlist	: varlist
		gettoken avar varlist	: varlist
		gettoken mvar cvars		: varlist
	}
	
	//get rid of possible duplicates in cvars
	local cvars : list uniq cvars
	
	
	//mvar takes on 0/1 values if mreg is logistic
	if ("`mreg'"=="logistic") {
		if ("`yreg'"=="logistic") & ("`casecontrol'"=="true") { 
			qui levelsof `mvar' if `yvar' == 0 & `touse', local(lom)
		}
		else {
			qui levelsof `mvar' if `touse', local(lom)
		}
		local lomallowed 0 1
		local lomx : list lom == lomallowed
		if `lomx' == 0 {
			display as error "Error: `mvar' must be coded as 0/1 (`mvar' is coded " /*
				*/ "as `lom' in the data)".
			error 198
		}
	}
	
	//avar takes on 0/1 values if avar is probably binary
	qui levelsof `avar' if `touse', local(loa)
	local nloa : word count `loa'
	if `nloa' == 2 {
		local loaallowed 0 1
		local loax : list loa == loaallowed
		if `loax' == 0 {
			local wrnngtxt `wrnngtxt' 7
		}
	}
	
	//casecontrol works only with yreg logistic
	if ("`casecontrol'"=="true") & !("`yreg'"=="logistic") {
		display as error "Error: option casecontrol can only be specified together with a " /*
			*/ "logistic regression model for the outcome: yreg(logistic)."
			error 198
	}
	
	//validate rare outcome
	if ("`yreg'"=="logistic") {
		if ("`casecontrol'"=="true") {
			local wrnngtxt `wrnngtxt' 2
		}
		else {
			qui logit_10 `yvar' if `touse'
			local prev = invlogit(_b[_cons])*100
			if `prev' > 10 {
				local wrnngtxt `wrnngtxt' 1
			}
		}
	}
	if ("`survoutcome'"=="true") {
		capture assert _d == 1 if _st == 1 & `touse'
		if _rc == 0 {
			local wrnngtxt `wrnngtxt' 3
		}
		if ("`yreg'"=="cox") {
			local wrnngtxt `wrnngtxt' 2
		}
	}

	//validate c and cvars
	validate_c if `touse', c(`c') cvars(`cvars') wrnngtxt(`wrnngtxt')
	tempname cmatrix
	mat `cmatrix' = r(cmatrix)
	
	//interaction
	local avartrunc = substr("`avar'", 1, 12)
	local mvartrunc = substr("`mvar'", 1, 12)
	
	local inter_var_names "_`avartrunc'X`mvartrunc'_000 _`mvartrunc'X`avartrunc'_111 _`avartruncX`mvartrunc'_001"
	local inter_var_names "`inter_var_names' _`avartrunc'X`mvartrunc'_010 _`avartrunc'X`mvartrunc'_100"
	local inter_var_names "`inter_var_names' _`mvartrunc'X`avartrunc'_001 _`mvartrunc'X`avartrunc'_010 _`mvartrunc'X`avartrunc'_100"

	foreach name of local inter_var_names {
		capture confirm new variable `name'
		if !_rc {
			local inter `name'
			continue, break
		}
	}
	//all the 8 suggested names have been used by existing variables, throw an error message
	if _rc {
		display as error "Error: The command needs to create an interaction variable " /*
			*/ "with one of the following names: `inter_var_names', " /*
			*/ "but these variables have already been defined."
		error 110
	}
	qui gen `inter' = `avar'*`mvar'
	
	// Compress a0 a1 m into one single vector
	tempname aam 
	mat `aam' = (`a0', `a1', `m')
	
	// Set up estimates' names
	if ("`yreg'"=="linear") {
		local namesr "te cde intref intmed pie"
		local namesf "p_cde p_intref p_intmed p_pie op_m op_ati op_e"
	}
	else if ("`yreg'"=="logistic") | ("`yreg'"=="cox") | ("`yreg'"=="aft") | /*
		*/ ("`yreg'"=="logbinomial") | ("`yreg'"=="poisson") | ("`yreg'"=="negbinomial") {
		local namesr "tereri ereri_cde ereri_intref ereri_intmed ereri_pie"
		local namesf "terira p_cde p_intref p_intmed p_pie op_m op_ati op_e"
	}
	
	if "`output'" == "reduced" {
		local names `namesr'
	}
	else if "`output'" == "full" {
		local names `namesr' `namesf'
	}
		
	local nnames : word count `names'

	// Count observations if touse
	qui count if `touse'
	local nobs = r(N)
	//==========================================================================	
	
	// Step 2===================================================================
	//print warnng text and summary of med4way
	if "`warning'"=="true" {
		local wrnngtxt : list sort wrnngtxt
		foreach w of local wrnngtxt {
			if `w' == 1 {
				display as error "Warning: this analysis assumes a rare outcome. " /*
					*/ "The outcome variable `yvar' has " /*
					*/  %2.0f `prev' `"% of "positive" outcomes."'
			}
			if `w' == 2 {
				display as error "Warning: this analysis assumes a rare outcome."
			}
			if `w' == 3 {
				display as error "Warning: no censored event times. Please, check that " /*
					*/ "the data was stset correctly."
			}
			if `w' == 4 {
				display as error "Warning: no covariates specified."
			}
			if `w' == 5 {
				display as error "Warning: fixed values for the covariates `cvar' " /*
					*/ "were not provided. All covariates are fixed at their respective mean."
			}
			if `w' == 6 {
				display as error "Warning: fixed values for the covariates `cvar' " /*
					*/ "were provided only for some variables. Covariates are fixed at " /*
					*/ "the provided values or at their respective mean."
			}
			if `w' == 7 {
				display as error "Warning: it looks like the treatment variable `avar' " /*
					*/ "is binary. Treatment variable must be coded as 0/1 (`avar' " /*
					*/ "is coded as `loa' in the data)."
			}
		}
	}
	
	display _n(1) as text "-> Summary" _n(1)
	display _col(4) as text "Outcome    (yvar):"  _col(24) as res "`=cond("`survoutcome'"=="true", "_t", "`yvar'")'"
	display _col(4) as text "Exposure   (avar):"  _col(24) as res "`avar'"
	display _col(4) as text "Mediator   (mvar):"  _col(24) as res "`mvar'"
	display _col(4) as text "Covariates (cvars):" _col(24) as res "`=cond(`nc'==0, "[ none ]", "`cvar'")'" _n(1)
	
	display _col(4) as text "Model for the outcome  (yreg):" _col(35) as res "`yreg'`=cond("`dist'"!="", ", `dist'", "")'"
	display _col(4) as text "Model for the mediator (mreg):" _col(35) as res "`mreg'" _n(1)
	
	display _col(4) as text "Referent exposure level (a0):" 			_col(46) as res "`a0'"
	display _col(4) as text "Actual exposure level   (a1):" 			_col(46) as res "`a1'"
	display _col(4) as text "Mediator level for the decomposition (m):" _col(46) as res "`m'"
	if `nc' != 0 {
		display _col(4) as text "Fixed values of the covariates (c):"   _col(46) as res "`cdisp'"
	}
	if "`bootstrap'"=="true" {
		display _n _col(4) as text "Bootstrap replications (reps):" 	_col(35) as res "`reps'"
	}
	//==========================================================================	

	// Step 3===================================================================	
	// Delta method (or only display yreg and mreg if nodeltamethod = true)
	med4way_engine if `touse', yvar(`yvar') avar(`avar') mvar(`mvar') /*
		*/ cvar(`cvar')  c(`cmatrix') nc(`nc') inter(`inter') /*
		*/ aam(`aam') yreg(`yreg') mreg(`mreg') dist(`dist') /*
		*/ casecontrol(`casecontrol') output(`output') /*
		*/ deltamethod(`deltamethod') bootstrap(false) `robust' /*
		*/ names(`names') nn(`nnames') level(`level') eststore(`eststore') /*
		*/ `yregoptions' `mregoptions'

	if "`deltamethod'" == "true" {
		display _n(2) as text "-> 4-way decomposition: delta method" _n(1)
		ereturn display, level(`level') `diopts'
	}
	//==========================================================================	

	// Step 4===================================================================	
	// Bootstrap (if desired)
	if "`bootstrap'" == "true" {
		if "`deltamethod'" == "true" { // if deltamethod, store the modelbased V for future ereturn
			tempname V_modelbased
			mat `V_modelbased' = e(V)
		} 
		
		display _n(2) as text "-> 4-way decomposition: bootstrap" _n(1)

		bootstrap _b, reps(`reps') level(`level') `seed' `saving' nodots /*
			*/ level(`level') `bca' noh nol `diopts': /*
			*/	med4way_engine if `touse', yvar(`yvar') avar(`avar') mvar(`mvar') /*
			*/	cvar(`cvar')  c(`cmatrix') nc(`nc') inter(`inter') /*
			*/  aam(`aam') yreg(`yreg') mreg(`mreg') dist(`dist') /*
			*/  casecontrol(`casecontrol') output(`output') /*
			*/  deltamethod(false) bootstrap(true) `robust' /*
			*/  names(`names') nn(`nnames') level(`level') eststore(false) /*
			*/  `yregoptions' `mregoptions'
	
		tempname b_bs repsm bias z0 se ci_normal ci_percentile ci_bc 
		matrix `b_bs' = e(b_bs)
		matrix `repsm' = e(reps)
		matrix `bias' = e(bias)
		matrix `z0' = e(z0)
		matrix `se' = e(se)
		matrix `ci_normal' = e(ci_normal)
		matrix `ci_percentile' = e(ci_percentile)
		matrix `ci_bc' = e(ci_bc)
		if "`bca'" == "bca" {
			tempname ci_bca accel
			matrix `ci_bca' = e(ci_bca)
			matrix `accel' = e(accel)
		}
		
		*estat bootstrap, noheader `bca'
	}
	//==========================================================================	

	// Step 5===================================================================
	// interaction is no longer needed
	drop `inter'
	//==========================================================================
	
	// Step 6===================================================================
	// print legend if requested
	if "`legend'"=="true" {
		local allnames "te cde intref intmed pie p_cde p_intref p_intmed p_pie op_m op_ati op_e"
		local allnames "`allnames' tereri ereri_cde ereri_intref ereri_intmed ereri_pie terira"
		
		local fullnames `""total effect" "controlled direct effect" "reference interaction" "mediated interaction" "pure indirect effect" "proportion controlled direct effect" "proportion reference interaction" "proportion mediated interaction" "proportion pure indirect effect" "overall proportion mediated" "overall proportion attributable to interaction" "overall proportion eliminated" "total excess relative risk" "excess relative risk due to controlled direct effect" "excess relative risk due to reference interaction" "excess relative risk due to mediated interaction" "excess relative risk due to pure indirect effect" "total effect risk ratio""'
		
		forval i = 1/`nnames' {
			local w : word `i' of `names'
			local p : list posof "`w'" in allnames
			local f : word `p' of `fullnames'
			local printlegend "`printlegend'" "`w'=`f'`=cond(`i'==`nnames', ".", "; ")'"
		}
		
		display
		display "{p 3 3 5 0}"
		display "`printlegend'"
		display "{p_end}"	
	}
	
	// Step 7===================================================================
	// ereturn stuff (if bootstrap, ereturns bootstrap)
	tempname b V
	matrix `b' = e(b)
	matrix `V' = e(V)
	ereturn post `b' `V', esample(`touse') obs(`nobs')
	
	ereturn local estimands = "`names'"
	
	ereturn local dist =  "`dist'"
	ereturn local yreg =  "`yreg'"
	ereturn local mreg =  "`mreg'"
	
	ereturn local yvar =  "`yvar'"
	ereturn local avar =  "`avar'"
	ereturn local mvar =  "`mvar'"
	ereturn local cvars = "`cvars'"
	
	ereturn scalar level = `level'
	ereturn scalar a0 = `a0'
	ereturn scalar a1 = `a1'
	ereturn scalar m = `m'

	if "`bootstrap'"=="true" { // bootstrap stuff
		ereturn scalar N_reps = `reps'
		ereturn local prefix = "bootstrap"
		ereturn local vcetype = "Bootstrap"
		ereturn local vce = "bootstrap"

		if "`bca'" == "bca" {
			ereturn matrix ci_bca = `ci_bca'
		}
		if "`deltamethod'" == "true" {
			ereturn matrix V_modelbased = `V_modelbased'
		}
		ereturn matrix ci_bc = `ci_bc'
		ereturn matrix ci_percentile = `ci_percentile'
		ereturn matrix ci_normal = `ci_normal'
		ereturn matrix se = `se'
		ereturn matrix z0 = `z0'
		ereturn matrix bias = `bias'
		ereturn matrix reps = `repsm'
		if "`bca'" == "bca" {
			ereturn matrix accel = `accel'	
		}
		ereturn matrix b_bs = `b_bs'	
	}
	
	if `nc' > 0	{
		ereturn matrix c = `cmatrix'
	}

	ereturn local med4way_version "2.3.0"
	ereturn local cmd "med4way"
	ereturn local cmdline "med4way `0'"
	//==========================================================================	
end med4way



/*********************
*
* 	 SUBROUTINES
*
**********************/

/*********************
* validate_c
**********************/
capture program drop validate_c
program define validate_c, rclass
	version 10.0
	syntax [if], [c(string) cvars(string)] wrnngtxt(numlist)
	//if c is missing, take the mean for all the variables in cvars
	//if c is not missing, is the number of elements in c = to the number of elements in cvars? If no, issue error.
	//if c is not missing and n = nc, replace the . with the variable's mean, if needed
	
	//[if] [in] marksample
	marksample touse
	
	local nc : word count `cvars'

	if "`cvars'"=="" {
		local wrnngtxt `wrnngtxt' 4
	}
	
	local n 0
	if "`c'"=="" {
		if "`cvars'"!="" {
			local wrnngtxt `wrnngtxt' 5
			
			foreach i of varlist `cvars' {
				su `i' if `touse', meanonly
				local ctemp `ctemp' `r(mean)'
				local cx: display %-8.4g `r(mean)'
				local cdisptemp `cdisptemp' `cx'
			}
			local c `ctemp' 			// actual c values to be used
			c_local cdisp `cdisptemp' 	// c values for display purposes only
		}
	}
	else {
		local n: word count `c'
		
		if "`cvars'"=="" {
// 			display as error "Warning: c values are ignored when no " /*
// 				*/ "covariates are included via cvars(varlist)."
			local c // return empty c
		}
		else {
			if `n'!=`nc' {
				display as error "Error: the number of c values (`n') " /*
					*/ "does not match the number of covariates (`nc')."
				error 198
			}
			
			else if `n'==`nc' {
				c_local cdisp `c'

				local dot : list posof "." in c // is there a . in c?
				if `dot' > 0 { // yes, there is. the following is needed only if there is a . in c. Otherwise, keep c as it is.
					local wrnngtxt `wrnngtxt' 6
					
					foreach i of numlist 1/`n' {
						local ci : word `i' of `c'
						local cvari : word `i' of `cvars'

						if "`ci'" == "." { // replace the . with the mean
							su `cvari' if `touse', meanonly
							local ctemp `ctemp' `r(mean)'
							local cx: display %-8.4g `r(mean)'
							local cdisptemp `cdisptemp' `cx'
						}
						else if "`ci'" != "." { // leave the value as provided by the user
							local ctemp `ctemp' `ci'
							//local cx: display %-8.4g `ci'
							local cdisptemp `cdisptemp' `ci'
						}
					}
					local c `ctemp' 			// actual c values to be used
					c_local cdisp `cdisptemp' 	// c values for display purposes only (rounded)
				}
			}
		}
	}

	c_local cvar `cvars'	// rename cvars to cvar to be consistent with yvar avar mvar
	c_local nc `nc' 		// total number of covariates
// 	c_local eretc "`c'"		// c used in med4way's ereturn
	c_local wrnngtxt `wrnngtxt' 
	
	mata: st_local("c_missing", strofreal(hasmissing(strtoreal(tokens(st_local("c")))))) // check that there are no missing values in c. It happens for example if c() contains strings by mistake.
	if (`c_missing' == 1) {
		di as error "Error: please, check the option c(). It must contain only numbers and/or dots (.)"
		error 198
	}
	tempname cmatrix // c needs to be a matrix to pass it on to mata -> dump c into cmatrix
	mata: st_matrix("`cmatrix'", strtoreal(tokens(st_local("c"))))
	matrix colnames `cmatrix' = `cvars'
	return mat cmatrix = `cmatrix'
end validate_c 


capture program drop Display
program define Display
	version 10.0
	syntax [, Level(cilevel) *]
	_get_diopts diopts options , `options'
	
	_prefix_display, nohead level(`level') `diopts'
end Display

