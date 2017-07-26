*! Hello, this is med4way.ado
*! v2.1.0 - 26jul2017

/* 
Previous versions:
v2.0.0 - 27mar2017
*/

capture program drop med4way
program define med4way, eclass
	version 10.0

	if replay() {
		if ("`e(cmd)'" != "med4way") error 301
		
		ereturn display
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
	*/ ROBUST /*
	*/ NOLEGEND /* undocumented
	*/ NOWARNing /* undocumented
	*/ level(cilevel) /*
	*/ BOOTstrap /*
	*/ reps(integer 1000) /*
	*/ seed(passthru) /*
	*/ SAving(passthru) /*
	*/ BCA ]
		
	//[if] [in] marksample
	marksample touse
	
	// Step 1===================================================================
	// Checks and parsing of program options
	local wrnngtxt 0 
		// initialize wrnngtxt local. The idea behind wrnngtxt is
		// that, in case of Error in Step 1, no Warnng msgs are issued. 
		// Plus, it is possible to customize the order of the warnngs

	//parse yreg
	gettoken yregx yreg : yreg, parse(",")
	gettoken comma yreg : yreg
	if length("`comma'") > 1 {
		local 0 = substr("`comma'",2,.) + "`opt'"
 		local comma = substr("`comma'", 1, 1)
	}
	gettoken dist opt : yreg
	
	local yregx = trim(`"`yregx'"')
    local opt = trim(`"`opt'"')
	if `"`yregx'"'!="" & `"`opt'"'=="" {
		if `"`comma'"'=="" | (`"`comma'"'=="," & `"`dist'"'=="") {
			local yreg `"`yregx'"'
			local dist ""   
		}
		
		if `"`comma'"'=="," & `"`dist'"'!="" {
			 local yreg `"`yregx'"'
			 local dist `"`dist'"'	
		}
	}
	
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
	if "`bootstrap'" !="" {
		local bootstrap true
	}
	else {
		local bootstrap false
	}
	
	//suppress delta method
	if "`nodeltamethod'" !="" {
		local deltamethod false
	}
	else {
		local deltamethod true
	} 

	//full output or reduced output
	if "`fulloutput'" != "" {
		local output full
	}
	else {
		local output reduced
	}

	//casecontrol or not
	if "`casecontrol'" !="" {
		local casecontrol true
	}
	else {
		local casecontrol false
	}
	
	//legend or not
	if "`nolegend'" !="" {
		local legend false
	}
	else {
		local legend true
	}
	if "`deltamethod'"=="false" & "`bootstrap'"=="false" { 
		// if deltamethod not requested and bootstrap not requested, don't print the legend (there is no legend to print!)
		local legend false 
	}

	//warnings or not
	if "`nowarning'" !="" {
		local warning false
	}
	else {
		local warning true
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
	
	//validate rare outcome
	if ("`yreg'"=="logistic") {
		if "`casecontrol'"=="true" {
			local wrnngtxt `wrnngtxt' 2
		}
		else {
			qui logit `yvar' if `touse'
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
	//all the 8 suggested names have been used by existing variables, give an error message
	if _rc {
		display as error "Error: The command needs to create an interaction variable " /*
			*/ "with one of the following names: `inter_var_names', " /*
			*/ "but these variables have already been defined."
		error 110
	}
	gen `inter' = `avar'*`mvar'
	
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
	//==========================================================================	
	
	// Step 2===================================================================
	//print warnng text and summary of med4way
	if "`warning'"=="true" {
		local wrnngtxt : list sort wrnngtxt
		foreach w of local wrnngtxt {
			if `w' == 1 {
				display as error "Warning: this analysis assumes a rare outcome. " /*
					*/ "The outcome variable `yvar' has " /*
					*/  %2.0f `prev' `"% of "positive" outcomes. Consider a logbinomial "' /*
					*/ "or a Poisson model (with robust option) for the outcome."
			}
			if `w' == 2 {
				display as error "Warning: this analysis assumes a rare outcome."
			}
			if `w' == 3 {
				display as error "Warning: no censored event times. Please, check that " /*
					*/ "the data was stset correctly."
			}
			if `w' == 4 {
				display as error "Warning: no covariates specified via cvars(varlist)."
			}
			if `w' == 5 {
				display as error "Warning: fixed values for the covariates `cvar' " /*
					*/ "were not provided. All covariates are fixed at their means."
			}
			if `w' == 6 {
				display as error "Warning: fixed values for the covariates `cvar' " /*
					*/ "were provided only for some variables. Covariates are fixed at " /*
					*/ "the provided values or at their mean."
			}
			if `w' == 7 {
				display as error "Warning: it looks like the treatment variable `avar' " /*
					*/ "is binary. Treatment variable must be coded as 0/1 (`avar' " /*
					*/ "is coded as `loa' in the data)."
			}
		}
	}
	
	display _n(1) as text "-> Summary" _n(1)
	display _col(4) as text "Outcome    (yvar):"  _col(24) as res "`=cond("`survoutcome'"=="true", "[ survival outcome ]", "`yvar'")'"
	display _col(4) as text "Treatment  (avar):"  _col(24) as res "`avar'"
	display _col(4) as text "Mediator   (mvar):"  _col(24) as res "`mvar'"
	display _col(4) as text "Covariates (cvars):" _col(24) as res "`=cond(`nc'==0, "[ none ]", "`cvar'")'" _n(1)
	
	display _col(4) as text "Model for the outcome  (yreg):" _col(35) as res "`yreg'`=cond("`dist'"!="", ", `dist'", "")'"
	display _col(4) as text "Model for the mediator (mreg):" _col(35) as res "`mreg'" _n(1)
	
	display _col(4) as text "Referent treatment level (a0):" 			_col(46) as res "`a0'"
	display _col(4) as text "Actual treatment level   (a1):" 			_col(46) as res "`a1'"
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
	m4w_engine if `touse', yvar(`yvar') avar(`avar') mvar(`mvar') /*
		*/ cvar(`cvar')  c(`cmatrix') nc(`nc') inter(`inter') /*
		*/ aam(`aam') yreg(`yreg') mreg(`mreg') dist(`dist') /*
		*/ casecontrol(`casecontrol') output(`output') /*
		*/ deltamethod(`deltamethod') bootstrap(false) `robust' /*
		*/ names(`names') nn(`nnames') level(`level')

	if "`deltamethod'" == "true" {
		display _n(2) as text "-> 4-way decomposition: delta method" _n(1)
		ereturn display, level(`level')
	}
	//==========================================================================	

	// Step 4===================================================================	
	// Bootstrap (if desired)
	if "`bootstrap'" == "true" {
		qui bootstrap _b, reps(`reps') level(`level') `seed' `saving' nodots /*
			*/ level(`level') `bca': /*
				*/ m4w_engine if `touse', /*
				*/ yvar(`yvar') avar(`avar') mvar(`mvar') /*
				*/ cvar(`cvar')  c(`cmatrix') nc(`nc') inter(`inter') /*
				*/ aam(`aam') yreg(`yreg') mreg(`mreg') /*
				*/ dist(`dist') casecontrol(`casecontrol') output(`output') /*
				*/ deltamethod(false)  bootstrap(true) `robust' /*
				*/ names(`names') nn(`nnames') level(`level')
	
		display _n(2) as text "-> 4-way decomposition: bootstrap" _n(1)

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
		
		estat bootstrap, noheader `bca'
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
	ereturn post `b' `V'
	
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

	if "`bootstrap'"=="true" {
		ereturn scalar N_reps = `reps'
		ereturn local prefix = "bootstrap"
		
		if "`bca'" == "bca" {
			ereturn matrix ci_bca = `ci_bca'
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
	
	tempname cmatrix // c needs to be a matrix to pass it on to mata -> dump c into cmatrix
	if `nc' > 0 {
		local s = 1
		foreach i of local c   {
			if `s++' == 1 mat `cmatrix' = `i'
			else mat `cmatrix' = (`cmatrix' , `i')
		}
	}
	else {
		mat `cmatrix' = .
	}
	
	return mat cmatrix = `cmatrix'
end validate_c 


/**********************
* m4w_engine
**********************/
capture program drop m4w_engine
program define m4w_engine, eclass
	version 10.0	

	syntax [if], [ yvar(string) level(cilevel) ] avar(string) mvar(string) /*
		*/ [ cvar(varlist numeric) c(string) ] aam(string) /*
		*/ yreg(string) mreg(string) inter(string) nc(real) [ dist(string) ] /*
		*/ casecontrol(string) output(string) bootstrap(string) /*
		*/ deltamethod(string) [ robust ] names(string) nn(real)
	
	//[if] [in] marksample
 	marksample touse
	
********************************************************************************
*********** Fit model for the outcome and for the mediator *********************
********************************************************************************
	local titley `"_n(2) as text "-> Model for the outcome""'
	local titlem `"_n(2) as text "-> Model for the mediator""'
	
	tempname Vy betay Vm betam
		
	// Each block is a combination of a different yreg/mreg. There are 7 yreg
	// and 2 mreg = 14 blocks in total
	
	// Block 1: yreg=linear, mreg=linear
	if (("`yreg'"=="linear") & ("`mreg'"=="linear")) {
		display `titley'
		regress `yvar' `avar' `mvar' `inter' `cvar' if `touse', level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
		
		display `titlem'
		regress `mvar' `avar' `cvar' if `touse', level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 1
	
	// Block 2: yreg=linear, mreg=logistic
	else if (("`yreg'"=="linear") & ("`mreg'"=="logistic")) {
		display `titley'
		regress `yvar' `avar' `mvar' `inter' `cvar' if `touse', level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)

		display `titlem'
		logit `mvar' `avar' `cvar' if `touse', level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 2
	
	// Block 3: yreg=logistic, mreg=linear
	else if (("`yreg'"=="logistic") & ("`mreg'"=="linear")) {
		display `titley'
		logit `yvar' `avar' `mvar' `inter' `cvar' if `touse', level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		if "`casecontrol'"=="true" {
			regressml `mvar' `avar' `cvar' if `yvar' == 0 & `touse', /*
				*/ onlybeta(`bootstrap') level(`level')
		}
		else if "`casecontrol'"=="false" {
			regressml `mvar' `avar' `cvar' if `touse', /*
				*/ onlybeta(`bootstrap') level(`level')
		}
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 3
	
	// Block 4: yreg=logistic, mreg=linear
	else if (("`yreg'"=="logistic") & ("`mreg'"=="logistic")) {
		display `titley'
		logit `yvar' `avar' `mvar' `inter' `cvar' if `touse', level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		if "`casecontrol'"=="true" {
			logit `mvar' `avar' `cvar' if `yvar' == 0 & `touse', level(`level')
		}
		else if "`casecontrol'"=="false" {
			logit `mvar' `avar' `cvar' if `touse', level(`level')
		}
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 4
	
	// Block 5: yreg=cox, mreg=linear
	else if (("`yreg'"=="cox") & ("`mreg'"=="linear")) {
		display `titley'
		stcox `avar' `mvar' `inter' `cvar' if `touse', level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		regressml `mvar' `avar' `cvar' if `touse', /*
			*/ onlybeta(`bootstrap') level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 5
	
	// Block 6: yreg=cox, mreg=logistic
	else if (("`yreg'"=="cox") & ("`mreg'"=="logistic")) {
		display `titley'
		stcox `avar' `mvar' `inter' `cvar' if `touse', nohr level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		logit `mvar' `avar' `cvar' if `touse', level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 6
	
	// Block 7: yreg=aft, mreg=linear
	else if (("`yreg'"=="aft") & ("`mreg'"=="linear")) {
		display `titley'
		streg `avar' `mvar' `inter' `cvar' if `touse', time dist(`dist') level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
		if "`dist'"=="weibull" { // drop ancillary parameters
			matrix `Vy' = `Vy'[1..colsof(`Vy')-1, 1..colsof(`Vy')-1]
			matrix `betay' = `betay'[1, 1..colsof(`betay')-1]
		}
		
		display `titlem'
		regressml `mvar' `avar' `cvar' if `touse', /*
			*/ onlybeta(`bootstrap') level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 7
	
	// Block 8: yreg=aft, mreg=logistic
	else if (("`yreg'"=="aft") & ("`mreg'"=="logistic")) {
		display `titley'
		streg `avar' `mvar' `inter' `cvar' if `touse', time dist(`dist') level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
		if "`dist'"=="weibull" { // drop ancillary parameters
			matrix `Vy' = `Vy'[1..colsof(`Vy')-1, 1..colsof(`Vy')-1]
			matrix `betay' = `betay'[1, 1..colsof(`betay')-1]
		}
	
		display `titlem'
		logit `mvar' `avar' `cvar' if `touse', level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 8
	
	// Block 9: yreg=logbinomial, mreg=linear
	else if (("`yreg'"=="logbinomial") & ("`mreg'"=="linear")) {
		display `titley'
		glm `yvar' `avar' `mvar' `inter' `cvar' if `touse', family(binomial) /*
			*/ link(log) ml level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		regressml `mvar' `avar' `cvar' if `touse', /*
			*/ onlybeta(`bootstrap') level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 9
	
	// Block 10: yreg=logbinomial, mreg=logistic
	else if (("`yreg'"=="logbinomial") & ("`mreg'"=="logistic")) {
		display `titley'
		glm `yvar' `avar' `mvar' `inter' `cvar' if `touse', family(binomial) /*
			*/ link(log) ml level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		logit `mvar' `avar' `cvar' if `touse', level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 10
	
	// Block 11: yreg=poisson, mreg=linear
	else if (("`yreg'"=="poisson") & ("`mreg'"=="linear")) {
		display `titley'
		poisson `yvar' `avar' `mvar' `inter' `cvar' if `touse', `robust' level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		regressml `mvar' `avar' `cvar' if `touse', /*
			*/ onlybeta(`bootstrap') level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 11
	
	// Block 12: yreg=poisson, mreg=logistic
	else if (("`yreg'"=="poisson") & ("`mreg'"=="logistic")) {
		display `titley'
		poisson `yvar' `avar' `mvar' `inter' `cvar' if `touse', `robust' level(`level')
		matrix `Vy' = e(V)
		matrix `betay' = e(b)
	
		display `titlem'
		logit `mvar' `avar' `cvar' if `touse', level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 12
	
	 // Block 13: yreg=negbinomial, mreg=linear
	else if (("`yreg'"=="negbinomial") & ("`mreg'"=="linear")) {
		display `titley'
		nbreg `yvar' `avar' `mvar' `inter' `cvar' if `touse', level(`level')
		matrix `Vy' = e(V)
		matrix `Vy' = `Vy'[1..colsof(`Vy')-1, 1..colsof(`Vy')-1] // drop ancillary parameters
		matrix `betay' = e(b)
		matrix `betay' = `betay'[1, 1..colsof(`betay')-1] // drop ancillary parameters
	
		display `titlem'
		regressml `mvar' `avar' `cvar' if `touse', /*
			*/ onlybeta(`bootstrap') level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 13
	
	// Block 14: yreg=negbinomial, mreg=logistic
	else if (("`yreg'"=="negbinomial") & ("`mreg'"=="logistic")) {
		display `titley'
		nbreg `yvar' `avar' `mvar' `inter' `cvar' if `touse', level(`level')
		matrix `Vy' = e(V)
		matrix `Vy' = `Vy'[1..colsof(`Vy')-1, 1..colsof(`Vy')-1] // drop ancillary parameters
		matrix `betay' = e(b)
		matrix `betay' = `betay'[1, 1..colsof(`betay')-1] // drop ancillary parameters
	
		display `titlem'
		logit `mvar' `avar' `cvar' if `touse', level(`level')
		matrix `Vm' = e(V)
		matrix `betam' = e(b)
	}	
	// End of Block 14
********************************************************************************
	
	
********************************************************************************
*********** 4-way decomposition with delta method SE if needed *****************
********************************************************************************
	// Step 1===================================================================
	// Prepare output empty matrices
	tempname bEstimates VEstimates

	matrix `bEstimates' = J(1, `nn', 0)
	matrix `VEstimates' = J(`nn', `nn', 0)
		
	matrix colnames `bEstimates' = `names'
	matrix colnames `VEstimates' = `names'
	matrix rownames `VEstimates' = `names'
	//==========================================================================

	// Step 2===================================================================
	// Estimation: delta method or bootstrap
	if (("`deltamethod'" == "true") & ("`bootstrap'" == "false")) {
		foreach j of local names {
			mata: m4w_deriv("`j'", "`betay'",  "`betam'", "`Vy'", "`Vm'", /*
				*/ "`c'", `nc', "`yreg'", "`mreg'", "`aam'", "false")
			matrix `bEstimates'[1, colnumb(`bEstimates', "`j'")] = `p_estimate'
			matrix `VEstimates'[colnumb(`VEstimates', "`j'"), /*
				*/ rownumb(`VEstimates', "`j'")] = `dm_variance'
		}
	}
	else if (("`deltamethod'" == "false") & ("`bootstrap'" == "true")) {
		foreach j of local names {
			mata: m4w_deriv("`j'", "`betay'",  "`betam'", "`Vy'", "`Vm'", /*
				*/ "`c'", `nc', "`yreg'", "`mreg'", "`aam'", "true")			
			matrix `bEstimates'[1, colnumb(`bEstimates', "`j'")] = `p_estimate'
		}
	}
	//==========================================================================
********************************************************************************
	  	
	ereturn post `bEstimates' `VEstimates'
	
	ereturn local names = "`names'"
	ereturn scalar nnames = `nn'
end m4w_engine


/**********************
* normalml_lf
**********************/
capture program drop normalml_lf
program normalml_lf
	version 10.0
	args lnfj mu sigma2
	*sigma parametrized as sqrt(sigma^2) to get directly the SE for sigma^2 
	quietly replace `lnfj' = lnnormalden($ML_y1, `mu', sqrt(`sigma2'))
end normalml_lf


/**********************
* regressml
**********************/
capture program drop regressml
program define regressml, eclass
	version 10.0
	syntax varlist(min=2 numeric) [if] [, onlybeta(string) level(cilevel)]
	// the idea behind the onlybeta option is that, when regressML is called
	// within a bootstrap, it's useless to maximize the likelihood to get the SE
	// of sigma2 (the reason I wrote this program). 
	// Might as well return a matrix of 0 for e(V) and save time.
	
	marksample touse
	
	if "`onlybeta'" == "" local onlybeta "false"
	
	gettoken dep indep: varlist
	
	qui _regress `dep' `indep' if `touse'
	tempname breg sigma2reg
	matrix `breg' = e(b)
	matrix coleq `breg' = "mu"
	matrix `sigma2reg' = e(rmse)^2 * (e(df_r) / e(N))
	matrix colnames `sigma2reg' = "sigma2:_cons"

	if "`onlybeta'"=="false" {
		qui ml model lf normalml_lf (mu: `dep' = `indep') (sigma2:) if `touse', /*
			*/ title("Linear regression (Maximum Likelihood)") //waldtest(0)
		qui ml init `breg' `sigma2reg'
		ml maximize, search(off) level(`level')
	}
	
	tempname bML VML
	if "`onlybeta'" =="false" {
		matrix `bML' = e(b)
		matrix `VML' = e(V)	
	}
	else if "`onlybeta'"== "true" {
		matrix `bML' = e(b), `sigma2reg'
		matrix `VML' = J(colsof(e(b))+1, colsof(e(b))+1, 0)
	}
	local c : colfullnames `bML'
	local c : subinstr local c "sigma2:_cons" "sigma2:sigma2"
	matrix colnames `bML' = `c'
	matrix colnames `VML' = `c'
	matrix rownames `VML' = `c' 
	local nm = e(N)
	
	ereturn post `bML' `VML', depname(`mvar') obs(`nm')
end regressml


/**********************
* m4w_deriv (mata)
**********************/
clear mata
mata: 
mata set matastrict on 

void function m4w_deriv(string scalar t, string scalar betay, string scalar betam, string scalar Vy, string scalar Vm, string scalar c, real scalar nc, string scalar yreg, string scalar mreg, string vector aam, string scalar boot) 
{	

	real vector betayx, betamx, b, aamx, cx
	real matrix Vyx, Vmx, V 
	real scalar p_estimate, dm_variance 
	pointer p
	transmorphic D, G

	betayx = st_matrix(betay)
	betamx = st_matrix(betam)

	Vyx = st_matrix(Vy)
	Vmx = st_matrix(Vm)
	
	b = betayx, betamx
	V = blockdiag(Vyx, Vmx)
	
	cx = st_matrix(c)
	aamx = st_matrix(aam)

	p = findexternal("m4w_formulas()") // pointer!
	
	D = deriv_init()
	deriv_init_evaluator(D, &m4w_eval())
	deriv_init_params(D, b)
	deriv_init_argument(D, 1, p)
	deriv_init_argument(D, 2, t)
	deriv_init_argument(D, 3, cx)
	deriv_init_argument(D, 4, nc)
	deriv_init_argument(D, 5, yreg)
	deriv_init_argument(D, 6, mreg)
	deriv_init_argument(D, 7, aamx)

	p_estimate = deriv(D, 0)
	st_local("p_estimate", strofreal(p_estimate))
	
	if (boot == "false") {
		G = deriv(D, 1)
		dm_variance = (G*V*G')
		st_local("dm_variance", strofreal(dm_variance))
	}
}


/**********************
* m4w_eval (mata)
**********************/
void m4w_eval(real vector b, pointer(real scalar function) scalar f, string scalar t, real vector c, real scalar nc, string scalar yreg, string scalar mreg, real vector aam, v)
{
	v = (*f)(b, t, c, nc, yreg, mreg, aam)
}


/**********************
* m4w_formulas (mata)
**********************/
real scalar m4w_formulas(real vector b, string scalar t, real vector c, real scalar nc, string scalar yreg, string scalar mreg, real vector aam)
{

	real scalar offsetcox, i, betaTc, a0, a1, m, a1Ma0, a12Ma02, a1Tm, a0Tm, pie, intmed, intref, cde, A, B, te, ereri_cde, ereri_intref, ereri_intmed, ereri_pie, tereri
	real vector theta, beta, bb
	

	//--------------theta------------------- | -------------beta----------------
	// avar mvar inter cvar1 ... cvarK _cons | avar cvar1 ... cvarK _cons sigma2
	//  1    2    3    .......nc......  nc+4 |  1   .......nc......  nc+2  nc+3

	// offset for cox to compensate the lack of _cons in theta
	if (yreg=="cox") { 
		offsetcox = -1
	}
	else {
		offsetcox = 0
	}

	// split b into beta and theta - this makes the formulas below easier to read
	// especially because offsetcox is no longer needed
	theta = b[1::4+nc+offsetcox] 		// outcome model
	beta  = b[5+nc+offsetcox::cols(b)] 	// mediator model

	if (nc>0) {
		bb = J(1, nc, .)
		for (i=1; i<=nc; i++) {
			bb[i] = beta[1+i]
		}
		betaTc = bb*c'
	}
	else {
		betaTc = 0
	}

	a0 = aam[1]
	a1 = aam[2]
	m  = aam[3]

	a1Ma0 = a1 - a0
	a12Ma02 = a1*a1 - a0*a0
	a1Tm = a1 * m
	a0Tm = a0 * m
	 
	if (yreg=="linear") {															// yreg linear

		if (mreg=="linear") {														// mreg linear
			pie = ((theta[2]*beta[1]+theta[3]*beta[1]*a0)*(a1Ma0))
			intmed = (theta[3]*beta[1]*(a1Ma0)*(a1Ma0))
			intref = (theta[3]*(beta[nc+2]+beta[1]*a0+betaTc-m)*(a1Ma0))
			cde = ((theta[1]+theta[3]*m)*(a1Ma0))
		}
		else if (mreg=="logistic") {												// mreg logistic
			A = (beta[nc+2]+beta[1]*a1+betaTc)
			B = (beta[nc+2]+beta[1]*a0+betaTc)
			
			pie = ((theta[2]+theta[3]*a0)*(exp(A)/(1+exp(A)) - exp(B)/(1+exp(B))))
			intmed = (theta[3]*(a1Ma0)*(exp(A)/(1+exp(A)) - exp(B)/(1+exp(B))))
			intref = (theta[3]*(a1Ma0)*(exp(B)/(1+exp(B)) - m))
			cde = ((theta[1]+theta[3]*m)*(a1Ma0))
		}
		te = (cde + pie + intmed + intref)
		
		if (t == "cde") {
			return(cde)
		}
		if (t == "intref") {
			return(intref)
		}
		if (t == "intmed") {
			return(intmed)
		}
		if (t == "pie") {
			return(pie)
		}
		if (t == "te") {
			return(te)
		}
		if (t == "p_cde") {
			return(cde/te)
		}
		if (t == "p_intref") {
			return(intref/te)
		}
		if (t == "p_intmed") {
			return(intmed/te)
		}
		if (t == "p_pie") {
			return(pie/te)
		}
		if (t == "op_m") {
			return((pie+intmed)/te)
		}
		if (t == "op_ati") {
			return((intref+intmed)/te)
		}
		if (t == "op_e") {
			return(1-(cde/te))
		}
	}

	else if  ((yreg=="logistic") | (yreg=="cox") | (yreg=="aft") | /*
			*/ (yreg=="logbinomial") | (yreg=="poisson") | (yreg=="negbinomial")) {	// yreg logistic
		
		if (mreg=="linear") {														// mreg linear
			ereri_cde = (exp(theta[1]*(a1Ma0)+theta[2]*m+theta[3]*a1Tm-(theta[2]+theta[3]*a0)*(beta[nc+2]+beta[1]*a0+betaTc)-0.5*(theta[2]+theta[3]*a0)*(theta[2]+theta[3]*a0)*beta[nc+3])-exp(theta[2]*m+theta[3]*a0Tm-(theta[2]+theta[3]*a0)*(beta[nc+2]+beta[1]*a0+betaTc)-0.5*(theta[2]+theta[3]*a0)*(theta[2]+theta[3]*a0)*beta[nc+3]))
			
			ereri_intref = (exp((theta[1]+theta[3]*(beta[nc+2]+beta[1]*a0+betaTc+theta[2]*beta[nc+3]))*(a1Ma0)+0.5*theta[3]*theta[3]*beta[nc+3]*(a12Ma02))-1-exp(theta[1]*(a1Ma0)+theta[2]*m+theta[3]*a1Tm-(theta[2]+theta[3]*a0)*(beta[nc+2]+beta[1]*a0+betaTc)-0.5*(theta[2]+theta[3]*a0)*(theta[2]+theta[3]*a0)*beta[nc+3])+exp(theta[2]*m+theta[3]*a0Tm-(theta[2]+theta[3]*a0)*(beta[nc+2]+beta[1]*a0+betaTc)-0.5*(theta[2]+theta[3]*a0)*(theta[2]+theta[3]*a0)*beta[nc+3]))
			
			ereri_intmed = (exp((theta[1]+theta[2]*beta[1]+theta[3]*(beta[nc+2]+beta[1]*a0+beta[1]*a1+betaTc+theta[2]*beta[nc+3]))*(a1Ma0)+0.5*theta[3]*theta[3]*beta[nc+3]*(a12Ma02))-exp((theta[2]*beta[1]+theta[3]*beta[1]*a0)*(a1Ma0))-exp((theta[1]+theta[3]*(beta[nc+2]+beta[1]*a0+betaTc+theta[2]*beta[nc+3]))*(a1Ma0)+0.5*theta[3]*theta[3]*beta[nc+3]*(a12Ma02))+1)
			
			ereri_pie = (exp((theta[2]*beta[1]+theta[3]*beta[1]*a0)*(a1Ma0))-1)
			
			tereri = ((exp((theta[1]+theta[3]*(beta[nc+2]+beta[1]*a0+betaTc+theta[2]*beta[nc+3]))*(a1Ma0)+0.5*theta[3]*theta[3]*beta[nc+3]*(a12Ma02))*exp((theta[2]*beta[1]+theta[3]*beta[1]*a1)*(a1Ma0)))-1)
		}
		
		else if (mreg=="logistic") {												// mreg logistic
			ereri_cde = (exp(theta[1]*(a1Ma0)+theta[2]*m+theta[3]*a1Tm)*(1+exp(beta[nc+2]+beta[1]*a0+betaTc))/(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0))-exp(theta[2]*m+theta[3]*a0Tm)*(1+exp(beta[nc+2]+beta[1]*a0+betaTc))/(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0)))
			
			ereri_intref = (exp(theta[1]*(a1Ma0))*(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a1))/(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0))-1-exp(theta[1]*(a1Ma0)+theta[2]*m+theta[3]*a1Tm)*(1+exp(beta[nc+2]+beta[1]*a0+betaTc))/(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0))+exp(theta[2]*m+theta[3]*a0Tm)*(1+exp(beta[nc+2]+beta[1]*a0+betaTc))/(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0)))
			
			ereri_intmed = (exp(theta[1]*(a1Ma0))*(1+exp(beta[nc+2]+beta[1]*a1+betaTc+theta[2]+theta[3]*a1))*(1+exp(beta[nc+2]+beta[1]*a0+betaTc))/((1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0))*(1+exp(beta[nc+2]+beta[1]*a1+betaTc)))-(1+exp(beta[nc+2]+beta[1]*a1+betaTc+theta[2]+theta[3]*a0))*(1+exp(beta[nc+2]+beta[1]*a0+betaTc))/((1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0))*(1+exp(beta[nc+2]+beta[1]*a1+betaTc)))-exp(theta[1]*(a1Ma0))*(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a1))/(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0))+1)
			
			ereri_pie = ((1+exp(beta[nc+2]+beta[1]*a0+betaTc))*(1+exp(beta[nc+2]+beta[1]*a1+betaTc+theta[2]+theta[3]*a0))/((1+exp(beta[nc+2]+beta[1]*a1+betaTc))*(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0)))-1)
			
			tereri = ((exp(theta[1]*a1)*(1+exp(beta[nc+2]+beta[1]*a0+betaTc))*(1+exp(beta[nc+2]+beta[1]*a1+betaTc+theta[2]+theta[3]*a1))/(exp(theta[1]*a0)*(1+exp(beta[nc+2]+beta[1]*a1+betaTc))*(1+exp(beta[nc+2]+beta[1]*a0+betaTc+theta[2]+theta[3]*a0))))-1)
		}

		if (t == "ereri_cde") {
			return(ereri_cde)
		}
		if (t == "ereri_intref") {
			return(ereri_intref)
		}
		if (t == "ereri_intmed") {
			return(ereri_intmed)
		}
		if (t == "ereri_pie") {
			return(ereri_pie)
		}
		if (t == "tereri") {
			return(tereri)
		}
		if (t == "terira") {
			return((tereri+1))
		}
		if (t == "p_cde") {
			return(ereri_cde/tereri)
		}
		if (t == "p_intref") {
			return(ereri_intref/tereri)
		}
		if (t == "p_intmed") {
			return(ereri_intmed/tereri)
		}
		if (t == "p_pie") {
			return(ereri_pie/tereri)
		}
		if (t == "op_m") {
			return((ereri_pie+ereri_intmed)/tereri)
		}
		if (t == "op_ati") {
			return((ereri_intref+ereri_intmed)/tereri)
		}
		if (t == "op_e") {
			return(1-(ereri_cde/tereri))
		}
	}
}
end mata
