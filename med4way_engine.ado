*! Hello, I am med4way_engine.ado
*! I am a program called by -med4way-
*! v2.1.0 - 26jul2017

capture program drop med4way_engine
program define med4way_engine, eclass
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
end med4way_engine


/*********************
*
* 	 SUBROUTINES
*
**********************/


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
		qui ml model lf med4way_normal_ll (mu: `dep' = `indep') (sigma2:) if `touse', /*
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
