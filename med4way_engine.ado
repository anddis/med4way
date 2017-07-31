*! Hello, I'm med4way_engine.ado
*! I'm a subprogram called by -med4way-

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
*********** 4-way decomposition with delta method / bootstrap SE ***************
********************************************************************************
	tempname bEstimates VEstimates

	// Estimation: delta method or bootstrap
	if (("`deltamethod'" == "false") & ("`bootstrap'" == "false")) {
		matrix `bEstimates' = J(1, `nn', 0)
		matrix `VEstimates' = J(`nn', `nn', 0)
	}
	else {
		mata: m4w_deriv(st_matrix("`betay'"),  st_matrix("`betam'"), /*
			*/ st_matrix("`Vy'"), st_matrix("`Vm'"), /*
			*/ st_matrix("`c'"), `nc', "`yreg'", "`mreg'", /*
			*/ st_matrix("`aam'"), "`bootstrap'", "`output'")			
	}
	
	matrix colnames `bEstimates' = `names'
	matrix colnames `VEstimates' = `names'
	matrix rownames `VEstimates' = `names'
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
	// the idea behind the onlybeta option is that, when regressml is called
	// within a bootstrap, it's useless to maximize the likelihood to get the SE
	// of sigma2 (the reason I wrote this program). 
	// Might as well return a matrix of 0 for e(V) and save time.
	
	marksample touse
	
	if "`onlybeta'" == "" local onlybeta "false"
	
	gettoken dep indep: varlist
	
	qui _regress `dep' `indep' if `touse'
	tempname sigma2reg
	matrix `sigma2reg' = e(rmse)^2 * (e(df_r) / e(N))
	matrix colnames `sigma2reg' = "sigma2:_cons"
	local nm = e(N)

	if "`onlybeta'"=="false" {
		tempname breg
		matrix `breg' = e(b)
		matrix coleq `breg' = "mu"
		
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
		
		local c : colfullnames e(b)
		matrix colnames `VML' = `c' "sigma2:_cons"
		matrix rownames `VML' = `c' "sigma2:_cons"
	} 
	
	ereturn post `bML' `VML', depname(`mvar') obs(`nm')
end regressml


/**********************
* m4w_deriv (mata)
**********************/
clear mata
set matastrict on 

mata: 
void function m4w_deriv(real vector betay, real vector betam, /*
					*/ real matrix Vy, real matrix Vm, real vector c, /* 
					*/ real scalar nc, string scalar yreg, string scalar mreg, /*
					*/ real vector aam, string scalar boot, string scalar output) {	

	real vector b, p_estimates
	real matrix V, dm_variance
	transmorphic D, G
	
	b = betay, betam
	V = blockdiag(Vy, Vm)
	
	D = deriv_init()
	deriv_init_evaluator(D, &m4w_formulas())
	deriv_init_evaluatortype(D, "t")
	deriv_init_params(D, b)
	deriv_init_argument(D, 1, c)
	deriv_init_argument(D, 2, nc)
	deriv_init_argument(D, 3, yreg)
	deriv_init_argument(D, 4, mreg)
	deriv_init_argument(D, 5, aam)
	deriv_init_argument(D, 6, output)

	m4w_formulas(b, c, nc, yreg, mreg, aam, output, p_estimates) // point estimates
	// note: it's much faster to call directly m4w_formulas than deriv(D, 0) to get p_estimates 
	
	if (boot == "false") {
		G = deriv(D, 1)
		dm_variance = G*V*G'	// variance covariance matrix
	}
	else if (boot == "true") {
		dm_variance = J(cols(p_estimates), cols(p_estimates), 0)
	}
	
	st_matrix(st_local("bEstimates"), p_estimates)
	st_matrix(st_local("VEstimates"), dm_variance)
}


/**********************
* m4w_formulas (mata)
**********************/
void m4w_formulas(real vector b, real vector c, /*
					*/ real scalar nc, string scalar yreg, string scalar mreg, /*
					*/ real vector aam, string scalar output, v) {

	real scalar offsetcox, betaTc, a0, a1, m, a1Ma0, a12Ma02, a1Tm, a0Tm, /*
			*/ pie, intmed, intref, cde, A, B, te, ereri_cde, ereri_intref, /*
			*/ ereri_intmed, ereri_pie, tereri
	real vector theta, beta, toreturn
	

	//--------------theta------------------- | -------------beta----------------
	// avar mvar inter cvar1 ... cvarK _cons | avar cvar1 ... cvarK _cons sigma2
	//  1    2    3    .......nc......  nc+4 |  1   .......nc......  nc+2  nc+3

	// offset to compensate the lack of _cons in theta if yreg == "cox"
	if (yreg=="cox") { 
		offsetcox = -1
	}
	else {
		offsetcox = 0
	}

	// split b into beta and theta - this makes the formulas below easier to read
	// especially because after this offsetcox is not needed anymore
	theta = b[|1 \ 4+nc+offsetcox|] 		// outcome model
	beta  = b[|5+nc+offsetcox \ .|] 		// mediator model

    if (nc>0) {
        betaTc = beta[|2 \ nc+1|]*c'		// linear combination covariates ("bcc" in vanderWeele)
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
		
		toreturn = (te, cde, intref, intmed, pie)
		if (output=="full") {
//			toreturn = te cde intref intmed pie p_cde p_intref p_intmed p_pie op_m op_ati op_e
//			always check that this order agrees with the one of the local macro `names'
			toreturn = (toreturn, cde/te, intref/te, intmed/te, pie/te, (pie+intmed)/te, (intref+intmed)/te, 1-(cde/te))
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

		toreturn = (tereri, ereri_cde, ereri_intref, ereri_intmed, ereri_pie)
		if (output=="full") {
// 			toreturn = tereri ereri_cde ereri_intref ereri_intmed ereri_pie terira p_cde p_intref p_intmed p_pie op_m op_ati op_e
//			always check that this order agrees with the one of the local macro `names'
			toreturn = (toreturn, (tereri+1), ereri_cde/tereri, ereri_intref/tereri, ereri_intmed/tereri, ereri_pie/tereri, (ereri_pie+ereri_intmed)/tereri, (ereri_intref+ereri_intmed)/tereri, 1-(ereri_cde/tereri))
		}
	}
	
	v = toreturn
}
end mata
