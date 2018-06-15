{smcl}
{* *! version 2.2.2 14jun2018}{...}

{cmd:help med4way}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{hi:med4way} {hline 2}}4-way decomposition using parametric regression models{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 13 2}
{cmd:med4way} [{it:yvar}] {it:avar} {it:mvar} [{it:cvars}] {ifin}, {opt a0(#)} {opt a1(#)} {opt m(#)} 
	{opt yreg(string)} {opt mreg(string)} [ {it:options} ]


{phang}
{it:yvar} is the variable name for the outcome. Note that {it:yvar} must not be specified when the model for the outcome is an Accelerated Failure Time or a Cox proportional hazards model. In these cases, you must {helpb stset} your data before running med4way.

{phang}
{it:avar} is the variable name for the treatment. If binary, it must be coded as 0/1.

{phang}
{it:mvar} is the variable name for the mediator. If binary, it must be coded as 0/1.

{phang}
{it:cvars} are the variable names for the covariates to be included in the model for the outcome and for the mediator.


{synoptset 31 tabbed}{...}
{synopthdr :options}
{synoptline}
{p2coldent: * {opt a0(#)}}referent treatment level{p_end}
{p2coldent: * {opt a1(#)}}actual treatment level{p_end}
{p2coldent: * {opt m(#)}}level of the mediator at which to compute the 4-way decomposition{p_end}
{p2coldent: * {opt yreg(string)}}form of the regression model for the outcome{p_end}
{p2coldent: * {opt mreg(string)}}form of the regression model for the mediator{p_end}

{synopt :{opth c(numlist)}}values of the covariates {it:cvars} at which to compute the 4-way decomposition{p_end}
{synopt :{opt casec:ontrol}}treat the data as coming from a case-control study{p_end}
{synopt :{opt full:output}}compute also the proportion of the Total Effect attributable to its 4 components and the proportions attributable to mediation, interaction, and either mediation or interaction or both{p_end}
{synopt :{opt nodeltam:ethod}}suppress the calculation of the standard errors using the delta method{p_end}
{synopt :{opt robust}}use the sandwich/robust estimator of variance when using a Poisson model for the outcome{p_end}

{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{it:display_options}}control formatting{p_end}
{synopt :{opt coef:legend}}display legend instead of statistics{p_end}

{synopt :{opt boot:strap}}use bootstrap to calculate confidence intervals{p_end}
{synopt :{opt reps(#)}}perform # bootstrap replications; default is {cmd: reps(1000)}{p_end}
{synopt :{opt seed(#)}}set random-number seed to #{p_end}
{synopt :{help prefix_saving_option:{bf:{ul:sa}ving(}{it:filename}{bf:, ...)}}}save results to filename; save statistics in double precision; save results to filename every # replications{p_end}
{synopt :{opt bca}}compute acceleration for BCa confidence intervals{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* are required.{p_end}


{title:Description}

{pstd}
{cmd:med4way} uses parametric regression models to estimate the components of the 4-way decomposition of the total effect of treatment {it:avar} on outcome {it:yvar} in the presence of mediator {it:mvar} with which the exposure may interact. This decomposition breaks down the total effect of the treatment on the outcome into components due to mediation alone, to interaction alone, to both mediation and interaction, and to neither mediation nor interaction. 

{pstd}
{cmd:med4way} provides standard errors and confidence intervals for the estimated components using the delta method (default) or the bootstrap.

{pstd}
Note: the 4-way decomposition holds without any assumptions about confounding. However, to interpret each of the components causally does require assumptions about confounding. See VanderWeele (2014) for a detailed exposition of those assumptions.


{title:Options}

{phang}
{opt a0(#)} specifies the referent level of the treatment.

{phang}
{opt a1(#)} specifies the actual level of the treatment.

{phang}
{opt m(#)} specifies the level of the mediator at which the 4-way decomposition is computed.

{phang}
{opt yreg(string)} specifies the form of the regression model for the outcome. The available forms are:{p_end}

{p 7 6 2}{opt lin:ear}: linear regression{p_end}
{p 7 6 2}{opt logi:stic}: logistic regression {p_end}
{p 7 6 2}{opt logb:inomial}: logbinomial regression (GLM with binomial distribution and log link function){p_end}
{p 7 6 2}{opt poi:sson}: Poisson regression{p_end}
{p 7 6 2}{opt negb:inomial}: negative binomial regression{p_end}
{p 7 6 2}{opt aft, {ul on}e{ul off}xponential}: Accelerated Failure Time (exponential survival distribution) ({helpb stset} required){p_end}
{p 7 6 2}{opt aft, {ul on}w{ul off}eibull}: Accelerated Failure Time (Weibull survival distribution) ({helpb stset} required){p_end}
{p 7 6 2}{opt cox}: Cox proportional hazards model ({helpb stset} required){p_end}

{phang}
{opt mreg(string)} specifies the form of the regression model for the mediator. The available forms are:

{p 7 6 2}{opt lin:ear}: linear regression{p_end}
{p 7 6 2}{opt logi:stic}: logistic regression regression{p_end}

{phang}
{opt c(string)} fixes the values of the covariates {it:cvars} at which to compute the 4-way decomposition. If {it:cvars} are specified but this option is omitted, {it:cvars} will be automatically fixed at their respective mean values. If this option is specified, the number of values of {opt c(numlist)} must correspond to the number of {it:cvars}. A dot (.) can be used to fix the value of a specific covariate to its mean. 
Example: the covariates specified are {it:cvar1 cvar2 cvar3} and the user wants to fix the value for {it:cvar2} to 6, while letting the values for {it:cvar1} and {it:cvar3} to be equal to their respective means. This can be achieved with the option {opt c(. 6 .)}.

{phang}
{opt casec:ontrol} specifies that the data comes from a case-control study (that is, sampling was done on the outcome).

{phang}
{opt full:output} specifies that, in addition to the 4 components of the total effect (controlled direct effect, reference interaction, mediated interaction, pure indirect effect), the following quantities are to be estimated: the proportions of the total effect due to each of the 4 components, the overall proportion mediated, the overall proportion due to interaction, and the overall proportion that would be eliminated if the mediator {it:mvar} were fixed to the value {opt m(#)}. 

{phang}
{opt nodeltam:ethod} suppresses the calculation of the standard errors of the estimated quantities using the delta method.

{phang}
{opt robust} uses the sandwich/robust estimator of the variance-covariance matrix when a Poisson regression model for the outcome is specified.

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

{phang}
{it:display_options}: {opt noci}, {opt nopv:alues}, {bf:cformat(}%fmt{bf:)}, {bf:pformat(}%fmt{bf:)}; see {helpb estimation options##coeflegend:[R] estimation options}.

{phang}
{opt coeflegend}; see
     {helpb estimation options##coeflegend:[R] estimation options}.

{phang}
{opt boot:strap} uses the bootstrap to calculate confidence intervals (CIs). By default, normal-based CIs are calculated. Bias-corrected and percentile CIs can be obtained with the post-estimation command {helpb estat bootstrap}.

{phang}
{opt reps(#)} specifies the number of bootstrap replications to be performed. The default is 1,000.

{phang}
{opt seed(#)} sets the random-number seed.

{phang}
{bf:{ul:sa}ving(}{it:filename} [{bf:,} {it:suboptions}]{bf:)} creates a Stata data file ({bf:.dta} file) consisting of, for each  quantity estimated by med4way, a variable containing the replicates. 

{pmore}
See {help prefix_saving_option} for details about {it:suboptions}.

{phang}
{opt bca} estimates the acceleration used to construct bias-corrected and accelerated confidence intervals.


{title:Examples}

{pstd}Download example datasets in the current working directory{p_end}
{phang2}{stata `"net get med4way, from("https://raw.githubusercontent.com/anddis/med4way/master/")"'}{p_end}

    {title:Binary outcome, binary mediator}

{pstd}Load simulated data{p_end}
{phang2}{stata use med4way_example_1.dta}{p_end}

{pstd}Logistic regression model for the outcome; Logistic regression model for the mediator; Delta method standard errors{p_end}
{phang2}{stata med4way y_bin treat m_bin cvar1 cvar2 cvar3, yreg(logistic) mreg(logistic) a0(0) a1(1) m(0) c(0 0 .)}

{pstd}Logistic regression model for the outcome; Logistic regression model for the mediator; Bootstrap confidence intervals (500 replicates){p_end}
{phang2}{stata med4way y_bin treat m_bin cvar1 cvar2 cvar3, yreg(logistic) mreg(logistic) a0(0) a1(1) m(0) c(0 0 .) bootstrap reps(500) seed(1234)}

{phang2}{stata estat bootstrap, all noheader}

    {title:Survival outcome, binary mediator}
    
{pstd}Load simulated data{p_end}
{phang2}{stata use med4way_example_2.dta}{p_end}

{pstd}Declare data to be survival-time data{p_end}
{phang2}{stata stset y_cens, failure(fail) noshow}

{pstd}Accelerated failure time regression model (exponential survival distribution) for the outcome; Logistic regression model for the mediator{p_end}
{phang2}{stata med4way treat m_bin cvar1, a0(0) a1(1) m(0) yreg(aft, e) mreg(logistic) c(1)}

{pstd}Show the legend of the coefficients{p_end}
{phang2}{stata med4way, coeflegend}

{pstd}Test whether the excess relative risk due to controlled direct effect (ereri_cde) is statistically different from the excess relative risk due to pure indirect effect (ereri_pie){p_end}
{phang2}{stata test _b[ereri_cde] = _b[ereri_pie]}{p_end}

{pstd}Given the 4 basic components of the total effect, additional derived quantities can be estimated with the post-estimation commands {helpb lincom} or {helpb nlcom}, as appropriate.
For example, to calculate the overall proportion mediated (op_m){p_end}
{phang2}{stata nlcom (_b[ereri_pie]+_b[ereri_intmed])/(_b[ereri_cde]+_b[ereri_intref]+_b[ereri_intmed]+_b[ereri_pie]), noheader}{p_end}

{pstd}Accelerated failure time regression model (exponential survival distribution) for the outcome; Logistic regression model for the mediator; Display full output{p_end}
{phang2}{stata med4way treat m_bin cvar1, a0(0) a1(1) m(0) yreg(aft, e) mreg(logistic) c(1) fulloutput}{p_end}


{title:Stored results}

{pstd}
{cmd:med4way} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(a0)}}referent treatment level{p_end}
{synopt:{cmd:e(a1)}}actual treatment level{p_end}
{synopt:{cmd:e(m)}}level of the mediator used for the 4-way decomposition{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{p2coldent: * {cmd:e(N_reps)}}number of requested bootstrap replications{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(cmd)}}{cmd:med4way}{p_end}
{synopt:{cmd:e(yvar)}}name of outcome variable{p_end}
{synopt:{cmd:e(avar)}}name of treatment variable{p_end}
{synopt:{cmd:e(mvar)}}name of mediator variable{p_end}
{synopt:{cmd:e(cvars)}}name of covariate variables{p_end}
{synopt:{cmd:e(yreg)}}form of the regression model for the outcome{p_end}
{synopt:{cmd:e(mreg)}}form of the regression model for the mediator{p_end}
{synopt:{cmd:e(estimands)}}acronyms of the estimated quantities{p_end}
{p2coldent: * {cmd:e(prefix)}}{cmd:bootstrap}{p_end}
{p2coldent: * {cmd:e(vce)}}{cmd:bootstrap}{p_end}
{p2coldent: * {cmd:e(vcetype)}}{cmd:Bootstrap}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(c)}}values of {it:cvars} at which the 4-way decomposition was computed{p_end}
{p2coldent: * {cmd:e(b_bs)}}bootstrap estimates{p_end}
{p2coldent: * {cmd:e(reps)}}number of nonmissing results{p_end}
{p2coldent: * {cmd:e(bias)}}estimated biases{p_end}
{p2coldent: * {cmd:e(se)}}estimated standard errors{p_end}
{p2coldent: * {cmd:e(z0)}}median biases{p_end}
{p2coldent: * {cmd:e(accel)}}estimated accelerations{p_end}
{p2coldent: * {cmd:e(ci_normal)}}normal-approximation CIs{p_end}
{p2coldent: * {cmd:e(ci_percentile)}}percentile CIs{p_end}
{p2coldent: * {cmd:e(ci_bc)}}bias-corrected CIs{p_end}
{p2coldent: * {cmd:e(ci_bca)}}bias-corrected and accelerated CIs{p_end}
{p2coldent: * {cmd:e(V_modelbased)}}model-based variance{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{p 4 6 2}* are stored only if the option {opt bootstrap} is specified.{p_end}


{title:References}

{phang}VanderWeele, T.J., 2014. A unification of mediation and interaction: a 4-way decomposition. Epidemiology (Cambridge, Mass.), 25(5), p.749.


{title:Authors}

{pstd}Andrea Discacciati [1], Andrea Bellavia [2,3], Linda Valeri [4,5]{p_end}

{pstd}[1] {it:Unit of Biostatistics, Karolinska Institutet, Stockholm, Sweden}{p_end}
{pstd}[2] {it:Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA}{p_end}
{pstd}[3] {it:Department of Biostatistics, Harvard T.H. Chan School of Public Health, Boston, MA}{p_end}
{pstd}[4] {it:Department of Psychiatry, Harvard Medical School, Boston, MA}{p_end}
{pstd}[5] {it:Psychiatric Biostatistics Laboratory, McLean Hospital, Belmont, MA}{p_end}
