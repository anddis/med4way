*! Hello, I'm med4way_normal_lf2.ado
*! I'm a subprogram called by -med4way-
cap program drop med4way_normal_lf2
program med4way_normal_lf2 
        version 11
        args todo b lnfj g1 g2 H
        tempvar mu 
		tempname sigma2 sigma
        mleval `mu' = `b', eq(1)
        mleval `sigma2' = `b', eq(2) scalar
        scalar `sigma' = sqrt(`sigma2')
        quietly {
                replace `lnfj' = ln( normalden($ML_y1,`mu',`sigma') )
                if (`todo'==0) exit

                tempvar z
                gen double `z' = ($ML_y1-`mu')/`sigma'
                replace `g1' = `z'/`sigma'
                replace `g2' = -1/(2*`sigma'^2) * (1 - `z'*`z')
                if (`todo'==1) exit

                tempname d11 d12 d22
                mlmatsum `lnfj' `d11' = -1/`sigma'^2      , eq(1)
                mlmatsum `lnfj' `d12' = -`z'/`sigma'^3  , eq(1,2)
                mlmatsum `lnfj' `d22' = 1/(2*`sigma'^4) * (1 - 2*`z'*`z')     , eq(2)
                matrix `H' = (`d11', `d12' \ `d12'', `d22')
        }
end
