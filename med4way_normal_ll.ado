*! Hello, I'm med4way_normal_ll.ado
*! I'm a subprogram called by -med4way-

capture program drop med4way_normal_ll
program med4way_normal_ll
	version 10.0
	args lnfj mu sigma2
	*sigma parametrized as sqrt(sigma^2) to get directly the SE for sigma^2 
	quietly replace `lnfj' = lnnormalden($ML_y1, `mu', sqrt(`sigma2'))
end med4way_normal_ll
