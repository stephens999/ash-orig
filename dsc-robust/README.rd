A dsc to assess the stability of the EM estimates from ash.

This isn't really a comparison - it is just using the dscr package as a convenient way to do the simulations.

input = (betahat, sebetahat)
method (ash.wrapper) - runs ash 10 times on each dataset from a random start and returns a vector of the likelihoods found.
score - just returns the likelihoods.

 
