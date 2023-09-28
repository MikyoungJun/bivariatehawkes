##### list of code provided to perform analysis presented in the paper with comments #####

### data processing ###

- data_process.R


### model fitting ###

- uni_model.R
- bi_model.R

comments: 

1. spatial distance: distance can be calculated with chordal distance or Euclidean distance. As long as the domain is small (as in the examples in the manuscript), it does not make much difference. If you use Euclidean distance, you need to conver the estimated triggering distance to the necessary unit (e.g. Km)
2. time points are scaled to (0,1) for numerical stability. However, the results should not change even if you do not scale the time points
3. the integral term in the likelihood function of Hawkes process is approximated numerically with Riemann sum. 
4. For each example -- univariate for Afghanistan and bivariate for Nigeria -- the code for the winning model is given as an example. Minimal modification in the triggering function is needed to fit the other versions of models.
5. In the loglikelihood function, parameters have been transformed to satisfy conditions such as positivity and the stability condition discussed in the Supplemental material for bivariate case.
6. For numerical optimization, nlm and optim were used interchangeably to ensure that we obtain stable numerical optima (with Hessians for standard errors).


### handling fitted results ###

- uni_fitted.R
- bi_fitted.R

comments:



