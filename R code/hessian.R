##### An example for calculating asymptotic standard errors based on Hessian matrix given by the optimization routine #####

### We show an example of M2-6 model for bivariate case
### Note that in the likelihood optimization, we used parameter transformation to ensure the numerical optimization routine stays in the correct parameter space. 
### However, calculation of Jacobian for some of the parameters is quite complex. Therefore, once we obtain parameter estimates, we use another R package to calculate 
### Hessian matrices as illustrated below.

## First, rewrite the (negative) logLiklelihood function in bi_model.R (or univ_model.R for univariate problem) WITHOUT parameter transformation 
## Let us call such function logLik1
## save the values of parameter estimates in "par" 

library(numDeriv)

Hes= hessian(func=logLik1, x=par) ## Hessian matrix of negative loglikelihood function evaluated at the parameter estimates

sqrt(diag(solve(Hes))) ## this should give you the asymptotic standard errors of parameter estimates 
