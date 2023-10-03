##### An example for calculating asymptotic standard errors based on Hessian matrix given by the optimization routine #####

### We show an example of M2-6 model for bivariate case
### Note that in the likelihood optimization, we used parameter transformation to ensure the numerical optimization routine stays in the correct parameter space. 
### However, calculation of Jacobian for some of the parameters is quite complex. Therefore, once we obtain parameter estimates, we use another R package to calculate 
### Hessian matrices as illustrated below.


library(numDeriv)

hessian(fun
