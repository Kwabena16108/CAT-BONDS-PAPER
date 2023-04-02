### The upper bound of the original integral is infinite, from 0 to +inf.
### We apply u=g(x) = 1/(1+x) tranfrom to get an integral with finite boundary, 
### from 0 to 1.  
### The first step is to estimate the E[f(x; theta)] on a grid of x values,
### where the expectation is taken under the posterior distribution of
### theta = (w1, w2, lambda1, lambda2), parameters contained in mixture 
### of two exponentials. In practice, we first generate our grid on interval [0,1] and then apply the 
### the inverse of g(x) to get back. 
### There are various ways of approximating an integral, such as Sympson's rule, 
### midpoint rule, and Trapezoidal rule. We opt for the second method which is 
### simple and straightforward tobe implemented and that gives reliable reults as the number of 
### subintervals is sufficiently large. 
### Define the limits of integration
a<- 0
b<- 1
###Number of subintervals 
N= 100
### The length of steps 
DeltaX= (b-a)/N
### Generate a sequence between 0 and 1 with step size of DeltaX:
h = seq(0, 1, by = DeltaX) 
### Generate midpoints. 
m= (h[2:length(h)]+h[1:(length(h)-1)])/2
### Generate a grid of x values corresponding to the midpoints stored in vecotor m
### uisng  the inverse of g(x) function. 
x_grid= (1-m)/m
### Lets estimate the E[f(x; theta)] uisng Gibbs sampling on a grid of "x_grid" values 
library(R2OpenBUGS)
library(coda)
library(lattice)
library(stats)
mixexp <- function() {
for (i in 1:m){
gr[i] ~dcat(p[])
w[i] ~ dexp(lambda[gr[i]])
}
p[1:2]~ddirich(alpha[])
for(k in 1:2){
lambda[k] ~ dgamma(tau[k],psi[k])
}
}
data_Bayes <- list(w=Data3, m=length(Data3),tau=c(0.01,0.01), psi=c(0.01,0.01), alpha=c(0.5,0.5))
inits <- function(){
  list(lambda=c(0.1,0.1))}
bayesout_mix<- bugs(data_Bayes, inits, model.file = mixexp,
               parameters = c("lambda", "p"), n.chains = 1, n.burnin=10000, n.iter = 40000,DIC=FALSE, codaPkg=TRUE,digits=3,debug=F)
out_mixexp <- read.bugs(bayesout_mix)
### As the output structure is a list, the follwoing function is used to covert 
### the output to a data frame, which is easier to work with. 
### The output is a data frame with 4 columns: lambda[1], lambda[2], p[1], p[2], where
### the two first are intensity parameters of exponential distributions, and the 
### last two are mixing weights of mixture model. 
out_mixexp = do.call(rbind.data.frame, out_mixexp)
### Define the density
dMix2Exp<-function(x,
                   w_exp1=NULL,
                   rate1 = NULL,
                   rate2 = NULL){
  w_exp2=1-w_exp1
  dexp1=dexp(x,rate = rate1)
  dexp2=dexp(x,rate = rate2)
  return(w_exp1*dexp1+w_exp2*dexp2)
}
Expect_f= numeric(length(x_grid))
for(i in 1:length(x_grid)){
Expect_f[i]= mean(dMix2Exp(x_grid[i],out_mixexp[,3], out_mixexp[,1], out_mixexp[,2]) )
}

obeject_func= function(param, point){

w_exp1 = param[1]

rate1 = param[2]

rate2 = param[3]


res = sum(Expect_f*DeltaX*log(dMix2Exp( (1-point)/point, w_exp1, rate1, rate2)))

return(-res)


}

optim(par=c(0.843, 0.00374, 0.000755), fn=obeject_func, method="L-BFGS-B", lower=c(0, 0, 0), upper=c(1, +Inf, +Inf), point=m  )


0.85266116 1.06341086 0.06296207





