# Data preparation-------------
path = "../rebuttal" # change path
library(lubridate)
library(tidyr)
# install.packages("nloptr", dependencies = TRUE)
library(nloptr)

interest = read.csv(file = paste0(path,"/modelling/Canada Interest Rate.csv"), header = T) # 1year, 3,6 months interest rate


###########################################################################
####   Finding initial values for CIR parameters by using OLS method   ####
###########################################################################

N = length(interest$X3MON)

#### The number of trading days is set to be 250 days per year. 

dt = 1/250

rate = interest$X3MON[1:(N-1)]

lagrate = interest$X3MON[2:N]

y = (lagrate - rate)/sqrt(dt*rate)

z1 = sqrt(dt/rate)

z2= sqrt(dt*rate)

data1 = data.frame(y = y, z1 = z1, z2 = z2)

fit = lm(y~-1+ z1 + z2, data = data1)

summary(fit)

# from Safarvesi, Domfeh & Chatterjee (2022) eqn 6.4-6.5, we find the relationship
# between the OLS estimates and parameters of the CIR model as follows:
# \beta_1 = \theta * m , \beta_2 = -theta. therefore,
# these are used as initial values for the optimization
theta_hat = -coef(fit)[2] # 0.9448066 
m_hat = coef(fit)[1] / theta_hat # 0.7087035 

Z = z1 * z2
#  beta_hat = ((X'X)^(-1))X'y
beta_hat = solve(t(Z)%*%Z)%*%t(Z)%*%y
# find sigma
euclidean <- function(a, b) sqrt(sum((a-b)^2))
sigma = (1 / (N-2)) * (euclidean(y, Z*beta_hat))^2 # 1.546442



# Finding the CIR model parameters via MLE -------------

CIRloglike = function(param , data , times , test=F , addsign =T ){
  
  # CIRloglike likelihood Function
  # param...parameters of the CIR model
  # dt ... time interval between the data points 
  # c ... multiplying term for the chi-square distribution 
  # df ... degree of freedom 
  # ncp ... non-centrality parameter 
  
  theta<- param[1]
  m<- param[2]
  sigma<- param[3]
  N<- length(data)
  
  if(test== T)
    dt= times
  else 
    dt = diff(times, 1)/250
  
  dt = as.numeric(dt)
  rate = data[1:(N-1)]
  lagrate = data[2:N]
  
  ncp = rate*(  (  4*theta*exp(-theta*dt)  )/( (sigma^2)*(1- exp(-theta*dt) )  )    )
  d = (4*theta*m)/(sigma^2)
  c = (  4*theta  )/( (sigma^2)*(1- exp(-theta*dt) )  ) 
  res = sum( dchisq(c*lagrate, df = d, ncp = ncp, log= TRUE) + log(c)  )
  
  if(addsign)
    return(-res)
  else 
    return(res)
}

MLE_CIR = optim(par=c(0.9448066, 0.7087035, 1.546442), fn=CIRloglike, method="L-BFGS-B",
                lower=c(0.0001, 0.0001, 0.0001), upper=c(+Inf, +Inf, +Inf),
                data=interest$X3MON, times= 1/250, test=T  )

MLE_CIR$par
# results [1] 0.9552238 0.7100216 1.1974436





# 2. Estimate with MLE  with constraints -----------------------
################################################################################
## Non-linear constraint optimization in R: We can do the above optimization #
## by considering contraints on parameters of the CIR model, i.e., theta > 0,  #
## m >0, sigma > 0, and 2*m*theta > sigma                                      #
################################################################################

## Objective function 
CIRloglike =function(param , data , times ){
  # CIRlolike log-likelihood Function to be maximized with respect to parameters
  # param...parameters of the CIR model
  # dt ... time interval between the data points 
  # c ... multiplying term for the chi-square distribution 
  # df ... degree of freedom 
  # ncp ... non-centrality parameter 
  
  theta<- param[1]
  m<- param[2]
  sigma<- param[3]
  N<- length(data)
  dt= times
  
  rate = data[1:(N-1)]/100
  lagrate = data[2:N]/100 
  
  ncp = rate*(  (  4*theta*exp(-theta*dt)  )/( (sigma^2)*(1- exp(-theta*dt) )  )    )
  d = (4*theta*m)/(sigma^2)
  c = (  4*theta  )/( (sigma^2)*(1- exp(-theta*dt) )  ) 
  res = sum( dchisq(c*lagrate, df = d, ncp = ncp, log= TRUE) + log(c)  )
  
  return(-res)
  
}

### Constraint function 
# Equality constraints

eval_g_ineq <- function(param)
{
  return ( (param[3])^2 - 2*param[1]*param[2] )
}

## Lower and upper bounds
lb <- c(0.001,0.001,0.001)
ub <- c(Inf,Inf,Inf)

## Initial values
x0 <- c(0.9448066, 0.7087035, 1.5464420)

## Set optimization options
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )

opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 160000,
              "local_opts" = local_opts,
              "print_level" = 0 )

res <- nloptr ( x0 = x0,
                eval_f = CIRloglike,
                lb = lb,
                ub = ub,
                eval_g_ineq = eval_g_ineq,
                opts = opts,
                data= interest$V80691344/100,
                times = 1/250
)

print(res)




###########################################################################
####        Computing the market price of interest rate risk          #####
###########################################################################

interest$Slope = ( (interest$X1YR) - (interest$X6MON) )/(180/250)

#### The MLE reults for parameters of the CIR model 
theta = 0.9552238
m =  0.7100216
sigma = 1.1974436


lambda_t = ((theta*m - 2*interest$Slope)/(interest$X3MON)) - theta # market price of risk
lambda_r = mean(lambda_t) #
# lambda_r : -0.5547602 # average market price of risk
# plot(interest$date, xlab = "", ylab = expression(lambda_t), lambda_t, type = "l")

# plot(interest$date, xlab = "", ylab = expression(lambda_t), lambda_t, type = "l")


###########################################################################
####        Computing the CIR risk-neutral Bond Yield         #####
###########################################################################


claim_amt = c(93980862,203573187,82438511,58207494,98789767,203784637,84601207,48266955,112214807,405363680,30273264,180974006,
              37739274,651461972,53756681,143107731,36090613,84401541,588539114,36781837,130601067,75331736,286491004,83436015,
              190295693,109419766,188648779,37859589,239306986,77273365,778381182,60095502,34004513,85415402,163555540,1719550417,
              992335953,134950418,28929260,248717795,40185538,216417898,132177018,573134416,38945825,53035617,158020158,83606360,
              273221882,132571995,218272961,35560468,3884408659,163731979,79632374,111151588,132945745,835223966,29211447,140578276,
              119570081,174425978,46653607,147798727,163938682,101706921,177644518,69793493,106730574,60693231,164691525,142759633,
              121663769,100681909,353522100,775291475,82083862,53310706,251714530,327836544,295925561,100576311,163783790,153859331,
              88008553,249893480,141914000,522868000,1306405000,219446000,17201000)/10e06 # claim amount size  (this is X_1, X_2,..X_n)

ann_to_inst <- function(r) {
  return (log1p(r))
}

inst_to_ann <- function(r) {
  return(expm1(r))
}


CirPriceYield <- function(r, tau, Param, priceyn=F) {
  # r .... r(t) sim. current value of short rate
  # tau .... (T-t), time to maturity
  # Params... vector holding the parameter of the CIR model
  h = sqrt(theta_star^2 + 2 * sigma2)
  B = 2 * (exp(h * tau) - 1) / (2 * h + (theta_star + h) * (exp(tau * h) -1))
  A = ((2*h*exp((h+theta_star)*tau/2))/(2*h+(h+theta_star)*(exp(h*tau)-1)))**(2*theta_star*m_star/sigma2)
  if (priceyn){
    if(tau==0) return(1) #price is par-value(1) at maturity
    else return(A * exp(-B * r))
  }
  else return((r * B - log(A)) / tau)
}

# this is based on our paper (section of CIR calibration)
theta_ = 0.9552238
m = 0.7100216
sigma2 = sqrt(1.1974436)
lambda_r = -0.5547602
theta_star = theta_ + lambda_r
m_star = (theta_ * m) / (theta_ + lambda_r)


r = ann_to_inst(m) # note that we use the long-run mean (m) as the starting short interest rate (i converted it to instantaneous because it was annualized)
t <- seq(from=0.0, to=2.0, length=100) # 2 years horizon
D <- seq(from=quantile(claim_amt, .25), to=quantile(claim_amt, .99), length.out = 100)
myGrid <- data.frame(expand.grid(t,D))
tau = myGrid$Var1

priceCIR <- c()
for (j in 1:length(tau)){
  priceCIR[j] <- CirPriceYield(r, tau[j], c(theta_star, m_star, sigma2), T)
}

###########################################################################
####        Get probability of trigger      #####
###########################################################################

# path = '/Users/dicksonnkwantabisa/Desktop/DBA Finance/PostDoc/rebuttal/modelling/cat_sim_20221031/'

lambda = as.matrix(read_csv(paste0(path,"lambda_2years.csv"))[,-1]) # assume that you have the intensity rates already.
shape =  2.846795406360715; scale = 3445794823.795233 # this is from fitting invgamma to the data in Python

Time = 2.0
prob_vec = c()
# start computing the probability of trigger of the myGrid
for (scen in 1:nrow(myGrid)) {
  
  D = myGrid[scen, 2]
  t = myGrid[scen, 1]
  
  Nft = array(NA, dim = c(nrow(lambda),length(lambda), nrow(myGrid)))
  
  # get an array of Nit for each time t
  for (i in 1:ncol(Nft)) {
    for (j in 1:nrow(lambda)) {
      Nft[j,i, scen] = (rpois(n=1, lambda = lambda[j,i]*t))
    }
  }
  # get an array of Sit for each time t
  Sft = array(NA, dim = c(nrow(lambda),length(lambda), nrow(myGrid)))
  
  for (i in 1:ncol(Sft)) {
    for (j in 1:length(kappa)) {
      Sft[j,i, scen] = sum(rgamma(n=Nft[j,i, scen], shape = shape, scale = scale))
    }
  }
  prob_vec[scen] = mean(Sft[,, scen] <= D * 10e06)
}


prob_thrshld = matrix(prob_vec, nrow = 100, ncol = 100); View(prob_thrshld)



###########################################################################
####        Plot the bond price evolution        #####
###########################################################################

# path = '/Users/dicksonnkwantabisa/Desktop/DBA Finance/PostDoc/rebuttal/modelling/cat_sim_20221031/'

# read in the pobability of trigger vector (multi-threshold) simulated in python for a 2 year horizon
# prob_vec = read.csv(paste0(path,"prob_vec_invgamma_multi_thrshld_2years_grp1.csv"))[,-1] # multi-threshold
# pistar = read.csv(paste0(path,"gammapistar_mlti_trshld_2years.csv"))$X0
# Calculate the payoff & final prices for a Coupon bond

# Payoff :  K + c, if Lτ < D (τ ≤ T )
#             K,   if Lτ > D then principal is protected

# Compute bond prices on a grid (so you ca nsee that the cirYield and Prob_vec were all calcuulatedon the myGrid)
# therfore you can just use Grid$t * Grid$D to get you the x, y axis

K = 100 ; c = 10 # coupon payment

payoff = ((K + c) * prob_vec) + (K * (1 - prob_vec))
bond_price = priceCIR * payoff 
# bond_price = priceCIR * payoff

library(lattice)


t <- seq(from=0., to=2.0, length=100)
D <- seq(from=quantile(claim_amt, .25), to=quantile(claim_amt, .99), length.out = 100)
Grid <- data.frame(expand.grid(t,D))
colnames(Grid) <- c("t","D")
bond_price <- unlist(bond_price)

wireframe(bond_price ~ Grid$t * Grid$D,
          scales = list(arrows = FALSE),
          xlab = list("Time to Maturity (years)", rot=-35),
          ylab = list("Threshold Level($10million)", rot=30),
          zlab=list("Value($)",rot=90),
          drape = TRUE, colorkey = TRUE,
          screen = list(z = -45, x = -60))


## Newly added 12/14/2022 -----

claim_amt = c(93980862,203573187,82438511,58207494,98789767,203784637,84601207,48266955,112214807,405363680,30273264,180974006,
              37739274,651461972,53756681,143107731,36090613,84401541,588539114,36781837,130601067,75331736,286491004,83436015,
              190295693,109419766,188648779,37859589,239306986,77273365,778381182,60095502,34004513,85415402,163555540,1719550417,
              992335953,134950418,28929260,248717795,40185538,216417898,132177018,573134416,38945825,53035617,158020158,83606360,
              273221882,132571995,218272961,35560468,3884408659,163731979,79632374,111151588,132945745,835223966,29211447,140578276,
              119570081,174425978,46653607,147798727,163938682,101706921,177644518,69793493,106730574,60693231,164691525,142759633,
              121663769,100681909,353522100,775291475,82083862,53310706,251714530,327836544,295925561,100576311,163783790,153859331,
              88008553,249893480,141914000,522868000,1306405000,219446000,17201000)/10e06 # claim amount size  (this is X_1, X_2,..X_n)

path = '/Users/dicksonnkwantabisa/Desktop/DBA Finance/PostDoc/rebuttal/modelling/cat_sim_20221231/'

prob_vec = read.csv(paste0(path,"prob_vec_invgamma_multi_thrshld_2years_grp3.csv"))[,-1] # multi-threshold
cir_rates = read.csv(paste0(path, "interest_rate_cir_vec_physical.csv"))$X0

K = 100 ; c = 10 # coupon payment

payoff = ((K + c) * prob_vec) + (K * (1 - prob_vec))

ann_to_inst <- function(r) {
  return (log1p(r))
}

inst_to_ann <- function(r) {
  return(expm1(r))
}


CirPriceYield <- function(r, tau, Param, priceyn=F) {
  # r .... r(t) sim. current value of short rate
  # tau .... (T-t), time to maturity
  # Params... vector holding the parameter of the CIR model
  h = sqrt(theta_star^2 + 2 * sigma2)
  B = 2 * (exp(h * tau) - 1) / (2 * h + (theta_star + h) * (exp(tau * h) -1))
  A = ((2*h*exp((h+theta_star)*tau/2))/(2*h+(h+theta_star)*(exp(h*tau)-1)))**(2*theta_star*m_star/sigma2)
  if (priceyn){
    if(tau==0) return(1) #price is par-value(1) at maturity
    else return(A * exp(-B * r))
  }
  else return((r * B - log(A)) / tau)
}

# this is based on our paper (section of CIR calibration)
theta_ = 0.9552238
m = 0.7100216
sigma2 = sqrt(1.1974436)
lambda_r = -0.5547602
theta_star = theta_ + lambda_r
m_star = (theta_ * m) / (theta_ + lambda_r)


MyCirPriceYield <- function(r, tau, Param, priceyn=F) {
  # r .... r(t) sim. current value of short rate
  # tau .... (T-t), time to maturity
  # Params... vector holding the parameter of the CIR model
  h = sqrt(a^2 + 2 * sigma2)
  B = 2 * (exp(h * tau) - 1) / (2 * h + (a + h) * (exp(tau * h) -1))
  A = ((2*h*exp((h+a)*tau/2))/(2*h+(h+a)*(exp(h*tau)-1)))**(2*a*theta_/sigma2)
  if (priceyn){
    if(tau==0) return(1) #price is par-value(1) at maturity
    else return(A * exp(-B * r))
  }
  else return((r * B - log(A)) / tau)
}

# from my model caliberation of Bayesian CIR 
a = 2.9734124627604723
b = 0.9237893398143259
sigma = sqrt(0.04012550059343504)

r = ann_to_inst(cir_rates) # note that  I am using the vector of interest rate here.
# r = ann_to_inst(m) # for your model
# r = ann_to_inst(a/b) # for my model
D <- seq(from=quantile(claim_amt, .25), to=quantile(claim_amt, .99), length.out = 100)
myGrid <- data.frame(expand.grid(t,D))
tau = myGrid$Var1


# choose one 

#either
priceCIR <- c()
for (j in 1:length(tau)){
  priceCIR[j] <- MyCirPriceYield(r[j], tau[j], c(theta_star, m_star, sigma2), T)
}

#or
priceCIR <- c()
for (j in 1:length(tau)){
  priceCIR[j] <- CirPriceYield(r[j], tau[j], c(theta_star, m_star, sigma2), T)
}


library(lattice)


colnames(myGrid) <- c("t","D")

bond_price = payoff * priceCIR

bond_price <- unlist(bond_price)

wireframe(bond_price ~ Grid$t * Grid$D,
          scales = list(arrows = FALSE),
          xlab = list("Time to Maturity (years)", rot=-35),
          ylab = list("Threshold Level($10million)", rot=30),
          zlab=list("Value($)",rot=90),
          drape = TRUE, colorkey = TRUE,
          screen = list(z = -45, x = -60))





