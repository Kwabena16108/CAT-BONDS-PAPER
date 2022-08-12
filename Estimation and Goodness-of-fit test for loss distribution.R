
#############################################################
##                                                         ##
##           Fitting Burr type XII uisng MLE method        ##
##                                                         ##
#############################################################

dburr= function(x, shape1, shape2, scale){

( (shape1*shape2)/(scale) )*( ( x/scale )^(shape1 - 1) )*( 1/( ( 1+  (x/scale)^shape1 )^(shape2+1) ) ) 

}

library(actuar)

library(fitdistrplus)

library(BMT)

N = 1000    ###  Number of bootstrap sample 

M = 230     ### Ransom sample size 

Fit_B = mpsedist(Data2$Adjusted_Losses, distr="burr", start = list(shape1 = 1, shape2 = 1, scale = 1))

shape1.ini = as.numeric(Fit_B$estimate[1])

shape2.ini = as.numeric(Fit_B$estimate[2])

scale.ini = as.numeric(Fit_B$estimate[3])

for(i in 1:N){

if(i==1){

lambda = c(shape1.ini, shape2.ini, scale.ini)

Y = try(gofstat(fitdist(Data2$Adjusted_Losses, "burr", start = list(shape1 = shape1.ini, shape2= shape2.ini, scale = scale.ini )), fitnames = c("Burr")), silent = TRUE)

R1.ks <- as.numeric(Y$ks)
 
R1.ad <- as.numeric(Y$ad)

R1.cvm <- as.numeric(Y$cvm)

R1.ks.vector = R1.ks 
 
R1.ad.vector = R1.ad 
 
R1.cvm.vector = R1.cvm

}else{
 
repeat{

y1.curent= rburr(M, shape1 = lambda[1], shape2 = lambda[2], scale = lambda[3] )

Fit_B = mpsedist(y1.curent, distr="burr", start = list(shape1 = 1, shape2 = 1, scale = 1))

shape1.ini = as.numeric(Fit_B$estimate[1])

shape2.ini = as.numeric(Fit_B$estimate[2])

scale.ini = as.numeric(Fit_B$estimate[3])

E <- try(fitdist(y1.curent, "burr", start = list(shape1 = shape1.ini+ abs(rnorm(1)), shape2= shape2.ini+ abs(rnorm(1)), scale = scale.ini + abs(rburr(1, 1, 1, 1)) ) ), silent = TRUE)

if(!inherits(E, "try-error") )break

}

R1.ks <- as.numeric(gofstat(E, fitnames=c("Burr"))$ks)

R1.cvm <- as.numeric(gofstat(E, fitnames=c("Burr"))$cvm)

R1.ks.vector = c(R1.ks.vector, R1.ks)
 
R1.ad.vector = c(R1.ad.vector,R1.ad )

R1.cvm.vector = c(R1.cvm.vector, R1.cvm)

lambda<- as.numeric(E$estimate)

}

print(i)
 
}


p.value.ks= length(which(R1.ks.vector[-1]>= R1.ks.vector[1] ))/(N-1)

p.value.ad= length(which(R1.ad.vector[-1]>= R1.ad.vector[1] ))/(N-1)
 
p.value.cvm= length(which(R1.cvm.vector[-1]>= R1.cvm.vector[1] ))/(N-1)

p.value.table = list(Name = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-von Mises"), sig.level = 0.05, p.values= c(p.value.ks, p.value.ad, p.value.cvm))

print(p.value.table)


### Output - Do not run the following code: 

$Name
[1] "Kolmogorov-Smirnov" "Anderson-Darling"   "Cramer-von Mises"  

$sig.level
[1] 0.05

$p.values
[1] 0.9409409 0.9489489 0.9559560

#### 


Fit_B = BMT::mpsedist(Data2$Adjusted_Losses, distr="burr", start = list(shape1 = 1, shape2 = 1, scale = 1)) ### we use this command line in order to get initial values for parameters. 

shape1.ini = as.numeric(Fit_B$estimate[1])

shape2.ini = as.numeric(Fit_B$estimate[2])

scale.ini = as.numeric(Fit_B$estimate[3])

Fit_B = fitdist(Data2$Adjusted_Losses, "burr", start = list(shape1 = shape1.ini, shape2= shape2.ini, scale = scale.ini ))

shape1.ini = as.numeric(Fit_B$estimate[1])

shape2.ini = as.numeric(Fit_B$estimate[2])

scale.ini = as.numeric(Fit_B$estimate[3])

ks.test(Data2$Adjusted_Losses, "pburr", shape1 =  shape1.ini , shape2 = shape2.ini  , scale = scale.ini  )

### Output - do not run the following 

       Asymptotic one-sample Kolmogorov-Smirnov test

data:  Data2$Adjusted_Losses
D = 0.080112, p-value = 0.1085
alternative hypothesis: two-sided


Parameters : 
           estimate   Std. Error
shape1 8.199442e-01 5.654303e-02
shape2 6.998999e-01 5.678553e-02
scale  2.847958e+07 2.097152e+03

### 

gofstat(Fit_B, fitnames = c("Burr"))

Goodness-of-fit statistics
                                   Burr
Kolmogorov-Smirnov statistic 0.08011238
Cramer-von Mises statistic   0.14938162
Anderson-Darling statistic   1.76218165

Goodness-of-fit criteria
                                   Burr
Akaike's Information Criterion 9217.523
Bayesian Information Criterion 9227.798



#############################################################
##                                                         ##
## Fitting Generalized Pareto distribution uisng MLE method##
##                                                         ##
#############################################################

#### Model estimation: Here, we provide different ways to fit a GPD in R. 


### First method to fit a Generalized Pareto distribution to the data. 

library(remotes)

remotes::install_url("https://cran.r-project.org/src/contrib/Archive/gPdtest/gPdtest_0.4.tar.gz")

library(gPdtest)

### "amle" method stands for the asympthotic maximum likelihood estimation method, if shape parameter is positive, otherewise "combined" method is used.  

gPdtest::gpd.fit(Data2$Adjusted_Losses, method = "amle" )

Parameter   estimate
shape       1.845493e+00
scale       5.100028e+07


##############################################################################

### Maximum-likelihood Fitting for the GPD Model

library(ismev)

the second argument is the thereshold level, which in our case is set to be zero.

ismev::gpd.fit(Data2$Adjusted_Losses , 0) 

$threshold
[1] 0

$nexc
[1] 227

$conv
[1] 0

$nllh
[1] 4601.922

$mle
[1] 2.595085e+07 2.201851e+00

### The first parameter is scale and the second is the shape parameter. 


########################################################################

library(evir)

### the second argument is the thereshold level. 

evir::gpd(Data2$Adjusted_Losses, 0)


$threshold
[1] 0

$p.less.thresh
[1] 0

$n.exceed
[1] 227

$method
[1] "ml"

$par.ests
          xi         beta 
2.200966e+00 2.588443e+07 

$par.ses
          xi         beta 
   0.1950104 2097.1520816 

$varcov
            [,1]          [,2]
[1,]  0.03802904 -1.140871e-01
[2,] -0.11408712  4.398047e+06

$information
[1] "observed"

$converged
[1] 0

$nllh.final
[1] 4601.922

attr(,"class")
[1] "gpd"

### Here, xi referes to the shape parameter and beta referes to the scale parameter. 

######################################################################################

## Fitting a GPD to Peaks Over a Threshold

library(POT)

POT::fitgpd(Data2$Adjusted_Losses, 0.0001, est = "mle")

Estimator: MLE 
 Deviance: 10016.35 
      AIC: 10020.35 

Varying Threshold: FALSE 

  Threshold Call: 1e-04 
    Number Above: 227 
Proportion Above: 1 

Estimates
    scale      shape  
2.309e+09  8.983e-01  

Standard Error Type: observed 

Standard Errors
    scale      shape  
1210.7913     0.1934  

Asymptotic Variance Covariance
       scale       shape     
scale   1.466e+06  -4.989e-02
shape  -4.989e-02   3.742e-02

Optimization Information
  Convergence: successful 
  Function Evaluations: 20 
  Gradient Evaluations: 10 

###############################################################

library(fitdistrplus) 

library(evir)

Fit_P= fitdist(Data2$Adjusted_Losses, "gpd", start=list(xi = 1, beta =1))

Fitting of the distribution ' gpd ' by maximum likelihood 
Parameters:
         estimate Std. Error
xi   2.200803e+00         NA
beta 2.589334e+07         NA


gofstat(Fit_P, fitnames = c("GPD"))
Goodness-of-fit statistics
                                    GPD
Kolmogorov-Smirnov statistic 0.04038288
Cramer-von Mises statistic   0.04399341
Anderson-Darling statistic   0.56223810

Goodness-of-fit criteria
                                    GPD
Akaike's Information Criterion 9207.845
Bayesian Information Criterion 9214.695

###############################################################

N = 1000

M = 230

for(i in 1:N){

if(i==1){

lambda = c(1, 1)

Y = try(gofstat(fitdist(Data2$Adjusted_Losses, "gpd", start = list(xi = lambda[1], beta= lambda[2] )), fitnames = c("GPD")), silent = TRUE)

R1.ks <- as.numeric(Y$ks)
 
R1.ad <- as.numeric(Y$ad)

R1.cvm <- as.numeric(Y$cvm)

R1.ks.vector = R1.ks 
 
R1.ad.vector = R1.ad 
 
R1.cvm.vector = R1.cvm

}else{
 
repeat{

y1.curent= rgpd(M, xi = lambda[1], beta = lambda[2] )

E <- try(fitdist(y1.curent, "gpd", start = list(xi = lambda[1], beta= lambda[2] ) ), silent = TRUE)

if(!inherits(E, "try-error") )break

}

R1.ks <- as.numeric(gofstat(E, fitnames=c("GPD"))$ks)

R1.ad <- as.numeric(gofstat(E, fitnames=c("GPD"))$ad)

R1.cvm <- as.numeric(gofstat(E, fitnames=c("GPD"))$cvm)

R1.ks.vector = c(R1.ks.vector, R1.ks)
 
R1.ad.vector = c(R1.ad.vector,R1.ad )

R1.cvm.vector = c(R1.cvm.vector, R1.cvm)

lambda<- as.numeric(E$estimate)

}

print(i)
 
}


p.value.ks= length(which(R1.ks.vector[-1]>= R1.ks.vector[1] ))/(N-1)

p.value.ad= length(which(R1.ad.vector[-1]>= R1.ad.vector[1] ))/(N-1)
 
p.value.cvm= length(which(R1.cvm.vector[-1]>= R1.cvm.vector[1] ))/(N-1)

p.value.table = list(Name = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-von Mises"), sig.level = 0.05, p.values= c(p.value.ks, p.value.ad, p.value.cvm))

print(p.value.table)


$Name
[1] "Kolmogorov-Smirnov" "Anderson-Darling"   "Cramer-von Mises"  

$sig.level
[1] 0.05

$p.values
[1] 0.5985986 0.2992993 0.7077077


#######################################################################
##                                                                   ##
## Fitting Generalized Extreme Value distribution uisng MLE method   ##
##                                                                   ##
#######################################################################


library(ismev)

ismev::gev.fit(Data2$Adjusted_Losses)
$conv
[1] 0

$nllh
[1] 4597.11

$mle
[1] 1.894765e+07 3.943892e+07 2.084837e+00

$se
[1]        NaN        NaN 0.01743185

Warning message:
In sqrt(diag(z$cov)) : NaNs produced

#### the MLE's for the location, scale and shape parameters, resp.

#################################################################

library(fExtremes)

gevFit(Data2$Adjusted_Losses, type = "mle")

Title:
 GEV Parameter Estimation 

Call:
 gevFit(x = Data2$Adjusted_Losses, type = "mle")

Estimation Type:
  gev mle 

Estimated Parameters:
          xi           mu         beta 
2.579438e+00 1.050361e+08 2.721944e+08 

Description
  Tue Aug  9 12:39:59 2022 

##################################################################

library(evir)

library(fitdistrplus)

N = 1000

M = 230

Fit_B = fExtremes::gevFit(Data2$Adjusted_Losses, type = "mle")

coef=slot(Fit_B, 'fit')$par.ests

shape.ini = as.numeric(coef['xi'])

loc.ini = as.numeric(coef['mu'])

scale.ini = as.numeric(coef['beta'])

for(i in 1:N){

if(i==1){

lambda = c(shape.ini, loc.ini, scale.ini)

Y = try(gofstat(fitdist(Data2$Adjusted_Losses, "gev", start = list(xi = lambda[1], mu= lambda[2], sigma= lambda[3] )), fitnames = c("GEV")), silent = TRUE)

R1.ks <- as.numeric(Y$ks)
 
R1.ad <- as.numeric(Y$ad)

R1.cvm <- as.numeric(Y$cvm)

R1.ks.vector = R1.ks 
 
R1.ad.vector = R1.ad 
 
R1.cvm.vector = R1.cvm

}else{
 
repeat{

y1.curent= rgev(M, xi = lambda[1], mu = lambda[2], sigma= lambda[3] )

E <- try(fitdist(y1.curent, "gev", start = list(xi = lambda[1], mu= lambda[2], sigma= lambda[3] ) ), silent = TRUE)

if(!inherits(E, "try-error") )break

}

R1.ks <- as.numeric(gofstat(E, fitnames=c("GEV"))$ks)

R1.ad <- as.numeric(gofstat(E, fitnames=c("GEV"))$ad)

R1.cvm <- as.numeric(gofstat(E, fitnames=c("GEV"))$cvm)

R1.ks.vector = c(R1.ks.vector, R1.ks)
 
R1.ad.vector = c(R1.ad.vector,R1.ad )

R1.cvm.vector = c(R1.cvm.vector, R1.cvm)

lambda<- as.numeric(E$estimate)

}

print(i)
 
}

p.value.ks= length(which(R1.ks.vector[-1]>= R1.ks.vector[1] ))/(N-1)

p.value.ad= length(which(R1.ad.vector[-1]>= R1.ad.vector[1] ))/(N-1)
 
p.value.cvm= length(which(R1.cvm.vector[-1]>= R1.cvm.vector[1] ))/(N-1)

p.value.table = list(Name = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-von Mises"), sig.level = 0.05, p.values= c(p.value.ks, p.value.ad, p.value.cvm))

print(p.value.table)


$Name
[1] "Kolmogorov-Smirnov" "Anderson-Darling"   "Cramer-von Mises"  

$sig.level
[1] 0.05

$p.values
[1] 0 0 0


Fit_GEV = fitdist( Data2$Adjusted_Losses, "gev", start = list(xi = shape.ini, mu= loc.ini, sigma= scale.ini ) )


Fitting of the distribution ' gev ' by maximum likelihood 
Parameters:
          estimate   Std. Error
xi    2.581427e+00  0.009519565
mu    1.050361e+08  4.693145410
sigma 2.723988e+08 12.299277627


gofstat(Fit_GEV, fitnames =c("GEV"))

Goodness-of-fit statistics
                                    GEV
Kolmogorov-Smirnov statistic  0.2974481
Cramer-von Mises statistic    9.2039662
Anderson-Darling statistic   51.7864068

Goodness-of-fit criteria
                                    GEV
Akaike's Information Criterion 9310.124
Bayesian Information Criterion 9320.398



R = gofstat(list(Fit_B, Fit_P, Fit_GEV), fitnames = c("Burr XII", "GPD", "GEV"))

Goodness-of-fit statistics
                                   Burr        GPD        GEV
Kolmogorov-Smirnov statistic 0.08011238 0.04038288  0.2974481
Cramer-von Mises statistic   0.14938164 0.04399341  9.2039662
Anderson-Darling statistic   1.76218179 0.56223810 51.7864068

Goodness-of-fit criteria
                                   Burr      GPD      GEV
Akaike's Information Criterion 9217.523 9207.845 9310.124
Bayesian Information Criterion 9227.798 9214.695 9320.398


par(mfrow = c(2, 2))
plot.legend <- c("Burr XII", "GPD", "GEV")

denscomp(list(Fit_B, Fit_P, Fit_GEV), legendtext = plot.legend)
qqcomp(list(Fit_B, Fit_P, Fit_GEV ), legendtext = plot.legend)
cdfcomp(list(Fit_B, Fit_P, Fit_GEV ), legendtext = plot.legend)
ppcomp(list(Fit_B, Fit_P, Fit_GEV ), legendtext = plot.legend)

library(dplyr)

Y=Data2%>%group_by(year)%>%summarize(freq = n(), aggr = sum(Adjusted_Losses))

plot(Y$year, Y$aggr, type = "o", xlab="Time (years)", ylab= "Adjusted earthquake losses (in 2020 Dollars)", tck = 0.01, log = "y", xaxt = "n", yaxt = "n")

axis(1, at = seq(1906, 2007, by= 10), tck = 0.01)

marks = c(0, 1e+6, 1e+7, 1e+8, 1e+9, 1e+10, 1e+11)

# formatting function
# https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific <- function(x){
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

scientific_10 <- function(x) {  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }

axis(2, at = marks, tck = 0.01, labels=scientific(marks), las= 1, cex.axis=0.75, mgp = c(1, 0.5, 0))

plot(Y$year, Y$freq, type = "o", xlab="Time (years)", ylab= "Number of earthquake events per year", tck = 0.01, xaxt = "n", las = 1)

axis(1, at = seq(1906, 2007, by= 10), tck = 0.01)


