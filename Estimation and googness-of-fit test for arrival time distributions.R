### @ The main goal here is to fit different distributions to arrival 
### @ times observed for earthquake data set at hand. 

### @ Estimation for parameters of Exponential, Gamma, Weibull, and
### @ a mixutre of two exponential distributions by using MLE method.

Data3 = as.numeric(diff( unique(Data2$Date)))

library(actuar)

library(fitdistrplus)

fw <- fitdist(Data3, "weibull")

fexp <- fitdist(Data3, "exp")

fg <- fitdist(Data3, "gamma", start = list(shape = 1, rate = 1), lower = c(0, 0))

### @ we talk about differnt ways to find the estimation of 
### @ a mixture of two exponential distributions 

############################################################## 
##                                                          ##
##       Estimation of a mixture of two exponential         ##
##                      distributions                       ##
##                                                          ##
##############################################################

### Arrival times real data associated with the catastrophe
### events observed during years 1906 - 2006 

Data3 = c(1290, 2062, 1, 1033,  173,  620, 1832,  1,  551, 2261,  952, 12,  259,  556,
847,  408,  136, 1026,  528,  45, 1108,  219,  635,  341,  32,  548,  303,  35,
223,  49,  502,   13,  879,  595,  25,  488,  576,  397,  269,  483,  865,  495,
743 ,  64,  701,  126,  120,  988,  358, 70,  101,  115,  7,   63,  104,  169,
736,   71,  108,   19,  160,  642,  163,  5,  8,  437,  54,  580,  114,  133,
485,  302,   64,  270,  180,  118,  343, 2078,  178,  613,  414,  217,   64)


### ------------/|/ Estimation via MLE method /|/----------###

# Define the quantile
qMix2Exp<- function(p,
                   w_exp1=NULL,
                   rate1 = NULL,
                   rate2 = NULL){
  w_exp2=1-w_exp1
  qexp1=qexp(p,rate = rate1)
  qexp2=qexp(p,rate = rate2)
  return(w_exp1*qexp1+w_exp2*qexp2)
}

# Define the distribution function
pMix2Exp<-function(q,
                   w_exp1=NULL,
                   rate1 = NULL,
                   rate2 = NULL){
  w_exp2=1-w_exp1
  pexp1=pexp(q,rate = rate1)
  pexp2=pexp(q,rate = rate2)
  return(w_exp1*pexp1+w_exp2*pexp2)
}

# Define the density
dMix2Exp<-function(x,
                   w_exp1=NULL,
                   rate1 = NULL,
                   rate2 = NULL){
  w_exp2=1-w_exp1
  dexp1=dexp(x,rate = rate1)
  dexp2=dexp(x,rate = rate2)
  return(w_exp1*dexp1+w_exp2*dexp2)
}

library(fitdistrplus)


fit_A <- fitdist(Data3, "Mix2Exp", start = list(w_exp1=0.87, rate1=0.002052493, rate2=0.02382071), lower=c(0,0,0), upper=c(1,2000,2000))

### ------------/|/ Estimation via vglm function /|/----------###


library(VGAM)

Data4 = data.frame(Y = Data3)

fit <- vglm(Y ~ 1, mix2exp, data = Data4, trace = TRUE)

Coef(fit)


###############################################################
##                                                           ##
##      Expectation-Maximization method for mixture of       ##
##             of two exponential distributions              ##
##                                                           ##
###############################################################



##------------\|\ k-means initialization technique \|\---------##


library("ggplot2")
library("dplyr")
library("reshape2")

Data3.kmeans <- kmeans(Data3, 2, nstart = 25)
Data3.kmeans.cluster <- Data3.kmeans$cluster

Data3.df <- data_frame(x = Data3, cluster = Data3.kmeans.cluster)

Data3.df %>%
  mutate(num = row_number()) %>%
  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  geom_point() +
  ylab("Values") +
  ylab("Data Point Number") +
  scale_color_discrete(name = "Cluster") +
  ggtitle("K-means Clustering")

##------------\|\ Initial value for lambda \|\---------------##

Data3.summary.df <- Data3.df %>%
  group_by(cluster) %>%
  summarize(lambda = n()/sum(x), size = n())

### Initial value for mixing weight 

Data3.summary.df <- Data3.summary.df %>%
  mutate(alpha = size / sum(size))

##--------\|\ Expectation Step of the EM Algorithm \|\-------##


e_step <- function(x, lambda.vector, alpha.vector){
  comp1.prod <- dexp(x, rate =lambda.vector[1])*alpha.vector[1]
  comp2.prod <- dexp(x, rate = lambda.vector[2])*alpha.vector[2]
  sum.of.comps <- comp1.prod + comp2.prod
  comp1.post <- comp1.prod/sum.of.comps
  comp2.post <- comp2.prod/sum.of.comps

  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)

  list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = cbind(comp1.post, comp2.post))
}

##--------\|\ Maximization Step of the EM Algorithm \|\-------##

m_step <- function(x, posterior.df){
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])

  comp1.lambda <- comp1.n*(1/sum(posterior.df[, 1]*x) )
  comp2.lambda <- comp2.n*(1/sum(posterior.df[, 2]*x) )

  comp1.alpha <- comp1.n/length(x)
  comp2.alpha <- comp2.n/length(x)

  list("lambda" = c(comp1.lambda, comp2.lambda),
       "alpha" = c(comp1.alpha, comp2.alpha))
}

##---------------\|\ Iteration process  \|\-------------------## 

for (i in 1:1000){
  if(i == 1){
    # Initialization
    e.step <- e_step(Data3, Data3.summary.df[["lambda"]], Data3.summary.df[["alpha"]])
    m.step <- m_step(Data3, e.step[["posterior.df"]])
    cur.loglik <- e.step[["loglik"]]
    loglik.vector <- e.step[["loglik"]]
  }else{
    # Repeat E and M steps till convergence
    e.step <- e_step(Data3, m.step[["lambda"]], m.step[["alpha"]])
    m.step <- m_step(Data3, e.step[["posterior.df"]])
    loglik.vector <- c(loglik.vector, e.step[["loglik"]])

    loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
    if(loglik.diff < 1e-6){
      break
    }else{
      cur.loglik <- e.step[["loglik"]]
    }
}

}

loglik.vector

m.step

plot_mix_comps <- function(x, lam, alpha) {
  alpha * dexp(x, rate=lam)
}

data.frame(x = Data3) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 50, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$lambda[1], m.step$alpha[1] ),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$lambda[2], m.step$alpha[2] ),
                colour = "blue", lwd = 1.5) +
  ylab("Density") +
  xlab("Values") +
  ggtitle("Final EMM Fit")


### Let's put all peices together 

library(ReIns)

library(Renext)

fit_A <- fitdist(Data3, "mixexp2", start = list(prob1= 0.1217421, rate1=0.002052493, rate2=0.02382071, delta = 0.02382071-0.002052493),lower=c(0,0,0, 0), upper=c(1,10000,10000, 10000) )

par(mfrow = c(2, 2))

plot.legend <- c("weibull", "exp", "gamma", "mix-exp")

denscomp(list(fw, fexp, fg, fit_A), legendtext = plot.legend)

qqcomp(list(fw, fexp, fg, fit_A ), legendtext = plot.legend)

cdfcomp(list(fw, fexp, fg, fit_A ), legendtext = plot.legend)

ppcomp(list(fw, fexp, fg, fit_A ), legendtext = plot.legend)


################################################################
##                                                            ##
## Non-parametric goodness-of-fit tests for exponential       ##
##           distribution uisng bootstrap method              ##
##                                                            ##
################################################################

N= 1000

M = 100

for(i in 1:N){

if(i==1){

lambda = fitdist(Data3, "exp" )$estimate

R1.ks <- try(gofstat(fitdist(Data3, "exp"), fitnames = c("exponential"))$ks, silent = TRUE)

R1.ad <- try( gofstat(fitdist(Data3, "exp"), fitnames = c("exponential"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(Data3, "exp"), fitnames = c("exponential"))$cvm, silent = TRUE)

R1.ks.vector = R1.ks 

R1.ad.vector = R1.ad 

R1.cvm.vector = R1.cvm

}else{
 
repeat{

y1.curent= rexp(M, rate = lambda )

R1.ks <- try(gofstat(fitdist(y1.curent, "exp", start = list(rate=length(y1.curent)/sum(y1.curent)) ), fitnames = c("exponential"))$ks, silent = TRUE)

R1.ad <- try( gofstat(fitdist(y1.curent, "exp", start = list(rate=length(y1.curent)/sum(y1.curent)) ), fitnames = c("exponential"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(y1.curent, "exp", start = list(rate=length(y1.curent)/sum(y1.curent)) ), fitnames = c("exponential"))$cvm, silent = TRUE)

if(!inherits(R1.ks, "try-error") && !inherits(R1.ad, "try-error") && !inherits(R1.cvm, "try-error") )break
}
R1.ks.vector = c(R1.ks.vector, R1.ks)
R1.ad.vector = c(R1.ad.vector,R1.ad )
R1.cvm.vector = c(R1.cvm.vector, R1.cvm)

lambda = fitdist(y1.curent, "exp", start = list(rate=length(y1.curent)/sum(y1.curent)) )$estimate


}

}


p.value.ks= length(which(R1.ks.vector[-1]>= R1.ks.vector[1] ))/(N-1)

p.value.ad= length(which(R1.ad.vector[-1]>= R1.ad.vector[1] ))/(N-1)

p.value.cvm= length(which(R1.cvm.vector[-1]>= R1.cvm.vector[1] ))/(N-1)

p.value.table = list(Name = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-von Mises"), sig.level = 0.05, p.values= c(p.value.ks, p.value.ad, p.value.cvm))

print(p.value.table)

################################################################
##                                                            ##
##      Non-parametric goodness-of-fit tests for Gamma        ##
##           distribution uisng bootstrap method              ##
##                                                            ##
################################################################

N = 1000

M = 100

for(i in 1:N){

if(i==1){

lambda = fitdist(Data3, "gamma", start = list(shape = 1, rate = 1), lower = c(0, 0))$estimate

R1.ks <- try(gofstat(fitdist(Data3, "gamma", lower=c(0, 0), start = list(shape = 1, rate = 1)), fitnames = c("Gamma"))$ks, silent = TRUE)

R1.ad <- try( gofstat(fitdist(Data3, "gamma", lower = c(0,0), start = list(shape = 1, rate = 1)), fitnames = c("Gamma"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(Data3, "gamma", lower = c(0, 0), start = list(shape = 1, rate = 1)), fitnames = c("Gamma"))$cvm, silent = TRUE)

R1.ks.vector = R1.ks 

R1.ad.vector = R1.ad 

R1.cvm.vector = R1.cvm

}else{
 
repeat{

y1.curent= rgamma(M, shape = lambda[1], rate = lambda[2] )

R1.ks <- try(gofstat(fitdist(y1.curent, "gamma", lower =c(0, 0), start = list(shape = 1, rate = 1) ), fitnames = c("Gamma"))$ks, silent = TRUE)

R1.ad <- try( gofstat(fitdist(y1.curent, "gamma", lower = c(0, 0), start = list(shape = 1, rate = 1) ), fitnames = c("Gamma"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(y1.curent, "gamma", lower= c(0, 0), start = list(shape = 1, rate = 1) ), fitnames = c("Gamma"))$cvm, silent = TRUE)

if(!inherits(R1.ks, "try-error") && !inherits(R1.ad, "try-error") && !inherits(R1.cvm, "try-error") )break
}
R1.ks.vector = c(R1.ks.vector, R1.ks)
R1.ad.vector = c(R1.ad.vector,R1.ad )
R1.cvm.vector = c(R1.cvm.vector, R1.cvm)

lambda = fitdist(y1.curent, "gamma", lower = c(0, 0), start = list(shape = 1, rate = 1) )$estimate


}

}


p.value.ks= length(which(R1.ks.vector[-1]>= R1.ks.vector[1] ))/(N-1)

p.value.ad= length(which(R1.ad.vector[-1]>= R1.ad.vector[1] ))/(N-1)

p.value.cvm= length(which(R1.cvm.vector[-1]>= R1.cvm.vector[1] ))/(N-1)

p.value.table = list(Name = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-von Mises"), sig.level = 0.05, p.values= c(p.value.ks, p.value.ad, p.value.cvm))

print(p.value.table)


################################################################
##                                                            ##
##      Non-parametric goodness-of-fit tests for Weibull      ##
##           distribution uisng bootstrap method              ##
##                                                            ##
################################################################

N = 1000

M = 100

for(i in 1:N){

if(i==1){

lambda = fitdist(Data3, "weibull", start = list(shape = 1, scale = 1), lower = c(0, 0))$estimate

R1.ks <- try(gofstat(fitdist(Data3, "weibull", lower=c(0, 0), start = list(shape = 1, scale = 1)), fitnames = c("Gamma"))$ks, silent = TRUE)

R1.ad <- try( gofstat(fitdist(Data3, "weibull", lower = c(0,0), start = list(shape = 1, scale = 1)), fitnames = c("Gamma"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(Data3, "weibull", lower = c(0, 0), start = list(shape = 1, scale = 1)), fitnames = c("Gamma"))$cvm, silent = TRUE)

R1.ks.vector = R1.ks 

R1.ad.vector = R1.ad 

R1.cvm.vector = R1.cvm

}else{
 
repeat{

y1.curent= rweibull(M, shape = lambda[1], scale = lambda[2] )

R1.ks <- try(gofstat(fitdist(y1.curent, "weibull", lower =c(0, 0), start = list(shape = 1, scale = 1) ), fitnames = c("Gamma"))$ks, silent = TRUE)

R1.ad <- try( gofstat(fitdist(y1.curent, "weibull", lower = c(0, 0), start = list(shape = 1, scale = 1) ), fitnames = c("Gamma"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(y1.curent, "weibull", lower= c(0, 0), start = list(shape = 1, scale = 1) ), fitnames = c("Gamma"))$cvm, silent = TRUE)

if(!inherits(R1.ks, "try-error") && !inherits(R1.ad, "try-error") && !inherits(R1.cvm, "try-error") )break
}
R1.ks.vector = c(R1.ks.vector, R1.ks)
R1.ad.vector = c(R1.ad.vector,R1.ad )
R1.cvm.vector = c(R1.cvm.vector, R1.cvm)

lambda = fitdist(y1.curent, "weibull", lower = c(0, 0), start = list(shape = 1, scale = 1) )$estimate

}

}

p.value.ks= length(which(R1.ks.vector[-1]>= R1.ks.vector[1] ))/(N-1)

p.value.ad= length(which(R1.ad.vector[-1]>= R1.ad.vector[1] ))/(N-1)

p.value.cvm= length(which(R1.cvm.vector[-1]>= R1.cvm.vector[1] ))/(N-1)

p.value.table = list(Name = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-von Mises"), sig.level = 0.05, p.values= c(p.value.ks, p.value.ad, p.value.cvm))

print(p.value.table)


################################################################
##                                                            ##
##    Non-parametric goodness-of-fit tests for a mixture of   ##
##    two exponential distributions uisng bootstrap method    ##
##                                                            ##
################################################################

library("ggplot2")
library("dplyr")
library("reshape2")

MED = function(data){

###### k-means initialization technique

Data3.kmeans <- kmeans(data, 2, nstart = 25)
Data3.kmeans.cluster <- Data3.kmeans$cluster

Data3.df <- data_frame(x = data, cluster = Data3.kmeans.cluster)

### Initial value for lambda 

Data3.summary.df <- Data3.df %>%
  group_by(cluster) %>%
  summarize(lambda = n()/sum(x), size = n())

### Initial value for mixing weight 

Data3.summary.df <- Data3.summary.df %>%
  mutate(alpha = size / sum(size))

# Expectation Step of the EM Algorithm

e_step <- function(x, lambda.vector, alpha.vector){
  comp1.prod <- dexp(x, rate =lambda.vector[1])*alpha.vector[1]
  comp2.prod <- dexp(x, rate = lambda.vector[2])*alpha.vector[2]
  sum.of.comps <- comp1.prod + comp2.prod
  comp1.post <- comp1.prod/sum.of.comps
  comp2.post <- comp2.prod/sum.of.comps

  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)

  list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = cbind(comp1.post, comp2.post))
}

# Maximization Step of the EM Algorithm

m_step <- function(x, posterior.df){
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])

  comp1.lambda <- comp1.n*(1/sum(posterior.df[, 1]*x) )
  comp2.lambda <- comp2.n*(1/sum(posterior.df[, 2]*x) )

  comp1.alpha <- comp1.n/length(x)
  comp2.alpha <- comp2.n/length(x)

  list("lambda" = c(comp1.lambda, comp2.lambda),
       "alpha" = c(comp1.alpha, comp2.alpha))
}

### Iteration process 

for (i in 1:1000){
  if(i == 1){
    # Initialization
    e.step <- e_step(data, Data3.summary.df[["lambda"]], Data3.summary.df[["alpha"]])
    m.step <- m_step(data, e.step[["posterior.df"]])
    cur.loglik <- e.step[["loglik"]]
    loglik.vector <- e.step[["loglik"]]
  }else{
    # Repeat E and M steps till convergence
    e.step <- e_step(data, m.step[["lambda"]], m.step[["alpha"]])
    m.step <- m_step(data, e.step[["posterior.df"]])
    loglik.vector <- c(loglik.vector, e.step[["loglik"]])

    loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
    if(loglik.diff < 1e-6){
      break
    }else{
      cur.loglik <- e.step[["loglik"]]
    }
}

}

Z = c(m.step[[2]][1], m.step[[1]][1],m.step[[1]][2])

return(Z)


}


N = 1000

M = 100

library(Renext)

for(i in 1:N){

if(i==1){

lambda = MED(Data3)

R1.ks <- try(gofstat(fitdist(Data3, "mixexp2", start = list(prob1= lambda[1], rate1=lambda[2], rate2=lambda[3], delta = abs(lambda[2] - lambda[3]) ) ), fitnames = c("Mixture"))$ks, silent = TRUE)

R1.ad <- try(gofstat(fitdist(Data3, "mixexp2", start = list(prob1= lambda[1], rate1=lambda[2], rate2=lambda[3], delta = abs(lambda[2] - lambda[3]) ) ), fitnames = c("Mixture"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(Data3, "mixexp2", start = list(prob1= lambda[1], rate1=lambda[2], rate2=lambda[3], delta = abs(lambda[2] - lambda[3]) ) ), fitnames = c("Mixture"))$cvm, silent = TRUE)

R1.ks.vector = R1.ks 

R1.ad.vector = R1.ad 

R1.cvm.vector = R1.cvm

}else{
 
repeat{

y1.curent= rmixexp2(M, prob1 = lambda[1], rate1 = lambda[2], rate2 = lambda[3]  )

R1.ks <- try(gofstat(fitdist(y1.curent, "mixexp2", start = list(prob1= lambda[1], rate1=lambda[2], rate2=lambda[3], delta = abs(lambda[2] - lambda[3]) ) ), fitnames = c("Mixture"))$ks, silent = TRUE)

R1.ad <- try(gofstat(fitdist(y1.curent, "mixexp2", start = list(prob1= lambda[1], rate1=lambda[2], rate2=lambda[3], delta = abs(lambda[2] - lambda[3]) ) ), fitnames = c("Mixture"))$ad, silent = TRUE)

R1.cvm <- try(gofstat(fitdist(y1.curent, "mixexp2", start = list(prob1= lambda[1], rate1=lambda[2], rate2=lambda[3], delta = abs(lambda[2] - lambda[3]) ) ), fitnames = c("Mixture"))$cvm, silent = TRUE)

if(!inherits(R1.ks, "try-error") && !inherits(R1.ad, "try-error") && !inherits(R1.cvm, "try-error") )break
}
R1.ks.vector = c(R1.ks.vector, R1.ks)
R1.ad.vector = c(R1.ad.vector,R1.ad )
R1.cvm.vector = c(R1.cvm.vector, R1.cvm)

lambda = MED(y1.curent)

}

}

p.value.ks= length(which(R1.ks.vector[-1]>= R1.ks.vector[1] ))/(N-1)

p.value.ad= length(which(R1.ad.vector[-1]>= R1.ad.vector[1] ))/(N-1)

p.value.cvm= length(which(R1.cvm.vector[-1]>= R1.cvm.vector[1] ))/(N-1)

p.value.table = list(Name = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-von Mises"), sig.level = 0.05, p.values= c(p.value.ks, p.value.ad, p.value.cvm))

print(p.value.table)

### To get the test statistics for all non-parametric test, we have 

R = gofstat(list(fw, fexp, fg, fit_A), fitnames = c("weibull", "exponential", "Gamma", "mix-exp"))

### Final outputs: 

### parameter estimations: 

Fitting of the distribution ' weibull ' by maximum likelihood 
Parameters:
         estimate  Std. Error
shape   0.8452841  0.07386646
scale 398.7077089 54.40132229

Fitting of the distribution ' exp ' by maximum likelihood 
Parameters:
        estimate   Std. Error
rate 0.002308249 0.0001862489

Fitting of the distribution ' gamma ' by maximum likelihood 
Parameters:
         estimate Std. Error
shape 0.759189325         NA
rate  0.001752204         NA


Fitting of the distribution ' mixexp2 ' by maximum likelihood 
Parameters:
         estimate Std. Error
prob1 0.878887045         NA
rate1 0.002052217         NA
rate2 0.023882659         NA

### estimated p-values by uisng bootstrap method for goodness-of-fit tests

Name
"Kolmogorov-Smirnov" "Anderson-Darling"   "Cramer-von Mises"  

Exponential: 
$p.values
[1] 0.13413413   0.08308308   0.16216216

Gamma:
$p.values
[1] 0.2652653  0.7387387   0.6756757

Weibull:
$p.values
[1] 0.1301301  0.5335335  0.4634635

Mix-exp:
$p.values
[1] 0.4244244 0.6196196 0.6556557

Goodness-of-fit statistics
                                weibull exponential      Gamma    mix-exp
Kolmogorov-Smirnov statistic 0.07703060  0.09375521 0.07021095 0.05936541
Cramer-von Mises statistic   0.05138475  0.14875820 0.04156689 0.03358249
Anderson-Darling statistic   0.31886075  1.17653905 0.25950797 0.26670488

Goodness-of-fit criteria
                                weibull exponential    Gamma  mix-exp
Akaike's Information Criterion 1173.815    1175.830 1173.201 1177.676
Bayesian Information Criterion 1178.653    1178.249 1178.039 1187.351
