### CAT bond Price under the assumption that the inter-arrival time 
### is exponentialy distributed. Moreover, the beta function is assumed to be 
### the one leading to the Esscher principle formula: 

library(actuar)

library(fitdistrplus)

tau <- seq(from=0.0, to=2.0, length=100) # 2 years horizon

D <- seq(5.75e+8, 6.93e+9, length.out = 100)

myGrid<- data.frame(expand.grid(tau, D))

colnames(myGrid)=list("tau" , "D")

T = 2

#### Zero-Coupon bond price under Risk-Neural measure Q

P<-function(tau, r, theta,m, sigma){
  
  gam<- sqrt((theta)^2+2*sigma^2)
  
  A<-((2*gam*exp((theta+gam)*(tau)/2))/((theta+gam)*(exp(gam*(tau))-1)+2*gam))^((2*m*theta)/(sigma^2)) 
  
  B<-(2*(exp(gam*(tau)))-1)/((theta+gam)*(exp(gam*(tau))-1)+2*gam)
  
  Price<- A*exp(-B*r)
  
  return(Price)
  
}

alpha= seq(-1, 1, by = (2)/20 )

r0<- log1p(0.011); theta<- 0.254; m<- 0.011; sigma<- 0.074; MP = -0.033
  
theta_star<- theta + MP; m_star<- (theta*m)/(theta +MP)

N=length(tau); dt <- (tau[length(tau)]- tau[1])/length(tau)

nsim= 100

X <- matrix(rnorm(N*nsim, mean=0,sd = sqrt(dt)), N, nsim)

R <- matrix(0,N,nsim);R[1, ]= r0

for(j in 1:nsim){for (i in 2:N){R[i,j]= max(0,R[i-1,j]+theta_star*(m_star-R[i-1,j])*dt+sigma*sqrt(R[i-1,j])*X[i,j])}}

R_new1 = R[dim(R)[1]:1, ]

R_new2 = c(R_new1)


V_CATbond1<-function(tau, r , T, D, alpha){ 

  theta<- 0.254
  
  m<- 0.011
  
  sigma<- 0.074

  MP = -0.033

  Z= 1

  q = 0.5

  lambda = 0.002 # intensity rate of event
  
  theta_star<- theta + MP
  
  m_star<- (theta*m)/(theta +MP)

return( P(tau, r, theta_star, m_star, sigma)*(Z-(Z-q*Z)*(1 - exp(-lambda*exp(alpha)*T*(1-pburr(D, shape1= 8.199442e-01 , shape2= 6.998999e-01, scale=2.847958e+07))  ) )  ) )

}

Data_new = matrix(numeric(nrow(myGrid)*length(alpha)), nrow = nrow(myGrid), ncol = length(alpha))

for(i in 1:nrow(myGrid)){
 for(j in 1:length(alpha)){

Data_new[i,j] = V_CATbond1(myGrid$tau[i] , R_new2[i], T, myGrid$D[i], alpha[j])

}

}

plot(alpha, Data_new[1, ], type = "o" )

# S1 = numeric(dim(Data5)[1])

# for(i in 1:length(S1)){

# Y= splinefun(Data_new[i, ], alpha, method = "hyman")

# S1[i] = Y(1) 

# }

# alpha_new1 = min(S1)

Price1 = V_CATbond1(myGrid$tau , R_new2, T, myGrid$D, alpha= 0.5)

library(lattice)

wireframe( Price1~ myGrid$tau*myGrid$D , scales = list(arrows=FALSE),
           xlab= list("T(Maturity time, Year)", rot=-35), ylab= list("D(Threshold level)", rot= 30) , zlab= list("V(Price, $)",
                                                                                                                               rot=90), drap = TRUE , colorkey = TRUE , screen=list(z=-45,x=-60))
#########