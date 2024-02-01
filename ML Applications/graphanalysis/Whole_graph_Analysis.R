### In R script
library(igraph)

graph <- read_graph("whole_graph.graphml", format = "graphml")

# Get the degree of each vertex
degrees <- degree(graph)

# Maximum degree
max_degree <- max(degrees)
cat("Maximum Degree:", max_degree, "\n")

# output: Maximum Degree: 753

# Minimum degree 

min_degree <- min(degrees)
cat("Minimum Degree:", min_degree, "\n")

# output: Minimum Degree: 2 

# Average degree
avg_degree <- mean(degrees)
cat("Average Degree:", avg_degree, "\n")

# output: Average Degree: 14.66667 


############### Fit power-law ###########################|

data <- degrees
data <- data[data>0]
g <- graph 

data.dist <- data.frame(k=0:max(data),p_k=degree_distribution(g))

library(ggplot2)
ggplot(data.dist) + geom_point(aes(x=k, y=p_k)) + theme_bw()+
labs(x = "Degree", y = "Degree Distribution")

#### Initial estimation

# This estimation calcualate Kmin, γ and makes CDF and determine minimum D. 

library(poweRlaw)
m_pl <- displ$new(data)
est_pl <- estimate_xmin(m_pl)

est_pl$xmin #k_min
# output:  6

est_pl$pars #gamma
# [1] 2.584408

est_pl$gof #D
# [1] 0.08087806

### Scanning whole range

data.s <- unique(data)

d_est <- data.frame(K_min=sort(data.s)[1:(length(data.s)-2)], gamma=rep(0,length(data.s)-2), D=rep(0,length(data.s)-2))

for (i in d_est$K_min){
  d_est[which(d_est$K_min == i),2] <- estimate_xmin(m_pl, xmins = i)$pars
  d_est[which(d_est$K_min == i),3] <- estimate_xmin(m_pl, xmins = i)$gof
}

K.min_D.min <- d_est[which.min(d_est$D), 1]

ggplot(data=d_est, aes(x=K_min, y=D)) + geom_line() + theme_bw() + 
  geom_vline(xintercept=K.min_D.min, colour="red") + annotate("text", x=K.min_D.min, y=max(d_est$D)/3*2, label=K.min_D.min)

ggplot(data=d_est, aes(x=K_min, y=gamma)) + geom_line() + theme_bw() + 
  geom_vline(xintercept=K.min_D.min, colour="red") + annotate("text", x=K.min_D.min, y=max(d_est$gamma)/3*2, label=K.min_D.min)

### And the fitted power-law on CDF curve.

m_pl$setXmin(est_pl)
plot.data <- plot(m_pl, draw = F)
fit.data <- lines(m_pl, draw = F)

ggplot(plot.data) + geom_point(aes(x=log(x), y=log(y))) + labs(x="log(k)", y="log(CDF)") + theme_bw() + 
  geom_line(data=fit.data, aes(x=log(x), y=log(y)), colour="darkred")  

## Investigate goodness of fit

bs_pl <- bootstrap_p(m_pl, no_of_sims=1000, threads=8, seed = 123)
#threads=core number of processor that used by function
#parallel::detectCores() determines how many cores in your computer

plot(bs_pl)

df_bs_pl <- bs_pl$bootstraps


ggplot(data=df_bs_pl, aes(pars)) + geom_histogram() + labs(x="gamma", y="frequency") + theme_bw()

ggplot(data=df_bs_pl, aes(xmin)) + geom_histogram() + labs(x="K_min", y="frequency") + theme_bw()

gamma_D.min <- d_est[which.min(d_est$D), 2]

ggplot(data=df_bs_pl, aes(x=xmin, y=pars)) + labs(x="K_min", y="gamma") + theme_bw() + 
  geom_point(shape=21, colour="black", fill="red", size=0.5, stroke=2, 
             position = position_jitter(), alpha=0.6) +
  geom_vline(xintercept=K.min_D.min, colour="blue") +
  geom_hline(yintercept=gamma_D.min, colour="blue") +
  annotate("text", x=K.min_D.min, y=min(df_bs_pl$pars), label=K.min_D.min, col="blue") +
  annotate("text", x=min(df_bs_pl$xmin), y=gamma_D.min, label=round(gamma_D.min, digits=2), col="blue")

D.min <- d_est[which.min(d_est$D), 3]

ggplot(data=df_bs_pl, aes(gof)) + geom_histogram() + labs(x="D", y="frequency") + geom_vline(xintercept=D.min, colour="red") + theme_bw()

bs_pl$p #p value
 0

### Fitting real distribution

#generate kmin & kmax pairs
pairs <- as.data.frame(t(combn(sort(data.s), 2)))
pairs$D <- rep(0, length(pairs$V1))
pairs$gamma <- rep(0, length(pairs$V1))

#scan D for all kmin-kmax pairs
for (i in 1:length(pairs$D)){
  m_pl$setXmin(pairs[i,1])
  pairs[i,3]<- estimate_xmin(m_pl, xmin = pairs[i,1], xmax = pairs[i,2], distance = "ks")$gof
  pairs[i,4]<- estimate_xmin(m_pl, xmin = pairs[i,1], xmax = pairs[i,2], distance = "ks")$pars
}

bs_pl_sat_cut <- bootstrap_p(m_pl, xmins = pairs[which.min(pairs$D), 1], xmax = pairs[which.min(pairs$D), 2], no_of_sims = 1000, threads = 8)

pairs[which.min(pairs$D), 1] #k_{sat}

43
pairs[which.min(pairs$D), 2] #k_{cut}

568

#in this range
pairs[which.min(pairs$D), 3] #D

0.07908323

pairs[which.min(pairs$D), 4] #gamma
2.017185

bs_pl_sat_cut$p #p-value

0.991

#ksat, kcut, D, γ, p-value (by bootstratpping)

pairs[which.min(pairs$D), 1] -> k_sat
pairs[which.min(pairs$D), 2] -> k_cut
pairs[which.min(pairs$D), 4] -> gamma

#powerlaw
m_pl = displ$new(data)
est_pl <- estimate_xmin(m_pl, xmins = k_sat, xmax = k_cut, distance = "ks")
m_pl$setXmin(est_pl)

#lognormal
m_ln = dislnorm$new(data)
est_ln <- estimate_xmin(m_ln)
m_ln$setXmin(est_ln)

#exponential
m_exp = disexp$new(data)
est_exp <- estimate_xmin(m_exp)
m_exp$setXmin(est_exp)

#poisson
m_poi = dispois$new(data)
est_poi <- estimate_xmin(m_poi)
m_poi$setXmin(est_poi)

m_pl$setXmin(est_pl)
plot.data <- plot(m_pl, draw = F)
fit.data1 <- lines(m_pl, draw = F)
fit.data2 <- lines(m_ln, draw = F)
fit.data3 <- lines(m_exp, draw = F)
fit.data4 <- lines(m_poi, draw = F)
#"Poisson" = "dotdash"
ggplot(plot.data) +
  geom_point(aes(x = log(x), y = log(y)), color = "black", size = 2) +
  labs(x = "log(k)", y = "log(CDF)") +
  theme_bw() +
  geom_line(data = fit.data1, aes(x = log(x), y = log(y), linetype = "Power-law"), color = "black", size = 1) +
  geom_line(data = fit.data2, aes(x = log(x), y = log(y), linetype = "Log-normal"), color = "black", size = 1) +
  geom_line(data = fit.data3, aes(x = log(x), y = log(y), linetype = "Exponential"), color = "black", size = 1)+ 
  #geom_line(data = fit.data4, aes(x = log(x), y = log(y), linetype = "Poisson"), color = "black", size = 1) +
  scale_linetype_manual(
    values = c("Power-law" = "solid", "Log-normal" = "twodash", "Exponential" = "dotted"),
    name = "Distribution"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.96),  # Adjust the position as needed
    legend.justification = c(1, 1),
    legend.box.just = "right",
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(color = "black", size = 1),
  )




