#function that loads the data of the Zackenberg networks
#inputs: year = 1996 or 1997, day = day of observation
#output: web = mutualsitic network
load_data <- function(year,day){
  name <- paste('matrices/',year,'_',day,'.txt',sep='')
  d <- read.table(file=name,header=TRUE)
  web <- as.matrix(d[,2:dim(d)[2]])
  rownames(web) <- d[,1]
  return(web)
}


#compute the interaction strength matrices according to the parameterization given in the methods section
#inputs: web = mutualistic network (binary matrix), gamma_avg = average mutualistic strength, 
#rho = mean field interspecific competition, delta = mutualistic trade-off
#outputs: alpha = interaction strength matrix (full), alphaA = competition among plants, betaP = competition among animals,
#gammaA = mutualistic effect of plants on animals, gammaP = mutualistic effect of animals on plants
interaction_matrix <- function(web,gamma_avg,rho,delta){
  SA <- nrow(web)
  SP <- ncol(web)
  alphaA <- matrix(rho,SA,SA) + (1-rho) * diag(rep(1,SA))
  alphaP <- matrix(rho,SP,SP) + (1-rho) * diag(rep(1,SP))
  gammaA <- diag(rowSums(web)^-delta) %*% web
  gammaP <- diag(colSums(web)^-delta) %*% t(web)
  f <- sum(gammaA[web == 1] + gammaP[t(web) == 1] ) / (2 * sum(web==1))
  gammaA <- gamma_avg/f * diag(rowSums(web)^-delta) %*% web
  gammaP <- gamma_avg/f * diag(colSums(web)^-delta) %*% t(web)
  alpha <- rbind(cbind(alphaA,-gammaA),cbind(-gammaP,alphaP))
  out <- list(alpha = alpha, alphaA = alphaA, alphaP = alphaP, gammaA = gammaA, gammaP = gammaP)
  return(out)
}

get_interaction_matrix <- function(A){
  gram <- t(A) %*% A
  for(i in 1:ncol(gram)){
    gram[,i] <- gram[,i]/sum(gram[,i])
  }
  diag(gram) <- 1
  return(gram)
}

#stability condition gamma_hat (average mutualistic strength at the stability threshold)
#inputs: web = mutualistic network (binary matrix),
#rho = mean field interspecific competition, delta = mutualistic trade-off
#output: gamma_hat = stability condition
gamma_hat <- function(web,rho,delta){
  f_eig <- function(gamma_avg,web,rho,delta){
    alpha <- interaction_matrix(web,gamma_avg,rho,delta)$alpha
    out <- (min(Re(eigen(alpha)$values)))^2
    out
  }
  out <- optimize(f_eig,c(0,1000),web = web, rho = rho, delta = delta)$minimum
  return(out)
}


#feasibility condition---Omega
#input: alpha = interaction strength matrix
#output: Omega = feasibility condition (on a log10 scale)
require(mvtnorm) #this function requires the librrary mvtnorm (Genz et al. 2009) for the
#numerical computation of the cumulative distribution of the multivariate normal distirbution
Omega <- function(alpha){
  S <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha)
  m <- matrix(0,S,1)
  a <- matrix(0,S,1)
  b <- matrix(Inf,S,1)  
  d <- pmvnorm(lower = rep(0,S), upper = rep(Inf,S), mean = rep(0,S), sigma = Sigma)
  out <- (d[1])^(1/S)
  return(out)
}