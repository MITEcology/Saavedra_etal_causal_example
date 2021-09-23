### load the abiotic factors for each season
load('abiotic/e96.RData')
mean_temperature_96 <- mean_temperature
load('abiotic/e97.RData')
mean_temperature_97 <- mean_temperature

### initialize variables for 1996 season
omega_96 <- 0 ## feasibility domain
size_96 <- 0 ## plants + pollinators
plants_96 <- 0 ## number of plants
pollinators_96 <- 0 ## number of pollinators
### obtain information for 1996 season
for (i in 1:24){ ## 24 days in this season
  web <- load_data(1996,i) ## load binary plant-animal network for each day
  t <- gamma_hat(web,rho=0,delta = 0.1) #compute the stability condition for rho=0 and delta=0.1
  m <- get_interaction_matrix(web) ## generate competition matrix for pollinators
  plants_96[i] <- ncol(web)
  pollinators_96[i] <- nrow(web)
  size_96[i] <- ncol(web) + nrow(web)
  
  o <- 0 ## this loop is for convergence of omega
  for (j in 1:30){
    o <- o + Omega(m) 
  }
  omega_96[i] <- o / 30
  
}

pollinators_var_96 <- abs(pollinators_96[-1]-pollinators_96[-24])/pollinators_96[-24]

ratio_96 <- plants_96/pollinators_96 ## relative abundance of plants


### transform variables into binary using the median as cut-off
bomega_96 <- omega_96
bomega_96[omega_96>median(omega_96)] <- 1
bomega_96[omega_96<=median(omega_96)] <- 0

btemp_96 <- mean_temperature_96
btemp_96[mean_temperature_96>median(mean_temperature_96)] <- 1
btemp_96[mean_temperature_96<=median(mean_temperature_96)] <- 0

bratio_96 <- ratio_96
bratio_96[ratio_96>median(ratio_96)] <- 1
bratio_96[ratio_96<=median(ratio_96)] <- 0

bpollinators_var_96 <- pollinators_var_96
bpollinators_var_96[pollinators_var_96>median(pollinators_var_96)] <- 1
bpollinators_var_96[pollinators_var_96<=median(pollinators_var_96)] <- 0
bpollinators_var_96 <- c(bpollinators_var_96,1)

### Do the same for 1997 season
omega_97 <- 0
size_97 <- 0
plants_97 <- 0
pollinators_97 <- 0
for (i in 1:26){
  web <- load_data(1997,i) 
  m <- get_interaction_matrix(web)
  plants_97[i] <- ncol(web)
  pollinators_97[i] <- nrow(web)
  size_97[i] <- ncol(web) + nrow(web)
  
  o <- 0
  for (j in 1:30){
    o <- o + Omega(m) 
  }
  omega_97[i] <- o / 30
  
}

pollinators_var_97 <- abs(pollinators_97[-1]-pollinators_97[-26])/pollinators_97[-26]

ratio_97 <- plants_97/pollinators_97


bomega_97 <- omega_97
bomega_97[omega_97>median(omega_97)] <- 1
bomega_97[omega_97<=median(omega_97)] <- 0

btemp_97 <- mean_temperature_97
btemp_97[mean_temperature_97>median(mean_temperature_97)] <- 1
btemp_97[mean_temperature_97<=median(mean_temperature_97)] <- 0

bratio_97 <- ratio_97
bratio_97[ratio_97>median(ratio_97)] <- 1
bratio_97[ratio_97<=median(ratio_97)] <- 0

bpollinators_var_97 <- pollinators_var_97
bpollinators_var_97[pollinators_var_97>median(pollinators_var_97)] <- 1
bpollinators_var_97[pollinators_var_97<=median(pollinators_var_97)] <- 0
bpollinators_var_97 <- c(bpollinators_var_97,0)


### merge the two seasons
bpollinators_var <- c(bpollinators_var_96,bpollinators_var_97)
bratio <- c(bratio_96,bratio_97)
bomega <- c(bomega_96,bomega_97)
btemp <- c(btemp_96,btemp_97)