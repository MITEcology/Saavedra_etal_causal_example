rm(list=ls())
library('pcalg')
source('toolbox.R') #load the toolbox
source('load_data.R') #load the data

## Following the text 
## variable A is named here as bratio
## variable B is named here as bpollinators_var
## variable C is named here as btemp
## variable D is named here as bomega (synthetic data)

### test unconditional d-separation between C (temp) and D (context/omega)
dat <- as.matrix(data.frame(bratio,bpollinators_var,btemp,bomega))
gSquareBin(3,4,c(NULL),dat,n.min = 20)

##ACE of C (temp) on A (relative plant abundance)
d <- mean(bratio[btemp==1],na.rm = T)-mean(bratio[btemp==0],na.rm = T)
d ## ACE
gSquareBin(3,1,c(NULL),dat,n.min = 20)


##ACE of C (temp) on B (pollinator variability) filtered by mediator of high relative plant abundance A
d1 <- mean(bpollinators_var[btemp==1&bratio==1&bomega==1],na.rm = T)*mean(bomega,na.rm = T)
d2 <- mean(bpollinators_var[btemp==1&bratio==1&bomega==0],na.rm = T)*(1-mean(bomega,na.rm = T))
d5 <- mean(bpollinators_var[btemp==0&bratio==1&bomega==1],na.rm = T)*mean(bomega,na.rm = T)
d6 <- mean(bpollinators_var[btemp==0&bratio==1&bomega==0],na.rm = T)*(1-mean(bomega,na.rm = T))
if (is.na(d1)) {d1 <- 0}
if (is.na(d2)) {d2 <- 0}
if (is.na(d5)) {d5 <- 0}
if (is.na(d6)) {d6 <- 0}
d <- d1 + d2 - d5 - d6
d ## ACE
dat2 <- as.matrix(data.frame(bpollinators_var[bratio==0],btemp[bratio==0],bomega[bratio==0]))
gSquareBin(1,2,c(3),dat2,n.min = 20)

##ACE of W (temp) on Y (pollinator variability) filtered by mediator of low relative plant abundance
d1 <- mean(bpollinators_var[btemp==1&bratio==0&bomega==1],na.rm = T)*mean(bomega,na.rm = T)
d2 <- mean(bpollinators_var[btemp==1&bratio==0&bomega==0],na.rm = T)*(1-mean(bomega,na.rm = T))
d5 <- mean(bpollinators_var[btemp==0&bratio==0&bomega==1],na.rm = T)*mean(bomega,na.rm = T)
d6 <- mean(bpollinators_var[btemp==0&bratio==0&bomega==0],na.rm = T)*(1-mean(bomega,na.rm = T))
if (is.na(d1)) {d1 <- 0}
if (is.na(d2)) {d2 <- 0}
if (is.na(d5)) {d5 <- 0}
if (is.na(d6)) {d6 <- 0}
d <- d1 + d2 - d5 - d6
d  ## ACE
dat2 <- as.matrix(data.frame(bpollinators_var[bratio==0],btemp[bratio==0],bomega[bratio==0]))
gSquareBin(1,2,c(3),dat2,n.min = 20)

##ACE of A (relative plant abundance) on B (pollinator variability)
d1 <- mean(bpollinators_var[bratio==1&btemp==1&bomega==1],na.rm = T)*mean(btemp,na.rm = T)*mean(bomega[btemp==1],na.rm = T)
d2 <- mean(bpollinators_var[bratio==1&btemp==1&bomega==0],na.rm = T)*mean(btemp,na.rm = T)*(1-mean(bomega[btemp==1],na.rm = T))
d3 <- mean(bpollinators_var[bratio==1&btemp==0&bomega==1],na.rm = T)*(1-mean(btemp,na.rm = T))*(mean(bomega[btemp==0],na.rm = T))
d4 <- mean(bpollinators_var[bratio==1&btemp==0&bomega==0],na.rm = T)*(1-mean(btemp,na.rm = T))*(1-mean(bomega[btemp==0],na.rm = T))
d5 <- mean(bpollinators_var[bratio==0&btemp==1&bomega==1],na.rm = T)*mean(btemp,na.rm = T)*mean(bomega[btemp==1],na.rm = T)
d6 <- mean(bpollinators_var[bratio==0&btemp==1&bomega==0],na.rm = T)*mean(btemp,na.rm = T)*(1-mean(bomega[btemp==1],na.rm = T))
d7 <- mean(bpollinators_var[bratio==0&btemp==0&bomega==1],na.rm = T)*(1-mean(btemp,na.rm = T))*(mean(bomega[btemp==0],na.rm = T))
d8 <- mean(bpollinators_var[bratio==0&btemp==0&bomega==0],na.rm = T)*(1-mean(btemp,na.rm = T))*(1-mean(bomega[btemp==0],na.rm = T))
if (is.na(d1)) {d1 <- 0}
if (is.na(d2)) {d2 <- 0}
if (is.na(d3)) {d3 <- 0}
if (is.na(d4)) {d4 <- 0}
if (is.na(d5)) {d5 <- 0}
if (is.na(d6)) {d6 <- 0}
if (is.na(d7)) {d7 <- 0}
if (is.na(d8)) {d8 <- 0}
d <- d1 + d2 + d3 + d4 - d5 - d6 - d7 -d8
d  ## ACE
gSquareBin(1,2,c(3,4),dat,n.min = 20)