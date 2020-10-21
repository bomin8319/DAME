#notes on prior search for tau_beta's
rm(list=ls())
load('~/Desktop/DAME_revised/UNdatafull.RData')
attach(UNdatafull)
dim(UNdatafull$X)
dimnames(UNdatafull$X)

library(FastGP)
library(mvtnorm)
library(fields)
library(reshape)
library(MCMCpack)
library(expm)
library(igraph)
library(coda)
library(ggplot2)
library(LaplacesDemon)
library(truncnorm)
library(gridExtra)
library(Rfast)
#construct GP
Time = 32
dist_ij = c()
 for (i in 1:Time) {
 for (j in 1:Time) {
 dist_ij = c(dist_ij, abs(i-j))
 }
 }
dist_ij = matrix(dist_ij, Time, Time)

uppertri = upper.tri(diag(97))

#Y 
s2 = 0.0250
X_scale = list()
for (p in 1:6) {
XX = sapply(1:Time, function(tp) sum(X[tp, , , p][uppertri]^2, na.rm = TRUE))
#X_scale[[p]] = diag(XX, Time, Time)/s2
X_scale[[p]] = s2* solve(diag(XX, Time, Time))
}
lapply(X_scale, function(x) sum(uppertri)*max(diag(x)))
#lapply(X_scale, function(x) median(diag(x)))

#(k, theta) search based on Wasserman (2000)
#2. Richardson and Green (1997)
for (p in 1:6) {
	print(summary(1/rgamma(10^4, shape = 2, rate = sum(uppertri)*(max(diag(X_scale[[p]]))))))
}

#Y*100
s2 = 250.0
X_scale = list()
for (p in 1:6) {
XX = sapply(1:Time, function(tp) sum(X[tp, , , p][uppertri]^2, na.rm = TRUE))
#X_scale[[p]] = diag(XX, Time, Time)/s2
X_scale[[p]] = s2* solve(diag(XX, Time, Time))
}
lapply(X_scale, function(x) sum(uppertri)*max(diag(x)))
#lapply(X_scale, function(x) median(diag(x)))
for (p in 1:6) {
	print(summary(1/rgamma(10^4, shape = 2, rate = sum(uppertri)*(max(diag(X_scale[[p]]))))))
}

source('~/Desktop/DAME_revised/DAME_pkg_revised/R/DAME.R', chdir = TRUE)



# not existing countries -> all missing values imputed using model (biased)
Time = 32
N = 97
avail1 = matrix(1, Time, N)
colnames(avail1) = colnames(UNdatafull$Y)
avail1[1:8, which(colnames(avail1) %in% c("ROK", "PRK"))] = 0 #North and South Korea did not joined UN voting until 1990
avail1[1:9, which(colnames(avail1) %in% c("RUS"))] = 0 #RUS X variables not existed until 1991
avail1[13:21, which(colnames(avail1) %in% c("IRQ"))] = 0 #IRQ under sanction

###rerun
UN3_test = DAME_MH_revised(Y[1:Time,,], X[1:Time,,,1:6], RE = c(), R = 2, gammapriors = c(2, 1, 0.025), avail = avail1, burn = 20000, nscan = 50000, odens = 50)
UN3_test2 = DAME_MH_revised(Y[1:Time,,]*100, X[1:Time,,,1:6], RE = c(), R = 2, gammapriors = c(2, 1, 250), avail = avail1, burn = 20000, nscan = 50000, odens = 50)
save(UN3_test, file = "/Users/bomin8319/Desktop/UN_full3_test.RData")
save(UN3_test2, file = "/Users/bomin8319/Desktop/UN_full3_test2.RData")

set.seed(1)
UN_newprior = DAME_MH_revised(Y[1:Time,,]*100, X[1:Time,,,1:6], RE = c("additive", "multiplicative"), R = 2, gammapriors = c(2, 1, 250), avail = avail1, burn = 50000, nscan = 100000, odens = 50)
