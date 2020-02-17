#full
rm(list=ls())
library(devtools)
setwd('/Users/bomin8319/Desktop/DAME_revision/DAME_code_revised/DAME_pkg_revised/R')
load_all()
load("/Users/bomin8319/Desktop/DAME_revision/DAME_code_revised/UNdatafull.RData")
attach(UNdatafull)
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
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# 97 country version
Time = 32
N = 97

# not existing countries -> all missing values imputed using model (biased)
avail1 = matrix(1, Time, N)
colnames(avail1) = colnames(UNdatafull$Y)
avail1[1:8, which(colnames(avail1) %in% c("ROK", "PRK"))] = 0 #North and South Korea did not joined UN voting until 1990
avail1[1:9, which(colnames(avail1) %in% c("RUS"))] = 0 #RUS X variables not existed until 1991
avail1[13:21, which(colnames(avail1) %in% c("IRQ"))] = 0 #IRQ under sanction

Degrees = vapply(1:Time, function(tp) {rowSums(Y[tp,,], na.rm = TRUE)}, rep(0, N))
corr = vapply(1:31, function(l) {cor(Degrees[1:(N*(Time - l))], Degrees[(1 + N*l):(N*Time)], use = "complete")}, 0)


load('~/Desktop/UN_full5.RData') #DAME

avail = avail1
  colnames(avail) = dimnames(Y)[[2]]
  meaningful_NA_rows = lapply(1:Time, function(tp) {
    which(avail[tp,]==0)
  })
  meaningful_NA = lapply(1:Time, function(tp) {
    pre = matrix(0, N, N)
    pre[meaningful_NA_rows[[tp]],] = NA
    pre[,meaningful_NA_rows[[tp]]] = NA
    which(is.na(pre)==TRUE)
  })
      
setwd('/Users/bomin8319/Desktop/GOF')
for (tp in 1:32) {
	mname = paste0("Year", tp, ".png")
    png(filename = mname)
    diag(Y[tp,,]) = 0
    Y[tp, , ][meaningful_NA[[tp]]] = 0
    	pureNA = which(rowSums(Y[tp,,], na.rm = TRUE)==0)
    	if (length(meaningful_NA_rows[[tp]]) > 0) {
    	pureNA = pureNA[-which(pureNA %in% meaningful_NA_rows[[tp]])]
    	}
    Y[tp, , ][is.na(Y[tp, , ])] = 0 
	if (tp < 10 | tp %in% 12:21) {
	par(mfrow = c(3,2))	
	hist(apply(UN$Degree[[tp]][,-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA)), 1, mean), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "1st moment mean")
	abline(v = mean(rowSums(Y[tp,,], na.rm = TRUE)[-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA))), col = "red")
	hist(apply(UN$Degree[[tp]][,-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA)), 1, sd), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "1st moment sd")
	abline(v = sd(rowSums(Y[tp,,], na.rm = TRUE)[-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA))), col = "red")
	hist(apply(UN$secondDegree[[tp]][,-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA)), 1, mean), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "2nd moment mean")
	abline(v = mean(rowSums(Y[tp,,] %*% Y[tp,,], na.rm = TRUE)[-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA))), col = "red")
	hist(apply(UN$secondDegree[[tp]][,-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA)), 1, sd), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "2nd moment sd")
	abline(v = sd(rowSums(Y[tp,,]%*% Y[tp,,], na.rm = TRUE)[-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA))), col = "red")
	hist(apply(UN$thirdDegree[[tp]][,-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA)), 1, mean), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "3rd moment mean")
	abline(v = mean(rowSums(Y[tp,,] %*% Y[tp,,] %*% Y[tp,,], na.rm = TRUE)[-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA))), col = "red")
	hist(apply(UN$thirdDegree[[tp]][,-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA)), 1, sd), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "3rd moment sd")
	abline(v = sd(rowSums(Y[tp,,]%*% Y[tp,,] %*% Y[tp,,], na.rm = TRUE)[-c(which(avail1[tp,]==0),pureNA)]/(sum(avail1[tp,])-length(pureNA))), col = "red")
	} else {
	par(mfrow = c(3,2))	
	hist(apply(UN$Degree[[tp]][,]/sum(avail1[tp,]), 1, mean), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "1st moment mean")
	abline(v = mean(rowSums(Y[tp,,], na.rm = TRUE)/sum(avail1[tp,])), col = "red")
	hist(apply(UN$Degree[[tp]][,]/sum(avail1[tp,]), 1, sd), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "1st moment sd")
	abline(v = sd(rowSums(Y[tp,,], na.rm = TRUE)/sum(avail1[tp,])), col = "red")
	hist(apply(UN$secondDegree[[tp]][,]/sum(avail1[tp,]), 1, mean), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "2nd moment mean")
	abline(v = mean(rowSums(Y[tp,,] %*% Y[tp,,], na.rm = TRUE)/sum(avail1[tp,])), col = "red")
	hist(apply(UN$secondDegree[[tp]][,]/sum(avail1[tp,]), 1, sd), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "2nd moment sd")
	abline(v = sd(rowSums(Y[tp,,]%*% Y[tp,,], na.rm = TRUE)/sum(avail1[tp,])), col = "red")
	hist(apply(UN$thirdDegree[[tp]][,]/sum(avail1[tp,]), 1, mean), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "3rd moment mean")
	abline(v = mean(rowSums(Y[tp,,] %*% Y[tp,,] %*% Y[tp,,], na.rm = TRUE)/sum(avail1[tp,])), col = "red")
	hist(apply(UN$thirdDegree[[tp]][,]/sum(avail1[tp,]), 1, sd), main = "", prob = FALSE, col = "lightblue", ylab = "", xlab = "3rd moment sd")
	abline(v = sd(rowSums(Y[tp,,]%*% Y[tp,,] %*% Y[tp,,], na.rm = TRUE)/sum(avail1[tp,])), col = "red")
	}
	dev.off()
}