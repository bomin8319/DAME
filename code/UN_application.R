#full
rm(list=ls())
load('~/Desktop/DAME_revised/UNdatafull.RData')
attach(UNdatafull)
source('~/Desktop/DAME_revised/code/DAME.R', chdir = TRUE)
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

#################

UN = DAME_MH_revised(Y[1:Time,,]*100, X[1:Time,,,1:6], RE = c("additive", "multiplicative"), R = 2, avail = avail1, burn = 50000, nscan = 100000, odens = 50)
#save(UN, file = "/Users/bomin8319/Desktop/UN_full_long.RData")
UN2 = DAME_MH_revised(Y[1:Time,,]*100, X[1:Time,,,1:6], RE = c("additive"), R = 2, avail = avail1, burn = 20000, nscan = 50000, odens = 50)
#save(UN2, file = "/Users/bomin8319/Desktop/UN_full2.RData")
UN3 = DAME_MH_revised(Y[1:Time,,]*100, X[1:Time,,,1:6], RE = c(), R = 2, avail = avail1, burn = 20000, nscan = 50000, odens = 50)
#save(UN3, file = "/Users/bomin8319/Desktop/UN_full3.RData")
UN5 = DAME_UU_fixed_revised(Y[1:Time,,]*100, X[1:Time,,,1:6], RE = c("multiplicative"), R = 2, avail = avail1, burn = 20000, nscan = 50000, odens = 50, kappas = rep(10, 6+1+2))
#save(UN5, file = "/Users/bomin8319/Desktop/UN_full5.RData")
UN4 = UN5

#########summary of results#######
## side by side plot

hi = factor(1983:2014)
hello = list()
hello2 = list()
hello3 = list()
hello4 = list()
hellonew = list()
mean = list()

Y = Y * 100
n2 = 1
for	(n in 1:N){
  print(n)
  hello[[n]] = matrix(NA, nrow = 0, ncol = 3)
  mean[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    diag(Y[d,, ])= 0
    Y[d, which(avail1[d,]==0), ] = 0
    Y[d, , which(avail1[d,]==0)] = 0
    Y[d, , ][which(is.na(Y[d, , ]))] = 0	
    hello[[n]]  = rbind(hello[[n]], cbind(UN$Degree[[d]][,n2], rep(hi[d], length(UN$Degree[[d]][,n2])), rep("DAME", length(UN$Degree[[d]][,n2]))))
    mean[[n]] = rbind(mean[[n]], c(rowSums(Y[d,,], na.rm=TRUE)[n2], hi[d], "DAME"))
  }
  colnames(hello[[n]]) = c("Degree", "Year", "Model")
  hello[[n]] = as.data.frame(hello[[n]])
  hello[[n]]$Year = factor(sort(as.numeric(hello[[n]]$Year)), labels = c(1983:2014))
  
  hello2[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello2[[n]]  = rbind(hello2[[n]], cbind(UN2$Degree[[d]][,n2], rep(hi[d], length(UN2$Degree[[d]][,n2])), rep("AE", length(UN2$Degree[[d]][,n2]))))
  }
  colnames(hello2[[n]]) = c("Degree", "Year", "Model")
  hello2[[n]] = as.data.frame(hello2[[n]])
  hello2[[n]]$Year = factor(sort(as.numeric(hello2[[n]]$Year)), labels = c(1983:2014))
  
  hello3[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello3[[n]]  = rbind(hello3[[n]], cbind(UN3$Degree[[d]][,n2], rep(hi[d], length(UN3$Degree[[d]][,n2])), rep("NO", length(UN3$Degree[[d]][,n2]))))
  }
  colnames(hello3[[n]]) = c("Degree", "Year", "Model")
  hello3[[n]] = as.data.frame(hello3[[n]])
  hello3[[n]]$Year = factor(sort(as.numeric(hello3[[n]]$Year)), labels = c(1983:2014))
  
  hello4[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello4[[n]]  = rbind(hello4[[n]], cbind(UN4$Degree[[d]][,n2], rep(hi[d], length(UN4$Degree[[d]][,n2])), rep("ME", length(UN4$Degree[[d]][,n2]))))
  }
  colnames(hello4[[n]]) = c("Degree", "Year", "Model")
  hello4[[n]] = as.data.frame(hello4[[n]])
  hello4[[n]]$Year = factor(sort(as.numeric(hello4[[n]]$Year)), labels = c(1983:2014))
  
  
  hellonew[[n]] = as.data.frame(rbind(hello[[n]], hello4[[n]], hello2[[n]], hello3[[n]]))
  colnames(mean[[n]]) =  c("Degree", "Year", "Model")
  mean[[n]] = as.data.frame(mean[[n]])
  n2 = n2 + 1
}

countryname = (rownames(UN$U[[Time]]))
countryname2 = (rownames(UN$U[[Time]]))


######degree statistics
setwd('/Users/bomin8319/Desktop/model_validation_new')
countrynames = ggplotColours(4)
p10 = list()
n2 = 1
for (n in 1:N){
  hellonew[[n]]$Degree = as.numeric(as.character(hellonew[[n]]$Degree))
  mean[[n]]$Degree = as.numeric(as.character(mean[[n]]$Degree))
  mean[[n]]$Year = hi
  p10[[n]] = ggplot(hellonew[[n]], aes(x = Year, y = Degree, color= Model, fill =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal() +scale_fill_manual(values = alpha(countrynames,0.5), name = countryname[n])+scale_color_manual(values = countrynames,name = countryname[n]) + geom_line(data = mean[[n]], color = "blue", size = 0.2, group = 1)+geom_point(data = mean[[n]], color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
  #annotate("text", 30, 10, size = 5, label = countryname2[n], colour = countrynames[n])
  n2 = n2 + 1
  
  mname = paste0(countryname2[n], n, ".png")
  print(p10[[n]])
  ggsave(filename = mname, width = 12, height = 6)
}


#secondDegree
hello = list()
hello2 = list()
hello3 = list()
hello4 = list()
hellonew = list()
mean = list()

n2 = 1
for	(n in 1:N){
  print(n)
  hello[[n]] = matrix(NA, nrow = 0, ncol = 3)
  mean[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    #diag(Y[d,, ])= 0
    #Y[d, which(avail1[d,]==0), ] = 0
    #Y[d, , which(avail1[d,]==0)] = 0
    #Y[d, , ][which(is.na(Y[d, , ]))] = 0	
    hello[[n]]  = rbind(hello[[n]], cbind(UN$secondDegree[[d]][,n2], rep(hi[d], length(UN$secondDegree[[d]][,n2])), rep("DAME", length(UN$secondDegree[[d]][,n2]))))
    mean[[n]] = rbind(mean[[n]], c(rowSums(Y[d,,] %*% Y[d,,], na.rm = TRUE)[n2], hi[d], "DAME"))
  }
  colnames(hello[[n]]) = c("SecondDegree", "Year", "Model")
  hello[[n]] = as.data.frame(hello[[n]])
  hello[[n]]$Year = factor(sort(as.numeric(hello[[n]]$Year)), labels = c(1983:2014))
  
  hello2[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello2[[n]]  = rbind(hello2[[n]], cbind(UN2$secondDegree[[d]][,n2], rep(hi[d], length(UN2$secondDegree[[d]][,n2])), rep("AE", length(UN2$secondDegree[[d]][,n2]))))
  }
  colnames(hello2[[n]]) = c("SecondDegree", "Year", "Model")
  hello2[[n]] = as.data.frame(hello2[[n]])
  hello2[[n]]$Year = factor(sort(as.numeric(hello2[[n]]$Year)), labels = c(1983:2014))
  
  hello3[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello3[[n]]  = rbind(hello3[[n]], cbind(UN3$secondDegree[[d]][,n2], rep(hi[d], length(UN3$secondDegree[[d]][,n2])), rep("NO", length(UN3$secondDegree[[d]][,n2]))))
  }
  colnames(hello3[[n]]) = c("SecondDegree", "Year", "Model")
  hello3[[n]] = as.data.frame(hello3[[n]])
  hello3[[n]]$Year = factor(sort(as.numeric(hello3[[n]]$Year)), labels = c(1983:2014))
  
  hello4[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello4[[n]]  = rbind(hello4[[n]], cbind(UN4$secondDegree[[d]][,n2], rep(hi[d], length(UN4$secondDegree[[d]][,n2])), rep("ME", length(UN4$secondDegree[[d]][,n2]))))
  }
  colnames(hello4[[n]]) = c("SecondDegree", "Year", "Model")
  hello4[[n]] = as.data.frame(hello4[[n]])
  hello4[[n]]$Year = factor(sort(as.numeric(hello4[[n]]$Year)), labels = c(1983:2014))
  
  
  hellonew[[n]] = as.data.frame(rbind(hello[[n]], hello4[[n]], hello2[[n]], hello3[[n]]))
  colnames(mean[[n]]) =  c("SecondDegree", "Year", "Model")
  mean[[n]] = as.data.frame(mean[[n]])
  n2 = n2 + 1
}

countryname = (rownames(UN$U[[Time]]))
countryname2 = (rownames(UN$U[[Time]]))

p10 = list()
n2 = 1
for (n in 1:N){
  hellonew[[n]]$SecondDegree = as.numeric(as.character(hellonew[[n]]$SecondDegree))
  mean[[n]]$SecondDegree = as.numeric(as.character(mean[[n]]$SecondDegree))
  mean[[n]]$Year = hi
  p10[[n]] = ggplot(hellonew[[n]], aes(x = Year, y = SecondDegree, color= Model, fill =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal() +scale_fill_manual(values = alpha(countrynames,0.5), name = countryname[n])+scale_color_manual(values = countrynames,name = countryname[n]) + geom_line(data = mean[[n]], color = "blue", size = 0.2, group = 1)+geom_point(data = mean[[n]], color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
  
  #annotate("text", 30, 10, size = 5, label = countryname2[n], colour = countrynames[n])
  n2 = n2 + 1
  
  mname = paste0(countryname2[n], n, "second.png")
  print(p10[[n]])
  ggsave(filename = mname, width = 12, height = 6)
}

#thirdDegree
hello = list()
hello2 = list()
hello3 = list()
hello4 = list()
hellonew = list()
mean = list()

n2 = 1
for	(n in 1:N){
  print(n)
  hello[[n]] = matrix(NA, nrow = 0, ncol = 3)
  mean[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    #diag(Y[d,, ])= 0
    #Y[d, which(avail1[d,]==0), ] = 0
    #Y[d, , which(avail1[d,]==0)] = 0
    #Y[d, , ][which(is.na(Y[d, , ]))] = 0	
    hello[[n]]  = rbind(hello[[n]], cbind(UN$thirdDegree[[d]][,n2], rep(hi[d], length(UN$thirdDegree[[d]][,n2])), rep("DAME", length(UN$thirdDegree[[d]][,n2]))))
    mean[[n]] = rbind(mean[[n]], c(rowSums(Y[d,,] %*% Y[d,,] %*% Y[d,,], na.rm = TRUE)[n2], hi[d], "DAME"))
  }
  colnames(hello[[n]]) = c("ThirdDegree", "Year", "Model")
  hello[[n]] = as.data.frame(hello[[n]])
  hello[[n]]$Year = factor(sort(as.numeric(hello[[n]]$Year)), labels = c(1983:2014))
  
  hello2[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello2[[n]]  = rbind(hello2[[n]], cbind(UN2$thirdDegree[[d]][,n2], rep(hi[d], length(UN2$thirdDegree[[d]][,n2])), rep("AE", length(UN2$thirdDegree[[d]][,n2]))))
  }
  colnames(hello2[[n]]) = c("ThirdDegree", "Year", "Model")
  hello2[[n]] = as.data.frame(hello2[[n]])
  hello2[[n]]$Year = factor(sort(as.numeric(hello2[[n]]$Year)), labels = c(1983:2014))
  
  hello3[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello3[[n]]  = rbind(hello3[[n]], cbind(UN3$thirdDegree[[d]][,n2], rep(hi[d], length(UN3$thirdDegree[[d]][,n2])), rep("NO", length(UN3$thirdDegree[[d]][,n2]))))
  }
  colnames(hello3[[n]]) = c("ThirdDegree", "Year", "Model")
  hello3[[n]] = as.data.frame(hello3[[n]])
  hello3[[n]]$Year = factor(sort(as.numeric(hello3[[n]]$Year)), labels = c(1983:2014))
  
  hello4[[n]] = matrix(NA, nrow = 0, ncol = 3)
  for (d in 1:Time){
    hello4[[n]]  = rbind(hello4[[n]], cbind(UN4$thirdDegree[[d]][,n2], rep(hi[d], length(UN4$thirdDegree[[d]][,n2])), rep("ME", length(UN4$thirdDegree[[d]][,n2]))))
  }
  colnames(hello4[[n]]) = c("ThirdDegree", "Year", "Model")
  hello4[[n]] = as.data.frame(hello4[[n]])
  hello4[[n]]$Year = factor(sort(as.numeric(hello4[[n]]$Year)), labels = c(1983:2014))
  
  
  hellonew[[n]] = as.data.frame(rbind(hello[[n]], hello4[[n]], hello2[[n]], hello3[[n]]))
  colnames(mean[[n]]) =  c("ThirdDegree", "Year", "Model")
  mean[[n]] = as.data.frame(mean[[n]])
  n2 = n2 + 1
}

countryname = (rownames(UN$U[[Time]]))
countryname2 = (rownames(UN$U[[Time]]))

p10 = list()
n2 = 1
for (n in 1:N){
  hellonew[[n]]$ThirdDegree = as.numeric(as.character(hellonew[[n]]$ThirdDegree))
  mean[[n]]$ThirdDegree = as.numeric(as.character(mean[[n]]$ThirdDegree))
  mean[[n]]$Year = hi
  p10[[n]] = ggplot(hellonew[[n]], aes(x = Year, y = ThirdDegree, color= Model, fill =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal() +scale_fill_manual(values = alpha(countrynames,0.5), name = countryname[n])+scale_color_manual(values = countrynames,name = countryname[n]) + geom_line(data = mean[[n]], color = "blue", size = 0.2, group = 1)+geom_point(data = mean[[n]], color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
  
  #annotate("text", 30, 10, size = 5, label = countryname2[n], colour = countrynames[n])
  n2 = n2 + 1
  
  mname = paste0(countryname2[n], n, "third.png")
  print(p10[[n]])
  ggsave(filename = mname, width = 12, height = 6)
}



##########estimated coefficients####
#beta
beta = lapply(1:Time, function(t){summary(mcmc(UN$BETA[[t]][1:2000,]))[[2]]})
betas = list()
for (i in 1:6) {betas[[i]] = sapply(1:Time, function(t){beta[[t]][i,]})}
betacols= ggplotColours(6)
plots = list()
i= 1
years = c(1983:2014)
data = data.frame(cbind(years,t(betas[[i]])))
 colnames(data)[4] = "beta"
 f = ggplot(data, aes(x = years))
plots[[1]]=f + geom_line(aes(y = beta), colour="black") + geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.1) + scale_x_continuous(breaks=number_ticks(8))+ylab("Intercept") + theme_minimal()

color =betacols[-1]
varname = c("log(distance)", "Polity", "Alliance", "Lower trade-to-GDP ratio", "Common Language")
for (i in 2:6){
years = c(1983:2014)
data = data.frame(cbind(years,t(betas[[i]])))
 colnames(data)[4] = "beta"
 f = ggplot(data, aes(x = years))
 
 plots[[i]] <- f + geom_line(aes(y = beta), colour="black") + geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.1) + ylab(varname[i-1]) + scale_x_continuous(breaks=number_ticks(8)) + geom_hline(yintercept = 0) + theme_minimal()

 }
 marrangeGrob(plots[1:6], nrow = 2, ncol = 3, top = NULL)

#D plots UPDATED
  meaningful_NA_rows = lapply(1:Time, function(tp) {
    which(avail1[tp,]==0)
  })
Dout = matrix(NA, nrow = 32, ncol = 0)
Dnew = UN$DPS
for (i in 1:2000) {
		UDUPM = list()
		for (t in 1:32) {
			U = UN$UPS[[t]][,c(2*i-1, i*2)]
			UDUPM[[t]] = U %*% diag(UN$DPS[[t]][i,]) %*% t(U)
		}
		 eULU = lapply(1:Time, function(tp) {
    	exclude = meaningful_NA_rows[[tp]]
    	if (length(exclude) > 0) {
     	 eigentp = eigen(UDUPM[[tp]][-exclude, -exclude])
    	} else {
      	eigentp = eigen(UDUPM[[tp]])
    	}
    	eigentp
  		})
  		eR = lapply(1:Time, function(tp) {
    	which(rank(-abs(eULU[[tp]]$val), ties.method = "first") <= 2)
  		})
 	 	L =  lapply(1:Time, function(tp){
    	eULU[[tp]]$val[eR[[tp]]]
  		})
  		for (t in 1:32) {
  			Dnew[[t]][i,] = L[[t]]
  		}
}	


D = lapply(1:Time, function(t){summary(mcmc(Dnew[[t]]))[[2]]})
Ds = list()
for (i in 1:2) {Ds[[i]] = sapply(1:Time, function(t){D[[t]][i,]})}
betacols= ggplotColours(2)
plots = list()
years = c(1983:2014)
data = data.frame(cbind(years,t(Ds[[i]])))
 colnames(data)[4] = "D"
 f = ggplot(data, aes(x = years))
varname = c("D_1", "D_2")
for (i in 1:2){
years = c(1983:2014)
data = data.frame(cbind(years,t(Ds[[i]])))
 colnames(data)[4] = "D"
 f = ggplot(data, aes(x = years))
 plots[[i]] <- f + geom_line(aes(y = D), colour="black" )+ geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.1) + ylab(varname[i]) + scale_x_continuous(breaks=number_ticks(8)) + geom_hline(yintercept = 0) + theme_minimal()
#Dout = cbind(Dout, data)
 }
 marrangeGrob(plots[1:2], nrow = 1, ncol = 2, top = NULL)




#thetaplot
colors = sort(rownames(UN$U[[Time]]))
colors[which(colors == "GFR")] = "GMY"
rownames(UN$U[[Time]])[31] = "GMY"
	thetanew = t(sapply(1:Time, function(t){colMeans(UN$theta[[t]])}))
	colnames(thetanew) = colnames(avail1)
	orders = sapply(1:N, function(n){which(colors[n]== rownames(UN$U[[Time]]))})
	thetanew = thetanew[,orders]
	data3 = data.frame(years = years, theta = thetanew)
	colnames(data3)[-1] = colors
	data3new = melt(data3, id = "years")
	colnames(data3new)[3] = "theta"
	f <- ggplot(data3new, aes(years, theta, colour = variable, label = variable))
	f + geom_line(size = 0.2)+
geom_text(data = data3new[data3new$years=="1983", ], aes(label = variable), check_overlap = F, hjust = 1.3, vjust = 1, size =3, show.legend = F )+
geom_text(data = data3new[data3new$years=="2014", ], aes(label = variable), check_overlap = F, hjust = -0.3, vjust = 1, size =3, show.legend = F )+ scale_x_continuous(breaks=number_ticks(8)) + scale_colour_discrete(name = "countries") + theme_minimal()


#thetaplot_reduced
ggcolors = ggplotColours(21)
colors = sort(rownames(UN$U[[Time]])[which(rownames(UN$U[[Time]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"))])
data4new = data3new[data3new$variable %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"),]

f <- ggplot(data4new, aes(years, theta, color = variable, label = variable))
f + geom_line(size = 0.5, aes(linetype =as.factor(c(sapply(c(6:1,5:1,5:1,5:1), function(x){rep(x, Time)})))))+geom_point(size =1, aes(shape =as.factor(c(sapply(c(1:6,1:5,1:5,1:5), function(x){rep(x, Time)})))))+ scale_x_continuous(breaks=number_ticks(8)) + theme_minimal() + 
  geom_text(data = data4new[data4new$years=="1983", ], aes(label = variable), check_overlap = T, hjust = 1.3, vjust = 1, size =3, show.legend = F )+
geom_text(data = data4new[data4new$years=="2014", ], aes(label = variable), check_overlap = T, hjust = -0.3, vjust = 1, size =3, show.legend = F )+
  scale_color_manual(values = ggplotColours(21), name = "countries")+scale_linetype(guide = FALSE)+scale_shape(guide = FALSE)+
  guides(colour = guide_legend(override.aes = list(shape = c(1:6,1:5,1:5,1:5), linetype = c(6:1,5:1,5:1,5:1)))) 

f <- ggplot(data4new, aes(years, theta, color = variable, label = variable))
f + geom_line( aes(linetype =as.factor(c(sapply(c(6:1,5:1,5:1,5:1), function(x){rep(x, Time)})))))+ scale_x_continuous(breaks=number_ticks(8)) + theme_minimal() + 
  geom_text(data = data4new[data4new$years=="1983", ], aes(label = variable), check_overlap = F, hjust = 1.3, vjust = 1, size =3, show.legend = F )+
  geom_text(data = data4new[data4new$years=="1988", ], aes(label = variable), check_overlap = F, hjust = 0, vjust = 1.2, size =3, show.legend = F )+
  geom_text(data = data4new[data4new$years=="2009", ], aes(label = variable), check_overlap = F, hjust = 0, vjust = 1.2, size =3, show.legend = F )+
  geom_text(data = data4new[data4new$years=="2014", ], aes(label = variable), check_overlap = F, hjust = -0.3, vjust = 1, size =3, show.legend = F )+
  scale_color_manual(values = ggplotColours(21), name = "countries")+scale_linetype(guide = FALSE)+scale_shape(guide = FALSE)+
  guides(colour = guide_legend(override.aes = list(linetype = c(6:1,5:1,5:1,5:1)))) 


#UD plots_full
setwd('/Users/bomin8319/Desktop/UDU')
years = c(1983:2014)

Xstar = matrix(0, nrow = N, ncol = 2)
rownames(Xstar) = rownames(UN$U[[Time]])
rownames(Xstar)[which(rownames(Xstar) == "GFR")] = "GMY"
Xstar[1, ]= c(0,1)
Xstar[90, ] = c(1,0)
plots= list()
data2 = list()
colors = rownames(Xstar)
for (t in 1:Time) {
UDmat = matrix(NA, N, 2)
rownames(UDmat) = colors
rownames(UN$U[[t]])[which(rownames(UN$U[[t]]) == "GFR")] = "GMY"
UDmat[which(rownames(UDmat) %in% names(UN$U[[t]][,1] * UN$D[[t]][1] )),1] =  UN$U[[t]][,1] * sqrt(UN$D[[t]][1]) 
UDmat[which(rownames(UDmat) %in% names(UN$U[[t]][,2] * UN$D[[t]][2] )),2] =  UN$U[[t]][,2] * sqrt(abs(UN$D[[t]][2]))
UDmat[which(rownames(UDmat) %in% names(UN$U[[t]][,2] * UN$D[[t]][2] )),] = procrustes(UDmat[which(rownames(UDmat) %in% rownames(UN$U[[t]])),], Xstar[which(rownames(Xstar) %in% rownames(UN$U[[t]])),])$X.new
data2[[t]] = data.frame(UDmat)
colnames(data2[[t]])[1:2] = c("r1", "r2")
Xstar = UDmat
Xstar[is.na(Xstar)] = rep(0, 2)
}

for (t in 1:Time) {
colors.t = colors[which(colors %in% rownames(UN$U[[t]]))]
p <- ggplot(data2[[t]], aes(x = r1, y = r2, label = rownames(data2[[t]])))

plots[[t]] = p+ geom_text(colour = "black", size = 5, show.legend = F, check_overlap = F)+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks=number_ticks(3)) + scale_y_continuous(breaks=number_ticks(3))

	if (UN$D[[t]][2] < 0) {
		plots[[t]] = plots[[t]] + geom_vline(colour = 'red', aes(xintercept = 0))
	}

mname = paste0("plot", t, "full.png")
print(plots[[t]])
ggsave(filename = mname)
}

for (t in c(4,8,12,16,20,24,28,Time)-2) {
colors.t = colors[which(colors %in% rownames(UN$U[[t]]))]
p <- ggplot(data2[[t]], aes(x = r1, y = r2, label = rownames(data2[[t]])))+ labs(x = "r = 1", y = "r = 2")

plots[[t]] = p+ geom_text(colour = "black", size = 2, show.legend = F, check_overlap =FALSE)+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))  + scale_x_continuous(breaks=number_ticks(3)) + scale_y_continuous(breaks=number_ticks(3))
	if (UN$D[[t]][2] < 0) {
		plots[[t]] = plots[[t]] + geom_vline(colour = 'red', aes(xintercept = 0))
	}
}
marrangeGrob(plots[c(4,20, 8, 24, 12, 28, 16,Time)-2], nrow = 2, ncol = 4, top = NULL)


#UD plots_reduced
data3 = list()
for (t in 1:Time) {
    data3[[t]] = data2[[t]][rownames(data2[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"),]
}
plots= list()
colors = rownames(data3[[1]])
for (t in 1:Time) {
    colors.t = colors[which(colors %in% rownames(UN$U[[t]]))]
    p <- ggplot(data3[[t]], aes(x = r1, y = r2,  label = rownames(data3[[t]])))
    
    plots[[t]] = p+ geom_text(colour = "black", size = 5, show.legend = F)+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks=number_ticks(3)) + scale_y_continuous(breaks=number_ticks(3))
    if (UN$D[[t]][2] < 0) {
		plots[[t]] = plots[[t]] + geom_vline(colour = 'red', aes(xintercept = 0))
	}
    mname = paste0("plot", t, "reduced.png")
   print(plots[[t]])
   ggsave(filename = mname)
}

for (t in c(4,8,12,16,20,24,28,Time)-2) {
colors.t = colors[which(colors %in% rownames(UN$U[[t]]))]
p <- ggplot(data3[[t]], aes(x = r1, y = r2, label = rownames(data3[[t]]))) + labs(x = "r = 1", y = "r = 2")

plots[[t]] = p+ geom_text(colour = "black", size = 3, show.legend = F)+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))  + scale_x_continuous(breaks=number_ticks(3)) + scale_y_continuous(breaks=number_ticks(3))
  if (UN$D[[t]][2] < 0) {
		plots[[t]] = plots[[t]] + geom_vline(colour = 'red', aes(xintercept = 0))
	}

}
marrangeGrob(plots[c(4,20, 8, 24, 12, 28, 16,Time)-2], nrow = 2, ncol = 4, top = NULL)





##Posterior predictive of degree
setwd('/Users/bomin8319/Desktop/external')
ggcolors = ggplotColours(3)
datacollapse = matrix(0, nrow = 84000, ncol = 3)
observedcollapse = matrix(0, 28, 3)
pp = list()
for (tp in 1:Time) {
  diag(Y[tp,,]) = 0
  Y[tp, which(avail1[tp, ]==0), ] = 0
  Y[tp, , which(avail1[tp, ]==0)] = 0
  Y[tp, , ][which(is.na(Y[tp, , ]))] = 0
  data = t(sapply(1:1000, function(r){tabulate(round(UN$Degree[[tp]][r,]/100), 95)}))[,-c(1:59, 88:95)]
  colnames(data) =  c(60:87)
  datamat = matrix(0, 1000, 28)
  colnames(datamat) = c(60:87)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat1 = data.frame(Proportion = c(datamat), Degree = factor(c(sapply(c(60:87), function(k){rep(k, 1000)}))), Model = as.factor("DAME"))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN4$Degree[[tp]][r,]/100), 95)}))[,-c(1:59, 88:95)]
  colnames(data) =  c(60:87)
  datamat = matrix(0, 1000, 28)
  colnames(datamat) = c(60:87)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat2 = rbind(datamat1, data.frame(Proportion = c(datamat), Degree = factor(c(sapply(c(60:87), function(k){rep(k, 1000)}))), Model = as.factor("ME")))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN2$Degree[[tp]][r,]/100), 95)}))[,-c(1:59, 88:95)]
  colnames(data) =  c(60:87)
  datamat = matrix(0, 1000,28)
  colnames(datamat) = c(60:87)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat3 = rbind(datamat2, data.frame(Proportion = c(datamat), Degree = factor(c(sapply(c(60:87), function(k){rep(k, 1000)}))), Model = as.factor("AE")))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN3$Degree[[tp]][r,]/100), 95)}))[,-c(1:59, 88:95)]
  colnames(data) =  c(60:87)
  datamat = matrix(0, 1000, 28)
  colnames(datamat) =c(60:87)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat4 = rbind(datamat3, data.frame(Proportion = c(datamat), Degree = factor(c(sapply(c(60:87), function(k){rep(k, 1000)}))), Model = as.factor("NO")))
  datamat4 = datamat3
  datamat4 = as.data.frame(datamat4)
  observedpp = as.numeric(tabulate(round(rowSums(Y[tp,,])), 95) / sum(tabulate(round(rowSums(Y[tp,,])), 95)))[-c(1:59, 88:95)]
  names(observedpp) = c(60:87)
  pvec = rep(0, 28)
  names(pvec) =c(60:87)
  pvec[which(names(pvec) %in% names(observedpp))] = observedpp
  datamat4$Degree = as.factor(as.numeric(as.character(datamat4$Degree))*100)
  observed = data.frame(Proportion = pvec, Degree = as.factor(c(60:87)*100), Model = "DAME")
  pp[[tp]] = ggplot(datamat4, aes(x = Degree, y = Proportion,fill = Model, color =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observed, color = "blue", size = 0.2, group = 1)+geom_point(data =observed, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
  datacollapse[,1] = datacollapse[,1] + datamat3[,1]
  observedcollapse[,1] = observedcollapse[,1] + observed[,1]
}
years = c(1983:2014)
for (t in 1:Time){
  mname = paste0(years[t], "overalldegree", ".png")
  print(pp[[t]])
  ggsave(filename = mname, width = 11, height = 6)
}
datacollapse = data.frame(Proportion = datacollapse[,1]/Time, Degree = as.factor(as.numeric(as.character(datamat3$Degree))*100), Model = datamat3$Model)
observedcollapse = data.frame(Proportion = observedcollapse[,1]/Time, Degree = observed$Degree, Model = observed$Model)
datacollapse2 = datacollapse[datacollapse$Model == "DAME",]
#pp = ggplot(datacollapse, aes(x = Degree, y = Proportion,fill = Model, color =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observedcollapse, color = "blue", size = 0.2, group = 1)+geom_point(data =observedcollapse, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
pp = ggplot(datacollapse2, aes(x = Degree, y = Proportion)) + geom_boxplot(colour = "black",outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observedcollapse, color = "blue", size = 0.2, group = 1)+geom_point(data =observedcollapse, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
pp


#seconddegree
ggcolors = ggplotColours(3)
datacollapse = matrix(0, nrow = 2*37500, ncol = 3)
observedcollapse = matrix(0, 25, 3)
pp = list()
for (tp in 1:Time) {
  diag(Y[tp,,]) = 0
  Y[tp, which(avail1[tp, ]==0), ] = 0
  Y[tp, , which(avail1[tp, ]==0)] = 0
  Y[tp, , ][which(is.na(Y[tp, , ]))] = 0
  data = t(sapply(1:1000, function(r){tabulate(round(UN$secondDegree[[tp]][r,]/(N*10000)), 77)}))[,-c(1:44, 70:77)]
  colnames(data) =  c(45:69)
  datamat = matrix(0, 1000, 25)
  colnames(datamat) = c(45:69)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat1 = data.frame(Proportion = c(datamat), secondDegree = factor(c(sapply(c(45:69), function(k){rep(k, 1000)}))), Model = as.factor("DAME"))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN4$secondDegree[[tp]][r,]/(N*10000)), 77)}))[,-c(1:44, 70:77)]
  colnames(data) =  c(45:69)
  datamat = matrix(0, 1000, 25)
  colnames(datamat) = c(45:69)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat2 = rbind(datamat1, data.frame(Proportion = c(datamat), secondDegree = factor(c(sapply(c(45:69), function(k){rep(k, 1000)}))), Model = as.factor("ME")))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN2$secondDegree[[tp]][r,]/(N*10000)), 77)}))[,-c(1:44, 70:77)]
  colnames(data) = c(45:69)
  datamat = matrix(0, 1000, 25)
  colnames(datamat) = c(45:69)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat3 = rbind(datamat2, data.frame(Proportion = c(datamat), secondDegree = factor(c(sapply(c(45:69), function(k){rep(k, 1000)}))), Model = as.factor("AE")))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN3$secondDegree[[tp]][r,]/(N*10000)), 77)}))[,-c(1:44, 70:77)]
  colnames(data) =  c(45:69)
  datamat = matrix(0, 1000, 25)
  colnames(datamat) =c(45:69)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat4 = rbind(datamat3, data.frame(Proportion = c(datamat), secondDegree = factor(c(sapply(c(45:69), function(k){rep(k, 1000)}))), Model = as.factor("NO")))
  datamat4 = datamat3
  datamat4 = as.data.frame(datamat4)
  datamat4$secondDegree = as.factor(as.numeric(as.character(datamat4$secondDegree))*10000)
  observedpp = as.numeric(tabulate(round(rowSums(Y[tp,,]%*% Y[tp,,]/(N))), 77) / sum(tabulate(round(rowSums(Y[tp,,] %*% Y[tp,,] /N)), 77)))[-c(1:44, 70:77)]
  names(observedpp) = c(45:69)
  pvec = rep(0,25)
  names(pvec) =c(45:69)
  pvec[which(names(pvec) %in% names(observedpp))] = observedpp
  observed = data.frame(Proportion = pvec, secondDegree = as.factor(c(45:69)*10000), Model = "DAME")
  datamat4$secondDegree =as.factor(paste0(substr(datamat4$secondDegree, 1, 3), ""))
  observed$secondDegree =as.factor(paste0(substr(observed$secondDegree, 1, 3), ""))
  pp[[tp]] = ggplot(datamat4, aes(x = secondDegree, y = Proportion,fill = Model, color = Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observed, color = "blue", size = 0.2, group = 1)+geom_point(data =observed, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
  datacollapse[,1] = datacollapse[,1] + datamat4[,1]
  observedcollapse[,1] = observedcollapse[,1] + observed[,1]
}

for (t in 1:Time){
  mname = paste0(years[t], "overallseconddegree", ".png")
  print(pp[[t]]+labs(x = "secondDegree (in thousand)"))
  ggsave(filename = mname, width = 11, height = 6)
}

datacollapse = data.frame(Proportion = datacollapse[,1]/Time, secondDegree = as.factor(as.numeric(as.character(datamat3$secondDegree))*10000), Model = datamat3$Model)
observedcollapse = data.frame(Proportion = observedcollapse[,1]/Time, secondDegree = observed$secondDegree, Model = observed$Model)
datacollapse2 = datacollapse[datacollapse$Model == "DAME",]
datacollapse$secondDegree =as.factor(paste0(substr(datacollapse$secondDegree, 1, 3), ""))
datacollapse2$secondDegree =as.factor(paste0(substr(datacollapse2$secondDegree, 1, 3), ""))
observedcollapse$secondDegree =as.factor(paste0(substr(observedcollapse$secondDegree, 1, 3), ""))
#pp = ggplot(datacollapse, aes(x = secondDegree, y = Proportion,fill = Model, color =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observedcollapse, color = "blue", size = 0.2, group = 1)+geom_point(data =observedcollapse, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
pp = ggplot(datacollapse2, aes(x = secondDegree, y = Proportion)) + geom_boxplot(colour = "black", outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observedcollapse, color = "blue", size = 0.2, group = 1)+geom_point(data =observedcollapse, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
pp + labs(x = "secondDegree (in thousand)")




#thirddegree
ggcolors = ggplotColours(3)
datacollapse = matrix(0, nrow = 2*40500, ncol = 3)
observedcollapse = matrix(0, 27, 3)

pp = list()
for (tp in 1:Time) {
  diag(Y[tp,,]) = 0
  Y[tp, which(avail1[tp, ]==0), ] = 0
  Y[tp, , which(avail1[tp, ]==0)] = 0
  Y[tp, , ][which(is.na(Y[tp, , ]))] = 0
  data = t(sapply(1:1000, function(r){tabulate(round(UN$thirdDegree[[tp]][r,]/(N*N*1000000)), 66)}))[,-c(1:31, 59:66)]
  colnames(data) = c(32:58)
  datamat = matrix(0, 1000, 27)
  colnames(datamat) = c(32:58)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat1 = data.frame(Proportion = c(datamat), thirdDegree = factor(c(sapply(c(32:58), function(k){rep(k, 1000)}))), Model = as.factor("DAME"))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN4$thirdDegree[[tp]][r,]/(N*N*1000000)), 66)}))[,-c(1:31, 59:66)]
  colnames(data) =  c(32:58)
  datamat = matrix(0, 1000, 27)
  colnames(datamat) = c(32:58)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat2 = rbind(datamat1, data.frame(Proportion = c(datamat), thirdDegree = factor(c(sapply(c(32:58), function(k){rep(k, 1000)}))), Model = as.factor("ME")))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN2$thirdDegree[[tp]][r,]/(N*N*1000000)), 66)}))[,-c(1:31, 59:66)]
  colnames(data) = c(32:58)
  datamat = matrix(0, 1000, 27)
  colnames(datamat) = c(32:58)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat3 = rbind(datamat2, data.frame(Proportion = c(datamat), thirdDegree = factor(c(sapply(c(32:58), function(k){rep(k, 1000)}))), Model = as.factor("AE")))
  
  data =  t(sapply(1:1000, function(r){tabulate(round(UN3$thirdDegree[[tp]][r,]/(N*N*1000000)), 66)}))[,-c(1:31, 59:66)]
  colnames(data) =  c(32:58)
  datamat = matrix(0, 1000, 27)
  colnames(datamat) =c(32:58)
  for (i in 1:1000) {
    datamat[i, which(colnames(datamat) %in% colnames(data))] = data[i,] / sum(data[i,])
  }
  datamat4 = rbind(datamat3, data.frame(Proportion = c(datamat), thirdDegree = factor(c(sapply(c(32:58), function(k){rep(k, 1000)}))), Model = as.factor("NO")))
  datamat4 = datamat3
  datamat4 = as.data.frame(datamat4)
  datamat4$thirdDegree = as.factor(as.numeric(as.character(datamat4$thirdDegree))*1000000)
  observedpp = as.numeric(tabulate(round(rowSums(Y[tp,,]%*% Y[tp,,]%*% Y[tp,,]/(N*N))), 66) / sum(tabulate(round(rowSums(Y[tp,,] %*% Y[tp,,]%*% Y[tp,,] /(N*N))), 66)))[-c(1:31, 59:66)]
  names(observedpp) = c(32:58)
  pvec = rep(0, 27)
  names(pvec) =c(32:58)
  pvec[which(names(pvec) %in% names(observedpp))] = observedpp
  observed = data.frame(Proportion = pvec, thirdDegree = as.factor(c(32:58)*1000000), Model = "DAME")
  datamat4$thirdDegree =as.factor(paste0(substr(datamat4$thirdDegree, 1, 3), ""))
  observed$thirdDegree =as.factor(paste0(substr(observed$thirdDegree, 1, 3), ""))
  
  pp[[tp]] = ggplot(datamat4, aes(x = thirdDegree, y = Proportion,fill = Model, color =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observed, color = "blue", size = 0.2, group = 1)+geom_point(data =observed, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
  datacollapse[,1] = datacollapse[,1] + datamat4[,1]
  observedcollapse[,1] = observedcollapse[,1] + observed[,1]
}
years = c(1983:2014)
for (t in 1:Time){
  mname = paste0(years[t], "overallthirddegree", ".png")
  print(pp[[t]]+labs(x = "thirdDegree (in million)"))
  ggsave(filename = mname, width = 11, height = 6)
}

datacollapse = data.frame(Proportion = datacollapse[,1]/Time, thirdDegree =as.factor(as.numeric(as.character(datamat3$thirdDegree))*1000000), Model = datamat3$Model)
datacollapse2 = datacollapse[datacollapse$Model == "DAME",]
observedcollapse = data.frame(Proportion = observedcollapse[,1]/Time, thirdDegree = observed$thirdDegree, Model = observed$Model)
datacollapse$thirdDegree =as.factor(paste0(substr(datacollapse$thirdDegree, 1, 3), ""))
datacollapse2$thirdDegree =as.factor(paste0(substr(datacollapse2$thirdDegree, 1, 3), ""))
observedcollapse$thirdDegree =as.factor(paste0(substr(observedcollapse$thirdDegree, 1, 3), ""))
#pp = ggplot(datacollapse, aes(x = thirdDegree, y = Proportion,fill = Model, color =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observedcollapse, color = "blue", size = 0.2, group = 1)+geom_point(data =observedcollapse, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
pp = ggplot(datacollapse2, aes(x = thirdDegree, y = Proportion)) + geom_boxplot(colour = "black", outlier.size = 0.5, position = position_dodge()) + theme_minimal()+scale_fill_manual(values = alpha(ggcolors, 0.5)) + geom_line(data = observedcollapse, color = "blue", size = 0.2, group = 1)+geom_point(data =observedcollapse, color = "blue", size =2, group = 1)+guides(colour = guide_legend(override.aes = list(shape = NA))) + theme(legend.title = element_blank())
pp+labs(x = "thirdDegree (in million)")
