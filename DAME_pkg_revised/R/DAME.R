#' @useDynLib DAME
#' @importFrom stats cor lm median rgamma rnorm runif
#' @importFrom fields Exponential
#' @importFrom FastGP rcppeigen_invert_matrix rcpp_rmvnorm rcppeigen_get_chol rcpp_log_dmvnorm
#' @importFrom matrixcalc is.positive.definite
#' @importFrom Matrix nearPD
#' @importFrom LaplacesDemon dhalfcauchy rhalfcauchy dinvgamma

NULL

#' @title rs2_fc2
#' @description Sample the variance in Normal error
#'
#' @param XB Time-length list of XB matrix
#' @param theta a T x n matrix of nodal effects
#' @param UDU Time-length list of UDU
#' @param Y relational array of relations
#' @param a inverse-gamma shape parameter
#' @param b inverse-gamma scale parameter
#'
#' @return One sample estimate of sigma^2 from inverse-gamma distribution
#'
#' @export
rs2_fc2 = function(XB, theta, UDU, Y, a, b) {
    Time = nrow(theta)
    N = ncol(theta)
    uppertri = upper.tri(diag(N))
    E = lapply(1:Time, function(tp) {
        Y[tp, , ] - XB[[tp]] - outer(theta[tp,], theta[tp,], "+") - UDU[[tp]]
    })
    sums = sum(vapply(E, function(x){sum(x[uppertri]^2, na.rm = TRUE)}, c(1)))
    missings = sum(vapply(E, function(x){sum(is.na(x[uppertri]))}, c(1)))
    samp = rgamma(1, a + (Time * N * (N-1) / 2 - missings) / 2, b + sums / 2)
    return(1 / samp)
}

#' @title rtaup_fc
#' @description Sample the variance in beta
#'
#' @param beta a p by 1 vector of coefficients
#' @param cinv inverse of Gaussian process covariance matrix
#' @param a inverse-gamma shape parameter
#' @param b inverse-gamma scale parameter
#'
#' @return One sample estimate of tau_p from inverse-gamma distribution
#'
#' @export
rtaup_fc = function(beta, cinv, a, b) {
    Time = nrow(beta)
    P = ncol(beta)
    samp = c()
    for (p in 1:P) {
        sums = t(beta[, p]) %*% cinv[[p]] %*% beta[, p]
        samp[p] = rgamma(1, a + Time / 2, b + sums / 2)
    }
    return(1 / samp)
}

#' @title rtaui_fc
#' @description Sample the variance in theta
#'
#' @param theta a Time by n vector of nodal effects
#' @param cinv inverse of Gaussian process covariance matrix
#' @param a inverse-gamma shape parameter
#' @param b inverse-gamma scale parameter
#'
#' @return One sample estimate of tau_i from inverse-gamma distribution
#'
#' @export
rtaui_fc = function(theta, cinv, a, b) {
    Time = nrow(theta)
    N = ncol(theta)
    sums = sum(vapply(1:N, function(i) {t(theta[, i]) %*% cinv %*% theta[, i]}, c(1)))
    samp = rgamma(1, a + Time * N / 2, b + sums / 2)
    return(1 / samp)
}


#' @title rtaur_fc
#' @description Sample the variance in d
#'
#' @param d current value of d
#' @param a inverse-gamma shape parameter
#' @param b inverse-gamma scale parameter
#'
#' @return One sample estimate of tau_u from inverse-gamma distribution
#'
#' @export
rtaur_fc = function(d, a, b) {
    Time = dim(d)[1]
    R = dim(d)[2]
    sums = vapply(1:R, function(r) {t(d[, r]) %*% diag(Time) %*% d[, r]}, c(1))
    samp = rgamma(R, a + Time * R / 2, b + sums / 2)
    return(1 / samp)
}

#' @title rtauu_fc
#' @description Sample the variance in d
#'
#' @param U current value of U
#' @param cinv inverse of Gaussian process covariance matrix
#' @param a inverse-gamma shape parameter
#' @param b inverse-gamma scale parameter
#'
#' @return One sample estimate of tau_u from inverse-gamma distribution
#'
#' @export
rtauu_fc = function(U, cinv, a, b) {
    Time = dim(U)[1]
    N = dim(U)[2]
    R = dim(U)[3]
    sums = rep(0, R)
    for (n in 1:N) {
    	if (sum(is.na(U[,n,]))==0) {
    		add = vapply(1:R, function(r) {t(U[, n, r]) %*% cinv[[r]] %*% U[, n, r]}, c(1))
    	} else {
    		delete = which(is.na(U[,n,1]))
    		add = vapply(1:R, function(r) {t(U[-delete, n, r]) %*% cinv[[r]][-delete, -delete] %*% U[-delete, n, r]}, c(1))
    	}
    	sums = sums + add
    }
    missings = sum(is.na(c(U[,,1]^2)))
    samp = rgamma(R, a + (N * Time - missings) / 2, b + sums / 2)
    return(1/samp)
}

#' @title rbeta_fc
#' @description Sample beta
#'
#' @param X Time x n x n x p covariate array
#' @param beta a p by 1 vector of coefficients
#' @param theta a T x n matrix of nodal effects
#' @param UDU Time-length list of UDU
#' @param Y relational array of relations
#' @param cinv inverse of Gaussian process covariance matrix
#' @param tau_p scale parameter of Gaussian process covariance
#' @param s2 variance of normal error
#'
#' @return One sample estimate of beta from normal distribution
#'
#' @export
rbeta_fc = function(X, beta, theta, UDU, Y, cinv, tau_p, s2) {
    Time = nrow(theta)
    P = ncol(beta)
    N = ncol(theta)
    uppertri = upper.tri(diag(N))
    Xp2 = c()
    E_p = c()
    beta = beta
    E = lapply(1:Time, function(tp) {
        Y[tp, , ] - outer(theta[tp,], theta[tp,], "+") - UDU[[tp]]
    })
    for (p in sample(1:P)){
        for (tp in 1:Time) {
            if (P == 1) {
                XB_p = matrix(0, N, N)
            } else {
                XB_p = Reduce('+', lapply((1:P)[-p], function(k) {
                    X[tp, , , k] * beta[tp, k]
                }))
            }
            Xp2[tp] = sum(X[tp, , , p][uppertri]^2, na.rm = TRUE)
            E_p[tp] = sum((E[[tp]] - XB_p)[uppertri] * X[tp, , , p][uppertri], na.rm = TRUE)
        }
        invvar = rcppeigen_invert_matrix(cinv[[p]] / tau_p[p] + diag(Xp2, Time, Time) / s2)
        mu = E_p / s2
        beta[, p] = rcpp_rmvnorm(1, invvar, invvar %*% mu)
    }
    return(beta)
}

#' @title rtheta_fc
#' @description Sample theta with structural NA
#'
#' @param XB Time-length list of XB matrix
#' @param theta a T x n matrix of nodal effects
#' @param UDU Time-length list of UDU
#' @param Y relational array of relations
#' @param cinv inverse of Gaussian process covariance matrix
#' @param tau_i scale parameter of Gaussian process covariance
#' @param s2 variance of normal error
#'
#' @return One sample estimate of theta from normal distribution
#'
#' @export
rtheta_fc = function(XB, theta, UDU, Y, cinv, tau_i, s2){
    Time = nrow(theta)
    N = ncol(theta)
    
    theta = theta
    E = lapply(1:Time, function(tp) {
        Y[tp, , ] - XB[[tp]] - UDU[[tp]]
    })
    
    for (i in sample(1:N)) {
        E_i = vapply(1:Time, function(tp) {
            sum(E[[tp]][i, -i] - theta[tp, -i], na.rm = TRUE)
        }, c(1))
        
        invvar = rcppeigen_invert_matrix(cinv / tau_i + (N-1) * diag(Time) / s2)
        mu = E_i / s2
        theta[, i] = rcpp_rmvnorm(1, invvar, invvar %*% mu)
    }
    return(theta)
}

#' @title rd_fc2_revised
#' @description Sample d with structural NA
#'
#' @param XB Time-length list of XB matrix
#' @param theta a T x n matrix of nodal effects
#' @param U current value of U
#' @param d current value of d
#' @param Y relational array of relations
#' @param s2 variance of normal error
#' @param meaningful_NA position of structural NA
#'
#' @return One sample estimate of d from normal distribution
#'
#' @export
rd_fc2_revised = function(XB, theta, U, d, Y, s2, meaningful_NA){
  Time = nrow(theta)
  R = ncol(d)
  N = ncol(theta)
  uppertri = upper.tri(diag(N))
  UUr2 = c()
  E_r2 = c()
  d = d
  E = lapply(1:Time, function(tp) {
    Y[tp, , ] - XB[[tp]] - outer(theta[tp,], theta[tp,], "+")
  })
  
  for (r in 1:R) {
    for (tp in 1:Time) {
      UUr = U[tp, , r] %*% t(U[tp, , r])
      if (R > 2) {
        UDU_r = U[tp, , -r] %*% diag(d[tp, -r]) %*% t(U[tp, , -r])
      } else {
        if (R == 2) {
          UDU_r = d[tp, -r] * U[tp, , -r] %*% t(U[tp, , -r])
        } else {
          UDU_r =  matrix(0, N, N)
        }
      }
      #UUr[meaningful_NA[[tp]]] = NA
      EU =  (E[[tp]] - UDU_r) * UUr
      UUr2[tp] = sum(UUr[uppertri]^2, na.rm = TRUE)
      E_r2[tp] = sum(EU[uppertri], na.rm = TRUE)
    }
    invvar = rcppeigen_invert_matrix(diag(Time) + diag(UUr2, Time, Time) / s2)
    mu = E_r2 / s2
    d[, r] = rcpp_rmvnorm(1, invvar, invvar %*% mu)
  }
  return(d)
}



#' @title ru_fc_NA_revised
#' @description Sample u with structural NA
#'
#' @param XB Time-length list of XB matrix
#' @param theta a T x n matrix of nodal effects
#' @param U current value of U
#' @param d current value of d
#' @param Y relational array of relations
#' @param cinv inverse of Gaussian process covariance matrix
#' @param tau_u Gaussian process covariance
#' @param s2 variance of normal errorU
#' @param meaningful_NA_rows position of structural NA
#'
#' @return One sample estimate of d from normal distribution
#'
#' @export
ru_fc_NA_revised = function(XB, theta, U, d, Y, cinv, tau_u, s2, meaningful_NA_rows){
    Time = nrow(theta)
    R = ncol(d)
    N = ncol(theta)
    U = U
    E = lapply(1:Time, function(tp) {
        Y[tp, , ] - XB[[tp]] - outer(theta[tp,], theta[tp,], "+")
    })
    	for (r in 1:R){
    		UDU_r = list()
    		for (tp in 1:Time) {
			if (R > 2) {
        		UDU_r[[tp]] = U[tp, , -r] %*% diag(d[tp, -r]) %*% t(U[tp, , -r])
      		} else {
        	if (R == 2) {
           UDU_r[[tp]] = d[tp, -r] * U[tp, , -r] %*% t(U[tp, , -r])
        	} else {
          		UDU_r[[tp]] = matrix(0, N, N)
        	}
      		}
      	}
    	for (i in rep(sample(1:N, N), 4)) {
    		A_nr = matrix(0, Time, Time)
    	 	b_nr = rep(0, Time)
    	 	for (tp in 1:Time) {
    	 		A_nr[tp, tp] = d[tp, r]^2 * sum(U[tp, -i, r]^2, na.rm = TRUE)
    	 		b_nr[tp] = d[tp, r]* sum((E[[tp]] - UDU_r[[tp]])[i, -i] * U[tp, -i, r], na.rm = TRUE)
		}	
		A = 1/2 * A_nr/s2 + cinv[[r]] / tau_u[r]
		A_inv = rcppeigen_invert_matrix(A)
		b = 1/2 * b_nr/s2
		U[, i, r] = rcpp_rmvnorm(1, A_inv, A_inv %*% b)	
    	}
    }
    return(U)
}


#' @title DAME_fix_revised
#' @description MCMC algorithm using Gibbs sampling for each variable with structural NA
#'
#' @param Y relational array of relations
#' @param X Time x n x n x p covariate array
#' @param RE random effect to be included ("additive" and/or "multiplicative")
#' @param R dimension of the multiplicative effects
#' @param dist standard Exponential or squared Exponential
#' @param gammapriors inverse-gamma shape and scale parameters for (s2, beta, theta, d)
#' @param avail Time x n matrix reperesenting the availibility of nodes (1 if avail, 0 if structural NA)
#' @param burn burn in for the Markov chain
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param odens output density for the Markov chain
#' @param kappas P+2 length vector of GP length parameters
#'
#' @return Final estimate of the parameters
#'
#' @export
DAME_fix_revised = function(Y, X, RE = c("additive", "multiplicative"), R = 2, dist = "Exponential", gammapriors = c(2, 1), avail = matrix(1, dim(Y)[1], dim(Y)[2]), burn = 1000, nscan = 5000, odens = 100, kappas) {
  Time = dim(Y)[1]
  N = dim(Y)[2]
  P = dim(X)[4]
  a = gammapriors[1]; b = gammapriors[2]
  # construct covariance matrix for each variable
  dist_ij = c()
  if (dist == "Exponential"){
    for (i in 1:Time) {
      for (j in 1:Time) {
        dist_ij = c(dist_ij, abs(i-j))
      }
    }
  } else {
    for (i in 1:Time) {
      for (j in 1:Time) {
        dist_ij = c(dist_ij, (i-j)^2)
      }
    }
  }
  dist_ij = matrix(dist_ij, nrow = Time, ncol = Time)
  cinv = lapply(1:(P+1+R), function(k) {
    rcppeigen_invert_matrix(Exponential(dist_ij, kappas[k]))
  })
  # select initial values
  beta = matrix(0, Time, P)
  d = matrix(0, Time, R)
  U = array(0, dim = c(Time, N, R))
  theta = matrix(0, Time, N)
  BETAPS = lapply(1:Time, function(tp) {
    matrix(0, nrow = nscan / odens, ncol = P)
  })
  UDUPS = lapply(1:Time, function(tp) {
    matrix(0, N, N)
  })
  YPSsum = lapply(1:Time, function(tp) {
    matrix(0, N, N)
  })
  s2PS = matrix(0, nrow = nscan / odens, ncol = 1)
  tauPS = matrix(0, nrow = nscan / odens, ncol = P + 1 + R)
  thetaPS = lapply(1:Time, function(tp) {
    matrix(0, nrow = nscan / odens, ncol = N)
  })
  DPS = lapply(1:Time, function(tp) {
    matrix(0, nrow = nscan / odens, ncol = R)
  })
  UPS = lapply(1:Time, function(tp) {
   c()
  })
  Degreestats = lapply(1:Time, function(tp) {
    matrix(0, nrow = 0, ncol = N)
  })
  secondDegreestats = lapply(1:Time, function(tp) {
    matrix(0, nrow = 0, ncol = N)
  })
  thirdDegreestats = lapply(1:Time, function(tp) {
    matrix(0, nrow = 0, ncol = N)
  })
  corrstats = matrix(0, nrow = nscan / odens, ncol = Time - 1)
  years = sample(1:Time, 1)
  # to begin with, use 0 for NA's except meaningful ones
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
  meaningful_NA_years = which(vapply(meaningful_NA_rows, function(i) {length(i)}, c(1)) > 0)
  # to begin with, use 0 for NA's
  for (p in 1:P) {
    X[, , , p][which(is.na(X[, , , p]))] = 0
    for (tp in 1:Time) {
      X[tp, , , p][meaningful_NA[[tp]]] = NA
      Y[tp, , ][meaningful_NA[[tp]]] = NA
    }
  }
  XB = lapply(1:Time, function(tp) {
    Reduce('+', lapply(1:P, function(p){
      X[tp, , , p] * beta[tp, p]
    }))
  })
  na.positions = lapply(1:Time, function(tp) {
    which(is.na(Y[tp, , ]))
  })
  for (tp in 1:Time) {
    mu = mean(Y[tp, , ], na.rm = TRUE)
    row = rowMeans(Y[tp, , ] - mu, na.rm = TRUE)
    row[is.na(row)] = 0
    YA = mu + outer(row, row, "+")
    Y[tp, , ][na.positions[[tp]]] = YA[na.positions[[tp]]]
    diag(Y[tp, , ]) = 0
    beta[tp, 1] = mu
    if ("additive" %in% RE) {
    	theta[tp, ] = row
    	}
    	if ("multiplicative" %in% RE) {
    eigenError = eigen(Y[tp, , ] - YA)
    eR = which(rank(-abs(eigenError$val), ties.method = "first") <= R)
    signs = eigenError$val[eR]
    d[tp, ] = ifelse(signs > 0, 1, -1)
    }
    Y[tp, , ][meaningful_NA[[tp]]] = NA
  }
  
  UDU = UDUPS
  uppertri = upper.tri(diag(N))
  tau_u = rep(1, R)

  # starting the Gibbs sampler
  for (iter in 1:(burn + nscan)) {
    if (iter %% 500 == 0) print(iter)
    s2 = rs2_fc2(XB, theta, UDU, Y, a, b)
    tau_p = rtaup_fc(beta, cinv[1:P], a, b)
    beta = rbeta_fc(X, beta, theta, UDU, Y, cinv[1:P], tau_p, s2)
    XB = lapply(1:Time, function(tp) {
      Reduce('+', lapply(1:P, function(p){
        X[tp, , , p] * beta[tp, p]
      }))
    })
    if ("additive" %in% RE) {
    tau_i = rtaui_fc(theta, cinv[[P+1]], a, b)
    theta = rtheta_fc(XB, theta, UDU, Y, cinv[[P+1]], tau_i, s2)
    for (tp in meaningful_NA_years) {
      theta[tp, meaningful_NA_rows[[tp]]] = 0
    }
    }
    if ("multiplicative" %in% RE) {
    if (iter > 0.2 * burn) {
    		tau_u = rtauu_fc(U, cinv[(P+2):(P+1+R)], a, b)
       	d = rd_fc2_revised(XB, theta, U, d, Y, s2, meaningful_NA)
   		U = ru_fc_NA_revised(XB, theta, U, d, Y, cinv[(P+2):(P+1+R)], tau_u, s2, meaningful_NA_rows)
   		for (tp in meaningful_NA_years) {
      		U[tp, meaningful_NA_rows[[tp]],] = NA
   		}
    	UDU = lapply(1:Time, function(tp) {
      	if (R <= 1) {
        	d[tp, ] * U[tp, , ] %*% t(U[tp, ,])
      	} else {
        	U[tp, , ] %*% diag(d[tp, ]) %*% t(U[tp, , ])
      	}
    	})
    }
    }
    if (iter > burn & (iter - burn) %% odens == 0) {
      Errormat = array(0, dim = c(Time, N, N))
      errors = rnorm(Time * N * (N-1) / 2, 0, sqrt(s2))
      for (tp in 1:Time) {
        Errormat[tp, , ][upper.tri(Errormat[tp, , ])] = errors[((tp-1)*N*(N-1)/2+1):(tp*N*(N-1)/2)]
        Errormat[tp, , ] = (Errormat[tp, , ] + t(Errormat[tp, , ]))
      }
      YPS = lapply(1:Time, function(tp) {
        YPSmat = Reduce('+', lapply(1:P, function(p) {X[tp, , , p] * beta[tp, p]})) +
          outer(theta[tp, ], theta[tp, ], "+") + UDU[[tp]] + Errormat[tp,,]
        diag(YPSmat) = 0
        YPSmat
      })
      
      # replace NA's with YPS
      for (tp in which(vapply(na.positions, function(i) {length(i)}, c(1)) > 0)) {
        YPS[[tp]][meaningful_NA[[tp]]] = NA
        Y[tp, , ][na.positions[[tp]]] = YPS[[tp]][na.positions[[tp]]]
        YPS[[tp]][na.positions[[tp]]] = 0
      }
      id = (iter - burn) / odens
      s2PS[id,] = s2
      tauPS[id, ] = c(tau_p, tau_i, tau_u)
      for (tp in 1:Time) {
        BETAPS[[tp]][id, ] = beta[tp, ]
        thetaPS[[tp]][id, ] = theta[tp, ]
        DPS[[tp]][id, ] = d[tp, ]
        UPS[[tp]] = cbind(UPS[[tp]], U[tp,,])
        UDUPS[[tp]] =  UDUPS[[tp]] + UDU[[tp]]
        YPSsum[[tp]] = YPSsum[[tp]] + YPS[[tp]]
        Degreestats[[tp]] = rbind(Degreestats[[tp]], rowSums(YPS[[tp]], na.rm = TRUE))
        secondDegreestats[[tp]] = rbind(secondDegreestats[[tp]], rowSums(YPS[[tp]] %*% YPS[[tp]], na.rm = TRUE))
        thirdDegreestats[[tp]] = rbind(thirdDegreestats[[tp]], rowSums(YPS[[tp]] %*% YPS[[tp]] %*% YPS[[tp]], na.rm = TRUE))
      }
      Degrees = c(vapply(1:Time, function(tp) { rowSums(YPS[[tp]], na.rm = TRUE)}, rep(0, N)))
      corrstats[id, ] = vapply(1:(Time-1), function(l) {cor(Degrees[1:(N*(Time - l))], Degrees[(1 + N*l):(N*Time)], use = "complete")}, 0)
    }
  }
  
  UDUPM = lapply(UDUPS, function(x) {
    x / length(s2PS)
  })
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
    which(rank(-abs(eULU[[tp]]$val), ties.method = "first") <= R)
  })
  U =  lapply(1:Time, function(tp) {
    Uest = eULU[[tp]]$vec[, seq(1, R, length = R), drop = FALSE]
    exclude = meaningful_NA_rows[[tp]]
    if (length(exclude) > 0) {
      rownames(Uest) = rownames(Y[1,,])[-exclude]
    } else {
      rownames(Uest) = rownames(Y[1,,])
    }
    Uest
  })
  L =  lapply(1:Time, function(tp){
    eULU[[tp]]$val[eR[[tp]]]
  })
  
  YPM = lapply(YPSsum, function(x) {
    x / length(s2PS)
  })
  final = list(YPM = YPM, BETA = BETAPS, theta = thetaPS, UDU = UDUPM,
               U = U, D = L, s2 = s2PS, tau = tauPS, DPS = DPS, UPS = UPS, corr = corrstats,
               Degree = Degreestats, secondDegree= secondDegreestats, thirdDegree = thirdDegreestats)
  return(final)
}


#' @title GPpost
#' @description log-posterior distribution of kappa in case of multiple dimensions
#'
#' @param kappa current value of kappa
#' @param parameter parameter which kappa is embedded
#' @param tau scale parameter of Gaussian process covariance
#' @param dist_ij matrix of squared distance between all timepoints
#' @param a inverse-gamma shape parameter
#' @param b inverse-gamma scale parameter
#'
#' @return log-posterior density of kappa
#'
#' @export
GPpost = function(kappa, parameter, tau, dist_ij, a, b) {
    if (kappa > 0 & tau > 0) {
        lprior1 = dhalfcauchy(kappa, 5, log = TRUE)
        lprior2 = sum(LaplacesDemon::dinvgamma(tau, a, b, log = TRUE))
        Time = nrow(parameter)
        covfcmat = Exponential(dist_ij, kappa, phi = tau)
        if (is.positive.definite(covfcmat) == FALSE) {covfcmat = as.matrix(nearPD(covfcmat)$mat)}
        dims = ncol(parameter)
        lpost = 0
        for (col in 1:dims) {
        		if (sum(is.na(parameter[,col])) == 0) {
            lpost = lpost + rcpp_log_dmvnorm(covfcmat, rep(0, Time), c(parameter[,col]), FALSE)
            }
        }
        sums = as.vector(lpost + lprior1 + lprior2)
    } else {
        sums = -Inf
    }
    return(sums)
}

#' @title DAME_MH_revised
#' @description MCMC algorithm using Gibbs sampling for each variable with structural NA
#'
#' @param Y relational array of relations
#' @param X Time x n x n x p covariate array
#' @param RE random effect to be included ("additive" and/or "multiplicative")
#' @param R dimension of the multiplicative effects
#' @param dist standard Exponential or squared Exponential
#' @param gammapriors inverse-gamma shape and scale parameters for (s2, beta, theta, d)
#' @param avail Time x n matrix reperesenting the availibility of nodes (1 if avail, 0 if structural NA)
#' @param burn burn in for the Markov chain
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param odens output density for the Markov chain
#'
#' @return Final estimate of the parameters
#'
#' @export
DAME_MH_revised = function(Y, X, RE = c("additive", "multiplicative"), R = 2, dist = "Exponential", gammapriors = c(2, 1), avail = matrix(1, dim(Y)[1], dim(Y)[2]), burn = 1000, nscan = 5000, odens = 100) {
  Time = dim(Y)[1]
  N = dim(Y)[2]
  P = dim(X)[4]
  a = gammapriors[1]; b = gammapriors[2]
  # construct covariance matrix for each variable
  dist_ij = c()
  if (dist == "Exponential"){
    for (i in 1:Time) {
      for (j in 1:Time) {
        dist_ij = c(dist_ij, abs(i-j))
      }
    }
  } else {
    for (i in 1:Time) {
      for (j in 1:Time) {
        dist_ij = c(dist_ij, (i-j)^2)
      }
    }
  }
  dist_ij = matrix(dist_ij, nrow = Time, ncol = Time)
  kappas = rhalfcauchy(P+1+R, scale = 5)
  #kappas = rep(0.001, P+1+R)
  cinv = lapply(1:(P+1+R), function(k) {
    rcppeigen_invert_matrix(Exponential(dist_ij, kappas[k]))
  })
  # select initial values
  beta = matrix(0, Time, P)
  d = matrix(1, Time, R)
  U = array(0, dim = c(Time, N, R))
  theta = matrix(0, Time, N)
  tau_p = rep(1, P)
  tau_i = 1
  tau_u = rep(N, R)

  BETAPS = lapply(1:Time, function(tp) {
    matrix(0, nrow = nscan / odens, ncol = P)
  })
  UDUPS = lapply(1:Time, function(tp) {
    matrix(0, N, N)
  })
  YPSsum = lapply(1:Time, function(tp) {
    matrix(0, N, N)
  })
  s2PS = matrix(0, nrow = nscan / odens, ncol = 1)
  tauPS = matrix(0, nrow = nscan / odens, ncol = P + 1 + R)
  kappaPS = matrix(0, nrow = nscan / odens, ncol = P + 1 + R)
  thetaPS = lapply(1:Time, function(tp) {
    matrix(0, nrow = nscan / odens, ncol = N)
  })
  DPS = lapply(1:Time, function(tp) {
    matrix(0, nrow = nscan / odens, ncol = R)
  })
  UPS = lapply(1:Time, function(tp) {
   c()
  })
  Degreestats = lapply(1:Time, function(tp) {
    matrix(0, nrow = 0, ncol = N)
  })
  secondDegreestats = lapply(1:Time, function(tp) {
    matrix(0, nrow = 0, ncol = N)
  })
  thirdDegreestats = lapply(1:Time, function(tp) {
    matrix(0, nrow = 0, ncol = N)
  })
  corrstats = matrix(0, nrow = nscan / odens, ncol = Time - 1)
  years = sample(1:Time, 1)
  # to begin with, use 0 for NA's except meaningful ones
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
  meaningful_NA_years = which(vapply(meaningful_NA_rows, function(i) {length(i)}, c(1)) > 0)
  # to begin with, use 0 for NA's
  for (p in 1:P) {
    X[, , , p][which(is.na(X[, , , p]))] = 0
    for (tp in 1:Time) {
      X[tp, , , p][meaningful_NA[[tp]]] = NA
      Y[tp, , ][meaningful_NA[[tp]]] = NA
    }
  }
  na.positions = lapply(1:Time, function(tp) {
    which(is.na(Y[tp, , ]))
  })
  for (tp in 1:Time) {
    mu = mean(Y[tp, , ], na.rm = TRUE)
    row = rowMeans(Y[tp, , ] - mu, na.rm = TRUE)
    row[is.na(row)] = 0
    YA = mu + outer(row, row, "+")
    Y[tp, , ][na.positions[[tp]]] = YA[na.positions[[tp]]]
    diag(Y[tp, , ]) = 0
    beta[tp, 1] = mu
    if ("additive" %in% RE) {
      theta[tp, ] = row}
    if ("multiplicative" %in% RE) {
      eigenError = eigen(Y[tp, , ] - YA)
      eR = which(rank(-abs(eigenError$val), ties.method = "first") <= R)
      signs = eigenError$val[eR]
      d[tp, ] = ifelse(signs > 0, 1, -1)
    }
    Y[tp, , ][meaningful_NA[[tp]]] = NA
  }
  XB = lapply(1:Time, function(tp) {
    Reduce('+', lapply(1:P, function(p){
      X[tp, , , p] * beta[tp, p]
    }))
  })
  UDU = UDUPS
  uppertri = upper.tri(diag(N))
  # starting the Gibbs sampler
  for (iter in 1:(burn + nscan)) {
    if (iter %% 500 == 0) print(iter)
    s2 = rs2_fc2(XB, theta, UDU, Y, a, b)
    #M-H of tau_p and kappa_p
    for (p in 1:P) {
      var.old = c(kappas[p], tau_p[p])
      var.new = exp(rcpp_rmvnorm(1, 0.2 * diag(2) , log(var.old)))
      u = log(runif(1))
      temp = GPpost(var.new[1], matrix(beta[,p], ncol = 1), var.new[2], dist_ij, a, b) -
             GPpost(var.old[1], matrix(beta[,p], ncol = 1), var.old[2], dist_ij, a, b)
      alpha = min(0, temp)
      if (u <= alpha) {
        kappas[p] = var.new[1]
        cinv[[p]] = rcppeigen_invert_matrix(Exponential(dist_ij, kappas[p]))
        tau_p[p] = var.new[2]
      }
    }
    beta = rbeta_fc(X, beta, theta, UDU, Y, cinv[1:P], tau_p, s2)
    XB = lapply(1:Time, function(tp) {
      Reduce('+', lapply(1:P, function(p){
        X[tp, , , p] * beta[tp, p]
      }))
    })
    if ("additive" %in% RE) {
      #M-H of tau_i and kappa_i
      var.old = c(kappas[P+1], tau_i)
      var.new = exp(rcpp_rmvnorm(1, 0.01 * diag(2) , log(var.old)))
      u = log(runif(1))
      temp = GPpost(var.new[1], theta, var.new[2], dist_ij, a, b) -
             GPpost(var.old[1], theta, var.old[2], dist_ij, a, b)
       alpha = min(0, temp)
      if (u <= alpha) {
        kappas[P+1] = var.new[1]
        cinv[[P+1]] = rcppeigen_invert_matrix(Exponential(dist_ij, kappas[P+1]))
        tau_i = var.new[2]
      }
      theta = rtheta_fc(XB, theta, UDU, Y, cinv[[P+1]], tau_i, s2)
      for (tp in meaningful_NA_years) {
        theta[tp, meaningful_NA_rows[[tp]]] = 0
      }
    }
    if ("multiplicative" %in% RE) {
    	if (iter < burn) { prop = 0.1 } else { prop = 0.025 }
      if (iter > 0.1 * burn) {
        #M-H of kappa_r
        for (r in 1:R) {
          var.old = c(kappas[P+1+r], tau_u[r])
          var.new = exp(rcpp_rmvnorm(1, prop * diag(2) , log(var.old)))
          u = log(runif(1))
          temp = GPpost(var.new[1], matrix(U[, , r], ncol = N), var.new[2], dist_ij, a, b) -
                 GPpost(var.old[1], matrix(U[, , r], ncol = N), var.old[2], dist_ij, a, b)
          alpha = min(0, temp)
          if (u <= alpha) {
            kappas[P+1+r] = var.new[1]
            cinv[[P+1+r]] = rcppeigen_invert_matrix(Exponential(dist_ij, kappas[P+1+r]))
            tau_u[r] = var.new[2]
          }
        }
	  U = ru_fc_NA_revised(XB, theta, U, d, Y, cinv[(P+2):(P+1+R)], tau_u, s2, meaningful_NA_rows)
      d = rd_fc2_revised(XB, theta, U, d, Y, s2, meaningful_NA)
      for (tp in meaningful_NA_years) {
        U[tp, meaningful_NA_rows[[tp]],] = NA
      }
      UDU = lapply(1:Time, function(tp) {
        if (R <= 1) {
          d[tp, ] * U[tp, , ] %*% t(U[tp, ,])
        } else {
          U[tp, , ] %*% diag(d[tp, ]) %*% t(U[tp, , ])
        }
      })
      }
    }
    if (iter > burn & (iter-burn) %% odens == 0) {
      Errormat = array(0, dim = c(Time, N, N))
      errors = rnorm(Time * N * (N-1) / 2, 0, sqrt(s2))
      for (tp in 1:Time) {
        Errormat[tp, , ][upper.tri(Errormat[tp, , ])] = errors[((tp-1)*N*(N-1)/2+1):(tp*N*(N-1)/2)]
        Errormat[tp, , ] = (Errormat[tp, , ] + t(Errormat[tp, , ]))
      }
      YPS = lapply(1:Time, function(tp) {
        YPSmat = Reduce('+', lapply(1:P, function(p) {X[tp, , , p] * beta[tp, p]})) +
          outer(theta[tp, ], theta[tp, ], "+") + UDU[[tp]] + Errormat[tp,,]
        diag(YPSmat) = 0
        YPSmat
      })
      
      # replace NA's with YPS
      for (tp in which(vapply(na.positions, function(i) {length(i)}, c(1)) > 0)) {
        YPS[[tp]][meaningful_NA[[tp]]] = NA
        Y[tp, , ][na.positions[[tp]]] = YPS[[tp]][na.positions[[tp]]]
        YPS[[tp]][na.positions[[tp]]] = 0
      }
      id = (iter - burn) / odens
      s2PS[id,] = s2
      kappaPS[id, ] = c(kappas)
      tauPS[id, ] = c(tau_p, tau_i, tau_u)
      for (tp in 1:Time) {
        BETAPS[[tp]][id, ] = beta[tp, ]
        thetaPS[[tp]][id, ] = theta[tp, ]
        DPS[[tp]][id, ] = d[tp, ]
        UPS[[tp]] = cbind(UPS[[tp]], U[tp,,])
        UDUPS[[tp]] =  UDUPS[[tp]] + UDU[[tp]]
        YPSsum[[tp]] = YPSsum[[tp]] + YPS[[tp]]
        Degreestats[[tp]] = rbind(Degreestats[[tp]], rowSums(YPS[[tp]], na.rm = TRUE))
        secondDegreestats[[tp]] = rbind(secondDegreestats[[tp]], rowSums(YPS[[tp]] %*% YPS[[tp]], na.rm = TRUE))
        thirdDegreestats[[tp]] = rbind(thirdDegreestats[[tp]], rowSums(YPS[[tp]] %*% YPS[[tp]] %*% YPS[[tp]], na.rm = TRUE))
      }
      Degrees = c(vapply(1:Time, function(tp) { rowSums(YPS[[tp]], na.rm = TRUE)}, rep(0, N)))
      corrstats[id, ] = vapply(1:(Time-1), function(l) {cor(Degrees[1:(N*(Time - l))], Degrees[(1 + N*l):(N*Time)], use = "complete")}, 0)
    }
  }
  UDUPM = lapply(UDUPS, function(x) {
    x / length(s2PS)
  })
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
    which(rank(-abs(eULU[[tp]]$val), ties.method = "first") <= R)
  })
  U =  lapply(1:Time, function(tp) {
    Uest = eULU[[tp]]$vec[, seq(1, R, length = R), drop = FALSE]
    exclude = meaningful_NA_rows[[tp]]
    if (length(exclude) > 0) {
      rownames(Uest) = rownames(Y[1,,])[-exclude]
    } else {
      rownames(Uest) = rownames(Y[1,,])
    }
    Uest
  })
  L =  lapply(1:Time, function(tp){
    eULU[[tp]]$val[eR[[tp]]]
  })
  
  YPM = lapply(YPSsum, function(x) {
    x / length(s2PS)
  })
  final = list(YPM = YPM, BETA = BETAPS, theta = thetaPS, UDU = UDUPM, 
               U = U, D = L, s2 = s2PS, tau = tauPS, kappas = kappaPS, DPS = DPS, UPS = UPS, corr = corrstats,
               Degree = Degreestats, secondDegree= secondDegreestats, thirdDegree = thirdDegreestats)
  return(final)
}


#' @title DAME_UU_fixed_revised
#' @description MCMC algorithm using Gibbs sampling for each variable with structural NA
#'
#' @param Y relational array of relations
#' @param X Time x n x n x p covariate array
#' @param RE random effect to be included ("additive" and/or "multiplicative)
#' @param R dimension of the multiplicative effects
#' @param dist standard Exponential or squared Exponential
#' @param gammapriors inverse-gamma shape and scale parameters for (s2, beta, theta, d)
#' @param avail = Time x n matrix reperesenting the availibility of nodes (1 if avail, 0 if structural NA)
#' @param burn burn in for the Markov chain
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param odens output density for the Markov chain
#' @param kappas P+2 length vector of GP length parameters
#'
#' @return Final estimate of the parameters
#'
#' @export
DAME_UU_fixed_revised = function(Y, X, RE = c("additive", "multiplicative"), R = 2, dist = "Exponential", gammapriors = c(2, 1), avail = matrix(1, dim(Y)[1], dim(Y)[2]), burn = 1000, nscan = 5000, odens = 100, kappas) {
    Time = dim(Y)[1]
    N = dim(Y)[2]
    P = dim(X)[4]
    a = gammapriors[1]; b = gammapriors[2]
    # construct covariance matrix for each variable
    dist_ij = c()
    if (dist == "Exponential"){
        for (i in 1:Time) {
            for (j in 1:Time) {
                dist_ij = c(dist_ij, abs(i-j))
            }
        }
    } else {
        for (i in 1:Time) {
            for (j in 1:Time) {
                dist_ij = c(dist_ij, (i-j)^2)
            }
        }
    }
    dist_ij = matrix(dist_ij, nrow = Time, ncol = Time)
    cinv = lapply(1:(P+1+R), function(k) {
        rcppeigen_invert_matrix(Exponential(dist_ij, kappas[k]))
    })
    # select initial values
    beta = matrix(0, Time, P)
    U = array(0, dim = c(Time, N, R))
    d = matrix(1, Time, R)
    theta = matrix(0, Time, N)
    BETAPS = lapply(1:Time, function(tp) {
        matrix(0, nrow = nscan / odens, ncol = P)
    })
    UDUPS = lapply(1:Time, function(tp) {
        matrix(0, N, N)
    })
    YPSsum = lapply(1:Time, function(tp) {
        matrix(0, N, N)
    })
    s2PS = matrix(0, nrow = nscan / odens, ncol = 1)
    tauPS = matrix(0, nrow = nscan / odens, ncol = P + 1 + R)
    thetaPS = lapply(1:Time, function(tp) {
        matrix(0, nrow = nscan / odens, ncol = N)
    })
    UPS = array(0, dim = c(Time, N, R))
    Degreestats = lapply(1:Time, function(tp) {
        matrix(0, nrow = 0, ncol = N)
    })
    secondDegreestats = lapply(1:Time, function(tp) {
        matrix(0, nrow = 0, ncol = N)
    })
    thirdDegreestats = lapply(1:Time, function(tp) {
        matrix(0, nrow = 0, ncol = N)
    })
    corrstats = matrix(0, nrow = nscan / odens, ncol = Time - 1)
    
    years = sample(1:Time, 1)
    # to begin with, use 0 for NA's except meaningful ones
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
    meaningful_NA_years = which(vapply(meaningful_NA_rows, function(i) {length(i)}, c(1)) > 0)
    # to begin with, use 0 for NA's
    for (p in 1:P) {
        X[, , , p][which(is.na(X[, , , p]))] = 0
        for (tp in 1:Time) {
            X[tp, , , p][meaningful_NA[[tp]]] = NA
            Y[tp, , ][meaningful_NA[[tp]]] = NA
        }
    }
    na.positions = lapply(1:Time, function(tp) {
        which(is.na(Y[tp, , ]))
    })
    for (tp in 1:Time) {
        mu = mean(Y[tp, , ], na.rm = TRUE)
        row = rowMeans(Y[tp, , ] - mu, na.rm = TRUE)
        row[is.na(row)] = 0
        YA = mu + outer(row, row, "+")
        Y[tp, , ][na.positions[[tp]]] = YA[na.positions[[tp]]]
        diag(Y[tp, , ]) = 0
        beta[tp, 1] = mu
        if ("additive" %in% RE) {
      	theta[tp, ] = row
      	}
        Y[tp, , ][meaningful_NA[[tp]]] = NA
    }
    XB = lapply(1:Time, function(tp) {
        Reduce('+', lapply(1:P, function(p){
            X[tp, , , p] * beta[tp, p]
        }))
    })
    UDU = UDUPS
    uppertri = upper.tri(diag(N))
    tau_p = rep(1, P)
  	tau_i = 1
  	tau_u = rep(N, R)
    # starting the Gibbs sampler
    for (iter in 1:(burn + nscan)) {
        if (iter %% 100 == 0) print(iter)
        s2 = rs2_fc2(XB, theta, UDU, Y, a, b)
        tau_p = rtaup_fc(beta, cinv[1:P], a, b)
        beta = rbeta_fc(X, beta, theta, UDU, Y, cinv[1:P], tau_p, s2)
        XB = lapply(1:Time, function(tp) {
            Reduce('+', lapply(1:P, function(p){
                X[tp, , , p] * beta[tp, p]
            }))
        })
        if ("additive" %in% RE) {
        tau_i = rtaui_fc(theta, cinv[[P+1]], a, b)
        theta = rtheta_fc(XB, theta, UDU, Y, cinv[[P+1]], tau_i, s2)
        for (tp in meaningful_NA_years) {
            theta[tp, meaningful_NA_rows[[tp]]] = 0
        }
        }
        if ("multiplicative" %in% RE) {
        if (iter > 0.1 * burn) {
        	 U = ru_fc_NA_revised(XB, theta, U, d, Y, cinv[(P+2):(P+1+R)], tau_u, s2, meaningful_NA_rows)
        	 tau_u = rtauu_fc(U, cinv[(P+2):(P+1+R)], a, b)
      for (tp in meaningful_NA_years) {
        U[tp, meaningful_NA_rows[[tp]],] = NA
      }
        UDU = lapply(1:Time, function(tp) {
            U[tp, , ] %*% t(U[tp, ,])
        })
        }
        }
        if (iter > burn & (iter - burn) %% odens == 0) {
            Errormat = array(0, dim = c(Time, N, N))
            errors = rnorm(Time * N * (N-1) / 2, 0, sqrt(s2))
            for (tp in 1:Time) {
                Errormat[tp, , ][upper.tri(Errormat[tp, , ])] = errors[((tp-1)*N*(N-1)/2+1):(tp*N*(N-1)/2)]
                Errormat[tp, , ] = (Errormat[tp, , ] + t(Errormat[tp, , ]))
            }
            YPS = lapply(1:Time, function(tp) {
                YPSmat = Reduce('+', lapply(1:P, function(p) {X[tp, , , p] * beta[tp, p]})) +
                outer(theta[tp, ], theta[tp, ], "+") + UDU[[tp]] + Errormat[tp,,]
                diag(YPSmat) = 0
                YPSmat
            })
            # replace NA's with YPS
            for (tp in which(vapply(na.positions, function(i) {length(i)}, c(1)) > 0)) {
                YPS[[tp]][meaningful_NA[[tp]]] = NA
                Y[tp, , ][na.positions[[tp]]] = YPS[[tp]][na.positions[[tp]]]
                YPS[[tp]][na.positions[[tp]]] = 0
            }
            id = (iter - burn) / odens
            s2PS[id,] = s2
            tauPS[id, ] = c(tau_p, tau_i, tau_u)
            UPS = UPS + U
            for (tp in 1:Time) {
                BETAPS[[tp]][id, ] = beta[tp, ]
                thetaPS[[tp]][id, ] = theta[tp, ]
                UDUPS[[tp]] =  UDUPS[[tp]] + UDU[[tp]]
                YPSsum[[tp]] = YPSsum[[tp]] + YPS[[tp]]
                Degreestats[[tp]] = rbind(Degreestats[[tp]], rowSums(YPS[[tp]], na.rm = TRUE))
                secondDegreestats[[tp]] = rbind(secondDegreestats[[tp]], rowSums(YPS[[tp]] %*% YPS[[tp]], na.rm = TRUE))
                thirdDegreestats[[tp]] = rbind(thirdDegreestats[[tp]], rowSums(YPS[[tp]] %*% YPS[[tp]] %*% YPS[[tp]], na.rm = TRUE))
            }
            Degrees = c(vapply(1:Time, function(tp) { rowSums(YPS[[tp]], na.rm = TRUE)}, rep(0, N)))
            corrstats[id, ] = vapply(1:(Time-1), function(l) {cor(Degrees[1:(N*(Time - l))], Degrees[(1 + N*l):(N*Time)], use = "complete")}, 0)
        }
    }
    UDUPM = lapply(UDUPS, function(x) {
        x / length(s2PS)
    })
    U = UPS / length(s2PS)
    YPM = lapply(YPSsum, function(x) {
        x / length(s2PS)
    })
    final = list(YPM = YPM, BETA = BETAPS, theta = thetaPS, UDU = UDUPM,
    U = U, s2 = s2PS, tau = tauPS, corr = corrstats,
    Degree = Degreestats, secondDegree= secondDegreestats, thirdDegree = thirdDegreestats)
    return(final)
}

