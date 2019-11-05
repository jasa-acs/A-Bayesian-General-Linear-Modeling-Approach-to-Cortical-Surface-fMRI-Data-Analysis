# #library(devtools)
# #install_github('muschellij2/gifti')
library(gifti)
library(matrixStats) #rowVars

data_dir <- '~/HCP/data/hpc/'
wb_cmd <- '~/workbench/bin_rh_linux64/wb_command'
maindir <- '~/BayesianGLM/' #contains the following subdirectories:
#EVs: contains task design matrix (task convolved with HRF) for each subject
#motion: contains motion covariates (6 plus their temporal derivatives) for each subject
#locations: contains coordinates of surface vertices and CIFTI/GIFTI files used in transformations, used to store SPDE objects
#timeseries: contains smoothed and resampled fMRI timeseries for classical and Bayesian GLMs, respectively
#prewhitening: contains AR coefficient estimates, and sparse prewhitening matrix for Bayesian GLM
#results_OLS: contains results from classical GLM
#results_sphere: contains results from Bayesian GLM
#code: additional functions

source(paste0(maindir,'code/BayGLMfun.R'))
setwd(maindir)


###########################################################
# LOAD IN SETUP OBJECTS (FROM BAYGLM2.R)
###########################################################

load(file='subjects.Rdata')
M <- length(subjects) 

load('locations/loc32K.Rdata') #ns.32K, isNA


###########################################################
# FIT CLASSICAL GLM TO EACH SUBJECT
###########################################################

for(task in c('MOTOR','GAMBLING')){

	if(task=='MOTOR') { ntime <- 284; K <- 6 }
	if(task=='GAMBLING') { ntime <- 253; K <- 3 }

	###########################################################
	# GET AR COEFFICIENTS
	###########################################################

	load(file=paste0('prewhitening/',task,'/coef_avg_32K.Rdata')) #AR_resid_avg, var_resid_avg
	p <- 6

	############################################################
	## SET UP PERMUTATION TEST 
	############################################################

	keep <- (p+2):(ntime-p-1) #remove first and last (p+1) time points after prewhitening
	nkeep <- length(keep) #length of timeseries after prewhitening
	load(file=paste0('results_OLS/',task,'/permutation/time_perm.Rdata'))

	############################################################
	## SUBJECT-LEVEL MODELS 
	############################################################

	#0=no smoothing, 1=smoothing
	for(sm in 0:1){

		if(sm==0) suff <- '.nosmooth'
		if(sm==1) suff <- '.smooth'

		for(mm in 1:M){

			print(mm)
			t0 <- Sys.time()

			###########################################################
			# READ IN DESIGN MATRIX
			###########################################################
		    
			# Build up covariates
			file.name <- paste("EVs/",task,'/',subjects[mm],"_RL.csv",sep="")
			Z <- as.matrix(read.csv(file.name, header=T))  ## task
			Z <- Z/max(Z) #standardize
			ntime <- dim(Z)[1]
			Z <- scale(Z, scale=FALSE) #center to remove baseline

			file.name <- paste("motion/",task,'/',subjects[mm],"_RL.txt",sep="")
			X <- read.table(file.name)  ## nuisance
			X <- X/max(X) #standardize
			trend1 <- (1:ntime)/ntime #linear trends
			trend2 <- trend1^2 #quadratic trends
			X <- as.matrix(cbind(X, trend1, trend2))
			X <- scale(X, scale=FALSE)

			# Regress intercept & nuisance parameters from task
			invxvx <- solve(t(X) %*% X)  #(X'X)^{-1}
	 		betaX <- invxvx %*% t(X) %*% Z #(X'X)^{-1} X'Z
	 		residZ <- Z - X %*% betaX

		 	#############################################
			# CLASSICAL GLM WITH PREWHITENING
			#############################################

			# Time Series
			file.name <- paste("timeseries/",task,'/',subjects[mm],suff,".csv",sep="")
		    data.all <- as.matrix(read.csv(file.name, header=F))
		    local_means <- matrix(rowMeans(data.all, na.rm=TRUE), nrow=nrow(data.all), ncol=ncol(data.all))
		    y.all <- t(100*(data.all - local_means)/local_means) #scale to units of pct local signal change AND CENTER 

			# Regress nuisance parameters
	 		betaX <- invxvx %*% t(X) %*% y.all #(X'X)^{-1} X'Z
	 		residY <- y.all - X %*% betaX
			y <- residY
			Z <- residZ


			beta <- stebeta <- sigma <- matrix(nrow=ns.32K*2, ncol=K)
			resid <- matrix(nrow=ns.32K*2, ncol=nkeep)
			AR1_resid <- rep(NA, ns.32K*2) #first AR coefficient after prewhitening

			ptest_time <- 0

			for(v in 1:(ns.32K*2)){

				print(v)
				y.v <- y[,v]
				if(isNA[v]) next

				#prewhiten data and design matrix
				Inv <- getInvCovAR(ar=AR_resid_avg[v,], ntime=ntime)
				Dinv <- diag(rep(1/sqrt(var_resid_avg[v]), ntime))
				Inv <- Dinv %*% Inv %*% Dinv
				sqrtInv <- getSqrtInv(Inv)
				y.v <- sqrtInv %*% y.v #apply weights
				X.v <- sqrtInv %*% Z #apply weights

				# Remove first and last 7 time points (p+1 just to be sure)
				y.v <- y.v[keep]
				X.v <- X.v[keep,]

				# Compute beta and standard errors
				invxvx <- solve(t(X.v) %*% X.v)  #(X'Inv X)^{-1}
				invxvx_diag_sqrt <- sqrt(diag(invxvx))
				beta.v <- invxvx %*% t(X.v) %*% y.v #(X' Inv X)^{-1} X' Inv y
				resid.v <- y.v - X.v %*% beta.v
				sigma.v <- sqrt(var(resid.v)*(nkeep-1)/(nkeep-K))
				stebeta.v <- c(sigma.v) * invxvx_diag_sqrt

				# Save values
				beta[v,] <- beta.v
				stebeta[v,] <- stebeta.v
				sigma[v] <- sigma.v

				# Check autocorrelation and SE of residuals
				if(sum(is.na(resid.v))==0) {
					ar_v <- ar(resid.v, method='yw')
					AR1_resid[v] <- ifelse(ar_v$order==0, 0, ar_v$ar[1])
				}
			
				# Permute time and recompute beta and standard errors
				t0v <- Sys.time()
				mat_blank <- matrix(nrow=100, ncol=K)
				beta_perm_v <- stebeta_perm_v <- mat_blank
				for(k in 1:100){
					perm_k <- time_perm[k,]
					y.v.perm <- y.v[perm_k]
					beta.v <- invxvx %*% t(X.v) %*% y.v.perm 
					resid.v <- y.v - X.v %*% beta.v
					sigma.v <- sqrt(var(resid.v)*(nkeep-1)/(nkeep-K))
					stebeta.v <- sigma.v * invxvx_diag_sqrt
					beta_perm_v[k,] <- beta.v[1:K]
					stebeta_perm_v[k,] <- stebeta.v[1:K]
				}
				fname_v <- paste0("results_OLS/",task,"/permutation/beta",suff,'.',subjects[mm],".",v,".RData")
				save(beta_perm_v, stebeta_perm_v, file=fname_v)	
				tv <- (Sys.time() - t0v)
				ptest_time <- ptest_time + as.numeric(tv)

			} #end loop across locations
			Sys.time() - t0

			# Save to data file
			fname <- paste0("results_OLS/",task,"/beta",suff,'.',subjects[mm],".RData")
			save(beta, stebeta, sigma, file=fname)	
			load(file=fname) #p_resid, var_resid, AR1_resid

			if(mm==1){ AR1_resid_all <- AR1_resid } else { AR1_resid_all <- cbind(AR1_resid_all,AR1_resid) }

		} #end loop across subjects

		# Visualize first AR coefficient after prewhitening (should be close to zero)
		AR1_resid_avg = rowMeans(AR1_resid_all)
		writeCIFTIs(AR1_resid_avg, paste0("results_OLS/",task,"/AR_resid_POST"), hemisphere='both')
			
		###################################################################
		# GROUP-LEVEL MODEL
		###################################################################

		t0 <- Sys.time()

		### SAVE SUBJECT-LEVEL RESULTS & COMPUTE GRP-LEVEL ESTIMATES (< 1 MIN)

		beta_all <- stebeta_all <- array(dim=c(ns.32K*2, K, M))
		for(mm in (1:M)){
			print(mm)
			fname <- paste0("results_OLS/",task,"/beta",suff,'.',subjects[mm],".RData")
			load(fname)
			beta_all[,,mm] <- beta
			stebeta_all[,,mm] <- stebeta
			# #Visualize subject-level estimates
			# if(mm==1){
			# 	names <- paste0("results_OLS/",task,"/beta",1:K,suff,'.',subjects[mm])
			# 	writeCIFTIs(beta, names, hemisphere='both')
			# }
		}

		# even though the data has been standardized to have residual 
		# variance close to 1 at all locations, there are noticable spatial
		# differences in stebeta.  this is actually because of changes in the 
		# design matrix X induced from prewhitening, since stebeta = sigma^2(X'X)^(-1)

		stebeta_avg <- apply(stebeta_all, c(1,2), mean, na.rm=TRUE)
		names_stebeta <- paste0("results_OLS/",task,"/stebeta",1:K,"_avg",suff)
		writeCIFTIs(stebeta_avg, names_stebeta, hemisphere='both')

		### FIT GROUP-LEVEL MODEL (WEIGHTS = 1/SIGMA^2)

		wts <- stebeta_all^(-2) #inverse weight each beta_i by its squared standard error
		sum_wts <- apply(wts, c(1,2), sum, na.rm=TRUE) #sum weights over subjects (denominator in beta_grp and varbeta_grp)
		beta_grp <- apply(beta_all*wts, c(1,2), sum, na.rm=TRUE) / sum_wts #beta_grp = weighted average of beta_i's
		resid <- beta_all - array(beta_grp, dim=dim(beta_all)) #compute residuals (beta_i - beta_grp)
		sigma2 <- apply(resid/stebeta_all, c(1,2), var, na.rm=TRUE) #variance of residuals after whitening/weighting
		varbeta_grp <- sigma2 / sum_wts #var(beta) = sigma^2 (X'X)^(-1) = sigma^2 / sum(weights)
		stebeta_grp <- varbeta_grp^(.5)
		tstat <- beta_grp / stebeta_grp  #compute t statistic = beta/ste(beta) (df=20)

		print(Sys.time() - t0)

		names_beta <- paste0("results_OLS/",task,"/beta",1:K,"_grp",suff)
		writeCIFTIs(beta_grp, names_beta, hemisphere='both')

		names_stebeta <- paste0("results_OLS/",task,"/stebeta",1:K,"_grp",suff)
		writeCIFTIs(stebeta_grp, names_stebeta, hemisphere='both')


		### PERFORM FDR CORRECTION (< 1 SEC)

		t0 <- Sys.time()

		pvals <- pt(tstat, df=M-1, lower.tail=FALSE) #one-sided test
		pvals_FDR <- array(dim=dim(pvals))
		for(k in 1:K){
			pvals_FDR[,k] <- p.adjust(pvals[,k], "BH")
		}
		active_FDR95 <- 1*(pvals_FDR < 0.05)
		active_FDR99 <- 1*(pvals_FDR < 0.01)

		print(Sys.time() - t0)


		### PERFORM PERMUTATION TEST FOR FWER CORRECTION

		t0 <- Sys.time()
		tstat_perm <- array(dim=c(ns.32K*2, 100, K))
		beta_perm <- array(dim=c(ns.32K*2, 100, K))
		for(v in 1:(ns.32K*2)){

			print(v)
			if(isNA[v]) next

			# 1. Read in subject-level beta estimates for each permutation (1-100)
			beta_all_perm_v <- array(dim=c(100, K, M))
			stebeta_all_perm_v <- array(dim=c(100, K, M))
			for(mm in 1:M){
				fname_v <- paste0("results_OLS/",task,"/permutation/beta",suff,'.',subjects[mm],".",v,".RData")
				load(file=fname_v)#beta_perm_v, stebeta_perm_v
				beta_all_perm_v[,,mm] <- beta_perm_v
				stebeta_all_perm_v[,,mm] <- sqrt(stebeta_perm_v^2 *(nkeep-1)/(nkeep-6))
			}

			# 2. Fit group-level model for each permutation & coefficient
			wts <- stebeta_all_perm_v^(-2)
			sum_wts <- apply(wts, c(1,2), sum, na.rm=TRUE) #sum across subjects
			beta_grp_perm_v <- apply(beta_all_perm_v*wts, c(1,2), sum, na.rm=TRUE) / sum_wts #weighted sum across subjects
			ind_v <- (1:(ns.32K*2) == v) #faster than indexing
			beta_perm[ind_v,,] <- beta_grp_perm_v
			resid <- beta_all_perm_v - array(beta_grp_perm_v, dim=dim(beta_all_perm_v)) #compute residuals (beta_i - beta_grp)
			sigma2 <- apply(resid/stebeta_all_perm_v, c(1,2), var, na.rm=TRUE) #variance of residuals after weighting (division by stebeta_all_perm_v is the weighting)
			varbeta_grp <- sigma2 / sum_wts
			stebeta_grp_perm_v <- varbeta_grp^(.5)
			tstat_perm_v <- beta_grp_perm_v / stebeta_grp_perm_v  #compute t statistic = beta/ste(beta) (df=20)
			tstat_perm[ind_v,,] <- tstat_perm_v
		}
		print(Sys.time() - t0)

		#3. Compute max t-statistic across image for each permutation and coefficient

		#find maximum t-statistic across all voxels
		tstat_max <- apply(tstat_perm, c(2,3), max, na.rm=TRUE)
		tstat_q95 <- apply(tstat_max, 2, quantile, p=0.95) 
		tstat_q99 <- apply(tstat_max, 2, quantile, p=0.99)
		active_FWER95 <- 1*(tstat > matrix(tstat_q95, nrow=ns.32K*2, ncol=6, byrow=TRUE))
		active_FWER99 <- 1*(tstat > matrix(tstat_q99, nrow=ns.32K*2, ncol=6, byrow=TRUE))

		# Visualize FWER + FDR
		active_combined95 <- active_FDR95 + active_FWER95
		active_combined99 <- active_FDR99 + active_FWER99
		active_combined95[active_combined95==0] <- 0.01
		active_combined99[active_combined99==0] <- 0.01
		names_combined95 <- paste0("results_OLS/",task,"/active",1:K,"_combo95",suff)
		names_combined99 <- paste0("results_OLS/",task,"/active",1:K,"_combo99",suff)
		writeCIFTIs(active_combined95, names_combined95, hemisphere='both')
		writeCIFTIs(active_combined99, names_combined99, hemisphere='both')

	} # end smoothing/no smoothing loop

} #end loop over tasks

