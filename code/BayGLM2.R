library(reshape2) #melt
library(ggplot2)
library(Matrix)
library(INLA)
INLA:::inla.dynload.workaround() #to avoid error on creating mesh

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
# GET LIST OF SUBJECTS
###########################################################

## 20 RANDOMLY SAMPLED HCP 500 SUBJECTS

subjects <- c("100307", "120111", "126628", "131924", "154936", "182840", "198451", "211417", "334635", "371843", "620434", "638049", "677968", "833148", "845458", "849971", "901442", "904044", "912447", "922854")
disks <- c("disk1", "disk1", "disk1", "disk1", "disk2", "disk3", "disk3", "disk3", "disk4", "disk4", "disk4", "disk4", "disk5", "disk5", "disk5", "disk5", "disk5", "disk5", "disk5", "disk5")
save(subjects, disks, file='subjects.Rdata')
M <- length(subjects)   

###########################################################
# GET COORDINATES OF VERTICES, BUILD MESH AND SPDE OBJECT (30 SEC)
###########################################################

t0 <- Sys.time()

## Index NA locations in 32K data (Classical GLM)
file.name <- paste("timeseries/MOTOR/",subjects[1],".nosmooth.csv",sep="")
data.all <- read.csv(file.name, header = F)
isNA <- is.na(data.all[,1]) #same with or without smoothing
ns.32K <- length(isNA)/2
save(ns.32K, isNA, file='locations/loc32K.Rdata')
#load(file='locations/loc32K.Rdata')

## Remove NA locations in 6K data (Bayesian GLM)
file.name <- paste("timeseries/MOTOR/",subjects[1],".6K.csv",sep="")
data.all <- read.csv(file.name, header = F)
ns.all <- nrow(data.all)
idx.obs <- which(data.all[,1] != 0)
idx.obs1 <- idx.obs[idx.obs <= ns.all/2] #half of locations in left hemisphere
idx.obs2 <- idx.obs[idx.obs > ns.all/2] #half of locations in right hemisphere
ns1 <- length(idx.obs1)
ns2 <- length(idx.obs2)
ns <- c(ns1, ns2)
save(ns.all, idx.obs, idx.obs1, idx.obs2, ns1, ns2, ns, file='locations/ind_obs.Rdata')
#load(file='locations/ind_obs.Rdata')

## Build Mesh and SPDE object
loc.all <- read.csv('locations/positions_sphere_resamp.csv', header=F)
loc.org1 <- as.matrix(loc.all[idx.obs1,]) #left
loc.org2 <- as.matrix(loc.all[idx.obs2,]) #right
loc.norm1 <- loc.org1/sqrt(rowSums(loc.org1^2)) #left
loc.norm2 <- loc.org2/sqrt(rowSums(loc.org2^2)) #right
mesh1 <- inla.mesh.2d(loc.norm1, cutoff = .022, max.edge = .2) #left
mesh2 <- inla.mesh.2d(loc.norm2, cutoff = .022, max.edge = .2) #right
save(loc.norm1, loc.norm2, mesh1, mesh2, file='locations/mesh.sphere.resamp6K.Rdata')
#load(file='locations/mesh.sphere.resamp6K.Rdata')

## Matrix to transform back to original locations
Amat1 <- inla.spde.make.A(mesh1, as.matrix(loc.norm1))
Amat2 <- inla.spde.make.A(mesh2, as.matrix(loc.norm2))
Amat <- list(Amat1, Amat2)
save(Amat, file='locations/Amat.Rdata')
#load(file='locations/Amat.Rdata')

## Build SPDE object
spatial1 <- mesh1$idx$loc #left
spatial2 <- mesh2$idx$loc #right
spde1 <- inla.spde2.matern(mesh1) #left
spde2 <- inla.spde2.matern(mesh2) #right
spatial <- list(spatial1, spatial2)
spde <- list(spde1, spde2)
save(spde, spatial, file='locations/spde.Rdata')
#load(file='locations/spde.Rdata')


###########################################################
# PREWHITENING
###########################################################

#1. FIT CLASSICAL GLM MODEL TO EACH SUBJECT
#2. FIT AR MODEL TO RESIDUALS
#3. AVERAGE AR COEFFICIENT ESTIMATES ACROSS SUBJECTS


p <- 6 #AR model order
fname_ts <- '_Atlas.dtseries.nii'

for(task in c('MOTOR','GAMBLING')){

t0 <- Sys.time()

	session <- paste0('tfMRI_',task,'_RL')

	#########################################
	### Estimate AR coefficients from residuals (1 HOUR)
	#########################################

	#6K DATA: 45 SEC / SUBJECT
	#32K DATA: 2 MIN / SUBJECT

	#loop over 6K, 32K data
	for(n in c(6,32)){

		print(paste0('Start ',n,'K Data'))

		for(mm in 1:M){

			print(mm)
			t0 <- Sys.time()

			#########################################
			### READ IN AND PROCESS COVARIATES

			file.name <- paste("EVs/",task,'/',subjects[mm],"_RL.csv",sep="")
			Z <- as.matrix(read.csv(file.name, header=T))  ## task
			Z <- Zorig <- Z/max(Z) #standardize
			ntime <- dim(Z)[1]
			Z <- scale(Z, scale=FALSE) #center to remove baseline/intercept

			file.name <- paste("motion/",task,'/',subjects[mm],"_RL.txt",sep="")
			X <- read.table(file.name)  ## nuisance
			X <- X/max(X) #standardize
			trend1 <- (1:ntime)/ntime #linear trends
			trend2 <- trend1^2 #quadratic trends
			X <- as.matrix(cbind(X, trend1, trend2))
			X <- scale(X, scale=FALSE) #center to remove baseline/intercept

			# Regress nuisance parameters from task
			invxvx <- solve(t(X) %*% X)  #(X'X)^{-1}
			betaX <- invxvx %*% t(X) %*% Z #(X'X)^{-1} X'Z
			residZ <- Z - X %*% betaX

			#########################################
			### PLOT HRFS BEFORE AND AFTER NUISANCE REGRESSION

			## Put data into long format

			Z.df <- as.data.frame(Zorig)
			residZ.df <- as.data.frame(residZ)
			Z.df$time <- residZ.df$time <- 1:ntime
			Z.long <- melt(Z.df, id.vars='time')
			residZ.long <- melt(residZ.df, id.vars='time')
			names(Z.long)[2] <- names(residZ.long)[2] <- 'Task'

			## Plot activation profiles

			if(task=='MOTOR') colors <- c('black', 'black','#e41a1c','black','#e41a1c','black')
			if(task=='GAMBLING') colors <- c('black','black','black')

			if(task=='MOTOR') lines <- c(1,2,2,3,3,4)
			if(task=='GAMBLING') lines <- c(1,2,3)

			pdf(paste0('~/Bayesian2D/plots/hrf_',task,'_subj',mm,'.pdf'), width=8, height=5)

			#before nuisance regression
			print(ggplot(data = Z.long, aes(x = time, y=value, colour=Task, linetype=Task)) + geom_line() +
			  scale_colour_manual('',values=colors) +
			  scale_linetype_manual('',values=lines) +
			  xlab("Time") + ylab("Expected Response") + ggtitle('Before Nuisance Regression') +
			  theme_bw() + theme(text=element_text(size = 20), 
			                     axis.text.y = element_text(angle = 90, hjust = 0.5),
			                     axis.title.y = element_text(margin=margin(0,10,0,0)),
			                     axis.title.x = element_text(margin=margin(0,0,0,10)),
			                     plot.title = element_text(hjust = 0.5),
			                     legend.position = 'bottom', legend.direction='horizontal',
			                     legend.background = element_rect(colour = "black"),
			                     panel.grid = element_blank()) + 
			  guides(linetype=guide_legend(nrow=1), colour=guide_legend(nrow=1)))

			#after nuisance regression
			print(ggplot(data = residZ.long, aes(x = time, y=value, colour=Task, linetype=Task)) + geom_line() +
			  scale_colour_manual('',values=colors) +
			  scale_linetype_manual('',values=lines) +
			  xlab("Time") + ylab("Expected Response") + ggtitle('After Nuisance Regression') +
			  theme_bw() + theme(text=element_text(size = 20), 
			                     axis.text.y = element_text(angle = 90, hjust = 0.5),
			                     axis.title.y = element_text(margin=margin(0,10,0,0)),
			                     axis.title.x = element_text(margin=margin(0,0,0,10)),
								 plot.title = element_text(hjust = 0.5),
			                     legend.position = 'bottom', legend.direction='horizontal',
			                     legend.background = element_rect(colour = "black"),
			                     panel.grid = element_blank()) + 
			  guides(linetype=guide_legend(nrow=1), colour=guide_legend(nrow=1)))

			dev.off()

			#########################################
			### READ IN AND PROCESS FMRI TIMESERIES 

			if(n==6)  file.name <- paste("timeseries/",task,'/',subjects[mm],".6K.csv",sep="")
			if(n==32) file.name <- paste("timeseries/",task,'/',subjects[mm],".smooth.csv",sep="")
	  	
			data.all <- as.matrix(read.csv(file.name, header=F))
		  	V.all <- nrow(data.all)
		    local_means <- matrix(rowMeans(data.all, na.rm=TRUE), nrow=nrow(data.all), ncol=ncol(data.all))
		    y.all <- t(100*(data.all - local_means)/local_means) #scale to units of pct local signal change AND CENTER

			# Regress out nuisance parameters
			betaX <- invxvx %*% t(X) %*% y.all #(X'X)^{-1} X'Z
			residY <- y.all - X %*% betaX

			#########################################
			### COMPUTE AND ANALYZE RESIDUALS

			# Get residuals
			y <- residY
			Z <- residZ
			Beta <- solve(t(Z) %*% Z) %*% t(Z) %*% y
			Resid <- y - Z %*% Beta

			# Fit AR 
			if(mm==1){
				p_resid <- matrix(NA, nrow=V.all, ncol=M)
				AR_resid <- array(NA, dim=c(V.all,p,M))
				var_resid <- matrix(NA, nrow=V.all, ncol=M)
			}
			for(v in 1:V.all){
				if(is.na(Resid[1,v])) next
				ar_v <- ar(Resid[,v], method='yw', aic=FALSE, order.max=6)
				p_resid[v,mm] <- ar(Resid[,v], method='yw')$order
				AR_resid[v,,mm] <- ar_v$ar
				var_resid[v,mm] <- ar_v$var.pred #should be *(nkeep-1)/(nkeep-6) to get SD=1 but current code should result in constant variance across brain
			}

			print(Sys.time()-t0)
		} #end loop across subjects

		# Save avg AR coefs and variance for prewhitening
		AR_resid_avg <- apply(AR_resid, c(1,2), mean, na.rm=TRUE)
		var_resid_avg = rowMeans(var_resid)
		save(AR_resid_avg, var_resid_avg, file=paste0('prewhitening/',task,'/coef_avg_',n,'K.Rdata'))
		load(file=paste0('prewhitening/',task,'/coef_avg_',n,'K.Rdata'))

		if(n==32){

			# Visualize AR properties and SE of residuals after prewhitening
			p_resid_avg = rowMeans(p_resid)
			AR1_resid_avg = rowMeans(AR_resid[,1,])
			SE_resid_avg <- sqrt(var_resid_avg)
			writeCIFTIs(p_resid_avg, paste0("results_OLS/",task,"/p_resid_PRE"), hemisphere='both')
			writeCIFTIs(AR1_resid_avg, paste0("results_OLS/",task,"/AR_resid_PRE"), hemisphere='both')
			writeCIFTIs(SE_resid_avg, paste0("results_OLS/",task,"/SE_resid_PRE"), hemisphere='both')

		}

		############################################################
		## CREATE SPARSE PREWHITENING MATRIX FOR 6K MESH LOCATIONS (30 MIN)
		############################################################

		if(n==6){

			coef_avg1 <- AR_resid_avg[idx.obs1,] #left
			coef_avg2 <- AR_resid_avg[idx.obs2,] #right
			var_avg1 <- var_resid_avg[idx.obs1]
			var_avg2 <- var_resid_avg[idx.obs2]

			ntime <- ncol(data.all)

			#identifies each off-diagonal in a TxT matrix
			offs <- row(diag(ntime)) - col(diag(ntime)) 

			#initialize block diagonal matrix (allow for one extra diagonal in the sqrt matrix)
			template <- bandSparse(ntime, k=0:(p+1), symm=TRUE)
			template_list1 <- rep(list(template), ns1) #left
			template_list2 <- rep(list(template), ns2) #right

			## LOOP OVER HEMISPHERES

			coef_avg <- list(coef_avg1, coef_avg2)
			var_avg <- list(var_avg1, var_avg2)
			template_list <- list(template_list1, template_list2)
			rows.rm <- c(list(c()), list(c()))
			# Compute sqrt inverse covariance matrix for each voxel's AR coefficients
			for(h in 1:2){

				t0 <- Sys.time()
				print(paste0('Hemisphere', h))
				for(v in 1:ns[h]){
					print(v)
					rows <- (1:ntime) + ntime*(v-1)
					rows.rm[[h]] <- c(rows.rm[[h]], rows[1:(p+1)], rows[(ntime-p):ntime]) # keep track of which rows to remove

					ar.v <- coef_avg[[h]][v,]
					var.v <- var_avg[[h]][v]
					if(is.na(ar.v[1])) {
						template_list[[h]][[v]] <- sparseMatrix(i=NULL, j=NULL, dims=c(ntime, ntime))
					} else {
						Inv.v <- getInvCovAR(ar.v, ntime) #inverse correlation (sort of)
						Dinv.v <- diag(rep(1/sqrt(var.v), ntime)) #inverse diagonal SD matrix
						Inv.v <- Dinv.v %*% Inv.v %*% Dinv.v #inverse covariance
						sqrtInv.v <- getSqrtInv(Inv.v) #find V^(1) and take sqrt
						#after the (p+1)th off-diagonal, values of sqrtInv.v very close to zero
						diags <- list(); length(diags) <- p+2
						for(k in 0:(p+1)){
							diags[[k+1]] <- sqrtInv.v[offs==k]
						}
						matband <- bandSparse(ntime, k=0:(p+1), symm=TRUE, diag=diags)
						template_list[[h]][[v]] <- matband
					}
				} #end loop across locations
			print(Sys.time()-t0) 
			} #end loop across hemispheres
			sqrtInv_all1 <- bdiag(template_list[[1]]) #~10 seconds
			sqrtInv_all2 <- bdiag(template_list[[2]]) #~10 seconds
			sqrtInv_all <- c(list(sqrtInv_all1), list(sqrtInv_all2))
			save(sqrtInv_all, rows.rm, file=paste0('prewhitening/',task,'/sqrtinvcovAR.6K.Rdata')) #~1 minute, 87MB
		}

		print(paste0('End ',n,'K Data'))

	} # end loop over 6K, 32K data
		
	############################################################
	## SET UP PERMUTATION TEST FOR CLASSICAL GLM
	############################################################

	keep <- (p+2):(ntime-p-1) #remove first and last (p+1) time points after prewhitening
	ntime2 <- length(keep) #length of timeseries after prewhitening
	time_perm <- matrix(nrow=100, ncol=ntime2)
	for(k in 1:100){
		time_perm[k,] <- sample(1:ntime2)
	}
	save(time_perm, file=paste0('results_OLS/',task,'/permutation/time_perm.Rdata'))

} #end loop over tasks



