#library(devtools)
#install_github('muschellij2/gifti')
library(gifti)
#install.packages('sp')
#install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/testing")
library(INLA)
INLA:::inla.dynload.workaround() #to avoid error on creating mesh
inla.setOption("pardiso.license", "~/pardiso.lic") #required for parallel INLA (much faster)
# inla.update(testing=T)
library(excursions)
library(fields)
library(expm) #sqrtm
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(MASS) #mvrnorm
library(parallel)

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
load('locations/ind_obs.Rdata') #ns.all, idx.obs, idx.obs1, idx.obs2, ns1, ns2, ns
load('locations/mesh.sphere.resamp6K.Rdata') #loc.norm1, loc.norm2, mesh1, mesh2
load('locations/Amat.Rdata') #Amat =list(Amat1, Amat2)
load('locations/spde.Rdata') #spde = list(spde1, spde2), spatial = list(spatial1, spatial2)
mesh <- list(mesh1, mesh2) 
loc.norm <- list(loc.norm1, loc.norm2)
idx.obs <- list(idx.obs1, idx.obs2)

###########################################################
# FOR VISUALIZATION
###########################################################

ts_6K <- paste0(maindir,'locations/ts.6K.dtseries.nii')
ts_32K <- paste0(maindir,'locations/ts.32K.dtseries.nii')
rsurface_32K <- paste0(maindir,'locations/Sphere.32K.R.surf.gii') 
lsurface_32K <- paste0(maindir,'locations/Sphere.32K.L.surf.gii')
rsurface_6K <- paste0(maindir,'locations/Sphere.6k.R.surf.gii')
lsurface_6K <- paste0(maindir,'locations/Sphere.6k.L.surf.gii')

thresholds <- c(0,0.5,1) #activation thresholds
U <- length(thresholds)

for(task in c('MOTOR','GAMBLING')){

	if(task=='MOTOR') { ntime <- 284; K <- 6 }
	if(task=='GAMBLING') { ntime <- 253; K <- 3 }

	###########################################################
	# GET SPARSE JOINT PREWHITENING MATRIX
	###########################################################

	p <- 6
	load(file=paste0('prewhitening/',task,'/sqrtinvcovAR.6K.Rdata')) #sqrtInv_all, rows.rm


	###########################################################
	###########################################################
	# RUN SUBJECT-LEVEL MODELS
	###########################################################
	###########################################################

	#for group models
	y_all_left <- c()
	y_all_right <- c()
	Xmat_all_left <- NULL
	Xmat_all_right <- NULL
	Xmat_list_left <- NULL
	Xmat_list_right <- NULL

	for(mm in 1:M){

		print(mm)
		t00 <- Sys.time()

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
		# Fit Bayesian GLM 
		#############################################

		# MOTOR TASK: 45 MIN / HEMISPHERE
		# GAMBLING TASK: 30 MIN / HEMISPHERE

		# Time Series
	    file.name <- paste("timeseries/",task,'/',subjects[mm],".6K.csv",sep="")
	    data.all <- as.matrix(read.csv(file.name, header=F))
	    local_means <- matrix(rowMeans(data.all, na.rm=TRUE), nrow=nrow(data.all), ncol=ncol(data.all))
	    y.all <- t(100*(data.all - local_means)/local_means) #scale to units of pct local signal change AND CENTER 

		# Regress nuisance parameters & save for group models
		betaX <- invxvx %*% t(X) %*% y.all #(X'X)^{-1} X'Z
		residY <- y.all - X %*% betaX
		y <- residY
		Z <- residZ
		save(y, file=paste0('timeseries/',task,'/y_',subjects[mm],'.Rdata'))
		save(Z, file=paste0('timeseries/',task,'/Z_',subjects[mm],'.Rdata'))
		#load(file=paste0('timeseries/',task,'/y_',subjects[mm],'.Rdata'))
		#load(file=paste0('timeseries/',task,'/Z_',subjects[mm],'.Rdata'))

	    # Loop through Hemispheres
	    for(h in 1:2){

	    	print(paste0('~~~~~~~~~~ HEMISPHERE ',h, ' ~~~~~~~~~~'))

		   	if(task=='MOTOR'){
				if(h==1) cols <- c(1,4:6) # Left: EVs 1, 4, 5, 6
				if(h==2) cols <- c(1:3,6) # Right: EVs 1, 2, 3, 6
			} else { cols <- 1:3 }
			nxz <- length(cols)

	    	y.h <- as.vector(y[,idx.obs[[h]]])
			ix <- 1:(ntime*ns[h])
			iy <- rep(1:ns[h], each = ntime)

			############################################
			# Set up design matrix
			############################################

			nxz <- length(cols)
			nmesh <- length(spatial[[h]])								

			for(k in cols){
				Z_k <- sparseMatrix(ix, iy, x=rep(Z[,k], ns[h])) #don't actually need to rep Z because it will be recycled
				if(k==cols[1]) A <- Z_k else A <- cbind(A, Z_k)
			}
						
			############################################
			# Pre-whiten y and A
			############################################

			A <- sqrtInv_all[[h]] %*% A #~10 seconds
			y.h <- as.vector(sqrtInv_all[[h]] %*% y.h) #<1 seconds

			A <- A[-rows.rm[[h]], ]
			y.h <- y.h[-rows.rm[[h]]]

	    	if(h==1) y_all_left <- c(y_all_left, y.h)
	    	if(h==2) y_all_right <- c(y_all_right, y.h)

	    	if(h==1) {
	    		Xmat_all_left <- rbind(Xmat_all_left, A) #for GE group model
	    		#Xmat_list_left <- c(Xmat_list_left, list(A)) #for FE group model
	    	}
	    	if(h==2) {
	    		Xmat_all_ringht <- rbind(Xmat_all_right, A) #for GE group model
				#Xmat_list_right <- c(Xmat_list_right, list(A)) #for FE group model
			}

			############################################
			# Construct formula and data
			############################################

			BOLD <- list(y=y.h)
			formula <- 'y ~ -1'
			bbeta_list <- vector('list',nxz)
			for(j in 1:nxz){
				bbeta_list[[j]] <- c(rep(NA, nmesh*(j-1)), spatial[[h]], rep(NA, (nxz-j)*nmesh))
				BOLD <- c(BOLD, list(bbeta_list[[j]]))
				formula_k <- paste0('f(bbeta',cols[j],', model=spde[[h]], hyper=list(theta=list(initial=c(-2,2))))')
				formula <- paste0(formula,' + ',formula_k)
			}
			names(BOLD)[-1] <- paste0('bbeta',cols)
			formula <- as.formula(formula)
						 
			############################################
			# Run model
			############################################

			t0 <- Sys.time()
			result <- inla(formula, data=BOLD, control.predictor=list(A=A), 
					       verbose = TRUE, keep = FALSE, num.threads = 4,
					       control.inla = list(strategy = "gaussian", int.strategy = "eb"),
					       control.family=list(hyper=list(prec=list(initial=0.1))),	
	    			  	   control.compute=list(config=TRUE)) #needed for excursions
			print(Sys.time()-t0)
           
			res.beta <- result$summary.random
			res.hyper <- result$summary.hyperpar
		    residual <- y - result$summary.fitted.values$mean[1:length(y)]     
			mu.tmp <- result$misc$theta.mode #for joint group model
			Q.tmp <- solve(result$misc$cov.intern) #for joint group model
			file.name <- paste0("results_sphere/",task,"/result.",subjects[mm],".sphere",h,".resamp6K.RData")
            save(res.beta, res.hyper, residual, mu.tmp, Q.tmp, file=file.name)
            load(file=file.name)


			############################################
			# Visualize estimates
			############################################

			betas.h <- matrix(nrow = length(res.beta[[1]]$mean), ncol = K)
			for(j in 1:nxz){
				k <- cols[j]
				betas.h[,k] <- res.beta[[j]]$mean
			}
			betas.h <- Amat[[h]] %*% betas.h #project from mesh back to original data locations

			if(h==1) betas.mm <- matrix(nrow=nrow(data.all), ncol=K)
			betas.mm[idx.obs[[h]],] <- as.matrix(betas.h)
		} #end loop across hemispheres

		#write CIFTI (6K)	
		names <- paste0("results_sphere/",task,"/beta",1:K,".",subjects[mm],"_6K")	
		writeCIFTIs(betas.mm, names, hemisphere='both', template=ts_6K, data_field = 'x6k')

		#convert to 32K
		fnames_6K <- paste0(maindir,names,'.dtseries.nii')
		fnames_32K <- gsub('6K','32K',fnames_6K)
	 	for(k in 1:K){
			print(k)
			fname_6K <- fnames_6K[k]
			fname_32K <- fnames_32K[k]
			cmd <- paste0(wb_cmd,' -cifti-resample ', fname_6K,' COLUMN ', ts_32K,' COLUMN BARYCENTRIC CUBIC ', 
					fname_32K, ' -left-spheres ',lsurface_6K,' ',lsurface_32K,' -right-spheres ',rsurface_6K,' ',rsurface_32K)
			system(cmd)
		}
	           
	}#end loop across subjects

	### PLOT HYPERPARAMETER ESTIMATES

	param_names <- list()
	if(task=='MOTOR'){
		param_names[[1]] <- c('prec','log_tau1','log_kappa1','log_tau4','log_kappa4','log_tau5','log_kappa5','log_tau6','log_kappa6')
		param_names[[2]] <- c('prec','log_tau1','log_kappa1','log_tau2','log_kappa2','log_tau3','log_kappa3','log_tau6','log_kappa6')
	} else if(task=='GAMBLING'){
		param_names[[1]] <- param_names[[2]] <- c('prec','log_tau1','log_kappa1','log_tau2','log_kappa2','log_tau3','log_kappa3')
	}

	for(h in 1:2){
		for(mm in (1:M)){
			print(mm)
			file.name <- paste0("results_sphere/",task,"/result.",subjects[mm],".sphere",h,".resamp6K.RData")
			load(file=file.name)
			params.mm.h <- res.hyper
			params.mm.h$subject <- mm
			params.mm.h$param <- param_names[[h]]
			if(mm==1) params <- params.mm.h else params <- rbind(params,params.mm.h)
		}
		names(params)[3:5] <- c('Q025','Q50','Q975')
		pdf(paste0('plots/',task,'/hyperparams',h,'.pdf'))
		print(ggplot(params, aes(x=param, y=mean, group=subject)) + geom_point() + 
		geom_linerange(aes(ymin=Q025, ymax=Q975)) + theme(legend.position='bottom'))
		dev.off()
	}


	###################################################################	
	###################################################################
	# GROUP-LEVEL MODELS
	###################################################################
	###################################################################


	###################################################################
	# FULLY BAYESIAN MODEL (10-13 HOURS)
	###################################################################

	#Note: Fully Bayes model can be run for all 20 subjects for gambling task
	#with ~280GB of memory.  For motor task, would need more memory or fewer subjects.

	if(task=='GAMBLING'){

		dir <- paste0('results_sphere/',task,'/pop_full')
		betas.all <- matrix(0, nrow=ns.all, ncol=K)
		probs.all <- array(0, dim=c(ns.all, K, U)) #last dimension is for different activation thresholds

		for(h in 1:2){

			if(task=='MOTOR'){
				if(h==1) cols <- c(1,4:6) # Left: EVs 1, 4, 5, 6
				if(h==2) cols <- c(1:3,6) # Right: EVs 1, 2, 3, 6
			} else { cols <- 1:3 }
			nxz <- length(cols)

			############################################
			# Group Effect (GE) Model
			############################################

			### Construct formula and data

			nmesh <- length(spatial[[h]])								

			if(h==1) BOLD <- list(y=y_all_left) else BOLD <- list(y=y_all_right)
			formula <- 'y ~ -1'
			bbeta_list <- vector('list',K)
			for(k in 1:K){
				bbeta_list[[k]] <- c(rep(NA, nmesh*(k-1)), spatial[[h]], rep(NA, (K-k)*nmesh))
				BOLD <- c(BOLD, list(bbeta_list[[k]]))
				formula_k <- paste0('f(bbeta',cols[k],', model=spde[[h]], hyper=list(theta=list(initial=c(-2,2))))')
				formula <- paste0(formula,' + ',formula_k)
			}
			names(BOLD)[-1] <- paste0('bbeta',cols)
			formula <- as.formula(formula)

			### Run INLA: 10-13 hrs for 20 subjects

			if(h==1) Xmat_all <- Xmat_all_left else Xmat_all <- Xmat_all_right

			inla.setOption("num.threads", 4)
			t0 <- Sys.time()
			result <- inla(formula, data=BOLD, family='gaussian',
					   control.predictor=list(A = Xmat_all, compute=FALSE), 
					   control.compute = list(config = TRUE),
					   control.family=list(hyper=list(prec=list(initial=-0.1))),	
					   control.inla = list(strategy = "gaussian", int.strategy = 'eb'), verbose=T, keep=F)
			print(Sys.time() - t0)

			### Visualize beta estimates

			res.beta <- result$summary.random
			#would need to change looping for motor task
			for(k in 1:K){ 
				betas.all[idx.obs[[h]],k] <- as.vector(Amat[[h]] %*% res.beta[[k]]$mean)
			}

			dir <- paste0('results_sphere/',task,'/pop_full')
			fname <- file.path(dir,'result.pop.sphere.RData')
			save(betas.all, file=fname)

			names_beta <- file.path(dir, paste0("beta",cols,"_6K"))
			writeCIFTIs(betas.all, names_beta, hemisphere='both', template=ts_6K, data_field = 'x6k')
			#convert to 32K
			fnames_6K <- paste0(maindir,names_beta,'.dtseries.nii')
			fnames_32K <- gsub('6K','32K',fnames_6K)
			for(k in 1:K){
				print(k)
				fname_6K <- fnames_6K[k]
				fname_32K <- fnames_32K[k]
				cmd <- paste0(wb_cmd,' -cifti-resample ', fname_6K,' COLUMN ', ts_32K,' COLUMN BARYCENTRIC CUBIC ', 
						fname_32K, ' -left-spheres ',lsurface_6K,' ',lsurface_32K,' -right-spheres ',rsurface_6K,' ',rsurface_32K)
				system(cmd)	
			}

			### Compute Excursion Functions:

			#0.0% = 6-7 min
			#0.5% = 5-6 min
			#1.0% = 4 min
			for(u in 1:U){
				print(paste0('threshold: ',thresholds[u], '%'))
				thr <- thresholds[u]
				for(k in cols){
					t0 <- Sys.time()
					print(beta.k <- paste0('bbeta',k))
					res.exc.k <- excursions.inla(result, name=beta.k, u=thr, type='>', method='QC')
					probs.all[idx.obs[[h]],k,u] <- as.vector(Amat[[h]]%*%res.exc.k$F)
					rm(res.exc.k)
					print(Sys.time()-t0)
					fname <- file.path(dir,'result.pop.sphere.RData')
					save(betas.all, probs.all, file=fname)
				}
			}
		}

		### Visualize activation maps

		for(u in 1:U){

			print(paste0('threshold: ',thresholds[u], '%'))

			#1. Write probs.all and probs.all0 to CIFTI 
			probs.all.u <- probs.all[,,u]
			names_probs <- file.path(dir,paste0("probs",1:K,"_thr",u,"_6K_RH"))
			probs.all.u[probs.all.u < 0.01] <- 0.01
			writeCIFTIs(probs.all.u, names_probs, hemisphere='both', template=ts_6K, data_field = 'x6k')
			fnames_6K <- paste0(maindir,names_probs,'.dtseries.nii')
			fnames_32K <- gsub('6K','32K',fnames_6K)

			#2. Resample to 32K and read into R
			probs.all.u.32K <- matrix(nrow=ns.32K*2, ncol=K)
			for(k in 1:K){
				print(k)
				fname_6K <- fnames_6K[k]
				fname_32K <- fnames_32K[k]
				cmd <- paste0(wb_cmd,' -cifti-resample ', fname_6K,' COLUMN ', ts_32K,' COLUMN BARYCENTRIC CUBIC ', fname_32K, ' -left-spheres ',lsurface_6K,' ',lsurface_32K,' -right-spheres ',rsurface_6K,' ',rsurface_32K)
				system(cmd)	
				probs.k <- readCIFTI_ts(fname_32K)
				probs.all.u.32K[,k] <- probs.k[,1]
			}

			#3. Threshold 32K PPMs at 0.95 and 0.99
			active_99.u <- 1*(probs.all.u.32K >= 0.99)

			#4. Combine activation thresholds
			if(u==1) { active_99 <- active_99.u } else { active_99 <- active_99 + active_99.u }
		}

		#5. Write to CIFTI
		active_99[active_99==0] <- 0.01
		names99 <- file.path(dir,paste0("active",1:K,"_99_32K"))
		writeCIFTIs(active_99, names99, hemisphere='both', template=ts_32K, data_field = 'x32k')

		############################################
		# Fixed Effect (FE) Model (VERY SLOW)
		############################################

		# for(h in 1:2){

		# 	beta1 <- rep(spatial[[h]], M)
		# 	beta2 <- rep(spatial[[h]], M)
		# 	beta3 <- rep(spatial[[h]], M)
		# 	nmesh <- length(spatial[[h]])								

		# 	bbeta1 <- c(beta1, rep(NA, nmesh*M), rep(NA, nmesh*M))
		# 	bbeta2 <- c(rep(NA, nmesh*M), beta2, rep(NA, nmesh*M))
		# 	bbeta3 <- c(rep(NA, nmesh*M), rep(NA, nmesh*M), beta3)
		# 	bbeta_list <- c(list(bbeta1), list(bbeta2), list(bbeta3))

		# 	rep1 <- c(rep(1:M, each = nmesh), rep(NA, nmesh*M), rep(NA, nmesh*M))
		# 	rep2 <- c(rep(NA, nmesh*M), rep(1:M, each = nmesh), rep(NA, nmesh*M))
		# 	rep3 <- c(rep(NA, nmesh*M), rep(NA, nmesh*M), rep(1:M, each = nmesh))

		# 	if(h==1) BOLD <- list(y=y_all_left) else BOLD <- list(y=y_all_right)
		# 	formula <- 'y ~ -1'
		# 	bbeta_list <- vector('list',K)
		# 	for(k in 1:K){
		# 		BOLD <- c(BOLD, list(bbeta_list[[k]]))
		# 		formula_k <- paste0('f(bbeta',cols[k],', model=spde[[h]], replicate=rep',k,', hyper=list(theta=list(initial=c(-2,2))))')
		# 		formula <- paste0(formula,' + ',formula_k)
		# 	}
		# 	names(BOLD)[-1] <- paste0('bbeta',cols)
		# 	formula <- as.formula(formula)

		# 	if(h==1) Xmat_all <- bdiag(Xmat_list_left) else Xmat_all <- bdiag(Xmat_list_right)

		# 	t0 <- Sys.time()
		# 	result <- inla(formula, data=BOLD, family='gaussian',
		# 			   control.predictor=list(A = Xmat_all, compute=FALSE), 
		# 			   control.compute = list(config = TRUE),
		# 			   control.family=list(hyper=list(prec=list(initial=0.1))),	
		# 			   control.inla = list(strategy = "gaussian", int.strategy = 'eb'), verbose = T)
		# 	print(Sys.time() - t0)

		# }

	}

	###################################################################
	###################################################################
	# JOINT BAYESIAN APPROACH (FASTER ALTERNATIVE TO FULLY BAYESIAN APPROACH)
	###################################################################
	###################################################################

	betas.all <- matrix(0, nrow=ns.all, ncol=K)
	probs.all <- array(0, dim=c(ns.all, K, U)) #last dimension is for different activation thresholds

	#30 min per hemisphere (with sampling in parallel)
	for(h in 1:2){
		t0h <- Sys.time()

   		if(task=='MOTOR'){
			if(h==1) cols <- c(1,4:6) # Left: EVs 1, 4, 5, 6
			if(h==2) cols <- c(1:3,6) # Right: EVs 1, 2, 3, 6
		} else { cols <- 1:3 }
		nxz <- length(cols)

		print('Collecting theta posteriors from subject models')
		#0.3-0.7 seconds per subject
		theta.sub <- NULL
		mu.theta.tmp <- Q.theta <- 0
		for(mm in 1:M){
			t0 <- Sys.time()
		  	file.name <- paste("results_sphere/",task,"/result.", subjects[mm],".sphere", h, ".resamp6K.RData", sep = "")
			load(file.name)
			# sum_m Q_m * mu_m
			mu.theta.tmp <- mu.theta.tmp + as.vector(Q.tmp%*%mu.tmp)
			# sum_m Q_m
			Q.theta <- Q.theta + Q.tmp
			theta.sub <- cbind(theta.sub, res.hyper$mode)
			rm(mu.tmp, Q.tmp)
		}
		#(sum_m Q_m)^(-1) * sum_m Q_m * mu_m
		mu.theta <- solve(Q.theta, mu.theta.tmp)

		print('Drawing samples from q(theta|y)')

		nsamp <- 50
		logwt <- rep(NA, nsamp)
		
		theta.tmp <- mvrnorm(nsamp, mu.theta, solve(Q.theta))
		for(i in 1:nsamp){ logwt[i] <- F.logwt(theta.tmp[i,], spde[[h]], mu.theta, Q.theta, M) }

		#weights to apply to each posterior sample of theta
		wt.tmp <- exp(logwt - max(logwt))
		wt <- wt.tmp/(sum(wt.tmp))

		print('Computing cross-products for each subject')

		Xcros.all <- Xycros.all <- vector("list", M)
		#9 seconds per subject
		for(mm in 1:M){

			print(mm)

			## Read response and design matrix after nuisance regession & centering
			load(file=paste0('timeseries/',task,'/y_',subjects[mm],'.Rdata')) #y
			load(file=paste0('timeseries/',task,'/Z_',subjects[mm],'.Rdata')) #Z
			y <- as.vector(y[,idx.obs[[h]]])
			Z <- Z[,cols]
			yZ <- list(y, Z)

			cross.mm <- compute.cross(yZ, Amat[[h]], sqrtInv_all[[h]], rows.rm[[h]])

		    Xcros.all[[mm]] <- cross.mm$Xcros
		    Xycros.all[[mm]] <- cross.mm$Xycros
		}

		print('Computing posterior quantities of beta for each theta')

		## Create index vectors 
		n.mesh <- mesh[[h]]$n
		ind_beta <- list()
		for(k in 1:nxz){
			ind_beta[[k]] <- 1:n.mesh + (k-1)*n.mesh
		}

		#get posterior quantities of beta, conditional on a value of theta
		no_cores <- min(detectCores() - 1, 25)
		cl <- makeCluster(no_cores)
		t0 <- Sys.time()
		#in sequence, 8 min per iteration for motor task,  5-6 min for gambling task (for 20 subjects and 3 activation thresholds)
		#in parallel, 21 min total for gambling task!
		#in parallel, 24 min for motor task!
		#with 50 iterations, we save 50*8 - 25 = 375 min (6.25 hours!)
		beta.post.samps <- parApply(cl, theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde[[h]], K=nxz, M, Xcros.all, Xycros.all, thresholds=thresholds, alpha=0.01, ind_beta=ind_beta)
		print(Sys.time() - t0)
		stopCluster(cl)

		#organize samples
		mu.tot <- matrix(nrow=nxz*n.mesh, ncol=nsamp)
		F.tot <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp)), K)), U) #for each activation threshold and task, a Vx50 matrix
		for(itheta in 1:nsamp){
			mu.tot[,itheta] <- beta.post.samps[[itheta]]$mu
			for(u in 1:U){
				for(k in 1:nxz){
					F.tot[[u]][[k]][,itheta] <- beta.post.samps[[itheta]]$F[[u]][,k]
				}
			}
		}

		print('Computing posterior quantities of beta, summing over theta')

		### Sum over samples using weights, combine hemispheres (< 1 sec)

		#posterior mean
		beta.pop <- as.vector(mu.tot%*%wt)
		for(k in 1:nxz){
			beta.pop.k <- beta.pop[ind_beta[[k]]]
			betas.all[idx.obs[[h]],cols[k]] <- as.vector(Amat[[h]]%*%beta.pop.k)
		}

		#posterior probabilities
		for(u in 1:U){
			for(k in 1:nxz){
				F.pop.uk <- as.vector(F.tot[[u]][[k]]%*%wt)	
				probs.all[idx.obs[[h]],cols[k],u] <- as.vector(Amat[[h]]%*%F.pop.uk)
			}
		}

		print(Sys.time() - t0h)
	} #end loop over hemispheres

	### VISUALIZE RESULTS

	dir <- paste0('results_sphere/',task,'/pop_joint/full_model')

	# BETA MAPS
	names_beta <- file.path(dir, paste0("beta",1:K,"_6K"))
	writeCIFTIs(betas.all, names_beta, hemisphere='both', template=ts_6K, data_field = 'x6k')
	#convert to 32K
	fnames_6K <- paste0(maindir,names_beta,'.dtseries.nii')
	fnames_32K <- gsub('6K','32K',fnames_6K)
	for(k in 1:K){
		print(k)
		fname_6K <- fnames_6K[k]
		fname_32K <- fnames_32K[k]
		cmd <- paste0(wb_cmd,' -cifti-resample ', fname_6K,' COLUMN ', ts_32K,' COLUMN BARYCENTRIC CUBIC ', 
				fname_32K, ' -left-spheres ',lsurface_6K,' ',lsurface_32K,' -right-spheres ',rsurface_6K,' ',rsurface_32K)
		system(cmd)	
	}

	# ACTIVATION MAPS
	for(u in 1:U){

		print(paste0('threshold: ',thresholds[u], '%'))

		#1. Write probs.all and probs.all0 to CIFTI 
		probs.all.u <- probs.all[,,u]
		names_probs <- file.path(dir,paste0("probs",1:K,"_thr",u,"_6K"))
		probs.all.u[probs.all.u < 0.01] <- 0.01
		writeCIFTIs(probs.all.u, names_probs, hemisphere='both', template=ts_6K, data_field = 'x6k')
		fnames_6K <- paste0(maindir,names_probs,'.dtseries.nii')
		fnames_32K <- gsub('6K','32K',fnames_6K)

		#2. Resample to 32K and read into R
		probs.all.u.32K <- matrix(nrow=ns.32K*2, ncol=K)
		for(k in 1:K){
			print(k)
			fname_6K <- fnames_6K[k]
			fname_32K <- fnames_32K[k]
			cmd <- paste0(wb_cmd,' -cifti-resample ', fname_6K,' COLUMN ', ts_32K,' COLUMN BARYCENTRIC CUBIC ', fname_32K, ' -left-spheres ',lsurface_6K,' ',lsurface_32K,' -right-spheres ',rsurface_6K,' ',rsurface_32K)
			system(cmd)	
			probs.k <- readCIFTI_ts(fname_32K)
			probs.all.u.32K[,k] <- probs.k[,1]
		}

		#3. Threshold 32K PPMs at 0.95 and 0.99
		active_99.u <- 1*(probs.all.u.32K >= 0.99)
		active_95.u <- 1*(probs.all.u.32K >= 0.95)

		#4. Combine activation thresholds
		if(u==1) {
			active_99 <- active_99.u 
			active_95 <- active_95.u 
		} else {
			active_99 <- active_99 + active_99.u
			active_95 <- active_95 + active_95.u
		}
	}
	active_99[active_99==0] <- 0.01
	active_95[active_95==0] <- 0.01

	#5. Write to CIFTI
	names99 <- file.path(dir,paste0("active",1:K,"_99_32K"))
	names95 <- file.path(dir,paste0("active",1:K,"_95_32K"))
	writeCIFTIs(active_99, names99, hemisphere='both', template=ts_32K, data_field = 'x32k')
	writeCIFTIs(active_95, names95, hemisphere='both', template=ts_32K, data_field = 'x32k')

	###################################################################
	###################################################################
	# TWO-LEVEL BAYESIAN APPROACH (NOT RECOMMENDED DUE TO OVERSMOOTHING)
	###################################################################
	###################################################################

	# COMBINE SUBJECT-LEVEL RESULTS AND RUN GROUP-LEVEL MODELS

	betas.all <- matrix(0, nrow=ns.all, ncol=K)
	probs.all <- array(0, dim=c(ns.all, K, U)) #last dimension is for different activation thresholds

	inla.setOption("num.threads", 16)

	for(h in 1:2){

		print(paste0('~~~~~~~~~~~~~~~~ HEMISPHERE ',h,' ~~~~~~~~~~~~~~~~'))

		t0h <- Sys.time()

    	if(task=='MOTOR'){
			if(h==1) cols <- c(1,4:6) # Left: EVs 1, 4, 5, 6
			if(h==2) cols <- c(1:3,6) # Right: EVs 1, 2, 3, 6
		} else { cols <- 1:3 }
		nxz <- length(cols)

		######################################################
		# COMBINE SUBJECT-LEVEL RESULTS
		######################################################

		print('Gathering Subject-Level Results')

		betas.tot <- matrix(nrow = ns[h]*M, ncol = K)
		sds.tot <- matrix(nrow = ns[h]*M, ncol = K)

		for(mm in 1:M){

			print(mm)
			file.name <- paste0("results_sphere/",task,"/result.",subjects[mm],".sphere",h,".resamp6K.RData")
			load(file.name)

			betas <- matrix(nrow = length(res.beta$bbeta1$mean), ncol = K)
			sds <- matrix(nrow = length(res.beta$bbeta1$mean), ncol = K)

			for(k in cols){
				colname <- paste0('bbeta',k)
				betas[,k] <- res.beta[[colname]]$mean
				sds[,k] <- res.beta[[colname]]$sd
			}

			rows <- (1:ns[h]) + (mm-1)*ns[h]
			for(k in 1:K){
				betas.tot[rows,k] <- as.vector(Amat[[h]]%*%betas[,k])
				sds.tot[rows,k] <- as.vector(Amat[[h]]%*%sds[,k])
			}	
		}

		######################################################
		# RUN GROUP-LEVEL MODELS
		######################################################

		mesh.tot <- inla.mesh.2d(loc.norm[[h]], max.edge = 0.07)
		spde.tot <- inla.spde2.matern(mesh.tot)
		node <- mesh.tot$idx$loc	
		bbeta <- rep(node, M) 
		Amat.tot <- inla.spde.make.A(mesh.tot, as.matrix(loc.norm[[h]]))

		beta.pop <- NULL

		result.tot <- list(); length(result.tot) <- K

		for(k in cols){

			print(k)
		    dat.inla <- list(y = betas.tot[,k], x = bbeta, z = 1:dim(sds.tot)[1], s = 1/(sds.tot[,k])^2)
		    formula <- y ~ -1 + f(x, model = spde.tot) + f(z, model = 'iid', hyper = list(theta=list(scale=s)))
		    #formula0 <- y ~ -1 + f(x, model = spde.tot)
		    #formula1 <- y ~ -1 + f(x, model = spde.tot) + f(z, model = 'iid')
		    
		    print(paste0('Fitting Model of Column ',k,' with INLA'))

		    #7-8 minutes
		    t0 <- Sys.time()
	    	result <- inla(formula, data=dat.inla, control.compute=list(config=TRUE)) 
	    	#result0 <- inla(formula0, data=dat.inla, control.compute=list(config=TRUE)) 
	    	#result1 <- inla(formula1, data=dat.inla, control.compute=list(config=TRUE)) 
	    	print(Sys.time()-t0)

	     	result.tot[[k]] <- result
			mu.post <- result$summary.random$x$mean
			mu.post <- as.vector(Amat.tot%*%mu.post) #resample to data size
			betas.all[idx.obs[[h]],k] <- mu.post
			# var.post <- (result$summary.random$x$sd)^2
			# var.post0 <- (result0$summary.random$x$sd)^2
			# var.post1 <- (result1$summary.random$x$sd)^2
			# mar.var.2level[idx.obs[[h]],i] <- as.vector(Amat.tot%*%var.post)
			# mar.var.2level0[idx.obs[[h]],i] <- as.vector(Amat.tot%*%var.post0)
			# mar.var.2level1[idx.obs[[h]],i] <- as.vector(Amat.tot%*%var.post1)

			print(paste0('Identifying Active Regions for Column ',k))

			#0.0% = 6-7 min
			#0.5% = 5-6 min
			#1.0% = 4 min
			for(u in 1:U){
				t0 <- Sys.time()
				thr <- thresholds[u]
				res.exc <- excursions.inla(result, name='x', u=thr, type='>', method='QC')
				
				#joint posterior probabilities
			  	F.post <- as.vector(Amat.tot%*%res.exc$F)
				probs.all[idx.obs[[h]],k,u] <- F.post
				print(Sys.time()-t0)
			}

		}

		print(Sys.time() - t0h)

		fname <- paste0('results_sphere/',task,'/pop_2level/original_approach/result.full.sphere',h,'.RData')
		save(result.tot, file = fname)
		fname <- paste0('results_sphere/',task,'/pop_2level/original_approach/result.pop.sphere',h,'.RData')
		save(beta.pop, mesh.tot, file=fname)
	}

	dir <- paste0('results_sphere/',task,'/pop_2level/original_approach')
	fname <- file.path(dir,'result.pop.sphere.RData')
	save(betas.all, probs.all, file=fname)
	#load(file=fname)

	#######################
	### VISUALIZE TWO-LEVEL MODEL RESULTS
	#######################

	# BETA MAPS 
	names_beta <- file.path(dir,paste0("beta",1:K,"_6K"))
	writeCIFTIs(betas.all, names_beta, hemisphere='both', template=ts_6K, data_field = 'x6k')
	#convert to 32K
	fnames_6K <- paste0(maindir,names_beta,'.dtseries.nii')
	fnames_32K <- gsub('6K','32K',fnames_6K)
	for(k in 1:K){
		print(k)
		fname_6K <- fnames_6K[k]
		fname_32K <- fnames_32K[k]
		cmd <- paste0(wb_cmd,' -cifti-resample ', fname_6K,' COLUMN ', ts_32K,' COLUMN BARYCENTRIC CUBIC ', 
				fname_32K, ' -left-spheres ',lsurface_6K,' ',lsurface_32K,' -right-spheres ',rsurface_6K,' ',rsurface_32K)
		system(cmd)	
	}

	# ACTIVATION MAPS 
	for(u in 1:U){
		#1. Write PPMs to 6K CIFTI 
		probs.all.u <- probs.all[,,u]
		names_probs <- file.path(dir,paste0("probs",1:K,"_thr",u,"_6K"))
		probs.all.u[probs.all.u < 0.01] <- 0.01
		writeCIFTIs(probs.all.u, names_probs, hemisphere='both', template=ts_6K, data_field = 'x6k')
		fnames_6K <- paste0(maindir,names_probs,'.dtseries.nii')
		fnames_32K <- gsub('6K','32K',fnames_6K)

		#2. Resample to 32K and read into R
		probs.all.u.32K <- matrix(nrow=ns.32K*2, ncol=K)
		for(k in 1:K){
			print(k)
			fname_6K <- fnames_6K[k]
			fname_32K <- fnames_32K[k]
			cmd <- paste0(wb_cmd,' -cifti-resample ', fname_6K,' COLUMN ', ts_32K,' COLUMN BARYCENTRIC CUBIC ', fname_32K, ' -left-spheres ',lsurface_6K,' ',lsurface_32K,' -right-spheres ',rsurface_6K,' ',rsurface_32K)
			system(cmd)	
			probs.k <- readCIFTI_ts(fname_32K)
			probs.all.u.32K[,k] <- probs.k[,1]
		}

		#3. Threshold 32K PPMs at 0.95 and 0.99
		active_99.u <- 1*(probs.all.u.32K >= 0.99)
		active_95.u <- 1*(probs.all.u.32K >= 0.95)

		#4. Combine both activation thresholds
		if(u==1) {
			active_99 <- active_99.u 
			active_95 <- active_95.u 
		} else {
			active_99 <- active_99 + active_99.u
			active_95 <- active_95 + active_95.u
		}
	}

	active_99[active_99==0] <- 0.01
	active_95[active_95==0] <- 0.01

	#5. Write to CIFTI
	names99 <- file.path(dir,paste0("active",1:K,"_99_32K"))
	names95 <- file.path(dir,paste0("active",1:K,"_95_32K"))
	writeCIFTIs(active_99, names99, hemisphere='both', template=ts_32K, data_field = 'x32k')
	writeCIFTIs(active_95, names95, hemisphere='both', template=ts_32K, data_field = 'x32k')
	
} #end loop over tasks (motor, gambling)


