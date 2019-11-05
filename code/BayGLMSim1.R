# THIS CODE PERFORMS SIMULATION 1 (SINGLE-SUBJECT)

setwd('code') #set this

library(INLA)
INLA:::inla.dynload.workaround() #to avoid error on creating mesh
inla.setOption(pardiso.license = '~/pardiso.lic')
# inla.update(testing=T)
#request a free PARDISO license at https://pardiso-project.org/r-inla/
#then activate it:
#library(INLA)
#inla.setOption(pardiso.license = "/path/to/pardiso.lic”)
#then check that it is working:
#inla.pardiso.check()

library(RColorBrewer) #brewer.pal
library(fields) #image.plot
library(matrixStats) #colVars
library(pROC) #roc
library(excursions) #excursions.inla

#set palette for fMRI images
pal <- brewer.pal(11,'Spectral')
pal <- c(pal[1],pal,pal[11])
rf <- colorRampPalette(rev(pal))   # make colors
r1 <- rf(64)

#set palette for beta maps
pal <- c('blue','turquoise','yellow','orange','red','darkred')
rf <- colorRampPalette(pal)   # make colors
r3 <- rf(64)

source('BayGLMfun.R') #getInvCovAR, getSqrtInv
source('functions_sim.R') #image.nan

########################################
### Setup
########################################

## Mask
mask <- as.matrix(read.table("Mask"))
mask2 <- as.matrix(read.table("Mask2"))
mask <- mask*mask2
ny <- dim(mask)[1]
nx <- dim(mask)[2]
mask_vec <- as.vector(mask)
xy.in <- which(mask==1, arr.ind=TRUE)
xy.in <- xy.in[,2:1]

zlim.mask <- c(0, 1)
obj.mask <- list(x=1:ny, y=1:nx, z=matrix(mask_vec, ny, nx)[ny:1,])
pdf('mask.pdf')
image.plot(obj.mask, xlab='', ylab='',  xlim=c(1,ny), ylim=c(1,nx), zlim=zlim.mask, axis.args = list(cex = 1.5), col=r1)
dev.off()

## True activation maps
q1 <- as.matrix(read.table('Q1_subj1'))
q2 <- as.matrix(read.table('Q2_subj1'))
q3 <- as.matrix(read.table('Q3_subj1'))

s1 <- 4*q1 + 4*q2
s2 <- 2*q2 + 2*q3
s1[s1>0] <- 1
s2[s2>0] <- 1 
s1 <- s1
s2 <- s2

#put NAs outside of the brain mask for s1 and s2
s1[mask==0] <- NA
s2[mask==0] <- NA

## Task activation profiles
z1 <- as.matrix(read.table('Z1.txt'))
z2 <- as.matrix(read.table('Z2.txt'))


##############################################################
##############################################################
### Visualize true beta maps and areas of activation
##############################################################
##############################################################

# TRUE ACTIVATION AMPLITUDES

Time <- length(z1)
X <- cbind(1, z1, z2)
datNoResid <- as.matrix(read.table('DatNoResid_subj1'))
Y <- as.matrix(t(datNoResid[mask_vec==1,]))
local_means <- matrix(colMeans(Y), nrow=nrow(Y), ncol=ncol(Y), byrow=TRUE)
Y <- (Y - local_means)/local_means
Beta_true <- solve(t(X) %*% X) %*% t(X) %*% Y

beta1 <- mask_vec; beta1[mask_vec > 0] <- Beta_true[2,]; beta1[mask_vec==0] <- NA
beta2 <- mask_vec; beta2[mask_vec > 0] <- Beta_true[3,]; beta2[mask_vec==0] <- NA
beta1 <- beta1*100; beta2 <- beta2*100

pdf('beta_true.pdf')
image.nan(matrix(beta1, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
image.nan(matrix(beta2, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
dev.off()

# AREAS OF ACTIVATION

pdf('active_true.pdf')
image.nan(matrix(s1, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white')
image.nan(matrix(s2, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white')
dev.off()


##############################################################
##############################################################
### Fit Models
##############################################################
##############################################################


#sm=0: no smoothing, sm=1: smoothing
for(sm in 1:0){

	suf <- ifelse(sm==1, '_sm', '_unsm')

	if(sm==0) dat <- as.matrix(read.table('DatAR_subj1'))	## AR(1) error, no smoothing
	if(sm==1) dat <- as.matrix(read.table('DatAR_sm_subj1'))	## AR(1) error, 6mm Gaussian smoothing

	########################################
	### Scale and prewhiten data

	Y <- as.matrix(t(dat[mask_vec==1,]))
	local_means <- matrix(colMeans(Y), nrow=nrow(Y), ncol=ncol(Y), byrow=TRUE)
	Y <- 100*(Y - local_means)/local_means
	N <- ncol(Y)

	Beta <- solve(t(X) %*% X) %*% t(X) %*% Y
	Resid <- Y - X %*% Beta

	# Fit AR(1) model to data at each voxel
	p <- 1
	coef <- sigma <- prec <- rep(NA, N)
	for(v in 1:N){
		ar_v <- ar(Resid[,v], method='yw', aic=FALSE, order.max=p)
		coef[v] <- ar_v$ar
		sigma[v] <- sqrt(ar_v$var.pred)
		prec[v] <- (ar_v$var.pred)^(-1)
	}

	# coef_img <- mask; coef_img[mask==0] <- NA; coef_img[mask==1] <- coef
	# sigma_img <- mask; sigma_img[mask==0] <- NA; sigma_img[mask==1] <- sigma
	# prec_img <- mask; prec_img[mask==0] <- NA; prec_img[mask==1] <- prec
	# pdf(paste0('coef_scale',suf,'.pdf'))
	# image.nan(coef_img,  zlim=c(0.1, 0.5), col=r1)
	# image.nan(sigma_img,  zlim=c(0.35, 0.45), col=r1)
	# image.nan(prec_img,  zlim=c(4.8, 8.9), col=r1)
	# dev.off()

	# Average AR coefficient across voxels
	coef <- mean(coef)
	Inv <- getInvCovAR(ar=coef, ntime=Time)
	sqrtInv <- getSqrtInv(Inv)

	##############################################################
	##############################################################
	### Classical GLM 
	##############################################################
	##############################################################

	tkeep <- 3:198 #remove first and last 2 time points after prewhitening
	Time2 <- length(tkeep)

	Y2 <- (sqrtInv %*% Y)[tkeep,]
	X2 <- (sqrtInv %*% X)[tkeep,]
	invxvx <- solve(t(X2) %*% X2)  #(X'Inv X)^{-1}
	Beta <- invxvx %*% t(X2) %*% Y2 #(X' Inv X)^{-1} X' Inv y
	Resid <- Y2 - X2 %*% Beta
	sigmaResid <- colVars(Resid)
	CovBeta <- diag(invxvx)
	SteBeta <- sqrt(matrix(CovBeta, nrow=3, ncol=1) %*% matrix(sigmaResid, nrow=1, ncol=N))
	tstat <- Beta / SteBeta  #compute t statistic = beta/ste(beta) (df=20)

	########################################
	### Visualize beta maps

	beta1 <- mask_vec; beta1[mask_vec > 0] <- Beta[2,]; beta1[mask_vec==0] <- NA
	beta2 <- mask_vec; beta2[mask_vec > 0] <- Beta[3,]; beta2[mask_vec==0] <- NA
	beta1 <- beta1
	beta2 <- beta2

	pdf(paste0('fit_beta_OLS',suf,'.pdf'))
	image.nan(matrix(beta1, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
	image.nan(matrix(beta2, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
	dev.off()

	########################################
	### Do permutation test to control FWER

	t0 <- Sys.time() #10 seconds
	M <- 1000
	maxt0 <- maxt1 <- maxt2 <- rep(0, M)
	for(perm in 1:M){
		
		#shuffle time in fMRI data
		print(perm)
		time.perm <- sample(1:Time2)
		Y2.perm <- Y2[time.perm,]

		#compute t-statistics 
		Beta.perm <- invxvx %*% t(X2) %*% Y2.perm #(X' Inv X)^{-1} X' Inv y
		Resid.perm <- Y2.perm - X2 %*% Beta.perm
		sigmaResid.perm <- colVars(Resid.perm)
		SteBeta.perm <- sqrt(matrix(CovBeta, nrow=3, ncol=1) %*% matrix(sigmaResid.perm, nrow=1, ncol=N))
		tstat.perm <- Beta.perm / SteBeta.perm  
		
		#save maximum test statistic for each beta
		maxt0[perm] <- max(tstat.perm[1,])
		maxt1[perm] <- max(tstat.perm[2,])
		maxt2[perm] <- max(tstat.perm[3,])
	}
	print(Sys.time() - t0)

	FWER1 <- mask_vec; FWER1[mask_vec > 0] <- tstat[2,]; FWER1[mask_vec==0] <- NA
	FWER2 <- mask_vec; FWER2[mask_vec > 0] <- tstat[3,]; FWER2[mask_vec==0] <- NA

	################################################
	### Perform FDR correction & visualize active regions

	# Perform FDR correction with Benjamini-Hochberg
	pvals <- pt(tstat, df=Time2-1, lower.tail=FALSE)
	pvals_FDR <- array(dim=dim(pvals))
	for(k in 1:3){
		pvals_FDR[k,] <- p.adjust(pvals[k,], "BH")
	}
	FDR1 <- mask_vec; FDR1[mask_vec > 0] <- pvals_FDR[2,]; FDR1[mask_vec==0] <- NA
	FDR2 <- mask_vec; FDR2[mask_vec > 0] <- pvals_FDR[3,]; FDR2[mask_vec==0] <- NA

	################################################
	### Compute ROC curve (same for voxel-wise FDR and FWER)

	roc1_OLS_sm <- roc(response=as.vector(s1)[mask_vec==1], predictor=pvals_FDR[2,])
	roc2_OLS_sm <- roc(response=as.vector(s2)[mask_vec==1], predictor=pvals_FDR[3,])
	if(sm==1){ #initialize roc1_OLS
		roc1_OLS <- list(roc1_OLS_sm)
		roc2_OLS <- list(roc2_OLS_sm)
	} else { #add to roc_OLS
		roc1_OLS <- c(roc1_OLS, list(roc1_OLS_sm))
		roc2_OLS <- c(roc2_OLS, list(roc2_OLS_sm))
	}

	################################################
	### Threshold and visualize

	for(cutoff in c(95,99)){

		alpha <- 1 - cutoff/100

		# Threshold FWER t-values
		thresh1 <- quantile(maxt1, cutoff/100)
		thresh2 <- quantile(maxt2, cutoff/100)
		FWER1_thr <- 1*(FWER1 > thresh1)
		FWER2_thr <- 1*(FWER2 > thresh2)

		# Threshold FDR q-values
		FDR1_thr <- 1*(FDR1 < (1-cutoff/100))
		FDR2_thr <- 1*(FDR2 < (1-cutoff/100))

		# Visualize FWER 
		main1 <- get_FPR_FNR(vals= FWER1_thr, mask, truth=s1)
		main2 <- get_FPR_FNR(vals= FWER2_thr, mask, truth=s2)
		pdf(paste0('FWER',cutoff,'_',suf,'.pdf'))
		image.nan(matrix(FWER1_thr, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main1, cex.main=2)
		image.nan(matrix(FWER2_thr, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main2, cex.main=2)
		dev.off()

		# Visualize FDR 
		main1 <- get_FPR_FNR(vals= FDR1_thr, mask, truth=s1)
		main2 <- get_FPR_FNR(vals= FDR2_thr, mask, truth=s2)
		pdf(paste0('FDR',cutoff,'_',suf,'.pdf'))
		image.nan(matrix(FDR1_thr, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main1, cex.main=2)
		image.nan(matrix(FDR2_thr, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main2, cex.main=2)
		dev.off()
	}

	##############################################################
	##############################################################
	### Bayesian GLM 
	##############################################################
	##############################################################

	if(sm==0){

		y <- as.vector(Y) #scaled but not prewhitened data
		ntot <- length(y)

		### Create triangulation and SPDE object
		boundary <- inla.nonconvex.hull(xy.in, resolution = 100) 
		mesh <- inla.mesh.2d(loc = xy.in, boundary = boundary, max.edge = c(2, 4))			
 		save(mesh, file = 'mesh_glm_simu.Rdata')
		# load(file = 'mesh_glm_simu.Rdata')
		pdf('mesh.pdf')
		plot(mesh)	
		dev.off()
		spde <- inla.spde2.matern(mesh)					


		### Create Sparse Z1, Z2 and XX
		ix <- 1:(Time*N)
		iy <- rep(1:N, each = Time)
		plot(z1, xlab='Time', ylab='Activation', type='l')
		plot(z2, xlab='Time', ylab='Activation', type='l')
		Z1 <- sparseMatrix(ix, iy, x = rep(z1, N))
		Z2 <- sparseMatrix(ix, iy, x = rep(z2, N))
		XX <- sparseMatrix(ix, iy, x = 1)
		DD <- Diagonal(n = Time*N, x = 1) # Design matrix for AR(1)
		A <- cbind(DD, XX, Z1, Z2)


		########################################
		### Fit INLA model

		spatial <- mesh$idx$loc	
		nmesh <- length(spatial)

		aalpha <- c(rep(1:Time, N), rep(NA, 3*nmesh))
		repli <- c(rep(spatial, each = Time), rep(NA, 3*nmesh))
		bbeta0 <- c(rep(NA, ntot), spatial, rep(NA, 2*nmesh)) #still need intercept since X not centered
		bbeta1 <- c(rep(NA, ntot), rep(NA, nmesh), spatial, rep(NA, nmesh))
		bbeta2 <- c(rep(NA, ntot), rep(NA, 2*nmesh), spatial)

		BOLD <- list(y=y, aalpha=aalpha, bbeta0=bbeta0, bbeta1=bbeta1, bbeta2=bbeta2)

		formula = y ~ -1 + f(aalpha, model='ar1', replicate=repli, compute = FALSE) + 
						   f(bbeta0, model=spde) + f(bbeta1, model=spde) +
						   f(bbeta2, model=spde) 
			
		inla.setOption("num.threads", 4)

		#different computation strategy
		k <- 1 #ccd only
		# for(k in 1:2){
			if(k==1) strategy <- 'ccd' 
			if(k==2) strategy <- 'eb' #46 minutes

			proj <- inla.mesh.projector(mesh, dims=c(nx,ny), xlim=c(1,nx), ylim=c(1,ny))

			#25 minutes for CCD, 11 min with iid errors (with initial error precision = 1)
			#20 minutes for EB, 9 min with iid errors (with initial error precision = 1)
			t0 <- Sys.time()
			result <- inla(formula, data=BOLD, family='gaussian',
				   control.family=list(hyper=list(prec=list(initial=12,fixed=TRUE))), 
				   control.predictor=list(A=A), 
				   control.compute = list(config = TRUE),
				   control.inla = list(strategy = "gaussian", int.strategy = strategy),
				   verbose=T, keep=F)
			print(Sys.time() - t0)

			fname <- paste0("result_",strategy,".Rdata")
			#save(result, file = fname)
			#load(file = fname)

			########################################
			### Plot estimates

			beta1 <- result$summary.random$bbeta1$mean
			beta2 <- result$summary.random$bbeta2$mean

			pdf(paste0("fit_beta_",strategy,".pdf"))
			obj.beta1 <- list(x=proj$y, y=proj$x, z=t(inla.mesh.project(proj, beta1)))
			obj.beta2 <- list(x=proj$y, y=proj$x, z=t(inla.mesh.project(proj, beta2)))
			obj.beta1$z[mask==0] <- obj.beta2$z[mask==0] <- NA #put NAs outside of brain
			obj.beta1$z <- obj.beta1$z[ny:1,]
			obj.beta2$z <- obj.beta2$z[ny:1,]
			image.nan2(obj.beta1, zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
			image.nan2(obj.beta2, zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
			dev.off()


			########################################
			### Plot regular and joint PPMs

			#10 minutes per task
			system.time(res.exc1 <- excursions.inla(result, name = 'bbeta1', u = 0, type = '>', method = 'QC'))
			system.time(res.exc2 <- excursions.inla(result, name = 'bbeta2', u = 0, type = '>', method = 'QC'))
			fname <- paste0('res_exc_',strategy,'.Rdata')
			save(res.exc1, res.exc2, file = fname)
			#load(fname) #res.exc1, res.exc2

			# Get regular PPM
			regPPM1 <- t(inla.mesh.project(proj,res.exc1$rho))
			regPPM2 <- t(inla.mesh.project(proj,res.exc2$rho))
			#regPPM1[mask[ny:1,]==0] <- NA #put NAs outside of brain
			#regPPM2[mask[ny:1,]==0] <- NA #put NAs outside of brain

			# Get joint PPM
			jointPPM1 <- t(inla.mesh.project(proj,res.exc1$F))
			jointPPM2 <- t(inla.mesh.project(proj,res.exc2$F))
			#jointPPM1[mask[ny:1,]==0] <- NA #put NAs outside of brain
			#jointPPM2[mask[ny:1,]==0] <- NA #put NAs outside of brain

			# Make histograms of joint and regular PPM by activated/not activated
			act1 <- as.vector(s1)[mask_vec==1]
			act2 <- as.vector(s1)[mask_vec==1]
			reg1 <- as.vector(regPPM1)[mask_vec==1]
			reg2 <- as.vector(regPPM2)[mask_vec==1]
			joint1 <- as.vector(jointPPM1)[mask_vec==1]
			joint2 <- as.vector(jointPPM2)[mask_vec==1]

			p_reg1 <- hist(reg1, breaks=seq(0,1,0.05), plot=FALSE) #all voxels
			p_reg2 <- hist(reg2, breaks=seq(0,1,0.05), plot=FALSE) #all voxels
			p_joint1 <- hist(joint1, breaks=seq(0,1,0.05), plot=FALSE) #all voxels
			p_joint2 <- hist(joint2, breaks=seq(0,1,0.05), plot=FALSE) #all voxels
			max_old1 <- max(max(p_reg1$counts), max(p_joint1$counts))
			max_old2 <- max(max(p_reg2$counts), max(p_joint2$counts))
			p_reg1$counts <- p_reg1$counts/max_old1*0.7
			p_reg2$counts <- p_reg2$counts/max_old2*0.7
			p_joint1$counts <- p_joint1$counts/max_old1*0.7
			p_joint2$counts <- p_joint2$counts/max_old2*0.7

			pdf('reg_vs_jointPPM.pdf', width=5, height=5.5)
			#Activation 1
			plot(reg1, joint1, pch='*', xlab = 'Marginal PPM', ylab = 'Joint PPM', xlim=c(0,1), ylim=c(0,1))#'Excursion Function Values')
			plot(p_reg1, col=adjustcolor('red', alpha.f=0.5), border=T, add=T)
			plot(p_joint1, col=adjustcolor('blue', alpha.f=0.5), border=T, add=T)
			points(reg1, joint1, pch='*')
			legend('topleft', legend=c('Marginal PPM','Joint PPM'), fill=c(adjustcolor('red', alpha.f=0.5),adjustcolor('blue', alpha.f=0.5)))
			#Activation 2
			plot(reg2, joint2, pch='*', xlab = 'Marginal PPM', ylab = 'Joint PPM', xlim=c(0,1), ylim=c(0,1))#'Excursion Function Values')
			plot(p_reg2, col=adjustcolor('red', alpha.f=0.5), border=T, add=T)
			plot(p_joint2, col=adjustcolor('blue', alpha.f=0.5), border=T, add=T)
			points(reg2, joint2, pch='*')
			legend('topleft', legend=c('Marginal PPM','Joint PPM'), fill=c(adjustcolor('red', alpha.f=0.5),adjustcolor('blue', alpha.f=0.5)))
			dev.off()


			### Threshold PPM

			for(alpha in c(0.01, 0.05)){

				print(paste0('alpha = ', alpha))
				level <- 100*(1-alpha)

				#Regular PPM
				beta1_set_reg <- 1*(res.exc1$rho > 1-alpha)
				beta2_set_reg <- 1*(res.exc2$rho > 1-alpha)
				beta1_set_reg <- t(inla.mesh.project(proj,beta1_set_reg))
				beta2_set_reg <- t(inla.mesh.project(proj,beta2_set_reg))
				beta1_set_reg[mask==0] <- NA
				beta2_set_reg[mask==0] <- NA

				#Joint PPM
				beta1_set_joint <- 1*(res.exc1$F > 1-alpha)
				beta2_set_joint <- 1*(res.exc2$F > 1-alpha)
				beta1_set_joint <- t(inla.mesh.project(proj,beta1_set_joint))
				beta2_set_joint <- t(inla.mesh.project(proj,beta2_set_joint))

				#Plot Joint PPM
				obj1.set <- list(x=proj$y, y=proj$x, z=beta1_set_joint)
				obj2.set <- list(x=proj$y, y=proj$x, z=beta2_set_joint)
				obj1.set$z[mask==0] <- obj2.set$z[mask==0] <- NA #put NAs outside of brain
				obj1.set$z <- obj1.set$z[ny:1,]
				obj2.set$z <- obj2.set$z[ny:1,]
				main1 <- get_FPR_FNR(vals= beta1_set_joint, mask, truth=s1)
				main2 <- get_FPR_FNR(vals= beta2_set_joint, mask, truth=s2)
				pdf(paste0('ppm_joint_',level,'_',strategy,'.pdf'))
				image.nan2(obj1.set, zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main1, cex.main=2)
				image.nan2(obj2.set, zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main2, cex.main=2)
				dev.off()

				#Plot Marginal PPM
				obj1.set <- list(x=proj$y, y=proj$x, z=beta1_set_reg)
				obj2.set <- list(x=proj$y, y=proj$x, z=beta2_set_reg)
				obj1.set$z[mask==0] <- obj2.set$z[mask==0] <- NA #put NAs outside of brain
				obj1.set$z <- obj1.set$z[ny:1,]
				obj2.set$z <- obj2.set$z[ny:1,]
				main1 <- get_FPR_FNR(vals= beta1_set_reg, mask, truth=s1)
				main2 <- get_FPR_FNR(vals= beta2_set_reg, mask, truth=s2)
				pdf(paste0('ppm_reg_',level,'_',strategy,'.pdf'))
				image.nan2(obj1.set, zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main1, cex.main=2)
				image.nan2(obj2.set, zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main2, cex.main=2)
				dev.off()
				
			}
			
			# BELOW HERE, CCD ONLY
			if(k==1){

				#make image of joint PPM 
				beta1_set <- res.exc1$F
				beta2_set <- res.exc2$F
				vals1_joint <- t(inla.mesh.project(proj,beta1_set))
				vals2_joint <- t(inla.mesh.project(proj,beta2_set))
				vals1_joint[mask==0] <- NA #put NAs outside of brain
				vals2_joint[mask==0] <- NA #put NAs outside of brain

				#ROC
				roc1_joint <- roc(response=as.vector(s1)[mask_vec==1], predictor=as.vector(vals1_joint)[mask_vec==1])
				roc2_joint <- roc(response=as.vector(s2)[mask_vec==1], predictor=as.vector(vals2_joint)[mask_vec==1])

				AUC1 <- round(c(roc1_joint$auc, roc1_OLS[[1]]$auc, roc1_OLS[[2]]$auc), 4)
				AUC2 <- round(c(roc2_joint$auc, roc2_OLS[[1]]$auc, roc2_OLS[[2]]$auc), 4)
			
				legend <- c('Bayesian GLM', 'Classical GLM (smoothing)', 'Classical GLM (no smoothing)')
				legend1 <- paste0(legend, ' (AUC = ', AUC1, ')')
				legend2 <- paste0(legend, ' (AUC = ', AUC2, ')')
			
				pdf('roc1.pdf', width=5, height=5)
				plot(roc1_joint)
				lines(roc1_OLS[[1]], lty=2)
				lines(roc1_OLS[[2]], lty=3)
				legend('bottomright', legend=legend1, lty=c(1,2,3), cex=0.75)
				dev.off()

				pdf('roc2.pdf', width=5, height=5)
				plot(roc2_joint)
				lines(roc2_OLS[[1]], lty=2)
				lines(roc2_OLS[[2]], lty=3)
				legend('bottomright', legend=legend2, lty=c(1,2,3), cex=0.75)
				dev.off()
			}
		#} #end CCD, EB
			
	} #end Bayesian GLM

} #end smoothing loop

	
















