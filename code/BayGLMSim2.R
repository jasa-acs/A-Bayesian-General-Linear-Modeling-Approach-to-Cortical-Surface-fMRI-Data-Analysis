# THIS CODE PERFORMS SIMULATION 2 (MULTI-SUBJECT)

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

library(RColorBrewer) #brewer.pal()
library(fields) #image.plot()
library(excursions)
library(matrixStats) #colVars
library(MASS) #mvrnorm
library(parallel) #detectCores()

#set palette for fMRI images
pal <- brewer.pal(11,'Spectral')
pal <- c(pal[1],pal,pal[11])
rf <- colorRampPalette(rev(pal))   # make colors
r1 <- rf(64)

#set palette for beta maps
pal <- c('blue','turquoise','yellow','orange','red','darkred')
rf <- colorRampPalette(pal)   # make colors
r3 <- rf(64)

source('BayGLMfun.R')
source('functions_sim.R')

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
N <- sum(mask_vec)
xy.in <- which(mask==1, arr.ind=TRUE)
xy.in <- xy.in[,2:1]

zlim.mask <- c(0, 1)
obj.mask <- list(x=1:ny, y=1:nx, z=matrix(mask_vec, ny, nx)[ny:1,])
# pdf('mask.pdf')
# image.plot(obj.mask, xlab='', ylab='',  xlim=c(1,ny), ylim=c(1,nx), zlim=zlim.mask, axis.args = list(cex = 1.5), col=r1)
# dev.off()

## Task activation profiles
z1 <- as.matrix(read.table('Z1.txt'))
z2 <- as.matrix(read.table('Z2.txt'))
z1 <- z1 - mean(z1)
z2 <- z2 - mean(z2)
Time <- length(z1)


## Create Sparse Z1 and Z2 (same for all subjects)
ix <- 1:(Time*N)
iy <- rep(1:N, each = Time)
Z1 <- sparseMatrix(ix, iy, x = rep(z1, N))
Z2 <- sparseMatrix(ix, iy, x = rep(z2, N))
Xmat <- cbind(Z1, Z2)


## Mesh and SPDE object
boundary <- inla.nonconvex.hull(xy.in, resolution = 100) 
mesh <- inla.mesh.2d(loc = xy.in, boundary = boundary, max.edge = c(2, 4))			
save(mesh, file = 'mesh_glm_simu.Rdata')
#load(file = 'mesh_glm_simu.Rdata')
pdf('mesh.pdf')
plot(mesh)	
dev.off()

spde <- inla.spde2.matern(mesh)					
proj <- inla.mesh.projector(mesh, dims=c(nx,ny), xlim=c(1,nx), ylim=c(1,ny))
Amat <- inla.spde.make.A(mesh, loc=xy.in)

spatial <- mesh$idx$loc	
nmesh <- length(spatial)
bbeta1 <- c(spatial, rep(NA, nmesh))
bbeta2 <- c(rep(NA, nmesh), spatial)

########################################
### Loop over subjects
########################################

M <- 10 #number of subjects to include in analysis 

beta1_all <- matrix(0, nx*ny, M)
beta2_all <- matrix(0, nx*ny, M)

s1_all <- s2_all <- matrix(0, ny, nx)

for(i in 1:M){

	print(i)

	## True activation maps
	q1 <- as.matrix(read.table(paste0('Q1_subj',i)))
	q2 <- as.matrix(read.table(paste0('Q2_subj',i)))
	q3 <- as.matrix(read.table(paste0('Q3_subj',i)))

	s1 <- 4*q1 + 4*q2
	s2 <- 2*q2 + 2*q3
	s1[s1>0] <- 1
	s2[s2>0] <- 1 

	#put NAs outside of the brain mask for s1 and s2
	s1[mask==0] <- NA
	s2[mask==0] <- NA
	s1_all <- s1_all + s1
	s2_all <- s2_all + s2

	##############################################################
	### Visualize true beta maps and areas of activation
	##############################################################

	# TRUE ACTIVATION AMPLITUDES

	X <- cbind(z1, z2)
	datNoResid <- as.matrix(read.table(paste0('Dat_noResid_subj',i)))
	Y <- as.matrix(t(datNoResid[mask_vec==1,]))
	local_means <- matrix(colMeans(Y), nrow=nrow(Y), ncol=ncol(Y), byrow=TRUE)
	Y <- (Y - local_means)/local_means
	Beta_true <- solve(t(X) %*% X) %*% t(X) %*% Y

	beta1 <- mask_vec; beta1[mask_vec > 0] <- Beta_true[1,]; beta1[mask_vec==0] <- NA
	beta2 <- mask_vec; beta2[mask_vec > 0] <- Beta_true[2,]; beta2[mask_vec==0] <- NA
	beta1 <- beta1*100; beta2 <- beta2*100
	beta1_all[,i] <- beta1
	beta2_all[,i] <- beta2

	pdf(paste0('beta_true_subj',i,'.pdf'))
	image.nan(matrix(beta1, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
	image.nan(matrix(beta2, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
	dev.off()

}

beta1_grp_true <- rowMeans(beta1_all)
beta2_grp_true <- rowMeans(beta2_all)

pdf('beta_true_grp.pdf')
image.nan(matrix(beta1_grp_true, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
image.nan(matrix(beta2_grp_true, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
dev.off()

active1_grp <- 1*(s1_all > 1)
active2_grp <- 1*(s2_all > 1)

pdf('active_true_grp.pdf')
image.nan(active1_grp, zlim=c(0, 1), col=c('lightgray','red'), na.color='white')
image.nan(active2_grp, zlim=c(0, 1), col=c('lightgray','red'), na.color='white')
dev.off()


##############################################################
### Fit Subject-Level Models (3 MIN/SUBJECT)
##############################################################

y_all <- c()
comptime <- rep(NA, M)
for(i in 1:M){

	print(paste0('~~~~~~~~~~~~~~~~~~~~~ SUBJECT ',i,' ~~~~~~~~~~~~~~~~~~~~~'))

	########################################
	### Read and scale data

	dat <- as.matrix(read.table(paste0('Dat_subj',i)))	
	Y <- as.matrix(t(dat[mask_vec==1,]))
	local_means <- matrix(colMeans(Y), nrow=nrow(Y), ncol=ncol(Y), byrow=TRUE) #center each voxel timecourse
	Y <- 100*(Y - local_means)/local_means
	y <- as.vector(Y) 
	y_all <- c(y_all, y)
	ntot <- length(y)

	########################################
	### Fit subject-level INLA model

	BOLD <- list(y=y, bbeta1=bbeta1, bbeta2=bbeta2)
	formula = y ~ -1 + f(bbeta1, model=spde) + f(bbeta2, model=spde) 

	#4 MINUTES
	inla.setOption("num.threads", 4)
	t0 <- Sys.time()
	result <- inla(formula, data=BOLD, family='gaussian',
		   control.predictor=list(A=Xmat), 
		   control.compute = list(config = TRUE),
		   control.family=list(hyper=list(prec=list(initial=1))),	
		   #control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
		   verbose=T, keep=F)
	print(comptime[i] <- Sys.time() - t0)

	########################################
	### Save stuff

	res.beta <- result$summary.random
	res.hyper <- result$summary.hyperpar
	mu.tmp <- result$misc$theta.mode #for joint group model
	Q.tmp <- solve(result$misc$cov.intern) #for joint group model
	save(res.beta, res.hyper, mu.tmp, Q.tmp, file = paste0("result_subj",i,".Rdata"))
	#load(file = fname)

	########################################
	### Plot estimates

	beta1 <- res.beta$bbeta1$mean
	beta2 <- res.beta$bbeta2$mean

	pdf(paste0("fit_beta_subj",i,".pdf"))
	obj.beta1 <- list(x=proj$y, y=proj$x, z=t(inla.mesh.project(proj, beta1)))
	obj.beta2 <- list(x=proj$y, y=proj$x, z=t(inla.mesh.project(proj, beta2)))
	obj.beta1$z[mask==0] <- obj.beta2$z[mask==0] <- NA #put NAs outside of brain
	obj.beta1$z <- obj.beta1$z[ny:1,]
	obj.beta2$z <- obj.beta2$z[ny:1,]
	image.nan2(obj.beta1, zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
	image.nan2(obj.beta2, zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)))
	dev.off()

} #end loop over subjects

#save(y_all, spde, spatial, nmesh, M, Xmat, Z1, Z2, file='group_stuff.Rdata')

##############################################################
##############################################################
### Fit Fully Bayesian Model
##############################################################
##############################################################

inla.setOption("num.threads", 4)

###############################
### GROUP EFFECT MODEL: 25-60 MINUTES
###############################

Xmat_all <- kronecker(matrix(1, M, 1), Xmat)

BOLD <- list(y = y_all, bbeta1 = bbeta1, bbeta2 = bbeta2)
formula <- y ~ -1 + f(bbeta1, model=spde) + f(bbeta2, model=spde) 
		
t0 <- Sys.time()
result <- inla(formula, data=BOLD, family='gaussian',
		   control.predictor=list(A = Xmat_all), 
		   control.compute = list(config = TRUE),
		   control.family=list(hyper=list(prec=list(initial=0.1))),	
		   control.inla = list(strategy = "gaussian", int.strategy = 'eb'), verbose=T, keep=F)
print(Sys.time() - t0)

### COMPUTE EXCURSION FUNCTION: 5 MIN PER TASK
system.time(res.exc1 <- excursions.inla(result, name='bbeta1', u=0, type='>', method='QC'))
system.time(res.exc2 <- excursions.inla(result, name='bbeta2', u=0, type='>', method='QC'))

### SAVE RESULTS
betas.all <- Amat %*% cbind(result$summary.random$bbeta1$mean, result$summary.random$bbeta2$mean)
probs.all <- Amat %*% cbind(res.exc1$F, res.exc2$F)
save(betas.all, probs.all, file='results_grp_full.Rdata')
#load(file='results_grp_full.Rdata')


###############################
### RANDOM EFFECTS MODEL: 5 HR WITH SEPARATE RE PRECISION
###############################

## Random slopes
r1 <- rep(spatial, M)
rr1 <- c(r1, r1, rep(NA, 2*nmesh))
rep1 <- c(rep(1:M, each = nmesh), rep(1:M, each = nmesh), rep(NA, 2*nmesh))
bbeta1 <- c(rep(NA, 2*length(r1)), spatial, rep(NA, nmesh))
bbeta2 <- c(rep(NA, 2*length(r1)), rep(NA, nmesh), spatial)

Xmat2 <- kronecker(matrix(1, M, 1), Xmat)
Xrr1 <- kronecker(matrix(1, M, M), Z1)
Xrr2 <- kronecker(matrix(1, M, M), Z2)
Xmat_all <- cbind(Xrr1, Xrr2, Xmat2)

BOLD <- list(y = y_all, bbeta1 = bbeta1, bbeta2 = bbeta2, rr1 = rr1)#, rr2 = rr2)
formula <- y ~ -1 + f(rr1, model = 'iid', replicate = rep1) + 
					#f(rr2, model = 'iid', replicate = rep2) + 
					f(bbeta1, model = spde) + 
					f(bbeta2, model = spde) 

t0 <- Sys.time()
result <- inla(formula, data=BOLD, family='gaussian',
		   control.predictor=list(A = Xmat_all), 
		   control.compute = list(config = TRUE),
		   control.family=list(hyper=list(prec=list(initial=0.1))),	
		   control.inla = list(strategy = "gaussian", int.strategy = 'eb'), verbose=T, keep=F)
print(Sys.time() - t0)

### COMPUTE EXCURSION FUNCTION: 5 MIN PER TASK
system.time(res.exc1 <- excursions.inla(result, name='bbeta1', u=0, type='>', method='QC'))
system.time(res.exc2 <- excursions.inla(result, name='bbeta2', u=0, type='>', method='QC'))

### SAVE RESULTS
betas.all <- Amat %*% cbind(result$summary.random$bbeta1$mean, result$summary.random$bbeta2$mean)
probs.all <- Amat %*% cbind(res.exc1$F, res.exc2$F)
save(betas.all, file='results_grp_RE.Rdata')
save(betas.all, probs.all, file='results_grp_RE.Rdata')
#load(file='results_grp_RE.Rdata')

###############################
### FIXED EFFECTS MODEL: 1.7 HOURS
###############################

beta1 <- rep(spatial, M)
beta2 <- rep(spatial, M)
bbeta1 <- c(beta1, rep(NA, length(beta2)))
bbeta2 <- c(rep(NA, length(beta1)), beta2)
rep1 <- c(rep(1:M, each = nmesh), rep(NA, length(beta2)))
rep2 <- c(rep(NA, length(beta1)), rep(1:M, each = nmesh))
Xmat_all <- kronecker(diag(1, M), Xmat)

## A matrix for linear combinations
A <- kronecker(t(rep(1, M)), Diagonal(x = 1/M, n = mesh$n))

BOLD <- list(y = y_all, bbeta1 = bbeta1, bbeta2 = bbeta2)
formula <- y ~ -1 + f(bbeta1, model = spde, replicate = rep1) + f(bbeta2, model = spde, replicate = rep2) 

t0 <- Sys.time()
result <- inla(formula, data=BOLD, family='gaussian',
		   control.predictor=list(A = Xmat_all, compute=FALSE), 
		   control.compute = list(config = TRUE),
		   control.family=list(hyper=list(prec=list(initial=0.1))),	
		   control.inla = list(strategy = "gaussian", int.strategy = 'eb'), verbose = T)
print(Sys.time() - t0)

### DRAW POSTERIOR SAMPLES FROM FULL BETA VECTOR (betas for all subjects)

## POSTERIOR PRECISION

#compute prior precision matrix of beta given hyperparameter estimates
Q_prior1 <- inla.spde2.precision(spde, theta=result$summary.hyperpar$mode[2:3]) #bbeta1
Q_prior2 <- inla.spde2.precision(spde, theta=result$summary.hyperpar$mode[4:5]) #bbeta2
Q_prior <- bdiag(c(rep(list(Q_prior1), M), rep(list(Q_prior2), M)))

#compute residual precision matrix (diagonal)
prec_err <- result$summary.hyperpar$mode[1]
dim_err <- nrow(Xmat_all)
Q_error <- Diagonal(dim_err, prec_err)

#compute posterior precision matrix of beta
Amat_all <- bdiag(rep(list(Amat), M*2))
Q_post <- Q_prior + t(Amat_all) %*% t(Xmat_all) %*% Q_error %*% Xmat_all %*% Amat_all

## POSTERIOR MEAN

mu <- result$misc$configs$config[[1]]$mean
#exactract indices for the parts we are interested in (the TWO latent components)
tags <- result$misc$configs$contents$tag #check to verify which elements correspond to latent fields
starts <- result$misc$configs$contents$start
lengths <- result$misc$configs$contents$length
latent.ind <- starts[3]:(starts[3]+sum(lengths[3:4])-1)
mu_post <- mu[latent.ind]

#compute posterior mean of linear combination and save
inds.bbeta1 <- rep(rep(1:0, each=mesh$n), M)
inds.bbeta2 <- rep(rep(0:1, each=mesh$n), M)
mu.bbeta1.LC <- as.vector(A %*% mu_post[inds.bbeta1==1])
mu.bbeta2.LC <- as.vector(A %*% mu_post[inds.bbeta2==1])
betas.all <- Amat %*% cbind(mu.bbeta1.LC, mu.bbeta2.LC)
save(betas.all, file='results_grp_FE.Rdata')


### DRAW SAMPLES: 2.5 MIN/BUCKET (CAN BE UP TO 12 MIN)

#do sampling in several buckets to avoid initializing full matrix of samples x latent locations (all subjects and tasks)
nsamp <- 1000 #samples per bucket
nbuckets <- 10 #number of buckets
bbeta1.LC.samps <- bbeta2.LC.samps <- NULL
for(b in 1:nbuckets){
	t0 <- Sys.time() 
	x <- inla.qsample(n=nsamp,Q_post,mu=mu_post, num.threads=4)
	bbeta1.LC.samps <- cbind(bbeta1.LC.samps, A %*% x[inds.bbeta1==1,]) #apply LC matrix A to beta1's
	bbeta2.LC.samps <- cbind(bbeta2.LC.samps, A %*% x[inds.bbeta2==1,]) #apply LC matrix A to beta1's
	rm(x)
	print(Sys.time() - t0)
}

### COMPUTE EXCURSION FUNCTION: 10 SEC PER TASK
system.time(res.exc1 <- excursions.mc(bbeta1.LC.samps, alpha=0.01, u=0, type='>'))
system.time(res.exc2 <- excursions.mc(bbeta2.LC.samps, alpha=0.01, u=0, type='>'))
probs1 <- res.exc1$F
probs2 <- res.exc2$F
probs.all <- Amat %*% cbind(res.exc1$F, res.exc2$F)
save(betas.all, probs.all, file='results_grp_FE.Rdata')

########################################
# VISUALIZE RESULTS
########################################


### VISUALIZE BETAS WITH MSE

beta1_grp <- mask_vec; beta1_grp[mask_vec > 0] <- betas.all[,1]
beta2_grp <- mask_vec; beta2_grp[mask_vec > 0] <- betas.all[,2]
beta1_grp[mask_vec==0] <- beta2_grp[mask_vec==0] <- NA

#MSE
MSE1 <- mean((beta1_grp[beta1_grp_true>0] - beta1_grp_true[beta1_grp_true>0])^2, na.rm=TRUE)
MSE2 <- mean((beta2_grp[beta2_grp_true>0] - beta2_grp_true[beta2_grp_true>0])^2, na.rm=TRUE)
main1 <- paste0('Active MSE: ', format(round(MSE1, 4), nsmall=4))
main2 <- paste0('Active MSE: ', format(round(MSE2, 4), nsmall=4))

pdf("fit_beta_grp_FE.pdf")
image.nan(matrix(beta1_grp, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)), main=main1, cex.main=2)
image.nan(matrix(beta2_grp, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)), main=main2, cex.main=2)
dev.off()


### VISUALIZE AREAS OF ACTIVATION WITH FPR/FNR
			
act1_grp <- mask_vec; act1_grp[mask_vec > 0] <- (probs.all[,1] >= .99); act1_grp[mask_vec==0] <- NA
act2_grp <- mask_vec; act2_grp[mask_vec > 0] <- (probs.all[,2] >= .99); act2_grp[mask_vec==0] <- NA
main1 <- get_FPR_FNR_vec(vals=act1_grp[mask_vec > 0], truth=active1_grp[mask==1])
main2 <- get_FPR_FNR_vec(vals=act2_grp[mask_vec > 0], truth=active2_grp[mask==1])

pdf("act99_grp_FE.pdf")
image.nan(matrix(act1_grp, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main1, cex.main=2)
image.nan(matrix(act2_grp, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main2, cex.main=2)
dev.off()


##############################################################
##############################################################
### Joint Group Model
##############################################################
##############################################################

#15 MINUTES
t0 <- Sys.time()

betas.all <- matrix(0, nrow=N, ncol=2)
probs.all <- matrix(0, nrow=N, ncol=2) 

### COLLECT THETA POSTERIORS FROM SUBJECT MODELS

theta.sub <- NULL
mu.theta.tmp <- Q.theta <- 0
for(i in 1:M){
	t0 <- Sys.time()
  	load(paste0("result_subj",i,".Rdata"))
	# sum_m Q_m * mu_m
	mu.theta.tmp <- mu.theta.tmp + as.vector(Q.tmp%*%mu.tmp)
	# sum_m Q_m
	Q.theta <- Q.theta + Q.tmp
	theta.sub <- cbind(theta.sub, res.hyper$mode)
	rm(mu.tmp, Q.tmp)
}
#(sum_m Q_m)^(-1) * sum_m Q_m * mu_m
mu.theta <- solve(Q.theta, mu.theta.tmp)


### DRAW SAMPLES FROM q(theta|y)

nsamp <- 50
logwt <- rep(NA, nsamp)
theta.tmp <- mvrnorm(nsamp, mu.theta, solve(Q.theta))
for(i in 1:nsamp){ logwt[i] <- F.logwt(theta.tmp[i,], spde, mu.theta, Q.theta, M) }
#weights to apply to each posterior sample of theta
wt.tmp <- exp(logwt - max(logwt))
wt <- wt.tmp/(sum(wt.tmp))


### COMPUTE NECESSARY CROSS-PRODUCTS FOR EACH SUBJECT

Z1 <- sparseMatrix(ix, iy, x = rep(z1, N)) %*% Amat
Z2 <- sparseMatrix(ix, iy, x = rep(z2, N)) %*% Amat
Xmat <- cbind(Z1, Z2)

Xcros <- crossprod(Xmat) #same for all subjects
Xcros.all <- rep(list(Xcros), M)
Xycros.all <- vector("list", M)
for(i in 1:M){

	print(i)
	dat <- as.matrix(read.table(paste0('Dat_subj',i)))	
	Y <- as.matrix(t(dat[mask_vec==1,]))
	local_means <- matrix(colMeans(Y), nrow=nrow(Y), ncol=ncol(Y), byrow=TRUE) #center each voxel timecourse
	Y <- 100*(Y - local_means)/local_means
	y <- as.vector(Y) 
	Xycros.all[[i]] <- crossprod(Xmat, y)
}


### COMPUTE POSTERIOR QUANTITITES OF BETA FOR EACH THETA

## Create index vectors 
n.mesh <- mesh$n
ind_beta <- list()
for(k in 1:2){ ind_beta[[k]] <- 1:n.mesh + (k-1)*n.mesh }

#get posterior quantities of beta, conditional on a value of theta
no_cores <- min(detectCores() - 1, 25)
cl <- makeCluster(no_cores)
beta.post.samps <- parApply(cl, theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, Xcros = Xcros.all, Xycros = Xycros.all, thresholds=0, alpha=0.01, ind_beta=ind_beta)
stopCluster(cl)

#organize samples
mu.tot <- matrix(nrow=2*n.mesh, ncol=nsamp)
F.tot <- rep(list(matrix(nrow=n.mesh, ncol=nsamp)), 2) #for each task, a Vx50 matrix
for(itheta in 1:nsamp){
	mu.tot[,itheta] <- beta.post.samps[[itheta]]$mu
	for(k in 1:2){ F.tot[[k]][,itheta] <- beta.post.samps[[itheta]]$F[[1]][,k] }
}


### COMPUTE POSTERIOR QUANTITITES OF BETA, SUMMING OVER THETA WITH WEIGHTS

#posterior mean
beta.pop <- as.vector(mu.tot%*%wt)
for(k in 1:2){
	beta.pop.k <- beta.pop[ind_beta[[k]]]
	betas.all[,k] <- as.vector(Amat%*%beta.pop.k)
}

#posterior probabilities
for(k in 1:2){
	F.pop.k <- as.vector(F.tot[[k]]%*%wt)	
	probs.all[,k] <- as.vector(Amat%*%F.pop.k)
}

Sys.time() - t0

#SAVE RESULTS

save(betas.all, probs.all, file='results_grp_joint.Rdata')
#load(file='results_grp_joint.Rdata')


### VISUALIZE BETAS WITH MSE

beta1_grp <- mask_vec; beta1_grp[mask_vec > 0] <- betas.all[,1]; beta1_grp[mask_vec==0] <- NA
beta2_grp <- mask_vec; beta2_grp[mask_vec > 0] <- betas.all[,2]; beta2_grp[mask_vec==0] <- NA

#MSE
MSE1 <- mean((beta1_grp[beta1_grp_true>0] - beta1_grp_true[beta1_grp_true>0])^2, na.rm=TRUE)
MSE2 <- mean((beta2_grp[beta2_grp_true>0] - beta2_grp_true[beta2_grp_true>0])^2, na.rm=TRUE)
main1 <- paste0('Active MSE: ', format(round(MSE1, 4), nsmall=4))
main2 <- paste0('Active MSE: ', format(round(MSE2, 4), nsmall=4))

pdf("fit_beta_grp_joint.pdf")
image.nan(matrix(beta1_grp, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)), main=main1, cex.main=2)
image.nan(matrix(beta2_grp, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)), main=main2, cex.main=2)
dev.off()


### VISUALIZE AREAS OF ACTIVATION WITH FPR/FNR

act1_grp <- mask_vec; act1_grp[mask_vec > 0] <- (probs.all[,1] >= .99); act1_grp[mask_vec==0] <- NA
act2_grp <- mask_vec; act2_grp[mask_vec > 0] <- (probs.all[,2] >= .99); act2_grp[mask_vec==0] <- NA
main1 <- get_FPR_FNR_vec(vals=act1_grp[mask_vec > 0], truth=active1_grp[mask==1])
main2 <- get_FPR_FNR_vec(vals=act2_grp[mask_vec > 0], truth=active2_grp[mask==1])

pdf("act99_grp_joint.pdf")
image.nan(matrix(act1_grp, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main1, cex.main=2)
image.nan(matrix(act2_grp, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main2, cex.main=2)
dev.off()


##############################################################
##############################################################
### Two-Level Group Model
##############################################################
##############################################################

#6 MINUTES
t0 <- Sys.time()

betas.all <- matrix(0, nrow=N, ncol=2)
probs.all <- matrix(0, nrow=N, ncol=2) 

inla.setOption("num.threads", 4)


### COMBINE SUBJECT-LEVEL RESULTS

betas.tot <- matrix(nrow = N*M, ncol = 2)
sds.tot <- matrix(nrow = N*M, ncol = 2)

for(i in 1:M){
  	load(paste0("result_subj",i,".Rdata"))

  	betas <- cbind(res.beta$bbeta1$mean, res.beta$bbeta2$mean)
  	sds <- cbind(res.beta$bbeta1$sd, res.beta$bbeta2$sd)

	rows <- (1:N) + (i-1)*N
	betas.tot[rows,] <- as.vector(Amat%*%betas)
	sds.tot[rows,] <- as.vector(Amat%*%sds)
}

### RUN GROUP-LEVEL MODEL FOR EACH TASK

bbeta <- rep(spatial, M) 

for(k in 1:2){

	print(k)
    dat.inla <- list(y = betas.tot[,k], x = bbeta, z = 1:dim(sds.tot)[1], s = 1/(sds.tot[,k])^2)
    formula <- y ~ -1 + f(x, model = spde) + f(z, model = 'iid', hyper = list(theta=list(scale=s)))
    
    print(paste0('Fitting Model of Column ',k,' with INLA'))

	print(system.time(result <- inla(formula, data=dat.inla, control.compute=list(config=TRUE), verbose=TRUE)))

	mu.post <- result$summary.random$x$mean
	mu.post <- as.vector(Amat%*%mu.post) #resample to data size
	betas.all[,k] <- mu.post

	print(system.time(res.exc <- excursions.inla(result, name='x', u=0, type='>', method='QC')))
				
	#joint posterior probabilities
  	F.post <- as.vector(Amat%*%res.exc$F)
	probs.all[,k] <- F.post

}

Sys.time() - t0


### SAVE RESULTS

save(betas.all, probs.all, file='results_grp_2level.Rdata')
#load(file='results_grp_2level.Rdata')


### VISUALIZE BETAS WITH MSE

beta1_grp <- mask_vec; beta1_grp[mask_vec > 0] <- betas.all[,1]; beta1_grp[mask_vec==0] <- NA
beta2_grp <- mask_vec; beta2_grp[mask_vec > 0] <- betas.all[,2]; beta2_grp[mask_vec==0] <- NA

#MSE
MSE1 <- mean((beta1_grp[beta1_grp_true>0] - beta1_grp_true[beta1_grp_true>0])^2, na.rm=TRUE)
MSE2 <- mean((beta2_grp[beta2_grp_true>0] - beta2_grp_true[beta2_grp_true>0])^2, na.rm=TRUE)
main1 <- paste0('Active MSE: ', format(round(MSE1, 4), nsmall=4))
main2 <- paste0('Active MSE: ', format(round(MSE2, 4), nsmall=4))

pdf("fit_beta_grp_2level.pdf")
image.nan(matrix(beta1_grp, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)), main=main1, cex.main=2)
image.nan(matrix(beta2_grp, ny, nx), zlim=c(0, 1.5), col=r3, na.color='white', axis.args=list(at=seq(0,1.5,0.3)), main=main2, cex.main=2)
dev.off()


### VISUALIZE AREAS OF ACTIVATION WITH FPR/FNR

act1_grp <- mask_vec; act1_grp[mask_vec > 0] <- (probs.all[,1] >= .99); act1_grp[mask_vec==0] <- NA
act2_grp <- mask_vec; act2_grp[mask_vec > 0] <- (probs.all[,2] >= .99); act2_grp[mask_vec==0] <- NA
main1 <- get_FPR_FNR_vec(vals=act1_grp[mask_vec > 0], truth=active1_grp[mask==1])
main2 <- get_FPR_FNR_vec(vals=act2_grp[mask_vec > 0], truth=active2_grp[mask==1])

pdf("act99_grp_2level.pdf")
image.nan(matrix(act1_grp, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main1, cex.main=2)
image.nan(matrix(act2_grp, ny, nx), zlim=c(0, 1), col=c('lightgray','red'), na.color='white', main=main2, cex.main=2)
dev.off()










