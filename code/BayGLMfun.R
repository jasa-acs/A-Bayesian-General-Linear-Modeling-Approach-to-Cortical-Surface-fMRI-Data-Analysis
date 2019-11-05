#######################################################################################
# FUNCTIONS FOR PREWHITENING
#######################################################################################

getInvCovAR <- function(ar, ntime){

	#ar = vector of p AR parameters
	#ntime = number of time points in timeseries

	# Compute inverse covariance matrix for AR process (up to a constant scaling factor)

	Inv0 <- diag(ntime)
	incr0 <- matrix(0, nrow=ntime, ncol=ntime)
	offs <- row(Inv0) - col(Inv0) #identifies each off-diagonal

	p <- length(ar)
	for(k in 1:p){
		incr <- incr0 #matrix of zeros
		incr[offs==k] <- -1*ar[k]
		Inv0 <- Inv0 + incr
	}
	Inv <- Inv0 %*% t(Inv0) 
	return(Inv)

}

getSqrtInv <- function(Inv){

	# compute square root of a symmetric matrix (such as inverse covariance mat) using SVD
	# Cov = UD^2U', then Inv = UD^(-2)U', and Inv_sqrt = UD^(-1)U' 
	ei <- eigen(Inv) 
	d2inv <- ei$values #diag(D^(-2))
	if(sum(d2inv<0)>0) print('negative eigenvalues')
	sqrtInv <- ei$vectors %*% diag(sqrt(d2inv)) %*% t(ei$vectors)
	return(sqrtInv)

}


#######################################################################################
# FUNCTIONS FOR READING/WRITING CIFTI FILES 
#######################################################################################


#make sure that gifti library is loaded
readCIFTI <- function(filename, wbcommand = '~/workbench/bin_rh_linux64/wb_command', scratchloc = '~/', tmpname = 'tmp', verbose=FALSE){

	#wbcommand = path to wb_command
	#scratchloc = location where temporary files should be written
	#tmpname = name of temporary files (should be unique if many jobs running in parallel)
	#verbose = if TRUE, print out progress

	require(gifti)
	tmpfile <- paste0(scratchloc,tmpname,'.gii')
	cmd <- paste0(wbcommand, ' -cifti-convert -to-gifti-ext ', filename, ' ', tmpfile)
	if(verbose) print('converting to GIFTI')
	system(cmd)
	if(verbose) print('reading GIFTI')
	cifti <- readgii(tmpfile)
	if(verbose) print('removing temporary files')
	file.remove(tmpfile)
	file.remove(paste0(tmpfile,'.data'))
	return(cifti$data$normal)
}

#NOTE: THESE FUNCTIONS CALL MATLAB FROM R. ALTERNATIVE IS TO USE MATLAB DIRECTLY OR USE THE READ FUNCTION ABOVE. ALSO CONSIDER THE NEW R CIFTI PACKAGE FOR READING.

readCIFTI_ts <- function(name, hemisphere='both', 
						toolboxloc = '~/matlab_toolboxes/cifti-matlab/',
						scratchloc = '~/',
						tsname = 'dtseries',
						V = 32492){

	# vector of name(s) of CIFTI dtseries file(s) to be read in
	# hemisphere = 'both', 'left', 'right', or 'all' (includes subcortical voxels)
	# toolboxloc = location and name of folder containing cifti-matlab or fieldtrip toolbox
	# scratchloc = location where CSV files will be written then deleted
	# tsname = name of field corresponding to time series
	# V = number of vertices in each hemisphere

	line0 <- paste0("cd '", getwd(), "'")
	line1 <- paste0("addpath '", toolboxloc, "'")
	line2 <- paste0("cifti = ft_read_cifti('",name,"');")  #reads in structure with cell 'indexmax'

	if(hemisphere=='all')  line3 <- paste0("ts = cifti.",tsname,";")
	if(hemisphere=='both') line3 <- paste0("ts = cifti.",tsname,"(1:",V*2,",:);")
	if(hemisphere=='left') line3 <- paste0("ts = cifti.",tsname,"(1:",V,",:);")
	if(hemisphere=='right') line3 <- paste0("ts = cifti.",tsname,"(",V+1,":",V*2,",:);")

	# Write timeseries to CSV
	fname <- file.path(scratchloc, 'tmp.csv')
	line4 <- paste0("csvwrite('",fname,"', ts);")

	# Run MATLAB script
	matlab_lines <- c(line0, line1, line2, line3, line4)
	fname_script <- file.path(scratchloc, 'myscript.m')
	writeLines(matlab_lines, con=fname_script)
	system(paste0("matlab -nodisplay -r \"run('",fname_script,"'); exit\""))
	file.remove(fname_script)

	# Read in written CSV file
	result <- read.csv(fname, header=F)	
	file.remove(fname)
	return(result)

}

writeCIFTIs <- function(table, names, hemisphere='both', 
						template = '~/HCP/data/groupICA/groupICA_3T_Q1-Q6related468_MSMsulc_d25.ica/melodic_IC_ftb.dlabel.nii', 
						data_field = 'indexmax',
						toolboxloc = '~/matlab_toolboxes/cifti-matlab/',
						scratchloc = '~/'){

	# table = VxQ matrix, each column will be written as a CIFTI file
	# table can be a vector if only a single CIFTI is needed
	# V = number of voxels in both or one hemisphere
	# Q = number of CIFTI files to be written
	# names = vector length Q, containing file names to be written (no extension)
	# hemisphere = 'both', 'left', 'right', or 'all' (includes subcortical voxels)
	# template = name of existing CIFTI file that can be used as a template
	# data_field = name of field in template that contains main data (e.g. parcellation labels, timeseries matrix)
	# toolboxloc = location and name of folder containing cifti-matlab or fieldtrip toolbox
	# scratchloc = location where CSV files will be written then deleted

	# this function writes a MATLAB script and executes it
	# must have MATLAB installed
	# must have the cifti-matlab or fieldtrip toolbox installed
	# currently can only write CIFTIs with a vector of values (e.g. not a timeseries or a set of label maps)

	# if hemisphere='both', first 32492 rows are left hemisphere and last 32492 rows are right hemisphere

	if(is.vector(table)) table <- matrix(table, nrow=length(table), ncol=1)
	if(!is.matrix(table)) table <- as.matrix(table)
	Q <- ncol(table)
	V <- nrow(table)

	# if(hemisphere=='all')  if(V != 96854) stop('Number of rows must equal 96854 with hemisphere="all"')
	# if(hemisphere=='both') if(V != 64984) stop('Number of rows must equal 64984 with hemisphere="both"')
	# if(hemisphere=='left') if(V != 32492) stop('Number of rows must equal 32492 with hemisphere="left"')
	# if(hemisphere=='right') if(V != 32492) stop('Number of rows must equal 32492 with hemisphere="right"')

	# Write table to CSV in scratch location
	tmp <- table
	tmp[is.na(tmp)] <- 0 #set NAs to zero
	fname <- file.path(scratchloc, 'tmp.csv')
	write.table(tmp, file=fname, row.names=FALSE, col.names=FALSE, sep=',')

	# Write MATLAB script
	line1 <- paste0("addpath '", toolboxloc, "'")
	line2 <- paste0("cifti = ft_read_cifti('",template,"');")  #reads in structure with cell 'indexmax'
	line3 <- paste0("indsL = (cifti.brainstructure == 1);", #L cortex only
					"indsR = (cifti.brainstructure == 2);", #R cortex only
					"indsLR = (cifti.brainstructure <= 2);", #L and R cortices only
					"nL = sum(indsL);",
					"nR = sum(indsR);",
					"nTOT = length(indsL);",
					paste0("inds_nan = isnan(cifti.",data_field,");")) #location of NaNs


	if(hemisphere=='all')  line4a <- paste0("if nTOT ~= ",V," error = 1; else error = 0; end")
	if(hemisphere=='both') line4a <- paste0("if (nL + nR) ~= ",V," error = 1; else error = 0; end")
	if(hemisphere=='left') line4a <- paste0("if nL ~= ",V," error = 1; else error = 0; end")
	if(hemisphere=='right') line4a <- paste0("if nR ~= ",V," error = 1; else error = 0; end")

	line4b <- paste0("if error == 1 disp('Incorrect number of rows for hemisphere = ",hemisphere," based on template brainstructure labels'); exit; end")

	line5 <- paste0("cd '", getwd(), "'")
	line6 <- paste0("fname = '",fname,"';",
					"coefs = csvread(fname);")

	line7 <- paste0("names = {'",paste(names, collapse="';'"),"'};")

	line8 <- paste0("for k=1:",Q)
	line9 <- "k"
	line10 <- "vals = coefs(:,k);"
	line11 <- paste0("cifti.",data_field," = cifti.",data_field,"*0;")

	if(hemisphere=='all') line12 <- paste0("cifti.",data_field," = vals;")
	if(hemisphere=='left') line12 <- paste0("cifti.",data_field,"(indsL) = vals;")
	if(hemisphere=='right') line12 <- paste0("cifti.",data_field,"(indsR) = vals;")
	if(hemisphere=='both') line12 <- paste0("cifti.",data_field,"(indsLR) = vals;")

    line13 <- paste0("cifti.",data_field,"(inds_nan) = 0/0;") #put in NaNs
	line14 <- paste0("ft_write_cifti(names{k}, cifti, 'parameter', '",data_field,"');")
	line15 <- "end"

	matlab_lines <- c(line1, line2, line3, line4a, line4b, line5, line6, line7, line8,
					  line9, line10, line11, line12, line13, line14, line15)
	fname_script <- file.path(scratchloc, 'myscript.m')
	writeLines(matlab_lines, con=fname_script)
	system(paste0("matlab -nodisplay -r \"run('",fname_script,"'); exit\""))
	file.remove(fname_script)

	file.remove(fname)
}



#######################################################################################
# FUNCTIONS FOR JOINT GROUP APPROACH
#######################################################################################

#for one posterior sample of theta, compute posterior mean and excursion probabilities for each latent field
beta.posterior.thetasamp <- function(theta, spde, Xcros, Xycros, thresholds=c(0,0.5,1), alpha=0.01, ind_beta){

	#theta - one sample of theta 

	require(excursions)
	require(INLA)
	INLA:::inla.dynload.workaround() 

	print('Constructing joint precision')
	prec.error <- exp(theta[1])
	K <- length(theta[-1])/2
	M <- length(Xcros)

	#construct prior precision matrix for beta, Q_theta,
	#for given sampled values of theta
	theta.beta <- list()
	Q.beta <- list()
	for(k in 1:K) { 
		theta.beta[[k]] <- theta[(2:3) + 2*(k-1)] #2:3, 4:5, ...
		Q.beta[[k]] <- inla.spde2.precision(spde, theta = theta.beta[[k]])
	}
	Q <- bdiag(Q.beta) 

	beta.samp.pop <- 0
	beta.mean.pop <- 0
	#~25 seconds per subject
	print('Looping over subjects')
	for(mm in 1:M){
		Xcros.mm <- Xcros[[mm]]
		Xycros.mm <- Xycros[[mm]]
		Q.m <- prec.error*Xcros.mm + Q 
		mu.m <- inla.qsolve(Q.m, prec.error*Xycros.mm) #20 sec
		beta.mean.pop <- beta.mean.pop + mu.m/M

		#draw samples from pi(beta_m|theta,y)

	  	#this only works when the pop-level quantity of interest is the average activation for each task
	  	#can use the linear combination matrix A to make more general
	  	#n=10: 20 sec, n=100: 90 sec
	    beta.samp.m <- inla.qsample(n = 100, Q = Q.m, mu = mu.m) #NKx100
	    beta.samp.pop <- beta.samp.pop + beta.samp.m/M #NKx100
	}
	mu.theta <- matrix(beta.mean.pop, ncol=1)

	#3.5-7 seconds per activation threshold
	print('Looping over activation thresholds')
	n.mesh <- spde$n.spde
	U <- length(thresholds)
	F.theta <- vector('list', U)
  	for(u in 1:U){
  		F.theta[[u]] <- matrix(nrow=n.mesh, ncol=K)
  		thr <- thresholds[u]
  		for(k in 1:K){
  			res_beta.theta.k <- excursions.mc(beta.samp.pop, u = thr, ind = ind_beta[[k]], type = '>', alpha = (1-alpha), verbose = FALSE)
  			F.theta[[u]][,k] <- res_beta.theta.k$F[ind_beta[[k]]]
  		}
  		F.theta[[u]][is.na(F.theta[[u]])] <- 0
  	}

  	result <- list(mu.theta, F.theta)
  	names(result) <- c('mu','F')
  	return(result)
}

F.logwt <- function(theta, spde, mu.theta, Q.theta, M){
  #theta - vector of hyperparameter values at which to compute the posterior log density
  #spde - spde object, determines prior precision matrix
  #mu.theta - posterior mean from combined subject-level models
  #Q.theta - posterior precision matrix from combined subject-level models
  #M - number of subjects
  a <- 1; b <- 5e-5
  n.spde <- (length(theta) - 1)/2
  mu.tmp <- spde$f$hyper$theta1$param[1:2]
  mu <- rep(mu.tmp, n.spde)
  Q.tmp <- matrix(spde$f$hyper$theta1$param[-(1:2)], 2, 2, byrow = TRUE)
  Q <- kronecker(diag(1, n.spde, n.spde), Q.tmp)
  
  ## Prior density
  pr.delta <- dgamma(exp(theta[1]), a, b, log = TRUE) #log prior density on residual precision
  pr.tk <- as.vector(-t(theta[-1] - mu)%*%Q%*%(theta[-1] - mu))/2 + log(det(Q))/2 - dim(Q)[1]*log(2*pi)/2 #joint log prior density on 2K spde parameters
  pr.theta <- pr.delta + pr.tk
  
  (1-M)*pr.theta 
}


compute.cross <- function(yZ, Amat, sqrtInv, rows.rm){

	y <- yZ[[1]]
	Z <- yZ[[2]]

	ntime <- dim(Z)[1]
	nvox <- length(y)/ntime
	ix <- 1:(ntime*nvox)
	iy <- rep(1:nvox, each = ntime)
	K <- ncol(Z)

	#construct design matrix Xmat
	for(k in 1:K){
		Z_col <- sparseMatrix(ix, iy, x=rep(Z[,k], nvox)) %*% Amat
		if(k==1) Xmat <- Z_col else Xmat <- cBind(Xmat, Z_col)
	}

	#prewhiten Xmat and y
    Xmat <- sqrtInv %*% Xmat #~10 seconds
	Xmat <- Xmat[-rows.rm, ]
	y <- as.vector(sqrtInv %*% y) #<1 seconds
	y <- y[-rows.rm]

	#compute and save cross products
    Xcros <- crossprod(Xmat)
    Xycros <- crossprod(Xmat, y)

    result <- c(list(Xcros), list(Xycros))
    names(result) <- c('Xcros','Xycros')
    return(result)
}




