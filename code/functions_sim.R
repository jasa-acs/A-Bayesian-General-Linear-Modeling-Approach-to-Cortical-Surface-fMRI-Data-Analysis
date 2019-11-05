################################################
### Functions to compute FPR and FNR
################################################

get_FPR_FNR <- function(vals, mask, truth){
	tmp <- truth[mask>0] - vals[mask>0]
	fpr <- round(100*length(which(tmp==-1))/length(tmp),2)
	fnr <- round(100*length(which(tmp==1))/length(tmp),2)
	print(paste0('FPR ',fpr,'%, FNR ',fnr,'%'))
}

get_FPR_FNR_vec <- function(vals, truth){
  #vals and truth are vectors, no need for a mask
  tmp <- truth - vals
  fpr <- round(100*length(which(tmp==-1))/length(tmp),2)
  fnr <- round(100*length(which(tmp==1))/length(tmp),2)
  print(paste0('FPR ',fpr,'%, FNR ',fnr,'%'))
}

#blue = false negative (truly activated, not detected)
#red = false positive (not activated, but detected)
img_FPR_FNR <- function(vals, mask, truth){
	img <- mask
	img[mask>=0] <- vals
	img_print <- mask
	img_print[mask==0] <- NA
	img_print[truth==1 & img==1] <- 2
	img_print[truth==1 & img==0] <- 3
	img_print[truth==0 & img==1] <- 4
	print(image(img_print, col=c('gray','black','blue','red')))
}

################################################
### Functions to visualize images
################################################

image.nan <- function(z, zlim, col, na.color='black', xlab='', axis.args=NULL, ...){

  zlim_orig <- zlim
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  #newz.below.outside <- zlim[1] - zstep # new z for values below zlim
  #newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA

  z[which(z<zlim[1])] <- zlim[1] # we affect newz.below.outside
  z[which(z>zlim[2])] <- zlim[2] # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na

  #zlim[1] <- zlim[1] - zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

  cols <- c(col, na.color) # we construct the new color range by including: na.color and na.outside

  nx <- ncol(z)
  ny <- nrow(z)
  obj.z <- list(x=1:ny, y=1:nx, z=z[ny:1,])

  axis.args <- c(as.list(axis.args), cex = 1.5)
  image.plot(obj.z, xlab=xlab, ylab='',  xlim=c(1,ny), ylim=c(1,nx), zlim=zlim, axis.args = axis.args, col=cols, xaxt="n", yaxt="n", ...)
  rect(47,0,60,55, col="white", border="white", xpd=TRUE)
  # Plot only the legend with original cutoffs
  image.plot(obj.z, xlab=xlab, ylab='',  xlim=c(1,ny), ylim=c(1,nx), zlim=zlim_orig, axis.args = axis.args, col=col, legend.only=TRUE, ...)
}

#for Bayesian GLM, where x and y might not be exactly 1:nx and 1:ny
image.nan2 <- function(zlist, zlim, col, na.color='black', axis.args=NULL, xlab='', ...){

  #zlist should contain three elements: x, y, and z (values)
  zlim_orig <- zlim

  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  #newz.below.outside <- zlim[1] - zstep # new z for values below zlim
  #newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA

  z <- zlist$z
  z[which(z<zlim[1])] <- zlim[1] # we affect newz.below.outside
  z[which(z>zlim[2])] <- zlim[2] # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  zlist$z <- z

  #zlim[1] <- zlim[1] - zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

  cols <- c(col, na.color) # we construct the new color range by including: na.color and na.outside

  axis.args <- c(as.list(axis.args), cex.axis = 1.5, cex.main=3)
  image.plot(zlist, xlab=xlab, ylab='',  xlim=c(1,ny), ylim=c(1,nx), zlim=zlim, axis.args = axis.args, col=cols, xaxt="n", yaxt="n", ...)
  rect(47,0,60,55, col="white", border="white", xpd=TRUE)
  # Plot only the legend with original cutoffs
  image.plot(zlist, xlab=xlab, ylab='',  xlim=c(1,ny), ylim=c(1,nx), zlim=zlim_orig, axis.args = axis.args, col=col, legend.only=TRUE, ...)
}
