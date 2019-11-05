% THIS CODE GENERATES SIMULATED FMRI DATA FOR SIMULATION 1 (SINGLE-SUBJECT) AND SIMULATION 2 (MULTI-SUBJECT)

addpath(genpath('~/matlab_toolboxes/spm12/'))
cd('code')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings

FOV = 192;
Nx = 46;
Ny = 55;
resx  = FOV/Nx;
resy  = FOV/Ny;
FWHM = [15 20 25];

V = spm_vol('grey.nii'); 
dat = spm_read_vols(V);

loc1 = [12,28];
loc2 = [36,28];
loc3 = [23,16];

intensity = [1,1,1];

% Paradigms (HRF of each activation)

Run1 = zeros(40,5);
Run1(1,:) = 1;
Run1 = reshape(Run1,200,1);

Run2 = zeros(40,5);
Run2(21,:) = 1;
Run2 = reshape(Run2,200,1);

h = spm_hrf(2);
h = h./max(h);

s1 = conv(Run1,h);
s2 = conv(Run2,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over Subjects

for(isubj = 1:20)

	isubj

	FWHM_i = FWHM + normrnd(0,1,[1,3]);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define kernel for areas of activations

	sig = FWHM_i./(2*sqrt(2*log(2)));
	sigx = sig./resx;
	sigy = sig./resy;

	len = -31:32;

	zx = normpdf(len,0,sigx(1));
	zy = normpdf(len,0,sigy(1));
	K1 = zx'*zy;
	K1 = K1./max(max(K1));

	zx = normpdf(len,0,sigx(2));
	zy = normpdf(len,0,sigy(2));
	K2 = zx'*zy;
	K2 = K2./max(max(K2));

	zx = normpdf(len,0,sigx(3));
	zy = normpdf(len,0,sigy(3));
	K3 = zx'*zy;
	K3 = K3./max(max(K3));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define areas of activations

	Q = (dat(:,:,45));
	Q = downsample(Q,2);
	Q = downsample(Q',2)';
	mask = 1*(Q>0);
	Q = Q + 0.1; %avoid problems calculating % local signal change with zero denominator

	loc1_i = loc1 + round(normrnd(0,.5,[1,2]));
	loc2_i = loc2 + round(normrnd(0,.5,[1,2]));
	loc3_i = loc3 + round(normrnd(0,.5,[1,2]));

	intensity_i = abs(normrnd(1,1,[1,3]));

	%Activation 1
	q1 = zeros(size(Q));
	q1(loc1_i(1),loc1_i(2)) = intensity_i(1); 
	q1 = conv2(q1,K1,'same');
	q1 = q1./max(max(q1));
	q1 = ((q1 > 0.2).*q1);

	%Activation 2
	q2 = zeros(size(Q));
	q2(loc2_i(1),loc2_i(2)) = intensity_i(2); 
	q2 = conv2(q2,K2,'same');
	q2 = q2./max(max(q2));
	q2 = ((q2 > 0.2).*q2);

	%Activation 3
	q3 = zeros(size(Q));
	q3(loc3_i(1),loc3_i(2)) = intensity_i(3); 
	q3 = conv2(q3,K3,'same');
	q3 = q3./max(max(q3));
	q3 = ((q3 > 0.2).*q3);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Generate timeseries data

	% Simulation 1: Single-subject, AR(1) residuals, with and without smoothing.  
	% 				Only generate data for first subject in loop.
	% Simulation 2: Multi-subject, independent residuals, no smoothing.

	DatNoResid = zeros(Nx,Ny,200);
	Dat = zeros(Nx,Ny,200); 

	%Simulation 1: generate AR residuals 
	if(isubj==1)  %simulation 1 (single subject) only
		DatAR = zeros(Nx,Ny,200); 
		E = zeros(Nx,Ny,200);
		phi = 0.3;
		E(:,:,1) = normrnd(0,2,Nx,Ny);
		for i=2:200
		    w=normrnd(0,2,Nx,Ny);
		    E(:,:,i) = phi*E(:,:,(i-1)) + w;
		end;
	end;

	%generate time series data for Simulations 1 & 2
	for i=1:200
	    w=normrnd(0,2,Nx,Ny);
	    DatNoResid(:,:,i) = 250*mask + 4*s1(i)*q1 + 2*s2(i)*q3 + 4*s1(i)*q2 + 2*s2(i)*q2;
	    Dat(:,:,i) = DatNoResid(:,:,i) + w; %simulation 2
	    if(isubj==1) DatAR(:,:,i) = DatNoResid(:,:,i) + E(:,:,i); end; %simulation 1
	end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Simulation 1: Smooth data with AR(1) errors for classical GLM

	if(isubj==1) %simulation 1 (single subject) only
		len = -32:32;
		FWHM0 = 6;
		sig0 = FWHM0./(2*sqrt(2*log(2)));
		sig0x = sig0./resx;
		sig0y = sig0./resy;

		zx0 = normpdf(len,0,sig0x(1));
		zy0 = normpdf(len,0,sig0y(1));
		K0 = zx0'*zy0;
		%K0 = K0./max(max(K0)); %this actually changes the scale of the data.  AUC of the kernel should equal 1.

		DatAR_sm = DatAR;
		for i=1:200,
			DatAR_sm(:,:,i) = conv2(DatAR_sm(:,:,i),K0,'same');
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Reshape & save data

	save(['Q1_subj',num2str(isubj)],'q1','-ascii')
	save(['Q2_subj',num2str(isubj)],'q2','-ascii')
	save(['Q3_subj',num2str(isubj)],'q3','-ascii')

	DatNoResid = reshape(DatNoResid,Nx*Ny,200);
	Dat = reshape(Dat,Nx*Ny,200);
	if(isubj==1) %simulation 1 (single subject) only
		DatAR = reshape(DatAR,Nx*Ny,200);
		DatAR_sm = reshape(DatAR_sm,Nx*Ny,200);
	end

	save(['Dat_noResid_subj',num2str(isubj)],'DatNoResid','-ascii')
	save(['Dat_subj',num2str(isubj)],'Dat','-ascii')
	if(isubj==1) %simulation 1 (single subject) only
		save(['DatAR_subj',num2str(isubj)],'DatAR','-ascii')
		save(['DatAR_sm_subj',num2str(isubj)],'DatAR','-ascii')
	end

end %loop over subjects

save('Mask2','mask','-ascii')
