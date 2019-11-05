maindir = '~/BayesianGLM/';
data_dir = '~/HCP/data/hpc/';
hcpdir = '~/HCP/data/hpc/';
wb_cmd = '~/workbench/bin_rh_linux64/wb_command';
toolbox_dir = '~/matlab_toolboxes/';

addpath(genpath(strcat(toolbox_dir,'cifti/'))) %ciftiopen, ciftisave, ciftisavereset
addpath(genpath(strcat(toolbox_dir,'spm8/')))
addpath(genpath(strcat(toolbox_dir,'gifti-1.6/')))
addpath(strcat(maindir,'code')) %CanonicalBasisSet.m

%this must be run after the other toolboxes above
addpath(genpath(strcat(toolbox_dir,'cifti-matlab/'))) %ft_read_cifti, ft_write_cifti


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET SUBJECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_dir)
subjects = zeros(0);
disks = zeros(0);
for idisk=1:5
    
   d = strcat('disk',num2str(idisk));
   s = dir(d);
   
   %get only directories
   s = struct2table(s);
   dirs = table2array(s(:,5));
   names = table2array(s(dirs==1,1));
   %remove . and ..
   names = names(3:end); 
   subjects = [subjects, names'];
   disks = [disks, repmat({d},[1 numel(names)])];
   
end

S = numel(subjects);
allsubjects = subjects;
alldisks = disks;
%samp = randsample(S, 20); %#1  50  65  81 143 216 249 280 324 336 405 409 417 463 467 468 490 491 494 496
samp = [1;  50;  65;  81; 143; 216; 249; 280; 324; 336; 405; 409; 417; 463; 467; 468; 490; 491; 494; 496];
mysubjects = subjects(samp); %"100307" "120111" "126628" "131924" "154936" "182840" "198451" "211417" "334635" "371843" "620434" "638049" "677968" "833148" "845458" "849971" "901442" "904044" "912447" "922854"
mydisks = disks(samp); %#"disk1" "disk1" "disk1" "disk1" "disk2" "disk3" "disk3" "disk3" "disk4" "disk4" "disk4" "disk4" "disk5" "disk5" "disk5" "disk5" "disk5" "disk5" "disk5" "disk5"
clear subjects disks
M = numel(mysubjects);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET COORDINATES OF GRAYORDINATES ON SPHERICAL SURFACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(strcat(maindir,'locations/'))
fname_giftiL = 'Sphere.6k.L.surf.gii'; %created in BayGLM0.sh
fname_giftiR = 'Sphere.6k.R.surf.gii'; %created in BayGLM0.sh

giftiL = gifti(fname_giftiL);
giftiR = gifti(fname_giftiR);

%get positions of each voxel
posR = giftiR.vertices;
posL = giftiL.vertices;

positions = [posL; posR];
V = size(positions,1); %11524 for 6K resampled data, 64984 for original 32K data
csvwrite('positions_sphere_resamp.csv', positions)
csvwrite('pos6K_right.csv', posR)
csvwrite('pos6K_left.csv', posL)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE DESIGN MATRIX (ALL SUBJECTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_dir)

for task = 1:2

    %For creating design matrix
    if(task==1) %Motor task
        EVs = {'cue','lf','lh','rf','rh','t'};
        len = 284;
        session = 'tfMRI_MOTOR_RL';
    else %Gambling task
        EVs = {'loss_event','neut_event','win_event'};
        len = 253;
        session = 'tfMRI_GAMBLING_RL';
    end
    K = numel(EVs); %Number of tasks
    TR = 0.72; %TR of data
    T = len*TR; % Total time of experiment in seconds 
    delta = 100;          % Downsample factor
    v0 = zeros(T*delta,1);% Define stick function
    h = CanonicalBasisSet(1/delta);
    ind = (TR*delta):(TR*delta):(T*delta); % Extract EVs in a function of TR

    for isubj=1:20

        isubj
        subj = mysubjects(isubj);
        disk = mysubjects(isubj);

        % Check whether Results directory exists
        dirs = ls(fullfile(alldisks{isubj},allsubjects{isubj},'MNINonLinear'));
        findResults = strfind(dirs, 'Results');
        if(numel(findResults)==0) continue; end

        % Check whether task session exists
        dirs = ls(fullfile(alldisks{isubj},allsubjects{isubj},'MNINonLinear/Results'));
        findsession = strfind(dirs, session);
        if(numel(findsession)==0) continue; end

        % CREATE DESIGN MATRIX
        idir = fullfile(alldisks{isubj},allsubjects{isubj},'MNINonLinear/Results',session,'EVs/');
        X = zeros(len,K);   
        for ev=1:K
              
            %first column is onset time, second column is duration
            itimes = dlmread(fullfile(idir,strcat(EVs{ev},'.txt')));
            ionsets = itimes(:,1);        
            idurations = itimes(:,2);
            
            %define stick functions
            iv = v0;
            for ii=1:numel(ionsets) 
                start_ii = round(ionsets(ii)*delta);
                end_ii = round((ionsets(ii) + idurations(ii))*delta);
                iv(start_ii:end_ii) = 1; 
            end

            is = conv(iv,h);
            is = is(ind);  % TR-specific EV
            X(:,ev) = is;  %add to design matrix

        end
        
        % Save design matrix
        if(task==1) fname = strcat('MOTOR/',subj,'_RL.csv'); end %original task
        if(task==2) fname = strcat('GAMBLING/',subj,'_RL.csv'); end %new task added in revisions
        writetable(array2table(X, 'VariableNames', EVs), fullfile('~/Bayesian2D/EVs/',fname{1}))


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% MOTION REGRESSORS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        if(task==1) fname = strcat('MOTOR/',subj,'_RL.txt'); end
        if(task==2) fname = strcat('GAMBLING/',subj,'_RL.txt'); end 
        idir = fullfile(alldisks{isubj},allsubjects{isubj},'MNINonLinear/Results',session);
        copyfile(fullfile(idir, 'Movement_Regressors.txt'), fullfile('~/Bayesian2D/motion/',fname{1}))    

    end %loop over subjects

end %loop over tasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FMRI TIMESERIES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FOR BAYESIAN GLM, FIRST RESAMPLE TO 6K 
% FOR CLASSICAL GLM, FIRST SMOOTH WITH 6MM FWHM 

V6K = 11524;
V32K = 64984;
fname_ts = '_Atlas.dtseries.nii';

sig = num2str(6/(2*sqrt(2*log(2)))); %convert 6mm FWHM to SD

for task = 1:2

    if(task==1) %MOTOR TASK
        session = 'tfMRI_MOTOR_RL'; 
        cd(strcat(maindir,'timeseries/MOTOR/'))
    else %GAMBLING TASK
        session = 'tfMRI_GAMBLING_RL'; 
        cd(strcat(maindir,'timeseries/GAMBLING/'))
    end 

    for isubj=1:20

        tic

        isubj
        subj = mysubjects(isubj);
        disk = mydisks(isubj);

        % Check whether Results directory exists
        dirs = ls(fullfile(hcpdir,mydisks{isubj},mysubjects{isubj},'MNINonLinear'));
        findResults = strfind(dirs, 'Results');
        if(numel(findResults)==0) continue; end

        % Check whether task session exists
        dirs = ls(fullfile(hcpdir,mydisks{isubj},mysubjects{isubj},'MNINonLinear/Results'));
        findsession = strfind(dirs, session);
        if(numel(findsession)==0) continue; end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %% PREPROCESSING (20 MIN)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        idir = fullfile(hcpdir,mydisks{isubj},mysubjects{isubj},'MNINonLinear/Results',session,'/');
        fname = strcat(idir,session, fname_ts);
        %names to write to in current directory
        fname_resamp = strcat(mysubjects{isubj},'.6K.dtseries.nii');
        fname_smooth = strcat(mysubjects{isubj},'.smooth.dtseries.nii');

        %Note: Helper files like Sphere.32K.R.surf.gii created previously in BayGLM0.sh

        % SMOOTHING (FOR CLASSICAL GLM)
        'Smoothing'
        rsurface = '../../locations/Sphere.32K.R.surf.gii'; 
        lsurface = '../../locations/Sphere.32K.L.surf.gii';
        cmd = [wb_cmd ' -cifti-smoothing ' fname ' ' sig ' ' sig ' COLUMN ' fname_smooth ' -left-surface ' lsurface ' -right-surface ' rsurface]; 
        system(cmd)

        % RESAMPLING (FOR BAYESIAN GLM)
        'Resampling'
        lsurface_6K = '../../locations/Sphere.6k.L.surf.gii';
        rsurface_6K = '../../locations/Sphere.6k.R.surf.gii';
        template = '../../locations/ts.6K.dtseries.nii';
        cmd = [wb_cmd ' -cifti-resample ' fname ' COLUMN ' template ' COLUMN BARYCENTRIC CUBIC ' fname_resamp ' -left-spheres '  lsurface ' ' lsurface_6K  ' -right-spheres ' rsurface ' ' rsurface_6K];
        system(cmd)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %% READ IN CIFTI FILES (3 SEC)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        'Reading'
        ts_nosmooth = ft_read_cifti(fname);
        dat_nosmooth = ts_nosmooth.dtseries(1:V32K,:);
        ts_smooth = ft_read_cifti(fname_smooth);
        dat_smooth = ts_smooth.smooth(1:V32K,:);
        ts_resamp = ft_read_cifti(fname_resamp);
        dat_resamp = ts_resamp.x6k(1:V6K,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %% WRITE TO CSV (40 SEC)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        'Writing'
        fname_csv_nosmooth = strcat(subj{1},'.nosmooth.csv');
        fname_csv_smooth = strcat(subj{1},'.smooth.csv');
        fname_csv_resamp = strcat(subj{1},'.6K.csv');
        csvwrite(fname_csv_nosmooth, dat_nosmooth);
        csvwrite(fname_csv_smooth, dat_smooth);
        csvwrite(fname_csv_resamp, dat_resamp);

    end %loop over subjects

end %loop over tasks


