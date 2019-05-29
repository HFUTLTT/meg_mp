%%%%%%%%%%%%%%% Max Voxel/Sensor - MP %%%%%%%%%%%%%%%
% This program computes time-frequency spectrum using MP algorithm for
% preprocessed and epoched MEG data:
% 1. Voxel with maximal gamma increase
% 2. Sensor with maximal gamma increase
% This implementation supports 2 types of MP dictionaries:
% 1. Dyadic.
% 2. Stochastic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment Setup

% Close and clear everything
clc;
clear all;
close all;

% Path Setup
restoredefaultpath;

% Code Paths

fieldtripfolder = 'fieldtrip/'; % Fieldtrip path 
%fieldtripfolder =  '/Applications/fieldtrip-20170515/';

addpath(fieldtripfolder);
ft_defaults; % Fieldtrip path setup
scriptPath = 'SpectralAnalysis';
addpath(scriptPath);

% File Paths
DataSource  = 'MaxSenVox'; % MEG recording source 'Max Sensor' and 'Max Voxel'
% Input File Paths
inputPath = strcat('SubjectData/', DataSource, '/'); % Source data
intermediatePathPrefix = 'MPIntermediateData_';
outputPathPrefix = 'MPOutputData_';

% Paths - Elena's Script
DataPath = 'FT_beamformer/';
SavePath = 'Test/';
megfolder =  '/max_vox_max_sens/';

%% MP Analysis Parameters
% Specified based on use-case by the user

start_time = clock; % [year month day hour minute seconds]

% Modes
MPDic = 'dyadic'; % Matching pursuit dicitonary: 'dyadic' or 'stochastic'
% Subtraction: Both lines below should be modified in an identical manner to relfect desired mode(s)
subtractAvg = [0]; % Subtract evoked response? 0:no, 1:yes e.g [0;1]
modeText = {'NoSubEvoked'}; % e.g. {'NoSubEvoked'; 'SubEvoked'}

% Options
runDecomposition = 1; % 0: Load previously decomposed data 1: Run decomposition from scratch
ClearPreviousOuput = 1; % Clear output (including intermediate data) from previous execution: 0: No, 1:Yes
DeleteIntermediateData = 0; % Delete intermediate MP data: 0:No, 1:Yes
MPmaxIterations = 500; % MP max number of iterations

% Data

% Subjects
% NT: 
SUBJ1 = [0101; 0102; 0103; 0104; 0105; 0136; 0137; 0138; 0140;0158; 0162; 0163; 0178; 0179; 0255; 0257; 0348; 0378; 0384];
% ASD: 
SUBJ2 = [0107; 0139; 0141; 0159; 0160; 0161; 0164; 0253; 0254; 0256; 0273; 0274; 0346; 0347; 0351; 0358; 0380; 0381;0382; 0383]; 

EVENTS = [2,4,8]; % List of event (i.e. condition) number as in data file names e.g. 0076_V1, 0076_V2, etc...

% Time ranges of interest (used for power change calculation and plots)
baseline_range = [-0.7 -0.1];
post_stimulus_range = [0.4 1];

% Special parameters for "Stochastic" dictionary only
norm_factor = 10^12; % Normalization factor for stochastic dictionary (magnetic field magnitude)
adaptiveDictionaryParam = 0.9; % ADP for subset selection (size = (1-ADP)*100%)
dictionarySize = 2500000; % Dictionary size

%% Automated Initialization

% Output folder setup
intermediateOuputPath = strcat(intermediatePathPrefix, MPDic, '_',DataSource, '/');
outputPath = strcat(outputPathPrefix, MPDic, '_',DataSource, '/');
MP_SetupOuputFolder(ClearPreviousOuput, outputPath); % MP intermediate data
MP_SetupOuputFolder(ClearPreviousOuput, intermediateOuputPath); % MP ouput data

% Set variable values
con_num = size(EVENTS,2); % Get the total number of events

%% Time Frequency Analysis - Elena's Script + MP

gridres = 6; % 6 mm grid
lambda='5%';
tapsmofrqnew = 5;
gammarange = [35 110];

filt = 'yes'; %'no'; %
hpfreq=gammarange(1)-5;
lpfreq=gammarange(2)+5;
% Select Subject Group
grp = 1; % 1:NT, Otherwise:ASD

% we defined subjects here:
SUBJ = [      '0076';   '0101';   '0102';   '0103';   '0104';   '0105';         '0106';   '0107';   
                        '0136';   '0137';   '0138';   '0139';   '0140';   '0141';      
                        '0158';   '0159';   '0160';   '0161';   '0162';   '0163';   '0164';   
                        '0178';   '0179';    
                        '0253';   '0254';   '0255';   '0256';   '0257';             '0259';                   
                        '0273';   '0274';   '0275';   '0276' ;  '0277';    
                                  '0346';   '0347';   '0348';             '0350';   '0351'; '0357'; '0358';
                        '0378';   '0380';   '0381';   '0382';   '0383';   '0384']; % 
                                                                                   
% Segmentations with different brain thresholds:
list05 = ['0380';  '0382';  '0383';  '0384';  '0385'; '0076';  '0102'; '0107'; '0139'; '0140'; '0141'; '0160'; '0161'; '0163'; '0178'; '0254'; '0259'; '0273'; '0274'; '0275';  '0277'; '0346'; '0347'; '0348'; '0350'; '0358' ]; %
list09 = [ '0103'; '0104'; '0106'; '0136'; '0137'; '0138'; '0158'; '0159';'0164'; '0179'; '0255'; '0257'; '0351'; '0357'; '0378'];
list098 = [ '0105'; '0162'; '0253'; '0256';'0381'] ;
list03 = ['0276'; '0101'];
list01 = ['0276']; %;

% load atlas, template grid, common for all subjects
subj = '0104';
%  load atlas
atlas = ft_read_atlas( strcat (fieldtripfolder, '/template/atlas/aal/ROI_MNI_V4.nii') ); 
atlas = ft_convert_units(atlas,'cm');% assure that atlas and template_grid are expressed in the %same units

%  template MRI
templatefile = strcat (fieldtripfolder, '/external/spm8/templates/T1.nii'); 
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'mni';

% load template grid. Use the same template gris, that was used for the
% warped grid construction!
load ( strcat(scriptPath, '/', DataPath,'0076/mri_linwarp_6mm_brthr0.5/0076_template_grid.mat') );

cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
cfg.inputcoord = 'mni';
mask          = ft_volumelookup(cfg, template_grid);

%chose only the 'brain tissue' 
template_grid.inside = false(template_grid.dim);
template_grid.inside(mask==1) = true;

% By doing this you get the label for every grid  :
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
atlas_template_grid = ft_sourceinterpolate(cfg, atlas, template_grid);
atlas_template_grid.tissue(isnan(atlas_template_grid.tissue)) =0;

ids      = find(atlas_template_grid.tissue);          % all interpolate regions
id       = atlas_template_grid.tissue(ids); %  all interpolate regions index
ROI      = atlas.tissuelabel(id);
occid1   = find(strncmpi(ROI,'Occipital',9));  %  indice
occid2   = find(strncmpi(ROI,'Calcarine',9));  %  indice
occid3   = find(strncmpi(ROI,'Cuneus',6));  %  indice
occid4   = find(strncmpi(ROI,'Lingual',7));  %  indice
% occid5   = find(strncmpi(ROI,'Precuneus', 9));  
occid    = sort([occid1, occid2, occid3, occid4]);        
OCC      = ROI(occid);  % label

mask_occ = zeros(size(template_grid.pos,1), 1); mask_occ(ids(occid)) = 1;

%% across subjects
foi_sensor   = [1 : 1:  gammarange(2)];
foi_source   = [25 : 1:  gammarange(2)]; % for source localisation with beamformers we  filter raw signal

MPData_All = {}; % List to save all MP ouput data
FsAll = {};

for s=1:size (SUBJ,1) % Loop over subjects 
    
    close all
    
    subj = SUBJ (s,:);
    
    VoxelSourceInfoFile = [scriptPath, '/', DataPath, subj, '/meg6mm_linwarp_tapsmofrq5_lcmv_norm/', subj, '_Lcmv_source_spectra.mat'];
    VoxelSensorInfoFile = [scriptPath, '/',DataPath, subj, '/sensors/', subj, '_sensors_res.mat'];

    %mkdir (strcat(SavePath, subj));
    %mkdir (strcat(SavePath, subj, megfolder));
    %savemegto = strcat(SavePath, subj, megfolder);
        
    if strmatch(subj,list05) 
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.5/');
    elseif strmatch(subj,list09) 
        mrifolder =strcat( '/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.9/');
    elseif strmatch(subj,list098) 
        mrifolder =strcat( '/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.98/');
    elseif strmatch(subj,list03) 
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.3/');
    elseif strmatch(subj,list01) 
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.1/');
    else
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.5/');
    end
    
    %% start PPTX report
    %exportToPPTX('new');
    % PPTXname  = strcat(savemegto, subj, '_lcmv_report');

    %% subject specific Loads   
    load ( strcat(scriptPath, '/', DataPath , subj, mrifolder, subj, '_grid_MNI_lf.mat') ); % load leadfield / % source model 
    load ( strcat(scriptPath, '/', DataPath , subj, mrifolder, subj, '_individ_hdm_vol.mat' )); % head model: individ_hdm_vol
    
    % load (strcat(DataPath , subj, megfolder, subj, '_DISCbroad_source_maxT.mat'));
    load (strcat(scriptPath, '/', DataPath , subj, '/meg6mm_linwarp_tapsmofrq5_individWF_maxchange/', subj, '_DISCbroad_source_maxT.mat')); % DICS results  
    
    load ( strcat(scriptPath, '/', DataPath , subj, mrifolder, subj, '_mri_orig_realigned.mat' )); % load individual MRI

    %% Load data epochs 
    %ep_fiff_file = strcat(DataPath, subj, '/epochs/', subj, '-noerror-lagcorrected-epo.fif');
    ep_fiff_file = strcat(inputPath,subj,'/epochs/', subj, '-noerror-lagcorrected-epo.fif');
    hdr = ft_read_header(ep_fiff_file);
    
    % find good epochs
    C = strsplit(hdr.orig.epochs.drop_log,'], ');
    find1=strcmp(C, '[[');
    find2=strcmp(C, '[');
    find3=strcmp(C, '[]]');
    ind = find(find1+find2+find3);
    % read events (saved by '/Users/mtw/MEG/Scripts/Karolinska/PY/Sensors/resample_raw.py' )
    % events = load (strcat (DataPath, subj, '/ICA_nonotch_crop/', subj, '_events.mat'));
    events = load ( strcat (inputPath,subj, '/ICA_nonotch_crop/', subj, '_events.mat') );
    events = events.events( ind, :);
    % Load raw data
    %fiff_file = strcat(DataPath, subj, '/ICA_nonotch_crop/', subj, '_rings_ICA_raw.fif');
    fiff_file = strcat(inputPath,subj,'/ICA_nonotch_crop/', subj, '_rings_ICA_raw.fif');
    hdrraw = ft_read_header(fiff_file);
    pre = -1.0* hdr.Fs ;
    post = 1.2* hdr.Fs ;
    trl=[];
    for i=1:size (events,1)
         trl(i, 1)=(events(i,1)+pre) ; 
         trl(i, 2)=(events(i,1)+post) ; 
         trl(i, 3)= -1.0*hdr.Fs ; % offset
         trl(i, 4) = events(i,3); % stimulus_value;
    end
    % extract data and epochs from the raw
    cfg = [];
    cfg.trl=trl;
    cfg.channel     = 'meg';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50 100]; % for power noise rejection
    cfg.demean = 'yes';
    cfg.dataset = fiff_file;
    cfg    = ft_definetrial(cfg);
    epochs_sensors = ft_preprocessing(cfg);

    if strcmp(filt,'yes')  % chech if you want to filter the data before localization
       cfg.lpfilter  = 'yes';
       cfg.lpfreq    = lpfreq;
       cfg.hpfilter  = 'yes';
       cfg.hpfreq    = hpfreq;
    end
    epochs_source = ft_preprocessing(cfg);

    %% Plot average
    cfg = [];
    cfg.channel=epochs_source.label;
    avg = ft_timelockanalysis(cfg,epochs_source);

    hh=figure;
    ttt = find(avg.time>-0.1 & avg.time<0.5);
    plot (avg.time(ttt), avg.avg(:, ttt));
    title ('Sensor averages');

    %exportToPPTX('addslide'); % slide with distributions
    %exportToPPTX('addpicture', hh, 'Position', [0.5,0.5,9,6]);

    %%  select epochs according to events
    % EVENTS = [2,4,8];
    EV = {};
    for v=1:size(EVENTS,2)
        EV{end+1} = find( events(:,3)== EVENTS(v) );
    end

    %% calculate covariance matrix for source localization
    cfg = [];
    cfg.channel           = epochs_source.label;
    cfg.covariance        = 'yes';
    %cfg.covariancewindow  =  [-0.9 -0.1; 0.3 1.1];%'all'; %[-0.9 -0.1]; %
    cfg.covariancewindow  =  'all'; % use the whole window
    cfg.vartrllength      = 2;
    total_avg = ft_timelockanalysis(cfg, epochs_source);   
    covariance.wind = cfg.covariancewindow;
    covariance.epochs = 'all conditions';
    
    %% perform source analysis, use only points inside the brain
    sourcemodel = grid_MNI_lf;
    sourcemodel.inside = reshape(template_grid.inside, [1, template_grid.dim(1)*template_grid.dim(2)*template_grid.dim(3)]);
    
    cfg=[];
    cfg.method='lcmv';
    cfg.grid = sourcemodel; %template_grid;%grid_MNI_lf;
    cfg.headmodel=individ_hdm_vol;
    cfg.lcmv.keepfilter='yes'; % keep filters in the output, which are later multiplied with the data
    cfg.lcmv.fixedori='yes'; % id 'yes' consider only the dominant orientation
    cfg.lcmv.lambda=lambda;
    %cfg.lcmv.projectmom = 'yes';
    cfg.reducerank = 2;
    cfg.normalize = 'yes';
    source_total=ft_sourceanalysis(cfg, total_avg);
    source_total.pos = template_grid.pos;
 
    %%  load information about sensor and source with maximal gamma response.
    sensorinfo = load(VoxelSensorInfoFile, 'ChName');
    sourceinfo = load(VoxelSourceInfoFile, 'MAX_VOXEL_OCC_NUMBER', 'NVOXELS_ID');
    %%

    for con =1:size(EVENTS,2) % Conditions       
        for mode=1:size(subtractAvg,1) % Mode
            
            %% select data according to condition (slow-1, medium-2, fast-3)
            cfg=[];
            cfg.trials = EV{con}; 
            [data_source] = ft_selectdata(cfg, epochs_source);
            [data_sensors] = ft_selectdata(cfg, epochs_sensors);

            %% Multiply filters with the data and organize into FieldTrip sensable data structure
            spatialfilter=cat(1,source_total.avg.filter{:});  % all voxels inside

            virtsens=[]; 
            for i=1:length(data_source.trial)
                virtsens.trial{i}=spatialfilter*data_source.trial{i};
            end
            virtsens.time=data_source.time;
            virtsens.fsample=data_source.fsample;
            virtsens.label= cellstr(string(find(sourcemodel.inside))); %(occid(isnotempt)))'
            virtsens.pos = template_grid.pos(find(sourcemodel.inside),:);   
            virtsens.grad =data_source.grad;
            % Now virtsens contains time cources at the brain sources ('voxels')


            %%  We already found brain voxel with maximal gamma response and the  24 voxels closest to the max-voxel
            % we  loaded this information from: '/Users/mtw/MEG/Karolinska/FT_beamformer/0101/meg6mm_linwarp_tapsmofrq5_individWF_maxchange/0101_Lcmv_source_spectra.mat'
            % figure, plot(virtsens.time{1}, virtsens.trial{33}(info.MAX_VOXEL_OCC_NUMBER{con},:))

            %% MP analysis on max virtual channel timecourse
            
            % Max voxel info
            maxVoxIdx = sourceinfo.MAX_VOXEL_OCC_NUMBER{con};
            maxVoxName = char(virtsens.label(maxVoxIdx));
            
            % MP Analysis
            
            % Max Voxel Analysis details
            Prefix = 'maxVoxel';        
            tag = strcat(Prefix, '_Subject_', subj, '_Condition_', num2str(EVENTS(con)), '_Channel_', num2str(maxVoxIdx), '_', maxVoxName, '_Mode_', modeText{mode} );
            % Display current analysis info:
            % MP dictionary, subject number, condition number, channel name, and mode             
            disp([ 'MP dictionary: ', MPDic, ' | MEG Signal Source: ', DataSource]);
            disp(tag);
           
            % Get channel data 
            channel_data = MP_GetChannelDataEpochs(virtsens, maxVoxIdx);

            if(subtractAvg(mode) == 1) % Subrtact evoked response
                % Subtract trial average from each trial
                AvgTrialData = mean(channel_data.data,1); % Calculate average of trial
                AvgTrialDataMatrix = repmat(AvgTrialData, size(channel_data.data,1), 1);
                channel_data.data = channel_data.data - AvgTrialDataMatrix; % Subtract
            end

            if( strcmp(MPDic,'stochastic') == 1 ) % For stochastic dictionary only
                % Handle femto tesla issue / normalize
                channel_data.data = channel_data.data*norm_factor;
                % Add special configurations
                MP_config.adaptiveDictionaryParam = adaptiveDictionaryParam; % Read about this value !!!!!!!!!!!!!!!!!!!!!
                MP_config.dictionarySize = dictionarySize; 
            end

            MP_config.dictionaryType = MPDic;
            MP_config.maxIterations = MPmaxIterations; % Max number of iterations for MP algorithm 
            MP_config.timeVals = channel_data.timeVals;
            MP_config.Fs = channel_data.Fs;
            MP_config.runDecomposition = runDecomposition;
            
            % Decomposition into atoms
            [gaborInfo, header, timeVals] = MP_Decomposition(channel_data, MP_config, intermediateOuputPath, tag);

            % Raw energy reconstruction
            MPoutput = MP_Reconstruction(gaborInfo, header, timeVals, MP_config.Fs, MP_config.dictionaryType);

            % Save raw MP output data
            filename = [outputPath, 'MPData_', tag, '.mat'];
            save (filename, 'MPoutput');  
            
            FsAll(end+1,:) = {MP_config.Fs, filename};
            
            %% MP analysis on  max sensor: 
           
            % Max sensor info
            maxChName = char(sensorinfo.ChName{con});
            maxChIdx = find(contains(data_sensors.label,maxChName));
            
            % MP Analysis
            
            % Max Voxel sensor details
            Prefix = 'maxSensor';        
            tag = strcat(Prefix, '_Subject_', subj, '_Condition_', num2str(EVENTS(con)), '_Channel_', num2str(maxChIdx), '_', maxChName, '_Mode_', modeText{mode} );
            % Display current analysis info:
            % MP dictionary, subject number, condition number, channel name, and mode             
            disp([ 'MP dictionary: ', MPDic, ' | MEG Signal Source: ', DataSource]);
            disp(tag);
            
            % Get channel data 
            channel_data = MP_GetChannelDataEpochs(data_sensors, maxChIdx);               
                                         
            if(subtractAvg(mode) == 1) % Subrtact evoked response
                % Subtract trial average from each trial
                AvgTrialData = mean(channel_data.data,1); % Calculate average of trial
                AvgTrialDataMatrix = repmat(AvgTrialData, size(channel_data.data,1), 1);
                channel_data.data = channel_data.data - AvgTrialDataMatrix; % Subtract
            end

            if( strcmp(MPDic,'stochastic') == 1 ) % For stochastic dictionary only
                % Handle femto tesla issue / normalize
                channel_data.data = channel_data.data*norm_factor;
                % Add special configurations
                MP_config.adaptiveDictionaryParam = adaptiveDictionaryParam; % Read about this value !!!!!!!!!!!!!!!!!!!!!
                MP_config.dictionarySize = dictionarySize; 
            end

            MP_config.dictionaryType = MPDic;
            MP_config.maxIterations = MPmaxIterations; % Max number of iterations for MP algorithm 
            MP_config.timeVals = channel_data.timeVals;
            MP_config.Fs = channel_data.Fs;
            MP_config.runDecomposition = runDecomposition;

            % Decomposition into atoms
            [gaborInfo, header, timeVals] = MP_Decomposition(channel_data, MP_config, intermediateOuputPath, tag);

            % Raw energy reconstruction
            MPoutput = MP_Reconstruction(gaborInfo, header, timeVals, MP_config.Fs, MP_config.dictionaryType);

            % Save raw MP output data
            filename = [outputPath, 'MPData_', tag, '.mat'];
            save (filename, 'MPoutput');  
            
            FsAll(end+1,:) = {MP_config.Fs, filename}; 
            
            % Delete intermediate data
            if ( DeleteIntermediateData == 1 && exist(intermediateOuputPath,'file') > 0 )
                rmdir(intermediateOuputPath, 's');
            end
        
        end % Mode
        
    end % Condition

    %exportToPPTX('saveandclose', PPTXname);
     
     
end % Subject

filename = [outputPath,'/FsAll.mat'];
save (filename, 'FsAll');

%%
end_time = clock; % [year month day hour minute seconds]
format bank
disp('start time: ');
start_time
disp('end time: ');
end_time
