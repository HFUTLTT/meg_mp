%%%%%%%%%%%%%%% LTC Data Analysis - MT %%%%%%%%%%%%%%%
% This program computes time-frequency spectrum using MT algorithm for preprocessed and epoched MEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Environment Setup

% Close and clear everything
clc;
clear all;
close all;

screensize = get( groot, 'Screensize' );

% Path Setup
restoredefaultpath;

% Code Paths
addpath('fieldtrip/');
ft_defaults; % Fieldtrip path setup
addpath('SpectralAnalysis');

DataSource  = 'LTC'; % MEG recording source 'LTC' or 'HTC'

% File path prefix
inputPath = 'SubjectData/MT/';
intermediatePathPrefix = 'MPIntermediateData_';
outputPathPrefix = 'MP_MT_Comp_Data_';
figurePathPrefix = 'MP_MT_Comp_Figures_';

%% MP Analysis Parameters
% Specified based on use-case by the user

% Modes
MPDic = 'dyadic'; % Matching pursuit dicitonary: 'dyadic' or 'stochastic'
MaxChannel = 1; % Maximum channel or multiple channels 0:multiple channels, 1:max channel
% Subtraction: Both lines below should be modified in an identical manner to relfect desired mode(s)
subtractAvg = [0]; % Subtract evoked response? 0:no, 1:yes e.g [0;1]
modeText = {'NoSubEvoked'}; % e.g. {'NoSubEvoked'; 'SubEvoked'}

% Options
runDecomposition = 1; % 0: Load previously decomposed data 1: Run decomposition from scratch
ClearPreviousOuput = 1; % Clear output (including intermediate data) from previous execution: 0: No, 1:Yes
DeleteIntermediateData = 0; % Delete intermediate MP data: 0:No, 1:Yes
SaveMPData = 1; % Save raw MP ouput data: 0:No, 1:Yes
SaveFigures = 1; % Save output figures: 0:No, 1:Yes
MPmaxIterations = 500 ; % Number of MP iterations

% Time ranges of interest (used for power change calculation and plots)
baseline_range = [-0.8 0];
post_stimulus_range = [0.4 1.2];

% Special parameters for "Stochastic" dictionary only
norm_factor = 10^12; % Normalization factor for stochastic dictionary (magnetic field magnitude)
adaptiveDictionaryParam = 0.9; % ADP for subset selection (size = (1-ADP)*100%)
dictionarySize = 2500000; % Dictionary size

%%
% Output folder setup
intermediateOuputPath = strcat(intermediatePathPrefix, MPDic, '_',DataSource, '/');
outputPath = strcat(outputPathPrefix, MPDic, '_',DataSource, '/');
figurePath = strcat(figurePathPrefix, MPDic, '_',DataSource, '/');
MP_SetupOuputFolder(ClearPreviousOuput, outputPath); % MP intermediate data
MP_SetupOuputFolder(ClearPreviousOuput, intermediateOuputPath); % MP ouput data
MP_SetupOuputFolder(ClearPreviousOuput, figurePath); % Figures
%%
% parameters for the multitaper analysis
tapsmofrq = 5;

% Here I define the broad gamma range in Hz
gammarange = [35 110];
gammamax = [45 90]; % low gamma often is suppresed at the fast velocity and EMG-artefact may dominate the response, 
                    % therefore the 'maximal gamma sensor' will be defined for the narrow
                    % gamma range

% filters
filt ='no';
if strcmp(filt,'yes')
    lpfreq=gammarange(2)+5;
    hpfreq=gammarange(1)-5;
end
% Some parameters for plots
X=[0.5, 3.5, 6.5; 0.5, 3.5, 6.5]; Y=[1, 1, 1; 4, 4, 4]; % positions for plot
LineColor{1} = 'b'; LineColor{2} = 'g'; LineColor{3} = 'r';

% Sensors used to search for the maximal gamma sensor. These are the posterior channels where the response is normally maximal. 
% This is done to exclude other channels, where the response may result
% from muscle contraction artefats
Ch = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
       'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};

SUBJ = ['0254'; '0259']; % Subjects
                            
%% 

for s=1:size(SUBJ,1)
    
    subj = SUBJ (s,:);
 
    %    
    % Load raw data
    %fiff_file = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_rings_ICA_raw.fif');
    %hdrraw = ft_read_header(fiff_file);
    %pre = -1.0* hdr.Fs ;
    %post = 1.2* hdr.Fs ;
    
    dataPrep = load(strcat(inputPath,subj,'.mat')); % Load epochs
    epochs = dataPrep.epochs;
    
    ev1 = find(dataPrep.events(:,3)==2);
    ev2 = find(dataPrep.events(:,3)==4);
    ev3 = find(dataPrep.events(:,3)==8);
    EV = {ev1, ev2, ev3};
      
    % Initiate the figure 
    %hhh=figure;
    %ax=gca; hold(ax,'on')

    %%
    for con=1:3 % for conditions 

        cfg = [];
        cfg.trials = EV{con}; % EV{con}';
        % Select the whole epochs for this condition
        [epochs_con] = ft_selectdata(cfg, epochs);
        % Select pre and poststimusus epochs for this condition: This is used to find sensor with max change in the gamma power
        cfg.latency = [-0.8 0.0];
        [dataPre] = ft_selectdata(cfg, epochs);
        cfg.latency = [0.4 1.2];   % here we excluded 1st 400 msec (when we have ERP) from analysis 
        [dataPost] = ft_selectdata(cfg, epochs);
        
        % Plot Event related potential (average relative to the stimulatiin start)
        %cfg = [];
        %cfg.channel=epochs.label;
        %avg = ft_timelockanalysis(cfg,epochs_con);
        %hh=figure;
        %ttt = find(avg.time>-0.1 & avg.time<0.5);
        %plot (avg.time(ttt), avg.avg(:, ttt));
        %title (['Sensor averages for condition...', num2str(con)]);
        
        % Do TMF        
        cfg = [];
        cfg.method    = 'mtmfft';
        cfg.output    ='pow'; % 'fourier'
        cfg.taper        = 'dpss';
        cfg.keeptrials = 'no';
        cfg.tapsmofrq = 5;
        cfg.foilim    = gammamax; %freq band to find max gamma power
        freqPre = ft_freqanalysis(cfg, dataPre);  % trials x Ch x freq
        freqPost = ft_freqanalysis(cfg, dataPost);

        % Numbers of posterior channels
        for j=1:length(Ch)
            [ch(j),x] = find(strcmp(freqPost.label,Ch{j}));
        end
        % b1 is the max ch
       
        % Find max Post/pre ratio in max channel, % DIMENSION: time x Ch x freq
        R =  squeeze(mean(freqPost.powspctrm(ch,: ),2))./squeeze(mean(freqPre.powspctrm(ch,: ),2));  
        [val,ind] =max(R);
        
        % MEG sensor with max gamma increase is: 
        t = ch(ind); % sensor number
        %ChName{con} = epochs.label{ch(ind)} % sensor name
        
        channel_data = MP_GetChannelDataEpochs(epochs_con, t);
        
        if( strcmp(MPDic,'stochastic') == 1 ) % For stochastic dictionary only
            % Handle femto tesla issue / normalize
            channel_data.data = channel_data.data*norm_factor;
            % Add special configurations
            MP_config.adaptiveDictionaryParam = adaptiveDictionaryParam; % Read about this value !!!!!!!!!!!!!!!!!!!!!
            MP_config.dictionarySize = dictionarySize; 
        end
        
        MP_config.dictionaryType = MPDic;
        MP_config.maxIterations = MPmaxIterations; % Max number of iterations for MP algorithm 
        MP_config.timeVals = epochs.time{1, 1};
        MP_config.Fs = epochs.fsample;
        MP_config.runDecomposition = runDecomposition;
               
        tag = strcat('Subject_', subj, '_Condition_', num2str(con), '_MaxChannel');
        currTitle = strcat('MP Analysis for Max Channel - Subject:', subj, '- Condition:', num2str(con));
        
        [gaborInfo, header, timeVals] = MP_Decomposition(channel_data, MP_config, intermediateOuputPath, tag);
    
        % Raw energy reconstruction
        MPoutput = MP_Reconstruction(gaborInfo, header, timeVals, MP_config.Fs, MP_config.dictionaryType);
          
        % Save raw MP output data
        filename = [outputPath, 'MPData_', tag, '.mat'];
        save (filename, 'MPoutput');  
 %%       
        MPdata = MPoutput;
        dEnDB = MP_CalculatePowerChangeDB(MPdata, baseline_range);

        figConfig.frequency = MPdata.frequency;
        figConfig.time = MPdata.time;               
        figConfig.freqLimsHz = [0 150]; % Frequencies to display (Hz).  
        figConfig.timeLimsS = [MPdata.time(1) MPdata.time(end)]; % Time interval (seconds). Stimulus onset is at 0.
        figConfig.cLims = [-26 -24]; % Colormap limits for spectrogram
        figConfig.cLimsDiff = [-6 6]; % Colormap limits for change in power
        figConfig.fontSizeLarge = 40;
        figConfig.fontSizeSmall = 20;
        figConfig.title = currTitle;

        figPC = MP_PlotPowerChange(dEnDB,figConfig);

        if(SaveFigures == 1)
            saveas(figPC, strcat(figurePath, 'Plot_PowerChange_', tag), 'jpeg');
        end
        
  %%      
    end  % end for conditions in sensor space    


end   % end for subj
 

%% Clear files if requested to do so
if ( SaveMPData == 0 && exist(outputPath,'file') > 0 )
    rmdir(outputPath, 's');
end

if ( SaveFigures == 0 && exist(figurePath,'file') > 0 )
    rmdir(figurePath, 's');
end
