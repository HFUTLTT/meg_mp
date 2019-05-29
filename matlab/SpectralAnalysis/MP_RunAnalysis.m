%%%%%%%%%%%%%%% HTC, LTC Data Analysis - MP %%%%%%%%%%%%%%%
% This program computes time-frequency spectrum using MP algorithm for preprocessed and epoched MEG data
% This implementation supports 2 types of MP dictionaries:
% 1. Dyadic.
% 2. Stochastic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment Setup

% Close and clear everything
clc;
clear all;
close all;

%%
% Path Setup
restoredefaultpath;

% Code Paths
addpath('fieldtrip/');
ft_defaults; % Fieldtrip path setup
addpath('SpectralAnalysis');

% File Paths
DataSource  = 'LTC'; % MEG recording source 'LTC' or 'HTC'
% Input File Paths
inputPath = strcat('SubjectData/', DataSource, '/'); % Source data
intermediatePathPrefix = 'MPIntermediateData_';
outputPathPrefix = 'MPOutputData_';
figurePathPrefix = 'MPFigures_';

%% MP Analysis Parameters
% Specified based on use-case by the user

% Modes
MPDic = 'stochastic'; % Matching pursuit dicitonary: 'dyadic' or 'stochastic'
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
MPmaxIterations = 5 ; % Number of MP iterations 500

% Data
SUBJ = ['0076']; % List of subjects e.g. ['0076'; '0177'; '0259']
EVENTS = [1,2]; % List of event (i.e. condition) number as in data file names e.g. 0076_V1, 0076_V2, etc...

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
figurePath = strcat(figurePathPrefix, MPDic, '_',DataSource, '/');
MP_SetupOuputFolder(ClearPreviousOuput, outputPath); % MP intermediate data
MP_SetupOuputFolder(ClearPreviousOuput, intermediateOuputPath); % MP ouput data
MP_SetupOuputFolder(ClearPreviousOuput, figurePath); % Figures

% Set variable values
con_num = size(EVENTS,2); % Get the total number of events

%% Analysis Process

for s=1:size(SUBJ,1) % Loop over subjects
    
    subj = SUBJ(s,:); % Get current subject number
        
    for con = 1:con_num % Loop over conditions for current subject
        
        currCon = EVENTS(con);
        % Load data for current condition and subject (dataPrep = preporocessed data)
        dataPrep = load( strcat( inputPath,subj,'_V', num2str(currCon) ) ); 

        % Get channel list based on data source and desired channel options 
        [ChNames, ChIndices] = MP_GetChannels(DataSource, MaxChannel, dataPrep);

        for channel = 1:size(ChIndices,2) % Loop over channels  
            
            currChIdx = ChIndices(channel); % Current channel index
            currChName = ChNames(channel, :);
            
            for mode=1:size(subtractAvg,1) % Loop over modes (evoked response subraction: no, yes)
                      
                % Current analysis details
                tag = strcat('Subject_', subj, '_Condition_', num2str(currCon), '_Channel_', num2str(currChIdx), '_', currChName, '_Mode_', modeText{mode} );
                % Display current analysis info:
                % MP dictionary, subject number, condition number, channel name, and mode             
                disp([ 'MP dictionary: ', MPDic, ' | MEG Signal Source: ', DataSource]);
                disp(tag);
                
                % Get channel data based on data source and desired channel options 
                channel_data = MP_GetChannelData(DataSource, MaxChannel, dataPrep, currChIdx);
             
                if(subtractAvg(mode) == 1) % Subrtact evoked response
                    % Subtract trial average from each trial
                    AvgTrialData = mean(channel_data.data,1); % Calculate average of trials
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
                
                % Perform MP analysis
                
                % Decomposition into atoms
                [gaborInfo, header, timeVals] = MP_Decomposition(channel_data, MP_config, intermediateOuputPath, tag);
     
                % Raw energy reconstruction
                MPoutput = MP_Reconstruction(gaborInfo, header, timeVals, MP_config.Fs, MP_config.dictionaryType);
                             
                % Save raw MP output data
                filename = [outputPath, 'MPData_', tag, '.mat'];
                save (filename, 'MPoutput');  
                
                % Delete intermediate data
                if ( DeleteIntermediateData == 1 && exist(intermediateOuputPath,'file') > 0 )
                    rmdir(intermediateOuputPath, 's');
                end
                
            end % Mode
            
        end % Channel
        
    end % Condition
 
end % Subject
    
%% Plots

for mode=1:size(subtractAvg,1)
    
    for s=1:size(SUBJ,1) % Loop over subjects
        
        subj = SUBJ(s,:);
     
        subj_channels = {}; % Loop over conditions
        for c=1:size(EVENTS,2)
            currC = EVENTS(c);          
            dataPrep = load( strcat( inputPath,subj,'_V', num2str(currC) ) ); 
            [ChNames, ChIndices] = MP_GetChannels(DataSource, MaxChannel, dataPrep);
            ChannelConst.ChNames = ChNames;
            ChannelConst.ChIndices = ChIndices;
            subj_channels{end+1} = ChannelConst;
        end
        
        % Assumption: 
        % All conditions have the same number of channels for the same subject
        for channel = 1:size(ChIndices,2) % Loop over channels
                       
            SpectraList = [];
                      
            for con=1:size(EVENTS,2) % Loop over conditions for current subject
             
                currCon = EVENTS(con);
                
                con_channels = subj_channels{con};
                
                currChIdx = con_channels.ChIndices(channel); % Current channel index
                currChName = con_channels.ChNames(channel, :);
                
                currDataTag = strcat('Subject_', subj, '_Condition_', num2str(currCon) , '_Channel_', num2str(currChIdx), '_', currChName, '_Mode_', modeText{mode});
                FigTitlePC = strcat('Subject:', subj, ' - Condition: ', num2str(currCon) , ' - Channel #', num2str(currChIdx), ': ', currChName, ' - Mode: ', modeText{mode});
                filePrefix = 'MPData_';
                
                MPdata = load( strcat(outputPath, filePrefix, currDataTag, '.mat' ) );
                MPdata = MPdata.MPoutput;
                
                dEnDB = MP_CalculatePowerChangeDB(MPdata, baseline_range);
                
                figConfig.frequency = MPdata.frequency;
                figConfig.time = MPdata.time;               
                figConfig.freqLimsHz = [0 150]; % Frequencies to display (Hz).  
                figConfig.timeLimsS = [MPdata.time(1) MPdata.time(end)]; % Time interval (seconds). Stimulus onset is at 0.
                figConfig.cLims = [-26 -24]; % Colormap limits for spectrogram
                figConfig.cLimsDiff = [-6 6]; % Colormap limits for change in power
                figConfig.fontSizeLarge = 40;
                figConfig.fontSizeSmall = 20;
                figConfig.title = FigTitlePC;
        
                figPC = MP_PlotPowerChange(dEnDB,figConfig);
                
                if(SaveFigures == 1)
                    savefig(figPC, strcat(figurePath, 'Plot_PowerChange_', currDataTag));
                end
               
                dEnRatio = MP_CalculatePowerChangeRatio(MPdata, baseline_range, post_stimulus_range);
                
                SpectraData.dEnRatio = dEnRatio;
                SpectraData.frequency = MPdata.frequency; 
                SpectraData.category = num2str(currCon);
           
                SpectraList = [SpectraList;SpectraData];
                
            end % Conditions
            
            FigTitleSP = strcat('Subject:', subj, ' - Channel #', num2str(currChIdx), ': ', currChName, ' - Mode: ', modeText{mode});
            currSpectraTag = strcat('Subject_', subj, '_Channel_', num2str(currChIdx), '_', currChName, '_Mode_', modeText{mode});
           
            figSpecConfig.freqLimsHz = [0 150]; % Frequencies to display (Hz).  
            figSpecConfig.timeLimsS = [MPdata.time(1) MPdata.time(end)]; % Time interval (seconds). Stimulus onset is at 0.
            figSpecConfig.fontSizeSmall = 14;
            figSpecConfig.title = FigTitleSP; 

            figSpectra = MP_PlotSpectra(SpectraList,figSpecConfig);
                        
            if(SaveFigures == 1)
                savefig(figSpectra, strcat(figurePath, 'Plot_Spectra_', currSpectraTag) );
            end
            
        end % Channels
        
    end % Subject
    
end % Mode

%% Clear files if requested to do so
if ( SaveMPData == 0 && exist(outputPath,'file') > 0 )
    rmdir(outputPath, 's');
end

if (SaveFigures == 0 && exist(figurePath,'file') > 0 )
    rmdir(figurePath, 's');
end
