%%%%%%%%%%%%%%% Max Voxel/Sensor - MP %%%%%%%%%%%%%%%
% This program calculates group averages for NT and ASD subjects and plots
% power change and spectra graphs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment Setup

% Close and clear everything
clear all;
close all;
clc;

%% Path and Folder Setup

restoredefaultpath;
addpath('SpectralAnalysis');

DataSource  = 'MaxSenVox';
MPDic = 'dyadic'; % Matching pursuit dicitonary: 'dyadic' or 'stochastic'
outputPathPrefix = 'MPOutputData_';
figurePathPrefix = 'MPFigures_';
outputPath = strcat(outputPathPrefix, MPDic, '_',DataSource, '/');
figurePath = strcat(figurePathPrefix, MPDic, '_',DataSource, '/');
ClearPreviousOuput = 1;
MP_SetupOuputFolder(ClearPreviousOuput, figurePath); % Figures

%% MP Analysis Parameters

% Subjects
% NT: 
SUBJ1 = ['0101'; '0102'; '0103'; '0104'; '0105'; '0136'; '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178'; '0179'; '0255'; '0257'; '0348'; '0378'; '0384'];
% ASD: 
SUBJ2 = ['0107'; '0139'; '0141'; '0159'; '0160'; '0161'; '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0346'; '0347'; '0351'; '0358'; '0380'; '0381'; '0382'; '0383']; 

EVENTS = [2,4,8]; % List of event (i.e. condition) number as in data file names e.g. 0076_V1, 0076_V2, etc...

% Analysis parameters
% They should be the same as in MP_RunAnalysis_Groups.m
% It is assumed that the parameters are identical for all data for comparison
subtractAvg = [0]; % Subtract trial average from data? 1:yes, 0:no e.g. [0;1]
modeText = {'NoSubEvoked'}; % For average subtraction above e.g. {'NoSubEvoked'; 'SubEvoked'}
baseline_range = [-0.7 -0.1];
post_stimulus_range = [0.4 1];


%% Plotting Parameters
plotSource = {'maxSensor','maxVoxel'};
SaveMPData = 1; % Save raw MP ouput data: 0:No, 1:Yes
SaveFigures = 1; % Save output figures: 0:No, 1:Yes

%% Get frequency value for data in question
% Make sure that all data is sampled with the same frequency rate
FsPath = strcat(outputPath, 'FsAll.mat');
FsList = load(FsPath);
FsList = FsList.FsAll;
FsUnique = unique(cell2mat(FsList(:,1)));
if size(FsUnique,1) > 1
    error('Fs is not the same for all data!');
    return;
end


%% Calculate group average and plot

for ps=1:size(plotSource,2)
    
    Prefix = plotSource{ps}; 
    
    for mode=1:size(subtractAvg,1)

        % Max Sensor

        
        curr_folder = outputPath;

        % Calculate group average NT
        Group = 'NT';
        Avg1 = MP_GetGroupAverage(curr_folder, SUBJ1, EVENTS, Prefix, modeText{mode}); 
        tag = strcat(Prefix, '_Group_', Group, '_Mode_', modeText{mode} );   
        filename = [outputPath, '/', 'MPGrpAvgData_', tag, '.mat'];
        save(filename, strcat('Avg1') );

        % Calculate group average ASD
        Group = 'ASD';
        Avg2 = MP_GetGroupAverage(curr_folder, SUBJ2, EVENTS, Prefix, modeText{mode});
        tag = strcat(Prefix, '_Group_', Group, '_Mode_', modeText{mode} ); 
        filename = [outputPath, '/', 'MPGrpAvgData_', tag, '.mat'];
        save (filename, strcat('Avg2') );

        % Get number of conditions
        con_num = size(EVENTS,2); % Should be the same for both groups

        % Assumption: Those are the same for all subjects
        figConfig.frequency = Avg1.frequency;
        figConfig.time = Avg1.time; 
        figConfig.timeLimsS = [Avg1.time(1) Avg1.time(end)]; % Time interval (seconds). Stimulus onset is at 0.

        if (ps==1) % Plot source: max sensor
            figConfig.freqLimsHz = [0 150]; % Frequencies to display (Hz).  
            figConfig.cLims = [-26 -24]; % Colormap limits for spectrogram
            figConfig.cLimsDiff = [-6 6]; % Colormap limits for change in power
        else  % Plot source: max voxel    
            figConfig.freqLimsHz = [0 150]; % Frequencies to display (Hz).   
            figConfig.cLims = [-15 -12];
            figConfig.cLimsDiff = [-4 4]; % Colormap limits for change in power
        end
        
        figConfig.fontSizeLarge = 40;
        figConfig.fontSizeSmall = 20;

        for c=1: con_num % Loop pver conditions

            SpectraList = [];

            % Plot Change in Power

            Group = 'NT';   
            con_data.meanEnergy = Avg1.Avg{1,c};
            con_data.time = Avg1.time;
            dEnDB = MP_CalculatePowerChangeDB(con_data, baseline_range);

            currDataTag = strcat(Prefix, '_Group_', Group, '_Mode_', modeText{mode}, '_Condition_', num2str(c) ); 
            FigTitlePC = strcat(Prefix, ' - Group:', Group, ' - Mode:', modeText{mode}, ' - Condition:', num2str(c) ); 
            figConfig.title = FigTitlePC; 

            figPC = MP_PlotPowerChange(dEnDB,figConfig);

            if(SaveFigures == 1)
                savefig(figPC, strcat(figurePath, 'Plot_PowerChange_', currDataTag));
            end

            % Data for spectra plots
            dEnRatio = MP_CalculatePowerChangeRatio(con_data, baseline_range, post_stimulus_range);              
            SpectraData.dEnRatio = dEnRatio;
            SpectraData.frequency = figConfig.frequency; 
            SpectraData.category = Group;
            SpectraList = [SpectraList;SpectraData];

            con_data.meanEnergy = Avg2.Avg{1,c};
            con_data.time = Avg2.time;
            Group = 'ASD'; 
            dEnDB = MP_CalculatePowerChangeDB(con_data, baseline_range);

            currDataTag = strcat(Prefix, '_Group_', Group, '_Mode_', modeText{mode}, '_Condition_', num2str(c) ); 
            FigTitlePC = strcat(Prefix, ' - Group:', Group, ' - Mode:', modeText{mode}, ' - Condition:', num2str(c) ); 
            figConfig.title = FigTitlePC;

            figPC = MP_PlotPowerChange(dEnDB,figConfig);

            if(SaveFigures == 1)
                savefig(figPC, strcat(figurePath, 'Plot_PowerChange_', currDataTag));
            end 

            % Data for spectra plots
            dEnRatio = MP_CalculatePowerChangeRatio(con_data, baseline_range, post_stimulus_range);              
            SpectraData.dEnRatio = dEnRatio;
            SpectraData.frequency = figConfig.frequency; 
            SpectraData.category = Group;
            SpectraList = [SpectraList;SpectraData];

            % Plot Spectra
            currSpectraTag = strcat(Prefix, '_Mode_', modeText{mode}, '_Condition_', num2str(c) );
            FigTitleSP = strcat(Prefix, ' - Mode:', modeText{mode}, ' - Condition:', num2str(c) );

            figSpecConfig.freqLimsHz = [0 150]; % Frequencies to display (Hz).  
            figSpecConfig.timeLimsS = [figConfig.time(1) figConfig.time(end)]; % Time interval (seconds). Stimulus onset is at 0.
            figSpecConfig.fontSizeSmall = 14;
            figSpecConfig.title = FigTitleSP; 

            figSpectra = MP_PlotSpectra(SpectraList,figSpecConfig);

            if(SaveFigures == 1)
                savefig(figSpectra, strcat(figurePath, 'Plot_Spectra_', currSpectraTag) );
            end

        end  % Condition

    end % Mode

end % Plot Source

%% Clear files if requested to do so
if ( SaveMPData == 0 && exist(outputPath,'file') > 0 )
    rmdir(outputPath, 's');
end

if ( SaveFigures == 0 && exist(figurePath,'file') > 0 )
    rmdir(figurePath, 's');
end


