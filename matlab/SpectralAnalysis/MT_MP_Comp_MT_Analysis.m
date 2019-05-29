%%%%%%%%%%%%%%% LTC Data Analysis - MT %%%%%%%%%%%%%%%
% This program computes time-frequency spectrum using MT algorithm for preprocessed and epoched MEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Environment Setup

% Close and clear everything
clc;
clear all;
close all;
%%
screensize = get( groot, 'Screensize' );

% Path Setup
restoredefaultpath;

% Code Paths
addpath('fieldtrip/');
ft_defaults; % Fieldtrip path setup
addpath('SpectralAnalysis');

inputPath = 'SubjectData/MT/';
outputPath = 'MT_Output';


ClearPreviousOuput = 1;
MP_SetupOuputFolder(ClearPreviousOuput, outputPath); % Figures

%%
% parameters for the multitaper analysis
tapsmofrq = 5;

% Here I define the broad gamma range in Hz
gammarange = [35 110];
gammamax = [45 90]; % low gamma often is suppressed at the fast velocity and EMG-artefact may dominate the response, 
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

SUBJ = [ '0254'; '0259']; % Subjects

foi = [5:1:150]; %freq band of interest
toi = [-1.0:0.01:1.2]; % Time points for window centers
                            
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
        [val,ind] =max(R)
        
        % MEG sensor with max gamma increase is: 
        t = ch(ind) % sensor number
        %ChName{con} = epochs.label{ch(ind)} % sensor name
        
        % Scale time window (t_ftimwin) with frequency.
        % Do TMF in max channel       
        cfg = [];
        cfg.method    = 'mtmconvol';
        cfg.output    ='pow'; % 'fourier'
        cfg.taper        = 'dpss';
        cfg.keeptrials = 'no';
        cfg.foi    = foi; %freq band of interest
        cfg.toi    = toi;    
        cfg.tapsmofrq  = cfg.foi*0.2; 
        cfg.t_ftimwin  = 5./cfg.foi;
        cfg.channel = t;
        cfg.pad = 'nextpow2';
        TMF = ft_freqanalysis(cfg, epochs_con);  % trials x Ch x freq
        % plot TMF
        cfg = [];
        cfg.baseline = [-0.8 0.0];
        cfg.baselinetype = 'db'; 	%  data = (data - meanVals) ./ meanVals;
        h=figure('Position', get(0, 'Screensize')); 
        ax = gca;
        ft_singleplotTFR(cfg, TMF);
        currTitle = strcat('MT Analysis for Max Channel - Subject:', subj, '- Condition:', num2str(con));
        title( strcat(currTitle, ' - Frequency-scaled Window:  0.4 sec'),  'FontSize', 20 );
        ax.XTick = [toi(1):0.1:toi(end)];
        ax.YTick = [foi(1):10:foi(end)];
        ax.FontSize = 20;
        ax.XLabel.String = "Time (s)";
        ax.YLabel.String = "Frequency (Hz)";
        c = colorbar;
        c.Label.String = "Change in Power (dB)";
        c.Label.FontSize = 20;
        currTag = strcat('Subject_', subj, '_MaxChannel_Condition_', num2str(con));
        saveas(h, strcat(outputPath,'/', 'MT_ScaledWindow_PowerChange_', currTag), 'jpeg');

        % Keep time window constant:  0.4 sec.
        cfg = [];
        cfg.method    = 'mtmconvol';
        cfg.output    ='pow'; % 'fourier'
        cfg.taper        = 'dpss';
        cfg.keeptrials = 'no';
        cfg.tapsmofrq = tapsmofrq;  
        cfg.foi    = foi; %freq band of interest
        cfg.toi    = toi;    
        cfg.t_ftimwin  = ones(1, length(cfg.foi))* 0.4 ;
        cfg.channel = t;
        cfg.pad = 'nextpow2';
        TMF = ft_freqanalysis(cfg, epochs_con);  % trials x Ch x freq
        % plot TMF
        cfg = [];
        cfg.baseline = [-0.8 0.0];
        cfg.baselinetype = 'db'; 	%  data = (data - meanVals) ./ meanVals;
        h=figure('Position', get(0, 'Screensize'));  
        ax = gca;
        ft_singleplotTFR(cfg, TMF);             
        currTag = strcat('Subject_', subj, '_MaxChannel_Condition_', num2str(con));
        currTitle = strcat('MT Analysis for Max Channel - Subject:', subj, '- Condition:', num2str(con));
        title( strcat(currTitle, ' - Fixed Window:  0.4 sec'),  'FontSize', 20 );
        ax.XTick = [toi(1):0.1:toi(end)];
        ax.YTick = [foi(1):10:foi(end)];
        ax.FontSize = 20;
        ax.XLabel.String = "Time (s)";
        ax.YLabel.String = "Frequency (Hz)";
        c = colorbar;
        c.Label.String = "Change in Power (dB)";
        c.Label.FontSize = 20;
        saveas(h, strcat(outputPath,'/', 'MT_FixedWindow_PowerChange_', currTag), 'jpeg');
        
    end  % end for conditions in sensor space    


end   % end for subj
 


