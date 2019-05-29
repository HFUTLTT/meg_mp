% Get channel data based on data source and desired channel options 
% This must be adjusted based on the dataset format/structure
function channel_data = MP_GetChannelData(DataSource, MaxChannel, dataPrep, ChIdx)

    if( DataSource == "LTC" ) 

        if(MaxChannel == 1)
            % Max Channel 
            data = dataPrep.MaxCh_data; % Load channel data            
            Fs = 1000; % Specify sampling frequency
            dt = 1/Fs; % Calculate delta t
            t_start = dataPrep.times_mt(1); % Start time
            t_end = dataPrep.times_mt(end); % End time
            timeVals = [t_start:dt:t_end]; % Generate time points based on the above
        else
            % Multiple Channels                 
            data = squeeze(dataPrep.data(:,ChIdx,:)); % Load channel data 
            Fs = dataPrep.SF; % Sampling frequency        
            timeVals = dataPrep.times; % Time points
        end

    elseif( DataSource == "HTC" )
        
        if(MaxChannel == 1 || MaxChannel == 0 ) % Max Channel or Multiple Channels
            data = squeeze(dataPrep.data(:,ChIdx,:)); % Load channel data
            Fs = dataPrep.SamplingFreqHz; % Sampling frequency
            timeVals = dataPrep.timepoints; % Time points
        end
        
    end
    
    channel_data.data = data;
    channel_data.timeVals = timeVals;
    channel_data.Fs = Fs;
    
end