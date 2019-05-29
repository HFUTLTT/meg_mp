% Get channel data based on data source and desired channel options 
% This must be adjusted based on the dataset format/structure
function channel_data = MP_GetChannelDataEpochs(epochs, ChIdx)

    % Reformat condition data to a matrix for MP
    % Target format: L X M X N
    % L = length of signal, M = number of trials and N = number of channels. 
    trialData = epochs.trial;

    data3D = cat(3,trialData{:});
    inputSignal = permute(data3D,[2 3 1]);

    % Get data from channel with maximum power increase
    tmp_data = inputSignal(:,:,ChIdx);
    
    % Final ouput
    channel_data.data = permute(tmp_data,[2 1]);   
    channel_data.Fs = epochs.fsample;
    channel_data.timeVals = epochs.time{1};
    
end