% Get channel list based on data source and desired channel options 
% This must be adjusted based on the dataset format/structure
function [ChNames, ChIndices] = MP_GetChannels(DataSource, MaxChannel, dataPrep)

    if( DataSource == "LTC" )   

        if(MaxChannel == 1)
            % Max Channel
            ChNames = [dataPrep.MaxCh]; % Channel name
            ChIndices = [find(ismember(dataPrep.ch_names ,dataPrep.MaxCh, 'rows'))]; % Channel Index
        else
            % Multiple Channels
            ChNames = dataPrep.CHnames;
            ChIndices = [1:size(dataPrep.CHnames,1)];
        end

    elseif( DataSource == "HTC" )

        if(MaxChannel == 1)
            % Max Channel
            ChIndices = [2];
            ChNames = [dataPrep.ch_names(ChIndices,:)];
        else
            % Get all channels
            ChIndices = [1:size(dataPrep.ch_names,1)]; % Channel indices
            ChNames = dataPrep.ch_names; % Channel names
        end

    end
    
end