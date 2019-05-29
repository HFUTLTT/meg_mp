% This function performs MP decomposition for a dataset of dimensions:
% number of trials x number of time points
function [gaborInfo, header, timeVals] = MP_Decomposition(channel_data, MP_config, intermediateOuputPath, tag)
        
    % Get data parameters
    dims = size(channel_data.data);
    L1 = dims(2); % Signal length
    
    % Setup Analysis Parameters    
    L = L1;
    data = channel_data.data;
    timeVals = MP_config.timeVals;
    gaborInfoFile = strcat(intermediateOuputPath, 'Gabor_',tag);
    
    if( strcmp(MP_config.dictionaryType,'dyadic') == 1 ) % Only for dyadic dictionary
        
        % A shorter signal is taken if the current length is not a power of 2
        % as this method only works with powers of 2
        p1 = nextpow2(L1);
        p2 = p1-1;
        if(2^p1 ~= L1)
            L = 2^p2;
            % Exact same action should done for data and time points
            data = channel_data.data(:,1:L);
            timeVals = MP_config.timeVals(1:L);          
        end
        
        if (MP_config.runDecomposition == 1)
            
            disp('Performing MP decomposition...');
            
            %Import data
            X = data';
            signalRange = [1 L]; % Signal length range
            importData(X, intermediateOuputPath, tag, signalRange, MP_config.Fs);
            
            % Perform Gabor decomposition
            runGabor(intermediateOuputPath, tag, L, MP_config.maxIterations);
            
        else
            disp('Using previously decomposed MP data...');
        end
        
        gaborInfo = getGaborData(intermediateOuputPath,tag,1);
        header = [];
        save(gaborInfoFile,'gaborInfo');
        
    elseif( strcmp(MP_config.dictionaryType,'stochastic') == 1 ) % Only for stochastic dictionary   

        if (MP_config.runDecomposition == 1)
            
            disp('Performing MP decomposition...');
            
            % Perform Gabor decomposition  
            [gaborInfo, header] = getStochasticDictionaryMP3p1(data, timeVals, MP_config.maxIterations, MP_config.adaptiveDictionaryParam, MP_config.dictionarySize);
            save(gaborInfoFile,'gaborInfo','header');
            
        else
            
            disp('Using previously decomposed MP data...');
            
            if ( exist(intermediateOuputPath,'file') > 0 )
                
                disp(['Opening saved file ' gaborInfoFile]);
                
                savedData = load(gaborInfoFile);
                
                header = savedData.header;
                gaborInfo = savedData.gaborInfo;
            end
            
        end
        
    end

end
             