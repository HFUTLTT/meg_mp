function meanE = getEnergyMPDyadic(gaborInfo,header,timeVals,freqVals)
    
    % Set parameters
    numTrials = size(gaborInfo,2);
    L = length(timeVals);
    wrap=0; % Time wrapping enabled = 1, disabled = 0 
    atomList=[]; % Empty list = use all atoms
    
    recEnAll=zeros(length(freqVals),L);  % Initialize energy matrix
    
    disp('Reconstructing energy...');
 
    for i=1:numTrials

        disp(['trial ' num2str(i) ' of ' num2str(numTrials)]);

        % Compute energy for current trial
        % GaborInfo includes the atom list with info e.g. frequency, time,
        % and phase for each atom used to recontruct the signal
        recEn = reconstructEnergyFromAtomsMPP(gaborInfo{i}.gaborData, L, wrap, atomList);
        recEnAll = recEnAll + recEn; % Sum energy over all trials

    end

    % Average energy over trials
    meanE = recEnAll/numTrials; 

end