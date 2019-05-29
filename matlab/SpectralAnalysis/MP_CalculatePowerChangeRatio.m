% Calculate change in power in DB
function dEnRatio = MP_CalculatePowerChangeRatio(MPdata, baseline_range, post_stimulus_range)

    % Caculate change in power 
    MeanEn = MPdata.meanEnergy; % For all trials
    
    % Calculate change in power from baseline in DB
    blL = find(MPdata.time>=baseline_range(1),1); % Lower bound index of baseline time
    blU = find(MPdata.time<baseline_range(2),1,'last'); % Upper bound index of baseline time   
    
    % Save spectra data for plotting later in the next step for this subject
    % Power change as ratio: [post-pre]/pre
    blEnAvg = mean(MeanEn(:,blL:blU),2); % Baseline average
    dEnRatio = (MeanEn-repmat(blEnAvg,1,size(MeanEn,2)))./blEnAvg; % Power change 
    psL = find(MPdata.time>=post_stimulus_range(1),1); % Lower bound index of post-stimulus time
    psU = find(MPdata.time<post_stimulus_range(2),1,'last'); % Upper bound index of post-stimulus time
    dEnPs = dEnRatio(:,psL:psU); % Change in power in post-stimulus range
    dEnRatio = mean(dEnPs,2); % Average change in power (over time) for post-stimulus range

end
             