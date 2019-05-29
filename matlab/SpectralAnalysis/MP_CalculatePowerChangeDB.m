% Calculate change in power in DB
function dEnDB = MP_CalculatePowerChangeDB(MPdata, baseline_range)

    % Caculate log mean for change in power plots
    logMeanEn = log10(MPdata.meanEnergy); % For all trials
    
    % Calculate change in power from baseline in DB
    blL = find(MPdata.time>=baseline_range(1),1); % Lower bound index of baseline time
    blU = find(MPdata.time<baseline_range(2),1,'last'); % Upper bound index of baseline time             
    baselineEn = mean(logMeanEn(:,blL:blU),2); % Average baseline energy (over time)
    dEnDB = 10*(logMeanEn-repmat(baselineEn,1,size(logMeanEn,2))); % Power change in dB

end
             