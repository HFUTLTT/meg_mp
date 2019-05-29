% This function performs MP energy reconstruction for decomposed signal
function MPoutput = MP_Reconstruction(gaborInfo, header, timeVals, Fs, dictionaryType)
    
    L = length(timeVals);
    frequency = 0: Fs/L : Fs/2; % Frequency: equal points from 0-Fs/2 (Hz)
    
    if( strcmp(dictionaryType,'dyadic') == 1 ) % Only for dyadic dictionary
        MeanEn = getEnergyMPDyadic(gaborInfo,header,timeVals, frequency);
    elseif( strcmp(dictionaryType,'stochastic') == 1 ) % Only for stochastic dictionary  
        MeanEn = getEnergyMP3p1(gaborInfo,header,timeVals, frequency);   
    end

    MPoutput.meanEnergy = MeanEn;
    MPoutput.frequency = frequency;
    MPoutput.time = timeVals;
    MPoutput.Fs = Fs;
  
end