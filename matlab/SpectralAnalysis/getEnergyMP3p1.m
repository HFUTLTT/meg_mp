function meanE = getEnergyMP3p1(gaborInfo,header,timeVals,freqVals)

numTrials = size(gaborInfo,1);
%Fs=round(1/(timeVals(2)-timeVals(1)));
N=length(timeVals);

%freqVals = 0:Fs/N:Fs/2-Fs/N;
%freqVals = 0:Fs/N:Fs/2;

%sumE = zeros(floor(N/2),N);
sumE = zeros(length(freqVals),N);

for i=1:numTrials
    disp(['Computing Energy for trial ' num2str(i) ' of ' num2str(numTrials)]);
    E = mp2tf(squeeze(gaborInfo(i,:,:)),header(i,:));
    sumE = sumE+E;
end

meanE = sumE/numTrials;

end