%Testing RMS calculation

myData = blinkData_csPlus_fullFiltered;
samplingInterval = 1000/samplingRate; %in ms
nSamples = floor(((csTime+traceTime)/samplingInterval)-2*phaseOffsetFactor); %in number of samples
nSamples_pretone = floor((preTime/samplingInterval)-2*phaseOffsetFactor); %in number of sammples
nTrials = size(myData,1);
samplingInterval = 1000/samplingRate; %in ms
nSamples_csNtrace = floor(((csTime+traceTime)/samplingInterval)-2*phaseOffsetFactor); %in number of samples
nSamples_pretone = floor((preTime/samplingInterval)-2*phaseOffsetFactor); %in number of sammples

nBins_pretone = floor(nSamples_pretone/nSamples_csNtrace) ; %number of bins

%Calculate RMS values
Y = myData(:,((phaseOffsetFactor):(nSamples_pretone-(phaseOffsetFactor))));
X = reshape(Y(:,(end-(nSamples_csNtrace*nBins_pretone)+1):end), nTrials, nSamples_csNtrace, nBins_pretone); %1st dimension = trials; 2nd dimension = samples in bin; 3rd dimension = bin index
Z = zeros(nTrials,1,nBins_pretone);

for trial = 1:nTrials
    %Pre-Tone Phase
    for bin = 1:nBins_pretone
        Z(trial,:,bin) = rms(X(trial,:,bin));
        preTone_rms = squeeze(mean(Z,3));
    end
    %CS Phase
    cs_rms(trial) = (rms(myData(trial,((nSamples_pretone+phaseOffsetFactor+1):(nSamples_pretone+phaseOffsetFactor+nSamples_csNtrace)))))'; %Transposed
end