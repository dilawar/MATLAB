% AUTHOR: Kambadur Ananthamurthy
% Eye-Blink Conditioning Analysis: Calculate RMS
% This function takes in multi-trial data.
% Typically, this function should be used with an eye-blink analysis code.

function [preTone_rms, cs_rms] = calculateRMS(myData, samplingRate, phaseOffsetFactor, traceTime, preTime, csTime)

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
        %disp(bin);
        Z(trial,1,bin) = rms(X(trial,:,bin),2);
    end
    %CS Phase
    cs_rms(trial) = rms(myData(trial,((nSamples_pretone+phaseOffsetFactor+1):(nSamples_pretone+phaseOffsetFactor+nSamples_csNtrace))));
end
preTone_rms = mean(Z,3);
end