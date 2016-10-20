% AUTHOR: Kambadur Ananthamurthy
% Eye-Blink Conditioning Analysis: Calculate Performance based on
% Thresholded Blinks
% This function takes in multi-trial data.
% Typically, this function should be used with an eye-blink analysis code.

function score = scorePerformance_thresholdedBlinks(myData, samplingRate, phaseOffsetFactor, traceTime, preTime, csTime)

samplingInterval = 1000/samplingRate; %in ms
nSamples_csNtrace = floor(((csTime+traceTime)/samplingInterval)-2*phaseOffsetFactor); %in number of samples
nSamples_pretone = floor((preTime/samplingInterval)-2*phaseOffsetFactor); %in number of sammples
nTrials = size(myData,1);
nBins_pretone = floor(nSamples_pretone/nSamples_csNtrace) ; %number of bins
continuous_blinkBinThreshold = nBins_pretone/2;

for trial = 1:nTrials
    myData_std(trial)= std(myData(trial,:));
    for sample = 1:size(myData,2)
        if myData(trial, sample)>= 2*myData_std(trial)
            thresholdedBlinks(trial,sample) = 1;
        else
            thresholdedBlinks(trial,sample) = 0;
        end
    end
end
clear trial
clear sample

%Score Pre-Tone Phase
preToneBlinks = thresholdedBlinks(:,((phaseOffsetFactor):(nSamples_pretone-(phaseOffsetFactor))));
Y=preToneBlinks(:, end-(nBins_pretone*nSamples_csNtrace)+1:end);
X=reshape(Y, size(thresholdedBlinks,1), nSamples_csNtrace, nBins_pretone); %1st dimension = trials; 2nd dimension = samples in bin; 3rd dimension = bin index
nPretoneBlinks = squeeze(sum(sum(X,2)>0,3));
%clipping
preToneBlinks(preToneBlinks > continuous_blinkBinThreshold) = continuous_blinkBinThreshold;
%Score CS phase
blinks_csPlus_cs = thresholdedBlinks(:,((nSamples_pretone+phaseOffsetFactor+1):(nSamples_pretone+phaseOffsetFactor+nSamples_csNtrace)));
nBlinks_csPlus = (sum(blinks_csPlus_cs,2)>0);
%Get Trial Score
csHitScore = (1 - (nPretoneBlinks/4)).*nBlinks_csPlus; %per CS+ trial
csMissScore = (nPretoneBlinks/4).*(1-nBlinks_csPlus); %per CS+ trial
score = mean(csHitScore - csMissScore); % per CS+ session

end