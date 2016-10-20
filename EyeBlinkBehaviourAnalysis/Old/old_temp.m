clear all
close all

mouse = 'K26';

baselineCorrection = 0;

samplingRate = 100; %Hz
preTime = 5000; %ms
csTime = 350; %ms
usTime = 50; %ms

controlSessionType = 3;
TECSessionType = 9;

if controlSessionType == 3
    traceTime = 250;
elseif controlSessionType == 5
    traceTime = 500;
end
if TECSessionType == 9
    traceTime = 250;
elseif TECSessionType == 11
    traceTime = 500;
end

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');
dataset = [mouse '_SessionType' num2str(controlSessionType) '_Session1'];
saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];

samplePoint = samplingRate/1000;
uscorei = strfind(dataset,'_');
sessionType = dataset((uscorei(1)+12):(uscorei(2)-1));

if str2double(sessionType) <=5
    %No-Puff Control
    probe = 0;
else
    %TEC
    probe = 1;
end

if baselineCorrection == 1
    csPlusTrials = csvread([saveFolder 'Mouse' dataset '_csPlus_dRR.csv']);
    csMinusTrials = csvread([saveFolder 'Mouse' dataset '_csMinus_dRR.csv']);
else
    csPlusTrials = csvread([saveFolder 'Mouse' dataset '_csPlus.csv']);
    csMinusTrials = csvread([saveFolder 'Mouse' dataset '_csMinus.csv']);
end

if probe == 1
    probeTrials = csvread([saveFolder 'Mouse' dataset '_probeTrials.csv']);
    csPlusProbeTrials = csPlusTrials(probeTrials,:);
    csPlusTrials(probeTrials,:) = [];
elseif probe == 0;
    csPlusProbeTrials = [];
end

%Preallocation
csPlus = zeros(size(csPlusTrials));
csPlusProbe = zeros(size(csPlusProbeTrials));
csMinus = zeros(size(csMinusTrials));
auc_csPlus = zeros(size(csPlusTrials,1),1);
auc_csPlusProbe = zeros(size(csPlusProbeTrials,1),1);
auc_csMinus = zeros(size(csMinusTrials,1),1);

%Rectify
rectifiedcsPlusTrials = abs(csPlusTrials);
rectifiedcsPlusProbeTrials = abs(csPlusProbeTrials);
rectifiedcsMinusTrials = abs(csMinusTrials);

%Set baseline to zero
if baselineCorrection == 1
    csPlus = rectifiedcsPlusTrials;
    csPlusProbe = rectifiedcsPlusProbeTrials;
    csMinus = rectifiedcsMinusTrials;
else
    for trial = 1:size(rectifiedcsPlusTrials,1)
        csPlus(trial,:) = rectifiedcsPlusTrials(trial,:)-median(rectifiedcsPlusTrials(trial,:));
    end
    for trial = 1:size(rectifiedcsPlusProbeTrials,1)
        csPlusProbe(trial,:) = rectifiedcsPlusProbeTrials(trial,:)-median(rectifiedcsPlusProbeTrials(trial,:));
    end
    for trial = 1:size(rectifiedcsMinusTrials,1)
        csMinus(trial,:) = rectifiedcsMinusTrials(trial,:)-median(rectifiedcsMinusTrials(trial,:));
    end
end
%Set all negative values to zero
csPlus(csPlus<0) = 0;
csPlusProbe(csPlusProbe<0) = 0;
csMinus(csMinus<0) = 0;

%Area under the curve for the csOnset to usOnset period
%Time points
csOnset = preTime*samplePoint;
usOnset = (preTime+csTime+traceTime)*samplePoint;
%usOffset = usOnset+(usTime*samplePoint);

for trial = 1:size(csPlus)
    auc_csPlus(trial) = trapz(csPlus(trial,csOnset:usOnset));
end
for trial = 1:size(csPlusProbe)
    auc_csPlusProbe(trial) = trapz(csPlusProbe(trial,csOnset:usOnset));
end
for trial = 1:size(csMinus)
    auc_csMinus(trial) = trapz(csMinus(trial,csOnset:usOnset));
end