% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Read saved raw data and analyse for TEC Performance

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions_MATLAB')

clear all
close all

saveFigures = 1;
baselineCorrection = 0;

samplingRate = 100; %Hz
preTime = 5000; %ms
csTime = 350; %ms
usTime = 50; %ms
samplePoint = samplingRate/1000;

%Time points
csOnset = preTime*samplePoint;
csOffset = (preTime+csTime)*samplePoint;
%usOnset = (preTime+csTime+traceTime)*samplePoint;
%usOffset = usOnset+(usTime*samplePoint);
criticalWindow = (csOnset:csOffset);
spontaneousWindow = (33:csOnset);

myString = 'Learners_250';

if strcmp(myString,'Learners_250')
    controlSessionType = 3;
    TECSessionType = 9;
    nSessions = 4;
    mice = [26, 27]; %learners 250 ms
elseif strcmp(myString,'NonLearners_250')
    controlSessionType = 3;
    TECSessionType = 9;
    nSessions = 4;
    mice = [28,32,33,34]; %non-learners 250 ms
elseif strcmp(myString,'NonLearners_500')
    controlSessionType = 5;
    TECSessionType = 11;
    nSessions = 7;
    mice = [29, 30, 31]; %non-learners 500 ms
end

nAnimals = length(mice);
nSessions = 4;
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

fontSize = 10;

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');
fontSize = 10;
for animal = 1:nAnimals
    mouse = ['K' num2str(mice(animal))];
    
    %Use No-Puff Controls to establish blink threshold
    dataset = [mouse '_SessionType' num2str(controlSessionType) '_Session1'];
    uscorei = strfind(dataset,'_');
    sessionType = dataset((uscorei(1)+12):(uscorei(2)-1));
    saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
    
    [auc_csPlus, auc_csPlusProbe, auc_csMinus] = areaUnderCurve(saveFolder, dataset, baselineCorrection, preTime, csTime, traceTime, usTime, samplingRate);
    
    % Blink Threshold (mean+2*stddev)
    blinkThreshold_csPlus = mean(auc_csPlus)+2*(std(auc_csPlus));
    blinkThreshold_csMinus = mean(auc_csMinus)+2*(std(auc_csMinus));
    
    %No-Puff Control Blink Probability
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
    
    %Rectify
    rectifiedcsPlusTrials = abs(csPlusTrials);
    rectifiedcsPlusProbeTrials = abs(csPlusProbeTrials);
    rectifiedcsMinusTrials = abs(csMinusTrials);
    
    %Preallocate
    csPlus = zeros(size(csPlusTrials));
    csPlusProbe = zeros(size(csPlusProbeTrials));
    csMinus = zeros(size(csMinusTrials));
    %Set baseline to zero
    for trial = 1:size(rectifiedcsPlusTrials,1)
        csPlus(trial,:) = rectifiedcsPlusTrials(trial,:)-median(rectifiedcsPlusTrials(trial,:));
    end
    for trial = 1:size(rectifiedcsPlusProbeTrials,1)
        csPlusProbe(trial,:) = rectifiedcsPlusProbeTrials(trial,:)-median(rectifiedcsPlusProbeTrials(trial,:));
    end
    for trial = 1:size(rectifiedcsMinusTrials,1)
        csMinus(trial,:) = rectifiedcsMinusTrials(trial,:)-median(rectifiedcsMinusTrials(trial,:));
    end
    
    %Set all negative values to zero
    csPlus(csPlus<0) = 0;
    csPlusProbe(csPlusProbe<0) = 0;
    csMinus(csMinus<0) = 0;
    
    nbins = floor(length(spontaneousWindow)/length(criticalWindow));
    
    %Preallocate
    auc_spont_csPlus = zeros(size(csPlusTrials,1),nbins);
    auc_spont_csPlusProbe = zeros(size(csPlusProbeTrials,1),nbins);
    auc_spont_csMinus = zeros(size(csMinusTrials,1),nbins);
    
    %Reshape
    for trial = 1:size(csPlus)
        A = reshape(csPlus(trial,spontaneousWindow),[length(criticalWindow),nbins]);
        auc_spont_csPlus(trial,:) = trapz(A);
    end
    for trial = 1:size(csPlusProbe)
        B = reshape(csPlusProbe(trial,spontaneousWindow),[length(criticalWindow),nbins]);
        auc_spont_csPlusProbe(trial,:) = trapz(B);
    end
    for trial = 1:size(csMinus)
        C = reshape(csMinus(trial,spontaneousWindow),[length(criticalWindow),nbins]);
        auc_spont_csMinus(trial,:) = trapz(C);
    end
    
    for trial = 1:size(auc_spont_csPlus,1)
        for bin = 1:nbins
            if auc_spont_csPlus(trial,bin) > blinkThreshold_csPlus
                X(trial,bin) = 1;
            else
                X(trial,bin) = 0;
            end
        end
    end
    
    for trial = 1:size(auc_spont_csPlusProbe,1)
        for bin = 1:nbins
            if auc_spont_csPlusProbe(trial,bin) > blinkThreshold_csPlus
                Y(trial,bin) = 1;
            else
                Y(trial,bin) = 0;
            end
        end
    end
    
    for trial = 1:size(auc_spont_csMinus,1)
        for bin = 1:nbins
            if auc_spont_csMinus(trial,bin) > blinkThreshold_csMinus
                Z(trial,bin) = 1;
            else
                Z(trial,bin) = 0;
            end
        end
    end
    
    prob_csPlus_npc(animal) = numel(X(X==1))/numel(X);
    %prob_csPlusProbe_npc(animal) = numel(Y(Y==1))/numel(Y);
    prob_csMinus_npc(animal) = numel(Z(Z==1))/numel(Z);
    
    for session = 1:nSessions
        %Estimate blink probability in TEC sessions
        dataset = [mouse '_SessionType' num2str(TECSessionType) '_Session' num2str(session)];
        uscorei = strfind(dataset,'_');
        sessionType = dataset((uscorei(1)+12):(uscorei(2)-1));
        saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
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
        
  
% beep;
disp('Done!');
%close all