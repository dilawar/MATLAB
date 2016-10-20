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

myString = 'NonLearners_500';

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
    
    [X, Y, Z] = crProbability(csPlusTrials, csPlusProbeTrials, csMinusTrials, sessionType, spontaneousWindow, criticalWindow, blinkThreshold_csPlus, blinkThreshold_csMinus);
    prob_csPlus_npc(animal) = numel(X(X==1))/numel(X);
    prob_csPlusProbe_npc(animal) = numel(Y(Y==1))/numel(Y);
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
        
        [X, Y, Z] = crProbability(csPlusTrials, csPlusProbeTrials, csMinusTrials, sessionType, spontaneousWindow, criticalWindow, blinkThreshold_csPlus, blinkThreshold_csMinus);
        prob_csPlus_tec(animal,session) = numel(X(X==1))/numel(X);
        prob_csPlusProbe_tec(animal,session) = numel(Y(Y==1))/numel(Y);
        prob_csMinus_tec(animal,session) = numel(Z(Z==1))/numel(Z);
    end
end

figure(2);
subplot(3,3,1);
plot(prob_csPlus_npc');
hold on
%errorbar(mean(perf_blinks_csPlus,1),std(perf_blinks_csPlus,1),'blue')
title('NPC CS+',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Animals',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Blink Probability ',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,2);
plot(prob_csPlusProbe_npc');
hold on
%errorbar(mean(perf_blinks_csPlusProbe,1),std(perf_blinks_csPlusProbe,1),'green')
title('NPC Probe',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Animals',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Blink Probability',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,3);
plot(prob_csMinus_npc');
hold on;
%errorbar(mean(perf_blinks_csMinus,1),std(perf_blinks_csMinus,1),'red')
title('NPC CS-',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Animals',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Blink Probability',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,4);
plot(prob_csPlus_tec');
hold on
%errorbar(mean(perf_blinks_csPlus,1),std(perf_blinks_csPlus,1),'blue')
title('TEC CS+',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Blink Probability ',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,5);
plot(prob_csPlusProbe_tec');
hold on
%errorbar(mean(perf_blinks_csPlusProbe,1),std(perf_blinks_csPlusProbe,1),'green')
title('TEC Probe',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Blink Probability',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,6);
plot(prob_csMinus_tec');
hold on;
%errorbar(mean(perf_blinks_csMinus,1),std(perf_blinks_csMinus,1),'red')
title('TEC CS-',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Blink Probability',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,7);
plot(((prob_csPlus_tec)*size(csPlusTrials,1))');
hold on
%errorbar(mean(perf_blinks_csPlus,1),std(perf_blinks_csPlus,1),'blue')
title('TEC CS+',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Expected Hit Trials/Total Trials (%) ',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,8);
plot(((prob_csPlusProbe_tec)*size(csPlusProbeTrials,1))');
hold on
%errorbar(mean(perf_blinks_csPlusProbe,1),std(perf_blinks_csPlusProbe,1),'green')
title('TEC Probe',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Expected Hit Trials/Total Trials (%) ',...
    'FontSize',fontSize,...
    'FontWeight','bold')

subplot(3,3,9);
plot(((prob_csMinus_tec)*size(csMinusTrials,1))');
hold on;
%errorbar(mean(perf_blinks_csMinus,1),std(perf_blinks_csMinus,1),'red')
title('TEC CS-',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Expected Hit Trials/Total Trials (%) ',...
    'FontSize',fontSize,...
    'FontWeight','bold')

print(['/Users/ananth/Desktop/Figures/BlinkProbability_' myString],'-djpeg');

disp('Done!');
%close all