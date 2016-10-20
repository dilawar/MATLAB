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
    saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
    
    [auc_csPlus, auc_csPlusProbe, auc_csMinus] = areaUnderCurve(saveFolder, dataset, baselineCorrection, preTime, csTime, traceTime, usTime, samplingRate);
    
    % Blink Threshold (mean+2*stddev)
    blinkThreshold_csPlus = mean(auc_csPlus)+2*(std(auc_csPlus));
    blinkThreshold_csMinus = mean(auc_csMinus)+2*(std(auc_csMinus));
    
    for session = 1:nSessions
        %Estimate blinks in TEC sessions
        dataset = [mouse '_SessionType' num2str(TECSessionType) '_Session' num2str(session)];
        saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
        [auc_csPlus, auc_csPlusProbe, auc_csMinus] = areaUnderCurve(saveFolder, dataset, baselineCorrection, preTime, csTime, traceTime, usTime, samplingRate);
        
        for trial = 1:length(auc_csPlus)
            if auc_csPlus(trial) > blinkThreshold_csPlus
                blinks_csPlus(trial) = 1;
            else
                blinks_csPlus(trial) = 0;
            end
        end
        
        for trial = 1:length(auc_csPlusProbe)
            if auc_csPlusProbe(trial) > blinkThreshold_csPlus
                blinks_csPlusProbe(trial) = 1;
            else
                blinks_csPlusProbe(trial) = 0;
            end
        end
        
        for trial = 1:length(auc_csMinus)
            if auc_csMinus(trial) > blinkThreshold_csMinus
                blinks_csMinus(trial) = 1;
            else
                blinks_csMinus(trial) = 0;
            end
        end
        
        save([saveFolder 'Mouse' dataset '_blinks_csPlus'], 'blinks_csPlus');
        save([saveFolder 'Mouse' dataset '_blinks_csPlusProbe'], 'blinks_csPlusProbe');
        save([saveFolder 'Mouse' dataset '_blinks_csMinus'], 'blinks_csMinus');
        
        %Performance
        perf_blinks_csPlus(animal,session) = length(find(blinks_csPlus))*100/length(blinks_csPlus);
        perf_blinks_csPlusProbe(animal,session) = length(find(blinks_csPlusProbe))*100/length(blinks_csPlusProbe);
        perf_blinks_csMinus(animal,session) = length(find(blinks_csMinus))*100/length(blinks_csMinus);
        
    end
end

figure(1);
subplot(1,3,1);
plot(perf_blinks_csPlus',...
    'blue',...
    'LineWidth', 1,...
    'MarkerSize', 10);
hold on
%errorbar(mean(perf_blinks_csPlus,1)',std(perf_blinks_csPlus,1)',...
%     'blue',...
%     'LineWidth', 3,...
%     'MarkerSize', 10);
title('CS+',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Hit Trials/Total CS+ Trials (%)',...
    'FontSize',fontSize,...
    'FontWeight','bold')
axis([1 length(perf_blinks_csPlus) 0 100])

subplot(1,3,2);
plot(perf_blinks_csPlusProbe',...
    'green',...
    'LineWidth', 1,...
    'MarkerSize', 10);
hold on
%errorbar(mean(perf_blinks_csPlusProbe,1)',std(perf_blinks_csPlusProbe,1)',...
%     'green',...
%     'LineWidth', 3,...
%     'MarkerSize', 10);
title('Probe',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Hit Trials/Total Probe Trials (%)',...
    'FontSize',fontSize,...
    'FontWeight','bold')
axis([1 length(perf_blinks_csPlusProbe) 0 100])

subplot(1,3,3);
plot(perf_blinks_csMinus',...
    'red',...
    'LineWidth', 1,...
    'MarkerSize', 10);
hold on;
%errorbar(mean(perf_blinks_csMinus,1)',std(perf_blinks_csMinus,1)',...
%     'red',...
%     'LineWidth', 3,...
%     'MarkerSize', 10);
title('CS-',...
    'FontSize',fontSize,...
    'FontWeight','bold')
xlabel('Sessions',...
    'FontSize',fontSize,...
    'FontWeight','bold')
ylabel('Hit Trials/Total CS- Trials (%)',...
    'FontSize',fontSize,...
    'FontWeight','bold')
axis([1 length(perf_blinks_csMinus) 0 100])

print(['/Users/ananth/Desktop/Figures/HitScore_' myString],'-djpeg');
beep;
disp('Done!');
%close all