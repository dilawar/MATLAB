% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Read saved raw data and analyse for TEC Performance
% This code additionally requires areaUnderCurve.m to work.

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions_MATLAB')

clear all
close all

saveFigures = 1;
baselineCorrection = 0;

samplingRate = 100; %Hz
preTime = 5000; %ms
csTime = 350; %ms
usTime = 50; %ms

traceTime = 250;
controlSessionType = 3;
TECSessionType = 9;
nSessions = 8;
mice = [39, 40, 41];

nAnimals = length(mice);

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/Analysis4Anamika/');
fontSize = 12;

for animal = 1:nAnimals
    mouse = ['K' num2str(mice(animal))];
    
    %Use No-Puff Controls to establish blink threshold
    dataset = [mouse '_SessionType' num2str(controlSessionType) '_Session1'];
    saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
    
    [auc_csPlus, auc_csPlusProbe, auc_csMinus] = areaUnderCurve(saveFolder, dataset, baselineCorrection, preTime, csTime, traceTime, usTime, samplingRate);
    
    % Blink Threshold (mean+2*stddev)
    blinkThreshold_csPlus = mean(auc_csPlus)+2*(std(auc_csPlus));
    blinkThreshold_csMinus = mean(auc_csMinus)+2*(std(auc_csMinus));
    
    %Estimate blinks in TEC sessions
    for session = 1:nSessions
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
%subplot(1,3,1);
subplot(1,2,1);
plot(perf_blinks_csPlus',...
    'LineWidth', 1,...
    'MarkerSize', 10);
hold on
errorbar(mean(perf_blinks_csPlus,1)',std(perf_blinks_csPlus,1)',...
    'blue',...
    'LineWidth', 3,...
    'MarkerSize', 10);
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
legend('K39', 'K40', 'K41','Mean');
xt = get(gca, 'XTick');
set(gca, 'FontSize', fontSize);
yt = get(gca, 'YTick');
set(gca, 'FontSize', fontSize);

%subplot(1,3,2);
subplot(1,2,2);
plot(perf_blinks_csPlusProbe',...
    'LineWidth', 1,...
    'MarkerSize', 10);
hold on
errorbar(mean(perf_blinks_csPlusProbe,1)',std(perf_blinks_csPlusProbe,1)',...
    'green',...
    'LineWidth', 3,...
    'MarkerSize', 10);
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
legend('K39', 'K40', 'K41', 'Mean');
xt = get(gca, 'XTick');
set(gca, 'FontSize', fontSize);
yt = get(gca, 'YTick');
set(gca, 'FontSize', fontSize);

% subplot(1,3,3);
% plot(perf_blinks_csMinus',...
%     'LineWidth', 1,...
%     'MarkerSize', 10);
% hold on;
% %errorbar(mean(perf_blinks_csMinus,1)',std(perf_blinks_csMinus,1)',...
% %     'red',...
% %     'LineWidth', 3,...
% %     'MarkerSize', 10);
% title('CS-',...
%     'FontSize',fontSize,...
%     'FontWeight','bold')
% xlabel('Sessions',...
%     'FontSize',fontSize,...
%     'FontWeight','bold')
% ylabel('Hit Trials/Total CS- Trials (%)',...
%     'FontSize',fontSize,...
%     'FontWeight','bold')
% axis([1 length(perf_blinks_csMinus) 0 100])
% legend('K39', 'K40', 'K41', 'Mean');
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', fontSize);
% yt = get(gca, 'YTick');
% set(gca, 'FontSize', fontSize);

print('/Users/ananth/Desktop/anamika_HitScore','-djpeg');
beep;
disp('Done!');
%close all