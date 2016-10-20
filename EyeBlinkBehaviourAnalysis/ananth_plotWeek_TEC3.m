% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Use this to plot the week's behaviour data

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions_MATLAB')

clear all
close all

saveFigures = 1;
baselineCorrection = 0;

samplingRate = 100; %Hz
samplePoint = samplingRate/1000;
preTime = 5000; %ms
csTime = 350; %ms
traceTime = 250; %ms
usTime = 50; %ms

week = 2;
mice = [39, 40, 41];
nAnimals = length(mice);

controlSessionType = 3;
TECSessionType = 9;
nSessions = 2;
startSession = 1;

% Specify if probe trials exist
if TECSessionType <= 5
    probe = 0;
else
    probe = 1;
    %probe = 0;
end

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');

fontSize = 12;

%Time Points
csOn = (preTime)*samplePoint+1;
csOff = (preTime+csTime)*samplePoint;
traceOn = (preTime+csTime)*samplePoint+1;
traceOff = (preTime+csTime+traceTime)*samplePoint;
usOn = (preTime+csTime+traceTime)*samplePoint+1;
usOff = (preTime+csTime+traceTime+usTime)*samplePoint;
criticalWindow = traceOn:(traceOff-1);

count2 = 1;

for animal = 1:nAnimals
    mouse = ['K' num2str(mice(animal))];
    
    %CONTROL NO-PUFF SESSION
    %Use No-Puff Controls to establish blink threshold
    dataset = [mouse '_SessionType' num2str(controlSessionType) '_Session1'];
    saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
    
    if baselineCorrection == 1
        csPlusTrials = csvread([saveFolder 'Mouse' dataset '_csPlus_dRR.csv']);
        csMinusTrials = csvread([saveFolder 'Mouse' dataset '_csMinus_dRR.csv']);
    else
        csPlusTrials = csvread([saveFolder 'Mouse' dataset '_csPlus.csv']);
        csMinusTrials = csvread([saveFolder 'Mouse' dataset '_csMinus.csv']);
    end
    csPlusProbeTrials = [];
    
    [auc_csPlus, auc_csPlusProbe, auc_csMinus] = areaUnderCurve(csPlusTrials, csPlusProbeTrials, csMinusTrials, baselineCorrection, criticalWindow);
    
    % Blink Threshold (mean+2*stddev)
    blinkThreshold_csPlus = mean(auc_csPlus)+2*(std(auc_csPlus));
    blinkThreshold_csMinus = mean(auc_csMinus)+2*(std(auc_csMinus));
    
    count = 1;
    %TEC SESSIONS
    %Estimate blinks in TEC sessions
    for session = startSession:nSessions
        dataset = [mouse '_SessionType' num2str(TECSessionType) '_Session' num2str(session)];
        disp(dataset);
        saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
        
        csPlusTrials = csvread([saveFolder 'Mouse' dataset '_csPlus.csv']);
        csMinusTrials = csvread([saveFolder 'Mouse' dataset '_csMinus.csv']);
        if probe == 1
            probeTrials = csvread([saveFolder 'Mouse' dataset '_probeTrials.csv']);
            csPlusProbeTrials = csPlusTrials(probeTrials,:);
        else
            csPlusProbeTrials = [];
        end
        
        [auc_csPlus, auc_csPlusProbe, auc_csMinus] = areaUnderCurve(csPlusTrials, csPlusProbeTrials, csMinusTrials, baselineCorrection, criticalWindow);
        
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
        
        figure(1);
        subplot(1,nSessions,session);
        mesh(csPlusTrials(find(blinks_csPlus),400:600));
        title([mouse ' ST' num2str(TECSessionType) ' S' num2str(session) ' Blue LED'],...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xlabel('Time/10 ms',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        %set(gca, 'FontSize', fontSize);
        ylabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        yt = get(gca, 'YTick');
        %set(gca, 'FontSize', fontSize);
        %caxis([0, 1024]);
        %axis([400 600]);
        colorbar;
        colormap(jet);
        count = count+1;
        
        
        %         figure(2);
        %         subplot(nAnimals,1,count2)
        %         if probe == 1
        %             plot(median(blinkData_probe,1), 'g')
        %             hold on
        %         end
        %         plot(median(blinkData_csPlus,1),'b');
        %
        %         title([mouse ' ST' num2str(TECSessionType) ' S' num2str(session)],...
        %             'FontSize',fontSize,...
        %             'FontWeight','bold');
        %         axis([0 1000 0 1024]);
        %         xlabel('Time/10 ms',...
        %             'FontSize',fontSize,...
        %             'FontWeight','bold');
        %         xt = get(gca, 'XTick');
        %         %set(gca, 'FontSize', fontSize);
        %         ylabel('Response',...
        %             'FontSize',fontSize,...
        %             'FontWeight','bold');
        %         yt = get(gca, 'YTick');
        %         %set(gca, 'FontSize', fontSize);
        %         if probe == 1
        %             legend('Median Probe Trials', 'Median Blue LED Trials');
        %         else
        %             legend('Median Blue LED Trials');
        %         end
        %         %         figure(2);
        %         %         print(['/Users/ananth/Desktop/' 'median_set'...
        %         %             num2str(set) '_SessionType_' num2str(TECSessionType)...
        %         %             '_Session' num2str(session)],...
        %         %             '-djpeg');
        %         count2 = count2+1;
    end
    if saveFigures ==1
        print(['/Users/ananth/Desktop/rawBlinks_MouseK' num2str(mouse) '_SessionType' num2str(TECSessionType) ...
            '_week' num2str(week) '_' num2str(session)],...
            '-djpeg');
    end
end

% if saveFigures == 1
%     figure(1);
%     print(['/Users/ananth/Desktop/rawBlinks_SessionType' num2str(TECSessionType) ...
%         '_week' num2str(week) '_' num2str(session)],...
%         '-djpeg');
%     %     figure(2);
%     %     print(['/Users/ananth/Desktop/median_SessionType' num2str(TECSessionType) ...
%     %         '_week' num2str(week) '_' num2str(session)],...
%     %         '-djpeg');
% end

%pause(5);
beep;
disp('Done!');