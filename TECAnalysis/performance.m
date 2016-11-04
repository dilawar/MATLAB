% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Performance analysis
% DEPENDENCIES - 1) Sort all trials as .mat files using sortingVideos.m
%                2) Find the image processing parameters using findTheEye.m
%                3) Perform FEC analysis

tic
clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions/')

% Operations (0 == Don't Perform; 1 == Perform)
saveData = 1;
scorePerformance = 1;
plotFigures = 1;

% Dataset details
sessionType = 11;
mice = [1 2 3 4 5];
%mice = 5;
if sessionType == 9
    nSessions = 12;
else
    nSessions = 3;
end

score = nan(length(mice), nSessions);

startSession = 1;
%startSession = nSessions;

fecDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
motionDirec = '/Users/ananth/Desktop/Work/Analysis/MotionAnalysis/';
saveDirec = '/Users/ananth/Desktop/Work/Analysis/PerformanceAnalysis/';

% Protocol details
preCSTime = 0.5; % in seconds
csTime = 0.05; % in seconds
usTime = 0.05; % in seconds

alpha = 0.05; % Significance level for kstest2
trialRejectThreshold = 0.1; % Fano's Factor based rejection

learningCutoff = 50; % in percent

percentDisqualified = nan(length(mice),nSessions);

fontSize = 12;
lineWidth = 3;
markerSize = 8;

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if scorePerformance == 1
            disp('Scoring Performance ...')
            
            % Load FEC data
            try
                load([fecDirec 'Mouse' mouseName '/' dataset '/fec.mat'])
            catch
                warning('Unable to find FEC data')
                continue
            end
            
            % Performance - kstest2 between the pre-CS and trace intervals
            if sessionType == 9
                traceTime = 0.25; % in seconds
            elseif sessionType == 11
                traceTime = 0.35; % in seconds
            end
            
            nTrials = size(fec,1);
            hitTrials = zeros(nTrials,1);
            fanoFactor = zeros(nTrials,1);
            nRejects = 0;
            
            doi = samplingRate*traceTime; % length of the duration of interest
            csON = samplingRate*preCSTime;
            csOFF = csON + samplingRate*csTime;
            
            for trial = 1:nTrials
                x1 = fec(trial,1:(csON-1)); % pre-CS period
                x2 = fec(trial,(csOFF+1):(csOFF+doi)); % trace period
                
                % Disqualifications (if any) based on Fano's factor
                fanoFactor(trial) = var(x1)/mean(x1);
                if fanoFactor(trial) > trialRejectThreshold
                    hitTrials(trial,1) = 0;
                    hitTrials(trial,2) = nan;
                    nRejects = nRejects+1;
                    disp(['Trial ' num2str(trial) ' rejected'])
                else
                    hitTrials(trial) = kstest2(x1, x2, 'Alpha', alpha);
                end
                
                if saveData == 1
                    saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
                    if ~isdir(saveFolder)
                        mkdir(saveFolder);
                    end
                    
                    % trialInfo
                    save([saveFolder 'Trial' num2str(trial) '.mat' ], ...
                        'fanoFactor',...
                        'hitTrials',...
                        'nRejects')
                end
            end
            score(mouse,session) = (sum(hitTrials(:,1))/(nTrials-nRejects))*100; % in percentage
            percentDisqualified(mouse,session) = (nRejects/nTrials)*100;
        end
    end
    
    if saveData == 1
        %     saveFolder = saveDirec;
        %     if ~isdir(saveFolder)
        %         mkdir(saveFolder);
        %     end
        
        % Save FEC curve
        save([saveDirec 'sessionPerformance.mat'],...
            'alpha',...
            'trialRejectThreshold', 'percentDisqualified')
    end
    
    if plotFigures == 1
        learningLine = ones(nSessions,1);
        figure(1)
        if sessionType == 9
            subplot(1,5,1:4)
            title('Performance - 250 ms ISI', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            ylabel('Hit Trials/Total Trials (%)',...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            %legend('M1', 'M2', 'M3', 'M4', 'M5', 'Learning') % Later, make this a cell array
        else
            subplot(1,5,5)
            title('350 ms ISI', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            %             legend('M1', 'M2', 'M3', 'M4', 'M5', 'Learning') % Later, make this a cell array
            %             ylabel('Performance (%)',...
            %             'FontSize', fontSize,...
            %             'FontWeight', 'bold')
        end
        hold on
        plot(score(mouse,:),'-*',...
            'LineWidth',lineWidth,...
            'MarkerSize',markerSize)
        xlabel('Sessions', ...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        axis([1 nSessions 0 100]);
        if mouse == length(mice)
            hold on
            plot(learningLine*learningCutoff,'--black')
            legend('M1', 'M2', 'M3', 'M4', 'M5', 'Learnt') % Later, make this a cell array
        end
        
        print('/Users/ananth/Desktop/figs/performance',...
            '-djpeg');
        
        figure(2)
        if sessionType == 9
            subplot(1,5,1:4)
            title('Disqualifications - 250 ms ISI', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            ylabel('Disqualified Trials/Total Trials (%)',...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            %legend('M1', 'M2', 'M3', 'M4', 'M5', 'Learning') % Later, make this a cell array
        else
            subplot(1,5,5)
            title('350 ms ISI', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            %             legend('M1', 'M2', 'M3', 'M4', 'M5', 'Learning') % Later, make this a cell array
            %             ylabel('Performance (%)',...
            %             'FontSize', fontSize,...
            %             'FontWeight', 'bold')
        end
        hold on
        plot(percentDisqualified(mouse,:), '-*', ...
            'LineWidth', lineWidth, ...
            'MarkerSize', markerSize)
        xlabel('Sessions', ...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        axis([1 nSessions 0 10]);
        set(gca,'YTick', [0 5 10])
        set(gca,'YTickLabel',[0 5 10])
        legend('M1', 'M2', 'M3', 'M4', 'M5') % Later, make this a cell array
        
        print('/Users/ananth/Desktop/figs/disqualifiedTrials',...
            '-djpeg');
    end
end
toc
disp('All done!')
