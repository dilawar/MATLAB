% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Performance analysis
% DEPENDENCIES - 1) Sort all trials as .mat files using sortingVideos.m
%                2) Find the image processing parameters using findTheEye.m
%                3) Perform FEC analysis


%%%%%%% to be edited!!! %%%%%%% - make sure saved file names are consistent

tic
clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions/')

% Operations (0 == Don't Perform; 1 == Perform)
saveData = 0;
scorePerformance = 1;
plotFigures = 1;

% Dataset details
sessionType = 11;
mice = [1 2 3 4 5];
%mice = 2;
%nSessions = 4;
if sessionType == 9
    nSessions = 12;
else
    nSessions = 3;
end

allScores = nan(length(mice), nSessions);

startSession = 1;
%startSession = nSessions;

fecDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
saveDirec = '/Users/ananth/Desktop/Work/Analysis/PerformanceAnalysis/';

% Protocol details
preCSTime = 0.5; % in seconds
csTime = 0.05; % in seconds
usTime = 0.05; % in seconds

alpha = 0.01; % Significance level for kstest2
trialRejectThreshold = 0.05; % Fano's Factor based rejection

learningCutoff = 50; % in percent
disqualificationCutoff = 50; % in percent

allDisqualified = nan(length(mice),nSessions);

fontSize = 20;
lineWidth = 3;
markerSize = 8;

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    score = nan(nSessions,1);
    disqualifications = nan(nSessions,1);
    
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
            fanoFactor = zeros(nTrials,2);
            hitTrials = [];
            rejectedTrials = [];
            
            doi = samplingRate*traceTime; % length of the duration of interest
            csON = samplingRate*preCSTime;
            csOFF = csON + samplingRate*csTime;
            
            for trial = 1:nTrials
                x1 = fec(trial,1:(csON-1)); % pre-CS period
                x2 = fec(trial,(csOFF+1):(csOFF+doi)); % trace period
                
                % Disqualifications (if any) based on Fano's factor
                fanoFactor(trial,1) = var(x1)/mean(x1);
                fanoFactor(trial,2) = var(x2)/mean(x2);
                
                if fanoFactor(trial,1) > trialRejectThreshold
                    reject = 1;
                    rejectedTrials = [rejectedTrials trial];
                    disp(['Trial ' num2str(trial) ' rejected'])
                    hit = 0;
                else
                    reject = 0;
                    hit = kstest2(x1, x2, 'Alpha', alpha);
                    if hit == 1
                        hitTrials = [hitTrials trial];
                    end
                end
                
                nHits = length(hitTrials);
                nRejects = length(rejectedTrials);
                
                if saveData == 1
                    saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
                    if ~isdir(saveFolder)
                        mkdir(saveFolder);
                    end
                    
                    % trialInfo
                    save([saveFolder 'Trial' num2str(trial) '.mat' ], ...
                        'fanoFactor', 'reject', ...
                        'hit')
                end
            end
            
            percentDisqualified = (nRejects/nTrials)*100;
            if percentDisqualified >= disqualificationCutoff
                sessionScore = nan;
            else
                sessionScore = (nHits/(nTrials-nRejects))*100; % in percentage
            end
        end
        if saveData == 1
            % session info
            save([saveFolder 'performance' num2str(session) '.mat'], ...
                'sessionScore', ...
                'alpha', 'hitTrials', ...
                'fanoFactor', 'trialRejectThreshold', 'rejectedTrials')
        end
        allScores(mouse,session) = sessionScore;
        allDisqualified(mouse,session) = percentDisqualified;
    end
    
    if saveData == 1
        
        % all scores
        save([saveDirec 'allScores.mat'], ...
            'allScores', ...
            'alpha', ...
            'trialRejectThreshold', 'allDisqualified')
    end
    
    if plotFigures == 1
        learningLine = ones(nSessions,1);
        disqualifyLine = ones(nSessions,1);
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
        plot(allScores(mouse,:),'-*',...
            'LineWidth',lineWidth,...
            'MarkerSize',markerSize)
        xlabel('Sessions', ...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        axis([1 nSessions -5 100]);
        set(gca,'FontSize', fontSize)
        if mouse == length(mice)
            hold on
            plot(learningLine*learningCutoff,'--black')
            %legend('M1', 'M2', 'M3', 'M4', 'M5', 'Learnt') % Later, make this a cell array
        end
        
        print('/Users/ananth/Desktop/figs/scores',...
            '-dpng');
        
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
        plot(allDisqualified(mouse,:), '-*', ...
            'LineWidth', lineWidth, ...
            'MarkerSize', markerSize)
        xlabel('Sessions', ...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        set(gca,'FontSize', fontSize)
        axis([1 nSessions 0 100]);
        set(gca,'YTick', [0 25 50 75 100])
        set(gca,'YTickLabel',[0 25 50 75 100])
        if mouse == length(mice)
            hold on
            plot(disqualifyLine*disqualificationCutoff,'--black')
            %legend('M1', 'M2', 'M3', 'M4', 'M5', 'Disqualified') % Later, make this a cell array
        end
        
        print('/Users/ananth/Desktop/figs/disqualifications',...
            '-dpng');
    end
end
toc
disp('All done!')
