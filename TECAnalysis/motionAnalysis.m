% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Read raw .csv files for the treadmill motion data
% EXPERIMENT SUMMARY - The position of the treadmill is read out every 10 ms
%                      by point detecting IR reflectance off of a circular
%                      pattern of black and white squares (0.5 cm side)
%                      along the periphery.

tic
clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions/')

%Please work with the data from all sessions of an animal, before proceeding to the next animal

clear all
%close all

%Operations (0 == Don't Perform; 1 == Perform)
saveData = 0;
doMotionAnalysis = 0;
plotFigures = 1;

rawDirec = '/Users/ananth/Desktop/Work/Behaviour/Motion/';
saveDirec = '/Users/ananth/Desktop/Work/Analysis/MotionAnalysis/';

%Dataset details
sessionType = 9;
%mice = [2 3 5];
mice = 2;
nSessions = 12;
nTrials = 61; % NOTE: The first trial is a dummy

startSession = nSessions;
%startSession = 1;
startTrial = 1;

nHeaderLines = 5;

dataSize = zeros(nTrials,1);

samplingRate = 100; % Hz
trialDuration = 1.5; % seconds
nSamples = floor(samplingRate*trialDuration);

distanceLC = 0.5; %cm
timeLC = 0.05; % seconds
threshold = 700;

win4avg = samplingRate*timeLC; % samples

fontSize = 20;
lineWidth = 3;
markerWidth = 7;

timeLine = zeros(150,1);
csLine = timeLine;
usLine = timeLine;

csLine(50:55,1) = 1;
if sessionType == 9
    usLine(80:85) = 1;
else
    usLine(90:95) = 1;
end

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if doMotionAnalysis == 1
            disp('Performing motion analysis ...')
            
            % Preallocation
            rawSession = zeros(nTrials, nSamples);
            motion = zeros(nTrials, (nSamples/win4avg));
            rawTrialLowEdge = zeros(nSamples,1);
            probeTrials = [];
            
            for trial = startTrial:nTrials
                disp(['Trial ' num2str(trial)])
                try
                    fileName = [rawDirec 'Mouse' mouseName '/' dataset '/Trial' num2str(trial) '.csv'];
                    rawData = csvread(fileName, nHeaderLines, 0);
                catch
                    disp(['[ERROR] ' dataset ' Trial ' num2str(trial) ' not found!'])
                    break
                end
                
                rawSession(trial,:) = rawData(1:nSamples,1); % has only motion data
                rawTrial = rawSession(trial,:)>threshold; % binarize
                %rawTrialDiff = diff(rawTrial); % differential (NOTE: n-1 element output)
                %rawTrialLowEdge(find(rawTrialDiff == -1) + 1) = 1;
                
                % Reshape to get the averaging window
                %rawTrialLowEdgeReshaped = reshape(rawTrialLowEdge, [win4avg nSamples/win4avg]);
                %motion(trial,:) = sum(rawTrialLowEdgeReshaped,1)*(distanceLC/timeLC);
                rawTrialReshaped = reshape(rawTrial, [win4avg nSamples/win4avg]);
                motion(trial,:) = mean(rawTrialReshaped,1)*(distanceLC/timeLC);
                
                % Find probe trials
                try
                    if csvread(fileName, nHeaderLines+1, 3, [ nHeaderLines+1 3 nHeaderLines+1 3] ) == 0   %reads only the row after nHeaderLines, column 4
                        probeTrials = [probeTrials trial];
                        %disp('Found a probe Trial!')
                        %disp(probeTrials)
                    else
                    end
                catch
                    warning('Unable to determine ProbeTrials');
                end
                disp('... done!')
            end
            
            % Get rid of the dummy trial
            motion(1,:) = [];
            rawSession(1,:) = [];
            probeTrials = probeTrials-1;
        else
            load([saveDirec 'Mouse' mouseName '/' dataset '/motion.mat']);
        end
        
        if plotFigures == 1
            figure(1)
            clf
            subplot(9,9,1:72)
            imagesc(motion)
            colormap(hot)
            if sessionType == 9
                title(['Treadmill Running - 250 ms ISI - ' ...
                    mouseName ' S' num2str(session)],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
            else
                title(['Treadmill Running - 350 ms ISI - ' ...
                    mouseName ' S' num2str(session)],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
            end
            %             title('Weak Learner - Treadmill Running - Session 12',...
            %                 'FontSize', fontSize,...
            %                 'FontWeight','bold')
            %             xlabel(['Time/' num2str(timeLC*1000) ' ms'], ...
            %                 'FontSize', fontSize,...
            %                 'FontWeight', 'bold')
            set(gca,'XTick',[])
            set(gca,'XTickLabel',[])
            ylabel('Trials', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            z = colorbar;
            ylabel(z,'Speed*s/cm',...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            set(gca,'FontSize', fontSize)
            
            subplot(9,9,73:80)
            plot(csLine,'-','LineWidth',lineWidth)
            hold on
            plot(usLine,'-','LineWidth',lineWidth)
            legend('CS','US')
            set(gca,'FontSize', fontSize-2)
            set(gca,'YTick',[0 1])
            set(gca, 'YTickLabel', {'Off' 'On'})
            xlabel('Time/ms', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            set(gca,'XTick', [10 30 50 ...
                70 90 110 ...
                130 150])
            set(gca,'XTickLabel', [100 300 500 ...
                700 900 1100 ...
                1300 1500])
            if sessionType == 9
                if mice(mouse) == 5
                    print('/Users/ananth/Desktop/figs/motion_weakLearner', ...
                        '-dpng');
                elseif mice(mouse) == 3
                    print('/Users/ananth/Desktop/figs/motion_nonLearner', ...
                        '-dpng');
                elseif mice(mouse) == 2
                    print('/Users/ananth/Desktop/figs/motion_strongLearner', ...
                        '-dpng');
                else
                    disp('[ERROR] Figures not saved!')
                end
            else
                print('/Users/ananth/Desktop/figs/motion_350', ...
                    '-dpng');
            end
            
            %             print(['/Users/ananth/Desktop/figs/motion_heatmap_' mouseName ...
            %                 '_ST' num2str(sessionType) ...
            %                 '_S' num2str(session)],...
            %                 '-djpeg');
            
            trialAvgMotion = mean(motion,2);
            trialAvgMotion_stddev = (std(motion'))'; %will fix this soon
            
            figure(2)
            clf
            subplot(9,1,1:8)
            plot(trialAvgMotion')
            shadedErrorBar([],trialAvgMotion,trialAvgMotion_stddev, ...
                {'red', 'LineWidth', lineWidth});
            ylabel('Average Speed *s/cm', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            xlabel('Trials', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            %colormap(hot)
            if sessionType == 9
                title(['Average Trial Speed - 250 ms ISI - ', ...
                    mouseName ' S' num2str(session)],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
            else
                title(['Average Trial Speed - 350 ms ISI - ', ...
                    mouseName ' S' num2str(session)],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
            end
            %axis([1 60 -1.5 10])
            set(gca,'XTick',[10 20 30 40 50 60])
            set(gca,'XTickLabel',[10 20 30 40 50 60])
            set(gca,'YTick',[0 5 10])
            set(gca,'YTickLabel',[0 5 10])
            %             z = colorbar;
            %             ylabel(z,'Speed (cm/s)',...
            %                 'FontSize', fontSize,...
            %                 'FontWeight', 'bold')
            %             set(gca,'FontSize', fontSize)
            view(90,-90) % flip plot
            set(gca, 'xdir', 'reverse'); %# Reverse the x-axis
            set(gca, 'FontSize', fontSize-2)
            if sessionType == 9
                if mice(mouse) == 5
                    print('/Users/ananth/Desktop/figs/avgMotion_weakLearner', ...
                        '-dpng');
                elseif mice(mouse) == 3
                    print('/Users/ananth/Desktop/figs/avgMotion_nonLearner', ...
                        '-dpng');
                elseif mice(mouse) == 2
                    print('/Users/ananth/Desktop/figs/avgMotion_strongLearner', ...
                        '-dpng');
                end
            else
                print('/Users/ananth/Desktop/figs/avgMotion_350', ...
                    '-dpng');
            end
            
            
            %             figure(3)
            %             imagesc(rawSession)
            %             colormap(hot)
            %             title(['Raw ' ...
            %                 mouseName ' ST' num2str(sessionType) ' S' num2str(session)],...
            %                 'FontSize', fontSize,...
            %                 'FontWeight', 'bold')
            %             xlabel('Samples', ...
            %                 'FontSize', fontSize,...
            %                 'FontWeight', 'bold')
            %             ylabel('Trials', ...
            %                 'FontSize', fontSize,...
            %                 'FontWeight', 'bold')
            %             z = colorbar;
            %             ylabel(z,'A.U.',...
            %                 'FontSize', fontSize,...
            %                 'FontWeight', 'bold')
            %             set(gca,'FontSize', fontSize)
            %
            %             print(['/Users/ananth/Desktop/figs/Motion/runningRaw_heatmap_' mouseName ...
            %                 '_ST' num2str(sessionType) ...
            %                 '_S' num2str(session)],...
            %                 '-djpeg');
        end
        
        if saveData == 1
            saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            
            % Save motion data
            save([saveFolder 'motion.mat' ], ...
                'rawSession', 'motion', 'probeTrials', ...
                'threshold', 'distanceLC', 'timeLC', ...
                'samplingRate', 'trialDuration')
        end
    end
    pause(0.5)
end
toc
beep
disp('All done!')