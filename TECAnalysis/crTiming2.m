% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - CR Timing - across sessions
% DEPENDENCIES - 1) Sort all trials as .mat files using sortingVideos.m
%                2) Find the image processing parameters using findTheEye.m
%                3) Perform FEC analysis
%                4) Perform Motion analysis (to get the probe Trials)

tic
clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions/')

% Operations (0 == Don't Perform; 1 == Perform)
saveData = 0;
estimateTiming = 1;
plotFigures = 1;

% Dataset details
sessionType = 11;
%mice = [1 2 3 4 5];
mice = 2;

nSessions = 3;

startSession = 1;
%startSession = nSessions;

fecDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
motionDirec = '/Users/ananth/Desktop/Work/Analysis/MotionAnalysis/';
scoreDirec = '/Users/ananth/Desktop/Work/Analysis/PerformanceAnalysis/';
saveDirec = scoreDirec;

% Protocol details
preCSTime = 0.5; % in seconds
csTime = 0.05; % in seconds
usTime = 0.05; % in seconds

fontSize = 12;
lineWidth = 3;
markerWidth = 7;
transparency = 0.5;
colours = {'b', 'r'};

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
        
        if estimateTiming == 1
            disp('Estimating Response Timing ...')
            % Load motion data
            try
                load([motionDirec 'Mouse' mouseName '/' ...
                    dataset '/motion.mat'])
            catch
                warning('Unable to motion data')
                continue
            end
            
            % Load FEC data
            try
                load([fecDirec 'Mouse' mouseName '/' ...
                    dataset '/fec.mat'])
            catch
                warning('Unable to find FEC data')
                continue
            end
            
            nTrials = size(fec,1);
            
            % Load scores
            crTrials = [];
            hitTrials = [];
            
            try
                load([scoreDirec 'Mouse' mouseName '/' ...
                    dataset '/Session' num2str(session) '.mat'])
            catch
                warning('Unable to find scores')
                continue
            end
            
            % Collect all hit Trials
            crTrials(:,:) = fec(hitTrials,:);
            
            stdDevCRAmp(mouse,:) = std(crTrials,1);
            meanCRAmp(mouse,:) = mean(crTrials,1);
            
            stdDevCRTime(mouse,:) = std(crTrials,1);
            meanCRTime(mouse,:) = mean(crTrials,1);
            
        end
        
        if saveData == 1
            saveFolder = [saveDirec 'Mouse' mouseName '/' ...
                dataset '/'];
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            
            save([saveFolder 'cr.mat'], ...
                'cr')
        end
    end
end
if plotFigures == 1
    lineProps1.col{1} = 'red';
    lineProps2.col{1} = 'blue';
    lineProps3.col{1} = 'green';
    %     lineProps4.col{1} = 'black';
    
    figure(1)
    subplot(6,1,1:5)
    hold on
    %             shadedErrorBar([],meanCR(mouse,:),stdDev(mouse,:), ...
    %                 {'LineWidth', lineWidth});
    mseb([],meanCR(1,:),stdDev(1,:), ...
        lineProps1, transparency);
    mseb([],meanCR(1,:),stdDev(mouse,:), ...
        lineProps2, transparency);
    axis([0 150 0 1.2])
    if sessionType == 9
    title(['CR Timing - 250 ms ISI - Session ' num2str(session)], ...
        'FontSize', fontSize,...
        'FontWeight', 'bold')
    else
        title(['CR Timing - 350 ms ISI - Session ' num2str(session)], ...
        'FontSize', fontSize,...
        'FontWeight', 'bold')
    end
    set(gca,'XTick', [10 20 30 ...
        40 50 60 ...
        70 80 90 ...
        100 110 120 ...
        130 140 150])
    set(gca,'XTickLabel',[]) 
    ylabel('FEC', ...
        'FontSize', fontSize,...
        'FontWeight', 'bold')
    set(gca,'YTick', [0 0.5 1])
    set(gca,'YTickLabel', [0 0.5 1])
    legend('mean M2 (errorbar: SD) Session 1','mean M5 (errorbar: SD) Session 2', 'mean M5 (errorbar: SD) Session 3')
    
    subplot(6,1,6)
    plot(csLine,'-','LineWidth',lineWidth)
    hold on
    plot(usLine,'-','LineWidth',lineWidth)
    set(gca)
    legend('CS','US')
    set(gca,'YTick',[0 1])
    set(gca, 'YTickLabel', {'Off' 'On'})
    xlabel('Time/ms', ...
        'FontSize', fontSize,...
        'FontWeight', 'bold')
    set(gca,'XTick', [10 20 30 ...
        40 50 60 ...
        70 80 90 ...
        100 110 120 ...
        130 140 150])
    set(gca,'XTickLabel', [100 200 300 ...
        400 500 600 ...
        700 800 900 ...
        1000 1100 1200 ...
        1300 1400 1500])
    
    print(['/Users/ananth/Desktop/figs/crTiming' ...
        '_ST' num2str(sessionType) ...
        '_S' num2str(session)], ...
        '-djpeg');
    
end
toc
disp('All done!')