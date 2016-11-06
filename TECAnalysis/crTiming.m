% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - CR Timing - between animals
% DEPENDENCIES - 1) Sort all trials as .mat files using sortingVideos.m
%                2) Find the image processing parameters using findTheEye.m
%                3) Perform FEC analysis
%                4) Perform Motion analysis (to get the probe Trials)

tic
clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions/')

% Operations (0 == Don't Perform; 1 == Perform)
saveData = 0;
estimateTiming = 1;
findCROn = 1;
findCRPeak = 0;
plotFigures = 0;

% Dataset details
%sessionType = 9;
%mice = [1 2 3 4 5];
mice = [2 5];
%mice = 2;
%nSessions = 12;

%startSession = 1;
%startSession = nSessions;

fecDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
motionDirec = '/Users/ananth/Desktop/Work/Analysis/MotionAnalysis/';
scoreDirec = '/Users/ananth/Desktop/Work/Analysis/PerformanceAnalysis/';
saveDirec = scoreDirec;

% Protocol details
preCSTime = 0.5; % in seconds
csTime = 0.05; % in seconds
usTime = 0.05; % in seconds

fontSize = 16;
lineWidth = 3;
markerWidth = 7;
transparency = 0.5;
colours = {'b', 'r'};

timeLine = zeros(150,1);
csLine = timeLine;
usLine = timeLine;
%
% csLine(50:55,1) = 1;
% if sessionType == 9
%     usLine(80:85) = 1;
% else
%     usLine(90:95) = 1;
% end

% % Finding CR Onset
% win = 5;
% traceOn = 56; %frame
% if sessionType == 9
%     traceOff = 79;
% else
%     traceOff = 89;
% end

allCROns = nan(length(mice),4,60);
animal = 2;
for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    for count = 1:4
        %for session = startSession:nSessions
        if count == 1;
            sessionType = 9;
            session = 12;
        else
            sessionType = 11;
            session = count-1;
        end
        csLine(50:55,1) = 1;
        if sessionType == 9
            usLine(80:85) = 1;
        else
            usLine(90:95) = 1;
        end
        
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
            
            stdDev(mouse,:) = std(crTrials,1);
            meanCR(mouse,:) = mean(crTrials,1);
            
            disp('... done!')
        end
        
        if findCROn == 1
            %Find CR Onset
            disp('Finding CR onsets...')
            % Finding CR Onset
            win = 5;
            traceOn = 56; %frame
            if sessionType == 9
                traceOff = 79;
            else
                traceOff = 89;
            end
            crOns = nan(size(crTrials,1),1);
            for trial = 1:length(crOns)
                h = 0;
                disp(['Trial ' num2str(trial)])
                for i = traceOn:traceOff
                    x1 = crTrials(trial,1:49);
                    x2 = crTrials(trial,i:i+win-1);
                    h = kstest2(x1,x2,'Alpha', 0.01);
                    if h == 1
                        crOns(trial) = i+ceil(win/2);
                        disp('Found a CR!')
                        break
                    end
                end
            end
            crOns(find(isnan(crOns))) = [];
            allCROns(mouse,count,1:length(crOns)) = crOns;
            %             save(['/Users/ananth/Desktop/crs/crOns_' ...
            %                 mouseName '_ST' num2str(sessionType) ...
            %                 '_S' num2str(session) '.mat'], ...
            %                 'crOns')
            disp('... done!')
        end
        
        if findCRPeak == 1
            %Find CR Peak
            disp('Finding CR peaks ...')
            crPeaks = nan(size(crTrials,1),1);
            for trial = 1:length(crPeaks)
                [val idx] = max(crTrials(trial,1:traceOff));
                crPeaks(trial) = idx;
            end
            %             save(['/Users/ananth/Desktop/crs/crPeaks_' ...
            %                 mouseName '_ST' num2str(sessionType) ...
            %                 '_S' num2str(session) '.mat'], ...
            %                 'crPeaks')
            disp('... done!')
        end
        %count = count+1;
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
group1 = [repmat({'250ms-S12'}, length(find(allCROns(animal,1,:))), 1); ...
    repmat({'350ms-S1'}, length(find(allCROns(animal,2,:))), 1); ...
    repmat({'350ms-S2'}, length(find(allCROns(animal,3,:))), 1); ...
    repmat({'350ms-S3'}, length(find(allCROns(animal,4,:))),1)];

figure(2)
z = boxplot([squeeze(allCROns(animal,1,:));...
    squeeze(allCROns(animal,2,:));...
    squeeze(allCROns(animal,3,:));...
    squeeze(allCROns(animal,4,:))],...
    group1);
set(z,'LineWidth',lineWidth)
set(gca,'FontSize',fontSize)
ylim([55 90])
if animal == 1;
    title('Strong Learner - M2',...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
else
    title('Weak Learner - M5',...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
end
ylabel('Time/ms (from CS onset)', ...
    'FontSize', fontSize, ...
    'FontWeight', 'bold')
set(gca,'YTick', [60 65 70 75 80 85 90])
set(gca,'YTickLabel', [100 150 200 250 300 350 400])

a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects
outliers = a(1:4);
medians = a(5:8);
boxes = a(9:12);
whiskers = a(21:28);
if animal == 1
    set(boxes, 'Color', 'red');
    set(whiskers, 'Color', 'red');
    set(medians, 'Color', 'black');
else
    set(boxes, 'Color', 'blue');
    set(whiskers, 'Color', 'blue');
    set(medians, 'Color', 'black');
end

% figure(3);
% boxplot([squeeze(allCROns(2,1,:));...
%     squeeze(allCROns(2,2,:));...
%     squeeze(allCROns(2,3,:));...
%     squeeze(allCROns(2,4,:))],...
%     group2)
% title('Weak Learner - M5',...
%     'FontSize', fontSize, ...
%     'FontWeight', 'bold')
% set(gca,'FontSize',fontSize)
% ylabel('Time from CS onset (ms)', ...
%     'FontSize', fontSize, ...
%     'FontWeight', 'bold')
% set(gca,'YTick', [60 65 70 75 80 85 90])
% set(gca,'YTickLabel', [100 150 200 250 300 350 400])
%
% b = get(get(gca,'children'),'children');   % Get the handles of all the objects
% t = get(b,'tag');   % List the names of all the objects
% box2 = b(9:12);
% set(box2, 'Color', 'blue');


if plotFigures == 1
    lineProps1.col{1} = 'red';
    lineProps2.col{1} = 'blue';
    %     lineProps3.col{1} = 'green';
    %     lineProps4.col{1} = 'black';
    
    figure(1)
    subplot(6,1,1:5)
    hold on
    %             shadedErrorBar([],meanCR(mouse,:),stdDev(mouse,:), ...
    %                 {'LineWidth', lineWidth});
    mseb([],meanCR(1,:),stdDev(1,:), ...
        lineProps1, transparency);
    mseb([],meanCR(2,:),stdDev(mouse,:), ...
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
    set(gca,'FontSize', fontSize)
    set(gca,'XTickLabel',[])
    ylabel('FEC', ...
        'FontSize', fontSize,...
        'FontWeight', 'bold')
    set(gca,'YTick', [0 0.5 1])
    set(gca,'YTickLabel', [0 0.5 1])
    legend('M2','M5')
    
    subplot(6,1,6)
    plot(csLine,'-','LineWidth',lineWidth)
    hold on
    plot(usLine,'-','LineWidth',lineWidth)
    set(gca)
    legend('CS','US')
    set(gca,'FontSize', fontSize)
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
    
    print(['/Users/ananth/Desktop/figs/crTiming' ...
        '_ST' num2str(sessionType) ...
        '_S' num2str(session)], ...
        '-djpeg');
end
toc
disp('All done!')