% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Read raw .csv files for the treadmill motion data

tic
clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Please work with the data from all sessions of an animal, before proceeding to the next animal

clear all
%close all

%Operations (0 == Don't Perform; 1 == Perform)
saveData = 1;
doMotionAnalysis = 1;
plotFigures = 1;

saveDirec = '/Users/ananth/Desktop/Work/Analysis/Motion/';
direc = '/Users/ananth/Desktop/Work/Behaviour/Motion/';

%Dataset details
sessionType = 9;
mice = [1 2 3 4 5];
%mice = 2;
nSessions = 12;
nTrials = 61; % NOTE: The first trial is a dummy

startSession = 1;
startTrial = 1;

nHeaderLines = 5;

fontSize = 12;

dataSize = zeros(nTrials,1);

samplingRate = 100; % Hz
trialDuration = 1.5; % seconds
nSamples = floor(samplingRate*trialDuration);

distanceLC = 0.5; %cm
timeLC = 0.05; % seconds
threshold = 700;

win4avg = samplingRate*timeLC; % samples

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if doMotionAnalysis == 1
            disp('Performing motion analysis ...')
            
            % Preallocation
            raw = zeros(nTrials, nSamples);
            motion = zeros(nTrials, (nSamples/win4avg));
            rawTrialLowEdge = zeros(nSamples,1);
            
            for trial = startTrial:nTrials
                disp(['Trial ' num2str(trial-1)])
                try
                    fileName = [direc 'Mouse' mouseName '/' dataset '/Trial' num2str(trial) '.csv'];
                    rawData = csvread(fileName, nHeaderLines, 0);
                catch
                    disp(['[ERROR] ' dataset ' Trial ' num2str(trial) ' not found!'])
                    break
                end
                
                raw(trial,:) = rawData(1:nSamples,1); % has only motion values
                
                rawTrial = raw(trial,:)>threshold; % binarize
                rawTrialDiff = diff(rawTrial); % differential (NOTE: n-1 element output)
                rawTrialLowEdge(find(rawTrialDiff == -1) + 1) = 1;
                
                % Reshape to get the averaging window
                rawTrialLowEdgeReshaped = reshape(rawTrialLowEdge, [win4avg nSamples/win4avg]);
                motion(trial,:) = sum(rawTrialLowEdgeReshaped,1)*(distanceLC/timeLC);
                disp('... done!')
            end
            % Get rid of the dummy trial
            motion(1,:) = [];
            
            if plotFigures == 1
                figure(1)
                imagesc(motion)
                colormap(jet)
                title(['Treadmill Running ' ...
                    mouseName ' ST' num2str(sessionType) ' S' num2str(session) ...
                    ' (' num2str(samplingRate) ' fps)'],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                xlabel(['Time/' num2str(timeLC*1000) ' ms'], ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                ylabel('Trials', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                z = colorbar;
                ylabel(z,'Speed (cm/s)',...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                
                print(['/Users/ananth/Desktop/figs/motion_heatmap_' mouseName ...
                    '_ST' num2str(sessionType) ...
                    '_S' num2str(session)],...
                    '-djpeg');
            end
            
            if saveData == 1
                saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
                if ~isdir(saveFolder)
                    mkdir(saveFolder);
                end
                
                % Save motion data
                save([saveFolder 'motion.mat' ], ...
                    'raw', 'motion', ...
                    'samplingRate', 'trialDuration')
            end
        end
    end
end
toc
disp('All done!')