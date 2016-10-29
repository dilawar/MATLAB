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
saveData = 0;
doMotionAnalysis = 1;
plotFigures = 0;

saveDirec = '/Users/ananth/Desktop/Work/Analysis/Motion/';
direc = '/Users/ananth/Desktop/Work/Behaviour/Motion/';

%Dataset details
sessionType = 9;
%mice = [1 2 3 4 5];
mice = 2;
nSessions = 1;
nTrials = 2; % NOTE: The first trial is a dummy

startSession = 1;
startTrial = 1;

nHeaderLines = 5;

fontSize = 12;

dataSize = zeros(nTrials,1);

samplingRate = 100; % Hz
trialDuration = 1.5; % seconds
nSamples = floor(samplingRate*trialDuration);

distanceLC = 1; %cm
timeLC = 0.03; % seconds
threshold = 700;

win4avg = samplingRate*timeLC; % samples

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if doMotionAnalysis == 1
            disp('Performing motion analysis ...')
            
            raw = zeros(nTrials, nSamples);
            %rawSig = zeros(nTrials,nSamples);
            %threshold = zeros(nTrials,1);
            rawTrial = zeros(nSamples,1);
            motion = zeros(nTrials, (nSamples/win4avg));
            rate = zeros(nTrials, (nSamples/win4avg));
            
            for trial = startTrial:nTrials
                disp(['Trial ' num2str(trial-1)])
                try
                    fileName = [direc 'Mouse' mouseName '/' dataset '/Trial' num2str(trial) '.csv'];
                    rawData = csvread(fileName, nHeaderLines, 0);
                catch
                    disp(['[ERROR] ' dataset ' Trial ' num2str(trial) ' not found!'])
                    break
                end
                raw(trial,:) = rawData(1:nSamples,1); %has only motion values
                
                rawTrial = raw(trial,:)>threshold;
                
                %Reshape to get the averaging window
                rawTrialReshaped = reshape(rawTrial, [win4avg nSamples/win4avg]);
                diffRawTrialReshaped = diff(rawTrialReshaped');
                %motion(trial,:) = (sum(rawTrialReshaped))*distanceLC;
                %rate(trial,:) = motion(trial,:)/timeLC;
                %velocity(trial,:) = 
                disp('... done!')
            end
            % Get rid of the dummy trial
        end
    end
end
toc
disp('All done!')