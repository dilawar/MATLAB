% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - CR Timing
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
plotFigures = 1;

% Dataset details
sessionType = 9;
mice = [1 2 3 4 5];
%mice = 5;
nSessions = 12;
score = nan(length(mice), nSessions);

fecDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
motionDirec = '/Users/ananth/Desktop/Work/Analysis/MotionAnalysis/';
performanceDirec = '/Users/ananth/Desktop/Work/Analysis/PerformanceAnalysis/';
saveDirec = performanceDirec;

% Protocol details
preCSTime = 0.5; % in seconds
csTime = 0.05; % in seconds
usTime = 0.05; % in seconds

fontSize = 12;

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if estimateTiming == 1
            disp('Estimating Response Timing ...')
            
            % Load FEC data
            try
                load([fecDirec 'Mouse' mouseName '/' dataset '/fec.mat'])
            catch
                warning('Unable to find FEC data')
                continue
            end
            
            nTrials = size(fec,1);
            
            for trials = 1:nTrials
                
            end
        end
        if saveData == 1
            saveFolder = saveDirec;
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            
            % trialInfo
            save([saveFolder 'Trial' num2str(trial) '.mat' ], ...
                'fanoFactor',...
                'hitTrials')
        end
    end
end
toc
disp('All done!')