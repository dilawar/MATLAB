% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Performance analysis
% DEPENDENCIES - 1) Sort all trials as .mat files using sortingVideos.m
%                2) Find the image processing parameters using findTheEye.m
%                3) Perform FEC analysis
%                4) Perform Motion analysis

tic
clear all
close all

% Operations (0 == Don't Perform; 1 == Perform)
saveData = 1;
scorePerformance = 1;
plotFigures = 1;

% Dataset details
sessionType = 9;
mice = [1 2 3 4 5];
nSessions = 12;

%startSession = 1;
startSession = nSessions;

fecDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
motionDirec = '/Users/ananth/Desktop/Work/Analysis/MotionAnalysis/';
saveDirec = '/Users/ananth/Desktop/Work/Analysis/PerformanceAnalysis/';

fontSize = 12;

for session = startSession:nSessions
    for mouse = 1:length(mice)
        mouseName = ['M' num2str(mice(mouse))];
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if scorePerformance == 1
            disp('Scoring Performance ...')
            
            % Load FEC data
            load([fecDirec 'Mouse' mouseName '/' dataset '/fec.mat'])
            
            % Load Motion data
            load([motionDirec 'Mouse' mouseName '/' dataset '/motion.mat'])
            
            
        end
    end
end
toc
disp('All done!')
