% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - To cycle through saved .mat files for trials and establish FEC
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
playVideos = 0;

%Dataset details
%mice = [1 2 3 4 5];
mice = 5;
sessionType = 9;
nSessions = 1;
nTrials = 61;

%Video details
samplingRate = 100; % in Frames Per Second (fps)
trialDuration = 2; % in seconds
nFrames = samplingRate*trialDuration; %per trial

%startSession = nSessions; %single sessions
startSession = 1;
startTrial = 2; % The first trial is only a dummy

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/';
direc = '/Users/ananth/Desktop/Work/Behaviour/VideoAnalysis/Videos/';

for mouse = 1:length(mice)
    mouseName = ['MouseM' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = [mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(dataset);
        saveFolder 
    end
end

