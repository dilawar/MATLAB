% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - To segment eye-blink behaviour videos into appropriate trials
% and save as .mat files (MATLAB-friendly)

tic
close all
clear all

% "raw" is defined as height x width x 3 x nFrames, and is saved for every trial

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
sortVideos = 1;
playVideos = 0;

singleVideo = 0; %Does the session have only one video or several?

%Dataset details
mice = [1 2 3 4 5];
sessionType = 9;
nSessions = 9;
nTrials = 61;

%Video details
height = 300;
width = 300;
samplingRate = 100; % in Frames Per Second (fps)
trialDuration = 2; % in seconds
nFrames = samplingRate*trialDuration; %per trial

startSession = nSessions; %single sessions
%startSession = 1;
startTrial = 2; % The first trial is only a dummy

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';
direc = '/Users/ananth/Desktop/Work/Behaviour/Videos/';

for mouse = 1:length(mice)
    
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(dataset);
        saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
        if ~isdir(saveFolder)
            mkdir(saveFolder);
        end
        if sortVideos == 1
            disp('Sorting Videos ...')
            if singleVideo == 1
                fileName = [direc 'Mouse' mouseName '/' dataset '/M' ...
                    num2str(mice(mouse)) '_ST' num2str(sessionType) ...
                    '_S' num2str(session) '.avi'];
                %disp(fileName);
                video = VideoReader(fileName);
                for trial = startTrial:nTrials
                    %beep;
                    disp(['Trial ' num2str(trial-1)]);
                    raw = read(video,[((trial-1)*nFrames+1) (((trial-1)*nFrames+1)+(nFrames-1))]);
                    save([saveFolder, dataset, '_Trial' num2str(trial-1)],'raw') % for every trial, independently
                    disp('... done')
                end
                clear video
                %disp('... done!')
            else
                for trial = startTrial:nTrials
                    %beep;
                    disp(['Trial ' num2str(trial-1)]);
                    fileName = [direc 'Mouse' mouseName '/' dataset '/Trial' num2str(trial) '.avi']; %post ffmpeg conversion
                    %disp(fileName);
                    try
                        video = VideoReader(fileName);
                        raw = read(video,[1 nFrames]);
                        save([saveFolder, dataset, '_Trial' num2str(trial-1)],'raw') % for every trial, independently
                        disp('... done')
                    catch
                        disp(['[ERROR] ' dataset ' Trial ' num2str(trial) ' not found!'])
                    end
                    clear video
                end
            end
        end
        
        if playVideos == 1
            disp('Playing Videos ...')
            for trial = startTrial:nTrials
                raw = load([saveFolder dataset '_Trial' num2str(trial)]);
                for frame = 45:55
                    figure(1);
                    imshow(raw.raw(:,:,:,frame));
                    title(['Trial ' num2str(trial-1) ' Frame ' num2str(frame)]);
                    pause(0.1);
                end
                disp(['Trial ' num2str(trial-1) '...done!'])
            end
        end
    end
end
disp('All Done!')
toc
beep
