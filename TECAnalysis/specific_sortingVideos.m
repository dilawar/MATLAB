% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - To segment eye-blink behaviour videos into appropriate trials
% and save as .mat files (lighter; MATLAB-friendly)

close all
clear all

% The variale raw is defined by the following dimensions: 276x276x3xnFrames
% and is saved for every trial

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
sortVideos = 0;
playVideos = 1;

singleVideo = 1; %Does the session have only one video or several?

%Dataset details
mice = [5];
sessionType = 9;
nSessions = 1;
startSession = 1;
nTrials = 61;

% Video details
height = 276;
width = height;
samplingRate = 100; % in Frames Per Second (fps)
trialDuration = 2; % in seconds
%nFrames = samplingRate*trialDuration; %per trial
nFrames = 200;

saveDirec = '/Users/ananth/Desktop/temp/';
direc = '/Users/ananth/Desktop/Work/Behaviour/Videos/';

%disp('Ready to begin!')
for mouse = 1:length(mice)
    
    mouseName = ['MouseM' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        raw = zeros(nTrials,nFrames,height,width,3);
        dataset = [mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(dataset);
        saveFolder = [saveDirec mouseName '/' dataset '/'];
        if ~isdir(saveFolder)
            mkdir(saveFolder);
        end
        if sortVideos == 1
            disp('Sorting Videos ...')
            if singleVideo == 1
                filename = [direc mouseName '/' dataset '/M' ...
                    num2str(mice(mouse)) '_ST' num2str(sessionType) ...
                    '_S' num2str(session) '.avi'];
                %disp(filename);
                video = VideoReader(filename);
                for trial = 1:nTrials
                    %beep;
                    disp(['Trial ' num2str(trial)]);
                    raw = read(video,[((trial-1)*nFrames+1) (((trial-1)*nFrames+1)+(nFrames-1))]);
                    save([saveFolder, dataset, '_Trial' num2str(trial)],'raw') % for every trial, independently
                    disp('... done')
                end
                clear video
            else
                for trial = 1:nTrials
                    %beep;
                    disp(['Trial ' num2str(trial)]);
                    filename = [direc mouseName '/' dataset '/Trial' num2str(trial) '.avi'];
                    %disp(filename);
                    try
                        video = VideoReader(filename);
                        raw = read(video,[1 nFrames]);
                        save([saveFolder, dataset, '_Trial' num2str(trial)],'raw') % for every trial, independently
                        disp('... done')
                    catch
                        disp(['Trial ' num2str(trial) ' not found!'])
                    end
                    clear video
                end
            end
        end
        if playVideos == 1
            disp('Playing Videos ...')
            for trial = 1:68
                raw = load([saveFolder dataset '_Trial' num2str(trial)]);
                for frame = 45:55
                    figure(1);
                    imshow(raw.raw(:,:,:,frame));
                    title(['Trial ' num2str(trial) ' Frame ' num2str(frame)]);
                    pause(0.1);
                end
                disp(['Trial ' num2str(trial) '...done!'])
            end
        end
    end
end
disp('All Done!')


