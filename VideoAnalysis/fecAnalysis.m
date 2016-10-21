% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Analysis code for TEC behaviour including FEC and Motion
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 0;
doFECAnalysis = 0;
doMotionAnalysis = 0;
playVideos = 0;

%Dataset details
%mice = [1 2 3 4 5];
mice = 5;
sessionType = 9;
nSessions = 5;
nTrials = 60; % NOTE: During sorting, the dummy trial was excluded

%Video details
samplingRate = 100; % in Frames Per Second (FPS)
trialDuration = 2; % in seconds
nFrames = samplingRate*trialDuration; %per trial
%height = 300;
%width = height;

startSession = nSessions; %single sessions
%startSession = 1;
startTrial = 1; % NOTE: During sorting, the dummy trial was excluded

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

for mouse = 1:length(mice)
    mouseName = ['MouseM' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = [mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        %Load Reference and Cropped Images
        refImage = load();
        croppedImage = load();
        
        if doFECAnalysis == 1
            disp('Performing FEC analysis ...')
            %Analyze every trial for FEC
            for trial = startTrial:nTrials
                disp(['Trial ' num2str(trial)])
                raw = load([direc mouseName '/' dataset '/' dataset '_Trial' num2str(trial)]);
                if playVideos == 1
                    for frame = 45:55
                        figure(1);
                        imshow(raw.raw(:,:,:,frame));
                        title(['M' num2str(mice(mouse)) ' ST' num2str(sessionType) ...
                            ' S' num2str(session) ' Trial ' num2str(trial-1) ...
                            'Frame ' num2str(frame)])
                        pause(0.1);
                    end
                end
                disp('... done!')
            end
            disp([dataset ' analyzed'])
        end
        
        if doMotionAnalysis == 1
            disp('Performing motion analysis ...')
            if playVideos == 1
                for frame = 45:55
                    figure(1);
                    imshow(raw.raw(:,:,:,frame));
                    title(['M' num2str(mice(mouse)) ' ST' num2str(sessionType) ...
                        ' S' num2str(session) ' Trial ' num2str(trial-1) ...
                        'Frame ' num2str(frame)])
                    pause(0.1)
                end
            end
        end
        
        if saveData == 1
            saveFolder = [saveDirec mouseName '/' dataset '/'];
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            
            %Save referenceimage
            
            %Save FEC curve
            
            %Save Motion data
            
            %Save figures
            
        end
    end
end
disp('All done!')
