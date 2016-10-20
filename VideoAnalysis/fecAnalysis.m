% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Analysis code for TEC behaviour including FEC and Motion
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 0;
cropRefImage = 1;
doFECAnalysis = 0;
doMotionAnalysis = 0;
playVideos = 0;

%Dataset details
%mice = [1 2 3 4 5];
mice = 5;
sessionType = 9;
nSessions = 1;
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
        
        if cropRefImage == 1
            %Load the reference image (first image in Trial 1)
            raw = load([direc mouseName '/' dataset '/' dataset '_Trial1']);
            referenceImage = raw.raw(:,:,:,1);
            figure(2)
            imagesc(referenceImage)
            colormap(jet)
            colorbar
            title(['Reference Image for ' ...
                'M' num2str(mice(mouse)) ' ST' num2str(sessionType) ' S' num2str(session)])
            
            %Crop image
            rectCrop = [100 35 70 50]; %[xmin ymin width height]
            croppedImage = imcrop(referenceImage,rectCrop);
            figure(3)
            imagesc(croppedImage)
            colormap(jet)
            colorbar
            title(['Cropped Image for ' ...
                'M' num2str(mice(mouse)) ' ST' num2str(sessionType) ' S' num2str(session)])
        else
            %Load Reference and Cropped Images
            referenceImage = load();
            croppedImage = load();
        end
        
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
