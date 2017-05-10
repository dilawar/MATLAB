% PURPOSE - To establish the best parameters to find the eye
% DEPENDENCIES - None.
% NOTE - sortingVideos.m is now obsolete

%clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 1;
playVideo = 0;

%Dataset details
sessionType = 9;
%mice = [1 2 3 4 5];
mice = 10;
nSessions = 7;

%Crop parameters
xmin1 = 215;
ymin1 = 140;
width1 = 200;
height1 = 120;
crop = [xmin1 ymin1 width1 height1]; %[xmin ymin width height] of refImage

%FEC parameters
xmin2 = 118;
ymin2 = 5;
width2 = 30; % Don't change
height2 = 120;% Don't change
fecROI = [xmin2 ymin2 width2 height2]; %[xmin ymin width height] of croppedImage

%Filters
%m = 2; %for median filter
%level = 0.01; %for binarization

% %Contrast adjustment parameters
% low_in = 0;
% high_in = 1;
% low_out = 0;
% high_out = 1;

nTrials = 1;
startSession = nSessions;
startTrial = 1;
startFrame = 1;

%Video details
samplingRate = 100; % in Frames Per Second (FPS)
trialDuration = 1.5; % in seconds
%nFrames = floor(samplingRate*trialDuration); %per trial
nFrames = 250;

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/ImageProcess/';
%rawDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';
rawDirec = '/Users/ananth/Desktop/Work/Behaviour/DATA/';

fontSize = 16;

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        %         dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        dataset = [mouseName '_' num2str(sessionType) '_' num2str(session)];
        disp(['Working on ' dataset])
        
        if playVideo == 1
            disp('Playing Video ...');
            for trial = startTrial:nTrials
                
                if trial <10
                    file = [rawDirec mouseName '/' dataset, ...
                        '/trial_00' num2str(trial) '.tif'];
                else
                    file = [rawDirec mouseName '/' dataset, ...
                        '/trial_0' num2str(trial) '.tif'];
                end
                
                for frame = startFrame:nFrames
                    %1 - Load the reference image (first image in Trial 1)
                    refImage = double(imread(file, frame));
                    
                    %2 - Crop image - for eye (absolute coordinates)
                    croppedImage = imcrop(refImage,crop);
                    
                    %3 - Crop again - for FEC (relative coordinates)
                    fecImage = imcrop(croppedImage,fecROI);
                    
                    %4 - Binarize
                    if frame == startFrame
                        threshold = prctile(reshape(fecImage,[((height2+1)*(width2+1)),1]),50);
                    end
                    binImage = fecImage > threshold; %binarize
                    
                    pause(0.1)
                    figure(2)
                    subplot(1,3,1)
                    imagesc(croppedImage)
                    colormap(gray)
                    z = colorbar;
                    ylabel(z,'Intensity (A.U.)', ...
                        'FontSize', fontSize,...
                        'FontWeight', 'bold')
                    title(['Eye - ' mouseName ...
                        ' ST' num2str(sessionType) ' S' num2str(session) ...
                        ' Trial ' num2str(trial) ...
                        ' Frame ' num2str(frame)], ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                    
                    subplot(1,3,2)
                    imagesc(fecImage)
                    colormap(gray)
                    z = colorbar;
                    ylabel(z,'Intensity (A.U.)', ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                    title(['Binarized Frame ' num2str(frame)], ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                    
                    subplot(1,3,3)
                    imagesc(binImage)
                    colormap(gray)
                    z = colorbar;
                    ylabel(z,'Intensity (A.U.)', ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                    title(['fecROI Frame ' num2str(frame)], ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                end
            end
        else
            %1 - Load the reference image (first image in Trial 1)
            file = [rawDirec mouseName '/' dataset '/trial_001.tif'];
            refImage = double(imread(file, 1));
            
            %2 - Crop image - for eye (absolute coordinates)
            croppedImage = imcrop(refImage,crop);
            
            %3 - Crop again - for FEC (relative coordinates)
            fecImage = imcrop(croppedImage,fecROI);
            
            %4 - Binarize
            %threshold = prctile(reshape(fecImage,[((height2+1)*(width2+1)),1]),50);
            fecImage_vector = reshape(fecImage,1,[]);
            threshold = prctile(fecImage_vector,50);
            binImage = fecImage > threshold; %binarize
            
            figure(1)
            subplot(2,2,1)
            imagesc(refImage)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            title('1. Original', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            set(gca,'FontSize',fontSize)
            
            subplot(2,2,2)
            imagesc(croppedImage)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            title('2. Cropped', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            set(gca,'FontSize',fontSize)
            
            subplot(2,2,3)
            imagesc(fecImage)
            colormap(gray)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            title('3. FEC ROI', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            set(gca,'FontSize',fontSize)
            
            subplot(2,2,4)
            imagesc(binImage);
            colormap(gray)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            title('4. Binarized', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            set(gca,'FontSize',fontSize)
            
        end
        
        if saveData == 1
            saveFolder = [saveDirec mouseName '/' dataset '/'];
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            %Save FEC curve
            save([saveFolder 'imageProcess.mat' ], ...
                'crop', 'fecROI', 'threshold')
        end
    end
end
disp('All done!')