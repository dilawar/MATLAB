% PURPOSE - To establish the best parameters to find the eye
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

%clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 1;
playVideo = 0;

%Dataset details
sessionType = 9;
%mice = [1 2 3 4 5];
mice = 5;
nSessions = 5;

%Cropping parameters
xmin = 127;
ymin = 49;
width = 39;
height = 39;
crop = [xmin ymin width height]; %[xmin ymin width height]
fecROI = 23:24;

%Filters
m = 4; %for median filter
%level = 0; %for binarization

%Contrast adjustment parameters
low_in = 0;
high_in = 0.05;
low_out = 0;
high_out = 1;

nTrials = 1;
startSession = nSessions;
startTrial = 1;
startFrame = 35;

%Video details
samplingRate = 100; % in Frames Per Second (FPS)
trialDuration = 1.5; % in seconds
%nFrames = floor(samplingRate*trialDuration); %per trial
nFrames = 110;

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/ImageProcess/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

fontSize = 12;
prctileVal = 5;

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if playVideo == 1
            disp('Playing Video ...');
            for trial = startTrial:nTrials
                %1 - Load the reference image (first image in Trial 1)
                raw = load([direc 'Mouse' mouseName '/' dataset, ...
                    '/' dataset '_Trial' num2str(trial)]);
                for frame = startFrame:nFrames
                    refImage = rgb2gray(raw.raw(:,:,:,frame));
                    
                    %2 - Adjust contrast
                    refImage2 = imadjust(refImage,[low_in; high_in],[low_out; high_out]);
                    
                    %3 - Median filter
                    refImage3 = medfilt2(refImage2,[m m]);
                    
                    %4 - Crop image
                    croppedImage = imcrop(refImage3,crop);
                    
                    pause(0.1)
                    figure(2)
                    subplot(1,3,1)
                    imagesc(croppedImage)
                    colormap(hot)
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
                    
                    %5 - Binarize
                    croppedImage2 = im2bw(croppedImage,level);
                    subplot(1,3,2)
                    imagesc(croppedImage2)
                    colormap(hot)
                    z = colorbar;
                    ylabel(z,'Intensity (A.U.)', ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                    title(['Binarized Frame ' num2str(frame)], ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                    
                    subplot(1,3,3)
                    imagesc(croppedImage2(:,fecROI))
                    colormap(hot)
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
            raw = load([direc 'Mouse' mouseName '/' dataset '/' dataset '_Trial1']);
            refImage = rgb2gray(raw.raw(:,:,:,1));
            
            %2 - Adjust contrast
            refImage2 = imadjust(refImage,[low_in; high_in],[low_out; high_out]);
            
            %3 - Median filter
            refImage3 = medfilt2(refImage2,[m m]);
            
            %4 - Crop image
            croppedImage = imcrop(refImage3,crop);
            %croppedImage_new = imcrop(refImage,crop);
            
            %5 - Binarize
            level = prctile(reshape(croppedImage,1,[]),prctileVal);
            croppedImage2 = croppedImage > level; %binarize
            
            figure(1)
            subplot(1,3,1)
            imagesc(refImage)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            title('Original', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            
            subplot(1,3,2)
            imagesc(refImage2)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            title('Contrast Adjusted', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            
            subplot(1,3,3)
            imagesc(refImage3)
            colormap(hot)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            title(['Median Filtered [' ...
                num2str(m) 'x' num2str(m) ']' ], ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            
            figure(3)
            subplot(1,3,1)
            imagesc(croppedImage);
            colormap(hot)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            title('Original (cropped)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            
            subplot(1,3,2)
            imagesc(croppedImage);
            colormap(hot)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            title(['Median Filtered [' ...
                num2str(m) 'x' num2str(m) ']' ],...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            
            subplot(1,3,3)
            imagesc(croppedImage2)
            colormap(hot)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            title('Binarized', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
        end
        if saveData == 1
            saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            
            %Save FEC curve
            save([saveFolder 'imageProcess.mat' ], ...
                'low_in', 'high_in', 'low_out', 'high_out',...
                'crop', 'fecROI', 'm', 'prctileVal'...
                'samplingRate','trialDuration')
        end
    end
end
disp('All done!')