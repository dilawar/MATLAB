% PURPOSE - Analysis code for TEC behaviour including FEC and Motion
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

%clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 0;
playVideo = 1;

%Contrast adjustment parameters
low_in = 0;
high_in = 0.2;
low_out = 0;
high_out = 1;

crop = [129 55 39 24]; %[xmin ymin width height]
fecROI = 15:30;
m = 5; %for median filter
level = 0.01; %for binarization

if playVideo == 1
    nTrials = 1;
    startTrial = 1;
    startFrame = 1;
    
    %Video details
    samplingRate = 100; % in Frames Per Second (FPS)
    trialDuration = 1; % in seconds
    %nFrames = 60;
    nFrames = samplingRate*trialDuration; %per trial
end

%Dataset details
mouseName = 'M2';
sessionType = 9;
session = 8;

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
disp(['Working on ' dataset])

fontSize = 12;

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
            
            figure(2)
            pause(0.1)
            subplot(1,2,1)
            imagesc(croppedImage)
            colormap(hot)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            title(['Eye ' mouseName ...
                ' ST' num2str(sessionType) ' S' num2str(session) ...
                ' Trial ' num2str(trial) ...
                ' Frame ' num2str(frame)], ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            
            %5 - Binarize
            croppedImage2 = im2bw(croppedImage,level);
            subplot(1,2,2)
            imagesc(croppedImage2(:,fecROI))
            colormap(hot)
            z = colorbar;
            ylabel(z,'Intensity (A.U.)', ...
                'FontSize', fontSize, ...
                'FontWeight', 'bold')
            title('Binary', ...
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
    
    %5 - Binarize
    croppedImage2 = im2bw(croppedImage,level);
    
    figure(1)
    subplot(1,2,1)
    imagesc(refImage2)
    z = colorbar;
    ylabel(z,'Intensity (A.U.)', ...
        'FontSize', fontSize,...
        'FontWeight', 'bold')
    title('Contrast Adjusted', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    
    subplot(1,2,2)
    imagesc(refImage3)
    colormap(hot)
    z = colorbar;
    ylabel(z,'Intensity (A.U.)', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    title('Median Filtered', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    
    figure(2)
    subplot(1,3,1)
    imagesc(croppedImage);
    colormap(hot)
    z = colorbar;
    ylabel(z,'Intensity (A.U.)', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    title('Normal', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    
    subplot(1,3,2)
    imagesc(croppedImage);
    colormap(hot)
    z = colorbar;
    ylabel(z,'Intensity (A.U.)', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    title('Filtered', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    
    subplot(1,3,3)
    imagesc(croppedImage2)
    colormap(hot)
    z = colorbar;
    ylabel(z,'Intensity (A.U.)', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
    title('Binary', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold')
end
disp('All done!')