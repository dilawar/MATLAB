% PURPOSE - To establish the best parameters to find the eye
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

%clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 0;
playVideo = 0;

crop = [123 54 39 29]; %[xmin ymin width height]
fecROI = (20:40);
m = 3; %for median filter
level = 0.025; %for binarization

%Contrast adjustment parameters
low_in = 0;
high_in = 0.9;
low_out = 0;
high_out = 1;

%Dataset details
sessionType = 9;
mice = 1;
nSessions = 3;
nTrials = 1;

startSession = nSessions;
startTrial = 1;
startFrame = 35;

%Video details
samplingRate = 100.5; % in Frames Per Second (FPS)
trialDuration = 1.5; % in seconds
%nFrames = floor(samplingRate*trialDuration); %per trial
nFrames = 110;

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

fontSize = 12;

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
                    
                    figure(2)
                    pause(0.05)
                    subplot(1,3,1)
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
                    subplot(1,3,2)
                    imagesc(croppedImage2)
                    colormap(hot)
                    z = colorbar;
                    ylabel(z,'Intensity (A.U.)', ...
                        'FontSize', fontSize, ...
                        'FontWeight', 'bold')
                    title(['Binary Frame ' num2str(frame)], ...
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
            croppedImage_new = imcrop(refImage,crop);
            
            %5 - Binarize
            croppedImage2 = im2bw(croppedImage,level);
            
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
            title('Median Filtered', ...
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
        if saveData == 1
            saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            
            %Save FEC curve
            save([saveFolder 'imageProcess.mat' ], ...
                'low_in', 'high_in', 'low_out', 'high_out',...
                'crop', 'fecROI', 'm', 'level',...
                'samplingRate','trialDuration')
        end
    end
end
disp('All done!')