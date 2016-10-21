% PURPOSE - Analysis code for TEC behaviour including FEC and Motion
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 0;

%Dataset details
mouseName = 'MouseM5';
sessionType = 9;
session = 5;

%Video details
% samplingRate = 100; % in Frames Per Second (FPS)
% trialDuration = 2; % in seconds
% nFrames = samplingRate*trialDuration; %per trial
%height = 300;
%width = height;

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

dataset = [mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
disp(['Working on ' dataset])

%Load the reference image (first image in Trial 1)
raw = load([direc mouseName '/' dataset '/' dataset '_Trial1']);
%refImage = rgb2gray(raw.raw(:,:,:,1));
refImage = rgb2gray(raw.raw(:,:,:,1));
figure(2)
subplot(1,2,1)
imagesc(refImage)
colormap(hot)
colorbar
title(['Reference Image for ' ...
    'M' num2str(mice(mouse)) ' ST' num2str(sessionType) ' S' num2str(session)])
%Equalize histogram
refImage_histeq = histeq(refImage);
subplot(1,2,2)
imagesc(refImage_histeq)
colormap(hot)
colorbar
title('Histogram Equalized')

%Crop image
%rectCrop = [100 35 50 30]; %[xmin ymin width height]
rectCrop = [125 55 50 30]; %[xmin ymin width height]
croppedImage = imcrop(refImage_histeq,rectCrop);

%Median filter
croppedImage_medfilt = medfilt2(croppedImage,[5 5]);

hist_crop = imhist(croppedImage);
hist_crop_medfilt = imhist(croppedImage_medfilt);
figure(3)
subplot(1,2,1)
plot(hist_crop)
%axis([0 60 0 500])
title(['Histogram of the Cropped Image for ' ...
    'M' num2str(mice(mouse)) ' ST' num2str(sessionType) ' S' num2str(session)])
xlabel('Value')
ylabel('Counts')
subplot(1,2,2)
plot(hist_crop_medfilt)
%axis([0 60 0 500])
title('Median Filtered')
xlabel('Value')
ylabel('Counts')

%Binarize
croppedImage_binary = im2bw(croppedImage_medfilt,0.5);
%croppedImage_binary = imbinarize(croppedImage_medfilt);

figure(4)
subplot(1,3,1)
imagesc(croppedImage)
colormap(hot)
colorbar
title(['Cropped Image for ' ...
    'M' num2str(mice(mouse)) ' ST' num2str(sessionType) ' S' num2str(session)])

subplot(1,3,2)
imagesc(croppedImage_medfilt)
colormap(hot)
colorbar
title('Median Filtered')

subplot(1,3,3)
imagesc(croppedImage_binary)
colormap(hot)
colorbar
title('Binarized')