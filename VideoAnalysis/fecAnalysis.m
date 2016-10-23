% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Analysis code for TEC behaviour including FEC and Motion
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

tic
%clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 0;
doFECAnalysis = 1;
doMotionAnalysis = 0;
plotFigures = 1;
playVideo = 0;

%Video details
samplingRate = 100; % in Frames Per Second (FPS)
trialDuration = 2; % in seconds
nFrames = samplingRate*trialDuration; %per trial
startFrame = 1;

crop = [120 55 50 30]; %[xmin ymin width height]
fecROI = 30:40; %to avoid the regions where the blue LED light flashes
m = 5; %for median filter
level = 0.5; %for binarization

%Dataset details
sessionType = 9;
%mice = [1 2 3 4 5];
mice = 5;
nSessions = 5;
nTrials = 60; % NOTE: During sorting, the dummy trial was excluded

startSession = nSessions; %single sessions
%startSession = 1;
startTrial = 1; % NOTE: During sorting, the dummy trial was excluded

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

fontSize = 12;

eyeClosure = zeros(nTrials,nFrames); %for every individual session
fec = zeros(nTrials,nFrames); %fractional; for every individual session
motion = zeros(nTrials,nFrames); %!! EDIT !!

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if doFECAnalysis == 1
            disp('Performing FEC analysis ...')
            %Analyze every trial for FEC
            for trial = startTrial:nTrials
                disp(['Trial ' num2str(trial)])
                raw = load([direc 'Mouse' mouseName '/' dataset '/' dataset '_Trial' num2str(trial)]);
                
                for frame = startFrame:nFrames
                    refImage = rgb2gray(raw.raw(:,:,:,frame));
                    
                    %Equalize histogram
                    refImage_histeq = histeq(refImage);
                    
                    %Median filter
                    refImage_medfilt = medfilt2(refImage_histeq,[m m]);
                    
                    %Crop image
                    croppedImage = imcrop(refImage_medfilt,crop);
                    
                    %Binarize
                    croppedImage_binary = im2bw(croppedImage,level);
                    %croppedImage_binary = imbinarize(croppedImage);
                    
                    eyeClosure(trial,frame) = length(find(~croppedImage_binary(:,fecROI)));
                    
                    if playVideo == 1
                        figure(2)
                        pause(0.1)
                        subplot(1,2,1)
                        imagesc(croppedImage)
                        colormap(hot)
                        z = colorbar;
                        ylabel(z,'Intensity (A.U.)', ...
                            'FontSize', fontSize,...
                            'FontWeight', 'bold')
                        title(['Cropped Image ' mouseName ...
                            ' ST' num2str(sessionType) ' S' num2str(session) ...
                            'Trial ' num2str(trial)], ...
                            'FontSize', fontSize, ...
                            'FontWeight', 'bold')
                        
                        subplot(1,2,2)
                        imagesc(croppedImage_binary)
                        colormap(hot)
                        z = colorbar;
                        ylabel(z,'Intensity (A.U.)', ...
                            'FontSize', fontSize,...
                            'FontWeight', 'bold')
                        title(['Binarized Frame ' num2str(frame)], ...
                            'FontSize', fontSize, ...
                            'FontWeight', 'bold')
                    end
                end
                
                if doMotionAnalysis == 1
                    disp('Performing motion analysis ...')
                end
                
                disp('... done')
            end
            
            eyeClosure_baseline = max(median(eyeClosure,2));
            fec = 1 - eyeClosure/eyeClosure_baseline;
            
            if saveData == 1
                saveFolder = [saveDirec mouseName '/' dataset '/'];
                if ~isdir(saveFolder)
                    mkdir(saveFolder);
                end
                
                %Save FEC curve
                save([saveFolder 'eyeClosure' ],'eyeClosure')
                save([saveFolder 'fec'],'fec')
                
                %Save Motion data
                save([saveFolder 'motion' ],'motion')
            end
            
            if plotFigures == 1
                figure(3)
                imagesc(fec)
                colormap(hot)
                colorbar
                %waterfall(fec)
                title([mouseName ' ST' num2str(sessionType) ' S' num2str(session)],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                xlabel('Frame', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                ylabel('Trial', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                z = colorbar;
                ylabel(z,'FEC', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                
                %                 figure(4)
                %                 plot(eyeClosure_baselines, 'LineWidth', 2)
                %                 title(dataset,...
                %                     'FontSize', fontSize,...
                %                     'FontWeight', 'bold')
                %                 xlabel('Frame', ...
                %                     'FontSize', fontSize,...
                %                     'FontWeight', 'bold')
                %                 ylabel('Trial', ...
                %                     'FontSize', fontSize,...
                %                     'FontWeight', 'bold')
            end
            disp([dataset ' analyzed'])
        end
    end
end
disp('All done!')
toc
