% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - FEC analysis
% DEPENDENCIES - 1) Sort all trials as .mat files using sortingVideos.m
%                2) Find the image processing parameters using findTheEye.m

tic
clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations (0 == Don't Perform; 1 == Perform)
saveData = 1;
doFECAnalysis = 1;
plotFigures = 1;
playVideo = 0;

%Dataset details
sessionType = 11;
mice = [2 5];
%mice = 5;
nSessions = 3;
nTrials = 60; % NOTE: During sorting, the dummy trial was excluded

startSession = nSessions; %single sessions
%startSession = 3;
startTrial = 1; % NOTE: During sorting, the dummy trial was excluded
startFrame = 1;

loadDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/ImageProcess/';
saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

fontSize = 12;

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if doFECAnalysis == 1
            disp('Performing FEC analysis ...')
            
            %load image processing parameters
            load([loadDirec 'Mouse' mouseName '/' dataset '/imageProcess.mat']);
            
            %Video details
            nFrames = floor(samplingRate*trialDuration); %per trial
            time = 1:(1*1000/samplingRate):nFrames*1000/samplingRate; % in ms
            
            %Preallocation
            eyeClosure = zeros(nTrials,nFrames); %for every individual session
            eyeClosure_baseline = zeros(nTrials,1);
            fec = zeros(nTrials,nFrames);
            
            %Analyze every trial for FEC
            for trial = startTrial:nTrials
                disp(['Trial ' num2str(trial)])
                %1 - Load the reference image (first image in Trial 1)
                raw = load([direc 'Mouse' mouseName '/' dataset, ...
                    '/' dataset '_Trial' num2str(trial)]);
                for frame = startFrame:nFrames
                    refImage = rgb2gray(raw.raw(:,:,:,frame));
                    
                    % 2 - Adjust contrast
                    refImage2 = imadjust(refImage,[low_in; high_in],[low_out; high_out]);
                    
                    %3 - Median filter
                    refImage3 = medfilt2(refImage2,[m m]);
                    
                    %4 - Crop image
                    croppedImage = imcrop(refImage3,crop);
                    
                    %5 - Binarize
                    croppedImage2 = im2bw(croppedImage,level);
                    
                    eyeClosure(trial,frame) = (length(find(~croppedImage2(:,fecROI))))/length(find(croppedImage(:,fecROI)));
                    
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
                        imagesc(croppedImage2)
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
                eyeClosure_baseline(trial) = max(eyeClosure(trial,:));
                fec(trial,:) = 1 - (eyeClosure(trial,:)/eyeClosure_baseline(trial));
                disp('... done')
            end
            
            if plotFigures == 1
                figure(4)
                clf
                for trial = 1:nTrials
                    hold on
                    plot(time,((fec(trial,:)*3+(1*(trial-1)))),...
                        'black', 'LineWidth', 2)
                end
                axis([0 max(time) 1 nTrials]);
                title(['FEC [0: OPEN; 1: CLOSED] - ' ...
                    mouseName ' ST' num2str(sessionType) ' S' num2str(session)],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                xlabel('Time/ms', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                set(gca,'xtick',[100 300 500 ...
                    700 900 1100 ...
                    1300 1500])
                ylabel('Trials', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                ylim([0 nTrials+3])
                set(gca,'YTickLabel',[1 60])
                set(gca,'YTick',[1 60])
                
                print(['/Users/ananth/Desktop/figs/fec_' mouseName ...
                    '_ST' num2str(sessionType) ...
                    '_S' num2str(session)],...
                    '-djpeg');
                
                figure(5)
                clf
                imagesc(fec)
                colormap(jet)
                title(['FEC [0: OPEN; 1: CLOSED] - ' ...
                    mouseName ' ST' num2str(sessionType) ' S' num2str(session) ...
                    ' (' num2str(samplingRate) ' fps)'],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                xlabel('Frame Number', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                set(gca,'XTick', [10 30 50 ...
                    70 90 110 130 150])
                set(gca,'XTickLabel', [10 30 50 ...
                    70 90 110 130 150])
                ylabel('Trials', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                z = colorbar;
                ylabel(z,'FEC',...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                
                print(['/Users/ananth/Desktop/figs/FEC/fec_heatmap_' mouseName ...
                    '_ST' num2str(sessionType) ...
                    '_S' num2str(session)],...
                    '-djpeg');
            end
            
            if saveData == 1
                saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
                if ~isdir(saveFolder)
                    mkdir(saveFolder);
                end
                
                %Save FEC curve
                save([saveFolder 'fec.mat' ], ...
                    'eyeClosure', 'fec', ...
                    'low_in', 'high_in', 'low_out', 'high_out',...
                    'crop', 'fecROI', 'm', 'level',...
                    'samplingRate','trialDuration')
            end
            disp([dataset ' analyzed'])
        end
    end
end
disp('All done!')
toc
beep
