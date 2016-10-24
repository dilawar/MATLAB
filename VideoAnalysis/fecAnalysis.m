% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - FEC analysis
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

tic
%clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 1;
doFECAnalysis = 1;
plotFigures = 1;
playVideo = 0;

%Dataset details
sessionType = 9;
%mice = [1 2 3 4 5];
mice = 2;
nSessions = 1;
nTrials = 60; % NOTE: During sorting, the dummy trial was excluded

startSession = nSessions; %single sessions
%startSession = 1;
startTrial = 1; % NOTE: During sorting, the dummy trial was excluded

%Video details
samplingRate = 118; % in Frames Per Second (FPS)
trialDuration = 1.5; % in seconds
nFrames = floor(samplingRate*trialDuration); %per trial
startFrame = 1;

%Contrast adjustment parameters
low_in = 0;
high_in = 1;
low_out = 0;
high_out = 1;

crop = [108 45 39 24]; %[xmin ymin width height]
fecROI = (22:37);
m = 3; %for median filter
level = 0.27; %for binarization

saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';
direc = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';

fontSize = 12;

eyeClosure = zeros(nTrials,nFrames); %for every individual session
eyeClosure_baseline = zeros(nTrials,1);
fec = zeros(nTrials,nFrames); %fractional; for every individual session

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
                    
                    eyeClosure(trial,frame) = length(find(~croppedImage2(:,fecROI)));
                    
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
                fec(trial,:) = (1 - eyeClosure(trial,:)/eyeClosure_baseline(trial));
                disp('... done')
            end
            
            %eyeClosure_baseline = max(max(eyeClosure));
            %fec = 1 - eyeClosure/eyeClosure_baseline;
            time = 1:(1*1000/samplingRate):nFrames*1000/samplingRate; % in ms
            if plotFigures == 1
                figure(3)
                for trial = 1:nTrials
                    plot(time,((fec(trial,:)*2+(1*(trial-1)))),...
                        'black', 'LineWidth', 2)
                    hold on
                end
                axis([0 max(time) 1 nTrials]);
                title([mouseName ' ST' num2str(sessionType) ' S' num2str(session)],...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                xlabel('Time/ms', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                set(gca,'xtick',[200 400 600 ...
                    800 1000 1200 1400 1600 ...
                    1800 2000])
                ylabel('Trials', ...
                    'FontSize', fontSize,...
                    'FontWeight', 'bold')
                ylim([0 nTrials+1])
                set(gca,'yticklabel',[])
                set(gca,'ytick',[1 60])
                
                print(['/Users/ananth/Desktop/figs/fec_' mouseName ...
                    '_ST' num2str(sessionType) ...
                    'S' num2str(session)],...
                    '-djpeg');
            end
            
            if saveData == 1
                saveFolder = [saveDirec mouseName '/' dataset '/'];
                if ~isdir(saveFolder)
                    mkdir(saveFolder);
                end
                
                %Save FEC curve
                save([saveFolder 'fec.mat' ], ...
                    'eyeClosure',...
                    'fec')
                
            end
            disp([dataset ' analyzed'])
        end
    end
end
disp('All done!')
toc
beep
