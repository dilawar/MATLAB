% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - FEC analysis
% DEPENDENCIES - 1) Sort all trials as .mat files using sortingVideos.m
%                2) Find the image processing parameters using findTheEye.m

tic
clear all
%close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

% Operations (0 == Don't Perform; 1 == Perform)
saveData = 0;
doFECAnalysis = 0;
plotFigures = 1;
playVideo = 0;

% Dataset details
sessionType = 9;
%mice = [2 3 5];
mice = 3;
nSessions = 12;
nTrials = 60; % NOTE: During sorting, the dummy trial was excluded

startSession = nSessions; %single sessions
%startSession = 5;
startTrial = 1; % NOTE: During sorting, the dummy trial was excluded
startFrame = 1;

imageProcessDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/ImageProcess/';
rawDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/Videos/';
motionDirec = '/Users/ananth/Desktop/Work/Analysis/MotionAnalysis/';
performanceDirec = '/Users/ananth/Desktop/Work/Analysis/PerformanceAnalysis/';
saveDirec = '/Users/ananth/Desktop/Work/Analysis/VideoAnalysis/FEC/';

fontSize = 20;
lineWidth = 3;
markerWidth = 7;

timeLine = zeros(150,1);
csLine = timeLine;
usLine = timeLine;

csLine(50:55,1) = 1;
if sessionType == 9
    usLine(80:85) = 1;
else
    usLine(90:95) = 1;
end

for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if doFECAnalysis == 1
            disp('Performing FEC analysis ...')
            
            % Load image processing parameters
            load([imageProcessDirec 'Mouse' mouseName '/' dataset '/imageProcess.mat'])
            
            % Video details
            nFrames = floor(samplingRate*trialDuration); %per trial
            time = 1:(1*1000/samplingRate):nFrames*1000/samplingRate; % in ms
            
            % Preallocation
            eyeClosure = zeros(nTrials,nFrames); %for every individual session
            eyeClosure_baseline = zeros(nTrials,1);
            fec = zeros(nTrials,nFrames);
            
            % Analyze every trial for FEC
            for trial = startTrial:nTrials
                disp(['Trial ' num2str(trial)])
                %1 - Load the reference image (first image in Trial 1)
                raw = load([rawDirec 'Mouse' mouseName '/' dataset, ...
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
                    croppedImage2 = im2bw(croppedImage,level); %binarize
                    
                    eyeClosure(trial,frame) = (length(find(~croppedImage2(:,fecROI))))/length(reshape(croppedImage(:,fecROI),1,[]));
                    
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
        else
            load([saveDirec 'Mouse' mouseName '/' dataset '/fec.mat']);
            load([motionDirec 'Mouse' mouseName '/' dataset '/motion.mat']);
            load([performanceDirec 'Mouse' mouseName '/' dataset '/Session' num2str(session) '.mat']);
        end
        
        if plotFigures == 1
            %                 figure(4)
            %                 clf
            %                 for trial = 1:nTrials
            %                     hold on
            %                     plot(time,((fec(trial,:)*3+(1*(trial-1)))),...
            %                         'black', 'LineWidth', 2)
            %                 end
            %                 axis([0 max(time) 1 nTrials]);
            %                 title(['FEC [0: Open; 1: Closed] - ' ...
            %                     mouseName ' ST' num2str(sessionType) ' S' num2str(session)],...
            %                     'FontSize', fontSize,...
            %                     'FontWeight', 'bold')
            %                 xlabel('Time/ms', ...
            %                     'FontSize', fontSize,...
            %                     'FontWeight', 'bold')
            %                 set(gca,'xtick',[100 200 300 ...
            %                     400 500 600 ...
            %                     700 800 900 ...
            %                     1000 1100 1200 ...
            %                     1300 1400 1500])
            %                 ylabel('Trials', ...
            %                     'FontSize', fontSize,...
            %                     'FontWeight', 'bold')
            %                 ylim([0 nTrials+3])
            %                 set(gca,'YTickLabel',[1 60])
            %                 set(gca,'YTick',[1 60])
            %                 set(gca,'FontSize', fontSize)
            %
            %                 print(['/Users/ananth/Desktop/figs/FEC/fec_' mouseName ...
            %                     '_ST' num2str(sessionType) ...
            %                     '_S' num2str(session)],...
            %                     '-djpeg');
            
            figure(5)
            clf
            subplot(9,9,1:72)
            %subplot(6,9,1:45)
            imagesc(fec)
            colormap(jet)
            if sessionType == 9
                title([' FEC - 250 ms ISI - ' mouseName ' S' num2str(session)], ...
                    'FontSize', fontSize, ...
                    'FontWeight', 'bold')
            else
                title([' FEC - 350 ms ISI - ' mouseName ' S' num2str(session)], ...
                    'FontSize', fontSize, ...
                    'FontWeight', 'bold')
            end
            %             title('Weak Learner - 250 ms ISI - Session 12', ...
            %                 'FontSize', fontSize,...
            %                 'FontWeight', 'bold')
            %                 xlabel('Time/ms', ...
            %                     'FontSize', fontSize,...
            %                     'FontWeight', 'bold')
            set(gca,'XTick', [])
            set(gca,'XTickLabel', [])
            ylabel('Trials', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            z = colorbar;
            ylabel(z,'FEC [0: Open; 1: Closed]',...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            set(gca,'FontSize', fontSize-2)
            
            subplot(9,9,73:80)
            %subplot(6,9,46:53)
            plot(csLine,'-','LineWidth',lineWidth)
            hold on
            plot(usLine,'-','LineWidth',lineWidth)
            set(gca)
            legend('CS','US')
            set(gca,'FontSize', fontSize-2)
            set(gca,'YTick',[0 1])
            set(gca, 'YTickLabel', {'Off' 'On'})
            xlabel('Time/ms', ...
                'FontSize', fontSize,...
                'FontWeight', 'bold')
            set(gca,'XTick', [10 30 50 ...
                70 90 110 ...
                130 150])
            set(gca,'XTickLabel', [100 300 500 ...
                700 900 1100 ...
                1300 1500])
            if sessionType == 9
                if mice(mouse) == 5
                    print('/Users/ananth/Desktop/figs/fec_weakLearner', ...
                        '-dpng');
                elseif mice(mouse) == 3
                    print('/Users/ananth/Desktop/figs/fec_nonLearner', ...
                        '-dpng');
                elseif mice(mouse) == 2
                    print('/Users/ananth/Desktop/figs/fec_strongLearner', ...
                        '-dpng');
                end
            else
                print('/Users/ananth/Desktop/figs/fec_350', ...
                    '-dpng');
            end
            
            %             print(['/Users/ananth/Desktop/figs/FEC/fec_heatmap_' ...
            %                 mouseName '_ST' num2str(sessionType) '_S' num2str(session)],...
            %                 '-djpeg');
            
            hitList = zeros(nTrials,1);
            hitList(hitTrials) = 1;
            
            probeList = zeros(nTrials,1);
            probeList(probeTrials) = 1;
            
            figure(6);
            clf
            subplot(9,2,1:2:15)
            imagesc(hitList)
            title(['Hits - ' mouseName ' S' num2str(session)], ...
                'FontSize', fontSize,...
                'FontWeight','bold')
            ylabel('Trials', ...
                'FontSize', fontSize,...
                'FontWeight','bold')
            set(gca,'FontSize', fontSize-2)
            set(gca,'XTick',[])
            set(gca,'XTickLabel',[])
            set(gca,'YTick',[10 20 30 40 50 60])
            set(gca, 'YTickLabel', [10 20 30 40 50 60])
            z = colorbar;
            set(z,'YTick', [0 1],...
                'FontSize', fontSize)
            
            subplot(9,2,2:2:16)
            imagesc(probeList)
            title(['Probes - ' mouseName ' S' num2str(session)], ...
                'FontSize', fontSize,...
                'FontWeight','bold')
            colormap(gray)
            set(gca,'FontSize', fontSize-2)
            set(gca,'XTick',[])
            set(gca,'XTickLabel',[])
            set(gca,'YTick',[10 20 30 40 50 60])
            set(gca, 'YTickLabel', [10 20 30 40 50 60])
            z = colorbar;
            set(z,'YTick', [0 1],...
                'FontSize', fontSize)
            if sessionType == 9
                if mice(mouse) == 5
                    print('/Users/ananth/Desktop/figs/extraFec_weakLeraner', ...
                        '-dpng');
                elseif mice(mouse) == 3
                    print('/Users/ananth/Desktop/figs/extraFec_nonLeraner', ...
                        '-dpng');
                elseif mice(mouse) == 2
                    print('/Users/ananth/Desktop/figs/extraFec_strongLeraner', ...
                        '-dpng');
                end
            else
                print('/Users/ananth/Desktop/figs/extraFec_350', ...
                    '-dpng');
            end
            
        end
        
        if saveData == 1
            saveFolder = [saveDirec 'Mouse' mouseName '/' dataset '/'];
            if ~isdir(saveFolder)
                mkdir(saveFolder);
            end
            
            % Save FEC curve
            save([saveFolder 'fec.mat' ], ...
                'eyeClosure', 'fec', ...
                'low_in', 'high_in', 'low_out', 'high_out',...
                'crop', 'fecROI', 'm', 'level',...
                'samplingRate','trialDuration')
        end
        disp([dataset ' analyzed'])
    end
end
toc
beep
disp('All done!')
