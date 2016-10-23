% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Treadmill motion analysis
% DEPENDENCIES - First sort all trials as .mat files using sortingVideos.m

tic
%clear all
close all

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions')

%Operations
saveData = 0;
doMotionAnalysis = 0;
plotFigures = 1;

%Video details
samplingRate = 100; % in Frames Per Second (FPS)
trialDuration = 2; % in seconds
nFrames = samplingRate*trialDuration; %per trial
startFrame = 1;


%Dataset details
sessionType = 9;
%mice = [1 2 3 4 5];
mice = 5;
nSessions = 5;
nTrials = 60; % NOTE: During sorting, the dummy trial was excluded

startSession = nSessions; %single sessions
%startSession = 1;
startTrial = 1; % NOTE: During sorting, the dummy trial was excluded

saveDirec = '/Users/ananth/Desktop/Work/Behaviour/Motion';
direc = '/Users/ananth/Desktop/Work/Analysis/Motion/';

fontSize = 12;

motion = zeros(nTrials,nFrames);
for mouse = 1:length(mice)
    mouseName = ['M' num2str(mice(mouse))];
    
    for session = startSession:nSessions
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        if doMotionAnalysis == 1
            disp('Performing motion analysis ...')
        end
        
        disp('... done')
    end
    
    if saveData == 1
        saveFolder = [saveDirec mouseName '/' dataset '/'];
        if ~isdir(saveFolder)
            mkdir(saveFolder);
        end
        
        %Save Motion data
        save([saveFolder 'motion' ],'motion')
    end
    
    if plotFigures == 1
        figure(3)
        imagesc(motion)
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
    end
    disp([dataset ' analyzed'])
end
disp('All done!')
toc
