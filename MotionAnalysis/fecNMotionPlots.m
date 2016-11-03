% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Plot FEC and motion data side by side

clear all
close all

%mice = [1 2 3 4 5];
mice = 5;
sessionType = 9;
nSessions = 1;
startSession = 1;
%saveDirec = '/Users/ananth/Desktop/Work/Analysis/Motion/';
direc = '/Users/ananth/Desktop/Work/Analysis/';

fontSize = 12;

for session = startSession:nSessions
    count = 1;
    for mouse = 1:length(mice)
        mouseName = ['M' num2str(mice(mouse))];
        dataset = ['Mouse' mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(['Working on ' dataset])
        
        try
            % Load fec data
            load([direc 'VideoAnalysis/FEC/Mouse' mouseName '/' dataset '/fec.mat' ])
        catch
            disp(['[ERROR] FEC ' dataset ' not found!'])
            continue
        end
        
        try
            % Load motion data
            load([direc 'MotionAnalysis/Mouse' mouseName '/' dataset '/motion.mat' ])
        catch
            disp(['[ERROR] Motion' dataset ' not found!'])
            continue
        end
        
        % Plot
        figure(1)
        subplot(length(mice),2,count)
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
        
        subplot(length(mice),2,count+1)
        imagesc(motion)
        colormap(hot)
        title(['Treadmill Running ' ...
            mouseName ' ST' num2str(sessionType) ' S' num2str(session)],...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        xlabel(['Time/' num2str(timeLC*1000) ' ms'], ...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        ylabel('Trials', ...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        z = colorbar;
        ylabel(z,'Speed (cm/s)',...
            'FontSize', fontSize,...
            'FontWeight', 'bold')
        
        count = count+2;
    end
end