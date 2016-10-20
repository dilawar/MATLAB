% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Custom plots

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions_MATLAB')

clear all
close all

saveFigures = 1;

month = 'April';
csType = 'csMinus';
Z=25;
nAnimals = 3;
totalTEC = 6;
ymin = -3000;
ymax = 700;

sessionType = 3; %NPC
sessionType = 9; %TEC

if sessionType == 3
    nSessions = 1;
else
    nSessions = totalTEC;
end
startSession = 1;

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');

fontSize = 12;

for animal = 1:nAnimals
    count = 1+(animal-1)*nSessions;
    for session = startSession:nSessions
        %close all
        mouse = ['K' num2str(animal+Z)];
        
        %session = 1;
        dataset = [mouse '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        %disp(dataset);
        saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
        
        %blinkData_csPlus = csvread([saveFolder 'Mouse' dataset '_csPlus.csv']);
        blinkData = csvread([saveFolder 'Mouse' dataset '_' csType '.csv']);
        %disp(count);
        
        figure(1);
        subplot(nAnimals,nSessions,count);
        [m,n] = size(blinkData);
        
        
        for trial = 1:m
            plot(blinkData(trial,400:600)-((trial-1)*100),'r');
            hold on
        end
        
        if sessionType == 3
            title([mouse ' NPC'],...
                'FontSize',fontSize,...
                'FontWeight','bold');
        else
            title([mouse ' TEC' ' S' num2str(session)],...
                'FontSize',fontSize,...
                'FontWeight','bold');
        end
        
        %axis([0 200 0 1024]);
        xlim([0 200]);
        ylim([ymin ymax]);
        ax = gca;
        %ax.XTick = [0 100 200 300 400];
        ax.XTick = [0 100 200];
        ax.YTick = [0 200];
        ax.XTickLabel = {'-1s','Tone', '+1s'};
        ax.YTickLabel = {'0 AU', '200 AU'};
        %         xlabel('Time (ms)',...
        %             'FontSize',fontSize,...
        %             'FontWeight','bold');
        %xt = get(gca, 'XTick');
        %set(gca, 'FontSize', fontSize);
            ylabel('Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
        yt = get(gca, 'YTick');
        set(gca, 'FontSize', fontSize);
        hold on
        count = count + 1;
    end
end
if saveFigures == 1
    figure(1);
    print(['/Users/ananth/Desktop/tone/' month '_sessionType' num2str(sessionType) '_' csType],...
        '-djpeg');
end

%pause(5);
beep;
disp('Done!');