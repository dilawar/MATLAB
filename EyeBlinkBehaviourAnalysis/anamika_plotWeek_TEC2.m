% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Use this to plot the week's behaviour data

addpath('/Users/ananth/Documents/MATLAB/CustomFunctions_MATLAB')

clear all
close all

saveFigures = 1;

%mySessionType = 'NPControl';
mySessionType = 'TEC';

week = 1;
nAnimals = 3;
nSessions = 9;
startSession = 9;

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/Analysis4Anamika/');

fontSize = 12;

for session = startSession:nSessions
    close all
    
    count = 1;
    count2 = 1;
    
    if strcmp(mySessionType,'TEC')
        sessionType = 9;
    else
        sessionType = 3;
    end
    Z=38;
    
    % Specify if probe trials exist
    if sessionType <= 5
        probe = 0;
    else
        probe = 1;
        %probe = 0;
    end
    
    for i = 1:nAnimals
        mouse = ['K' num2str(i+Z)];
        
        
        %session = 1;
        dataset = [mouse '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        disp(dataset);
        saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
        
        blinkData_csPlus = csvread([saveFolder 'Mouse' dataset '_csPlus.csv']);
        blinkData_csMinus = csvread([saveFolder 'Mouse' dataset '_csMinus.csv']);
        if probe == 1
            probeTrials = csvread([saveFolder 'Mouse' dataset '_probeTrials.csv']);
            blinkData_probe = blinkData_csPlus(probeTrials,:);
        end
        
        figure(1);
        subplot(nAnimals,1,count);
        imagesc(blinkData_csPlus);
        title(['Anamika ' mouse ' ST' num2str(sessionType) ' S' num2str(session) ' Tone'],...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xlabel('Time/10 ms',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        %set(gca, 'FontSize', fontSize);
        ylabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        yt = get(gca, 'YTick');
        %set(gca, 'FontSize', fontSize);
        caxis([0, 1024]);
        colorbar;
        colormap(jet);
        
        count = count+1;
        
        figure(2);
        subplot(nAnimals,1,count2)
        if probe == 1
            plot(median(blinkData_probe,1), 'g')
            hold on
        end
        plot(median(blinkData_csPlus,1),'b');
        title(['Anamika ' mouse ' ST' num2str(sessionType) ' S' num2str(session) ' Tone'],...
            'FontSize',fontSize,...
            'FontWeight','bold');
        axis([0 1000 0 1024]);
        xlabel('Time/10 ms',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        %set(gca, 'FontSize', fontSize);
        ylabel('Response',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        yt = get(gca, 'YTick');
        %set(gca, 'FontSize', fontSize);
        if probe == 1
            legend('Median Probe Trials', 'Median Tone Trials');
        else
            legend('Median Tone Trials');
        end
        %         figure(2);
        %         print(['/Users/ananth/Desktop/' 'median_set'...
        %             num2str(set) '_SessionType_' num2str(sessionType)...
        %             '_Session' num2str(session)],...
        %             '-djpeg');
        count2 = count2+1;
    end
    if saveFigures == 1
        figure(1);
        print(['/Users/ananth/Desktop/anamika_rawBlinks_' mySessionType ...
            '_week' num2str(week) '_' num2str(session)],...
            '-djpeg');
        figure(2);
        print(['/Users/ananth/Desktop/anamika_median_' mySessionType ...
            '_week' num2str(week) '_' num2str(session)],...
            '-djpeg');
    end
    
end
%pause(5);
beep;
disp('Done!');