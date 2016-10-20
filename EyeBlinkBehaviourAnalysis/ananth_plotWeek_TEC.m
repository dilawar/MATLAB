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
nSessions = 1;
startSession = 1;
multipleSets = 0;

if multipleSets == 1
    nSets = 2;
else
    nSets = 1;
end

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');

fontSize = 12;

for session = startSession:nSessions
    close all
    
    for set = 1:nSets
        if set ~= 1
            count = 1;
            count2 = 1;
            if strcmp(mySessionType,'TEC')
                sessionType = 9;
            else
                sessionType = 3;
            end
            Z=38;
        else
            count = 1;
            count2 = 1;
            
            if strcmp(mySessionType,'TEC')
                sessionType = 9;
            else
                sessionType = 3;
            end
            Z=38;
        end
        
        % Specify if probe trials exist
        if sessionType <= 5
            probe = 0;
        else
            probe = 1;
            %probe = 0;
        end
        
        for i = 1:nAnimals
            if set == 1
                mouse = ['K' num2str(i+Z)];
            else
                mouse = ['K' num2str(i+Z)];
            end
            
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
            subplot(nAnimals,nSets*2,count);
            imagesc(blinkData_csPlus);
            title([mouse ' ST' num2str(sessionType) ' S' num2str(session) ' CS+'],...
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
            
            subplot(nAnimals,nSets*2,count+1);
            imagesc(blinkData_csMinus);
            title([mouse ' ST' num2str(sessionType) ' S' num2str(session) ' CS-'],...
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
            %         print(['/Users/ananth/Desktop/' 'rawBlinks_set'...
            %             num2str(set) '_SessionType_' num2str(sessionType)...
            %             '_Session' num2str(session)],...
            %             '-djpeg');
            count = count+nSets*2;
            
            figure(2);
            subplot(nAnimals,nSets,count2)
            if probe == 1
                plot(median(blinkData_probe,1), 'g')
                hold on
            end
            plot(median(blinkData_csPlus,1),'b');
            hold on;
            plot(median(blinkData_csMinus,1),'r');
            
            title([mouse ' ST' num2str(sessionType) ' S' num2str(session)],...
                'FontSize',fontSize,...
                'FontWeight','bold');
            axis([0 1000 0 1024]);
            xlabel('Time/10 ms',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xt = get(gca, 'XTick');
            %set(gca, 'FontSize', fontSize);
            ylabel('Blink (fold change)',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            yt = get(gca, 'YTick');
            %set(gca, 'FontSize', fontSize);
            if probe == 1
                legend('Median Probe Trials', 'Median CS+ Trials','Median CS- Trials');
            else
                legend('Median CS+ Trials','Median CS- Trials');
            end
            %         figure(2);
            %         print(['/Users/ananth/Desktop/' 'median_set'...
            %             num2str(set) '_SessionType_' num2str(sessionType)...
            %             '_Session' num2str(session)],...
            %             '-djpeg');
            count2 = count2+nSets;
        end
    end
    if saveFigures == 1
        figure(1);
        print(['/Users/ananth/Desktop/rawBlinks_' mySessionType ...
            '_week' num2str(week) '_' num2str(session)],...
            '-djpeg');
        figure(2);
        print(['/Users/ananth/Desktop/median_' mySessionType ...
            '_week' num2str(week) '_' num2str(session)],...
            '-djpeg');
    end
end
%pause(5);
beep;
disp('Done!');