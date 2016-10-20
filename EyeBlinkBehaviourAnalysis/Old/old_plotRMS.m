%Plots - RMS Values
clear all
close all

sessionType = 0;
nSessions = 1;
nAnimals = 4;
saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');

markerSize = 10;
fontSize = 15;
lineWidth = 2;
minTrials2Plot = 13;

for animal = 1:nAnimals
    mouseNumber = animal+8;
    mouse = ['MouseK' num2str(mouseNumber)];
    
    for currentSession = 1:nSessions
        dataset = [mouse, '_SessionType', num2str(sessionType), '_Session' num2str(currentSession)];
        disp(dataset);
        saveFolder = [saveDirec mouse '/' dataset '/'];
        
        cs_csPlus_rms(currentSession,1:minTrials2Plot) = csvread([saveFolder, dataset, '_cs_csPlus_rms.csv'], 0, 0, [0 0 (minTrials2Plot-1) 0])';
        preTone_csPlus_rms(currentSession,1:minTrials2Plot) = csvread([saveFolder, dataset, '_preTone_csPlus_rms.csv'], 0, 0, [0 0 0 (minTrials2Plot-1)]);
        cs_csMinus_rms(currentSession,1:minTrials2Plot) = csvread([saveFolder, dataset, '_cs_csMinus_rms.csv'], 0, 0, [0 0 (minTrials2Plot-1) 0])';
        preTone_csMinus_rms(currentSession,1:minTrials2Plot) = csvread([saveFolder, dataset, '_preTone_csMinus_rms.csv'], 0, 0, [0 0 0 (minTrials2Plot-1)]);
        
    end
    
    %Plots
    figure(animal)
    for i = (nSessions):-1:1
        subplot(3,2,1);
        hold on;
        plot(preTone_csPlus_rms(i,1:minTrials2Plot)+((i-1)*10), 'black',...
            'LineWidth', lineWidth,...
            'MarkerSize', markerSize);
        axis([1 minTrials2Plot 0 50]);
%         title([mouse, ': Pre-Tone Phase (CS+) RMS Values'], ...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        xlabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize)
%         ylabel('RMS blinks (fold change)',...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        yt = get(gca, 'YTick');
        set(gca,'YTick',[0 20]);
        set(gca, 'FontSize', fontSize)
        %hold on;
        subplot(3,2,2);
        hold on;
        plot(preTone_csMinus_rms(i,1:minTrials2Plot)+((i-1)*10), 'black',...
            'LineWidth', lineWidth,...
            'MarkerSize', markerSize);
        axis([1 minTrials2Plot 0 50]);
%         title([mouse, ': Pre-Tone Phase (CS-) RMS Values'], ...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        xlabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize)
%         ylabel('RMS blinks (fold change)',...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        yt = get(gca, 'YTick');
        set(gca,'YTick',[0 20]);
        set(gca, 'FontSize', fontSize)
        %hold on;
        subplot(3,2,3);
        hold on;
        plot(cs_csPlus_rms(i,1:minTrials2Plot)+((i-1)*10), 'black',...
            'LineWidth', lineWidth,...
            'MarkerSize', markerSize);
        axis([1 minTrials2Plot 0 50]);
%         title([mouse, ': CS and Trace Phases (CS+) RMS Values'], ...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        xlabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize)
%         ylabel('RMS blinks (fold change)',...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        yt = get(gca, 'YTick');
        set(gca,'YTick',[0 20]);
        set(gca, 'FontSize', fontSize)
        %hold on;
        subplot(3,2,4);
        hold on;
        plot(cs_csMinus_rms(i,1:minTrials2Plot)+((i-1)*10), 'black',...
            'LineWidth', lineWidth,...
            'MarkerSize', markerSize);
        axis([1 minTrials2Plot 0 50]);
%         title([mouse, ': CS and Trace Phases (CS-) RMS Values'], ...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        xlabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize)
%         ylabel('RMS blinks (fold change)',...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        yt = get(gca, 'YTick');
        set(gca,'YTick',[0 20]);
        set(gca, 'FontSize', fontSize)
        
        
        %Differences
        subplot(3,2,5);
        hold on;
        plot((cs_csPlus_rms(i,1:minTrials2Plot)+((i-1)*10))-preTone_csPlus_rms(i,1:minTrials2Plot)+((i-1)*10), 'red',...
            'LineWidth', lineWidth,...
            'MarkerSize', markerSize);
        axis([1 minTrials2Plot -30 50]);
%         title([mouse, ': CS and Trace - Pre-Tone (CS+) RMS Values'], ...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        xlabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize)
%         ylabel('RMS blinks (fold change)',...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        yt = get(gca, 'YTick');
        set(gca,'YTick',[0 30]);
        set(gca, 'FontSize', fontSize)
        %hold on;
        subplot(3,2,6);
        hold on;
        plot((cs_csMinus_rms(i,1:minTrials2Plot)+((i-1)*10))-preTone_csMinus_rms(i,1:minTrials2Plot)+((i-1)*10), 'red',...
            'LineWidth', lineWidth,...
            'MarkerSize', markerSize);
        axis([1 minTrials2Plot -30 50]);
%         title([mouse, ': CS and Trace - Pre-Tone (CS-) RMS Values'], ...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        xlabel('Trials',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize)
%         ylabel('RMS blinks (fold change)',...
%             'FontSize',fontSize,...
%             'FontWeight','bold');
        yt = get(gca, 'YTick');
        set(gca,'YTick',[0 30]);
        set(gca, 'FontSize', fontSize)
    end
    print(['/Users/ananth/Desktop/', mouse, 'TEB'], '-dpng')
end