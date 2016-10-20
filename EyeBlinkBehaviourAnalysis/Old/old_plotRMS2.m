%Plots - RMS Values
clear all
close all

sessionType = 0;
nSessions = 1;
nAnimals = 3;
saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');

markerSize = 10;
fontSize = 15;
lineWidth = 2;
minTrials2Plot = 13;

cs_csPlus_rms = NaN((nSessions),minTrials2Plot);
cs_csMinus_rms = NaN((nSessions),minTrials2Plot);
preTone_csPlus_rms = NaN((nSessions),minTrials2Plot);
preTone_csMinus_rms = NaN((nSessions),minTrials2Plot);

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
    subplot(3,2,1);
    clims = [0 10];
    imagesc(preTone_csPlus_rms, clims)
    colorbar;
    title([mouse ' CS+'],...
        'FontSize', (fontSize+2),...
        'FontWeight', 'bold');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontSize)
    ylabel('Sessions',...
        'FontSize',fontSize,...
        'FontWeight','bold');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', fontSize)
    
    subplot(3,2,2);
    clims = [0 10];
    imagesc(preTone_csMinus_rms, clims)
    colorbar;
    title([mouse ' CS-'],...
        'FontSize', (fontSize+2),...
        'FontWeight', 'bold');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontSize)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', fontSize)
    
    subplot(3,2,3);
    clims = [0 10];
    imagesc(cs_csPlus_rms,clims)
    colorbar;
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontSize)
    ylabel('Sessions',...
        'FontSize',fontSize,...
        'FontWeight','bold');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', fontSize)
    
    subplot(3,2,4);
    clims = [0 10];
    imagesc(cs_csMinus_rms,clims)
    colorbar;
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontSize)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', fontSize)
    
    %Differences
    subplot(3,2,5);
    clims = [-5 5];
    imagesc(cs_csPlus_rms-preTone_csPlus_rms, clims)
    colorbar;
    xlabel('Trials',...
        'FontSize',fontSize,...
        'FontWeight','bold');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontSize)
    ylabel('Sessions',...
        'FontSize',fontSize,...
        'FontWeight','bold');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', fontSize)
    
    subplot(3,2,6);
    clims = [-5 5];
    imagesc(cs_csMinus_rms-preTone_csMinus_rms,clims)
    colorbar;
    xlabel('Trials',...
        'FontSize',fontSize,...
        'FontWeight','bold');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontSize)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', fontSize)
    
    print(['/Users/ananth/Desktop/', mouse, 'Session', num2str(sessionType) 'TEB'], '-dpng')
end
