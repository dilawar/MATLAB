%Plots - Performance Values
clear all
close all

nSessions = 12;
nAnimals = 3;
saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');

csPlus_fullScore = zeros(nAnimals,nSessions+1); % Including No-Puff Controls
csMinus_fullScore = zeros(nAnimals,nSessions+1); % Including No-Puff Controls

for animal = 1:nAnimals
    mouseNumber = animal+8;
    mouse = ['MouseK' num2str(mouseNumber)];
    for sessionType_ind = 1:2;
        sessionType = sessionType_ind-1;
        if (sessionType == 0)
            currentSession = 1;
            dataset = [mouse, '_SessionType', num2str(sessionType), '_Session' num2str(currentSession)];
            disp(dataset);
            saveFolder = [saveDirec mouse '/' dataset '/'];
            
            csPlus_fullScore(animal,currentSession) = csvread([saveFolder, dataset, '_csPlus_fullScore.csv'])*100;
            csMinus_fullScore(animal,currentSession) = csvread([saveFolder, dataset, '_csMinus_fullScore.csv'])*100;
        else
            for currentSession = 2:nSessions
                dataset = [mouse, '_SessionType', num2str(sessionType), '_Session' num2str(currentSession)];
                disp(dataset);
                saveFolder = [saveDirec mouse '/' dataset '/'];
                csPlus_fullScore(animal,currentSession) = csvread([saveFolder, dataset, '_csPlus_fullScore.csv'])*100;
                csMinus_fullScore(animal, currentSession) = csvread([saveFolder, dataset, '_csMinus_fullScore.csv'])*100;
            end
        end
    end
    
    %All animals
    subplot(1,2,1)
    hold on
    plot(csPlus_fullScore(animal,:), 'blue',...
        'LineWidth', 1,...
        'MarkerSize', 10);
    title('CS+ Trials', ...
        'FontSize',20,...
        'FontWeight','bold');
    axis([0 (nSessions+2) -50 50]);
    xlabel('Sessions',...
        'FontSize',20,...
        'FontWeight','bold');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 20)
    ylabel('Performance/%',...
        'FontSize',20,...
        'FontWeight','bold');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 20)
    
    subplot(1,2,2)
    hold on
    plot(csMinus_fullScore(animal,:), 'red',...
        'LineWidth', 1,...
        'MarkerSize', 10);
    axis([0 (nSessions+2) -50 50]);
    title('CS- Trials', ...
        'FontSize',20,...
        'FontWeight','bold');
    xlabel('Sessions',...
        'FontSize',20,...
        'FontWeight','bold');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 20)
    ylabel('Performance/%',...
        'FontSize',20,...
        'FontWeight','bold');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 20)
    
end

%Mean

subplot(1,2,1)
hold on
errorbar(mean(csPlus_fullScore,1),std(csPlus_fullScore,1), 'blue',...
    'LineWidth', 3,...
    'MarkerSize', 10);
title('CS+ Trials', ...
    'FontSize',20,...
    'FontWeight','bold');
axis([0 (nSessions+2) -50 50]);
xlabel('Sessions',...
    'FontSize',20,...
    'FontWeight','bold');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
ylabel('Performance/%',...
    'FontSize',20,...
    'FontWeight','bold');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)

subplot(1,2,2)
hold on
errorbar(mean(csMinus_fullScore,1),std(csPlus_fullScore,1), 'red',...
    'LineWidth', 3,...
    'MarkerSize', 10);
axis([0 (nSessions+2) -50 50]);
title('CS- Trials', ...
    'FontSize',20,...
    'FontWeight','bold');
xlabel('Sessions',...
    'FontSize',20,...
    'FontWeight','bold');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
ylabel('Performance/%',...
    'FontSize',20,...
    'FontWeight','bold');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)

%print(['/Users/ananth/Desktop/', mouse, 'Performance' ], '-dpng')