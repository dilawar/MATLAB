close all
clear all

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');
nAnimals = 4;
nSessions = 3;

%first pass to get the length of the vectors
for animal = 5:1:8
    for session = 1:3
        mouseName = ['MouseK' num2str(animal)];
        dataset = [mouseName '_SessionType1_Session' num2str(session)];
        dataSize(animal, session, :) = length(mean(csvread([saveDirec, mouseName, '/', dataset, '/', dataset, '_csPlus_filtered.csv']),1));
    end
end

csPlus = zeros(nAnimals, nSessions, max(max(dataSize)));
csMinus = zeros(nAnimals, nSessions, max(max(dataSize)));

for animal = 5:1:8
    for session = 1:3
        mouseName = ['MouseK' num2str(animal)];
        dataset = [mouseName '_SessionType1_Session' num2str(session)];
        csPlus (animal, session, :) = mean(csvread([saveDirec, mouseName, '/', dataset, '/', dataset, '_csPlus_filtered.csv']),1);
        csMinus (animal, session, :) = mean(csvread([saveDirec, mouseName, '/', dataset, '/', dataset, '_csMinus_filtered.csv']),1);
    end
end
% Get rid of some
%csPlus(:, (1:13),:) = [];
%csMinus(:, (1:13),:) = [];

% %PLOTS
% for i = 1:4
%     figure(i)
%     
%     for ii = 1:4
%         plot(squeeze(squeeze(csPlus(i, ii, :))), 'blue');
%         hold on;
%         plot(squeeze(squeeze(csMinus(i, ii, :))), 'red');
%         hold on;
%     end
% end
close all
%1. To see differences between animals
figure(1);
for i = 1:4
    subplot(4,1,i)
    plot(squeeze(squeeze(csPlus(i,:,:)))); % Here each subplot represents a different animal
end

figure(2);
for i = 1:4
    subplot(4,1,i)
    plot(squeeze(squeeze(csPlus(:,i,:)))); % Here each subplot represents a different session
end