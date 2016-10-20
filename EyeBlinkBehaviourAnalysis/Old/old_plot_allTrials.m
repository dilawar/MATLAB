close all
clear all

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');
animal = 5;
sessionType = 1;
session = 4;

mouseName = ['MouseK' num2str(animal)];
dataset = [mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
csPlus = csvread([saveDirec, mouseName, '/', dataset, '/', dataset, '_csPlus_filtered.csv']);
csMinus = csvread([saveDirec, mouseName, '/', dataset, '/', dataset, '_csMinus_filtered.csv']);

%figures
%Session 1
figure(1); %CS+
imagesc(csPlus);
title([mouseName ' SessionType ' num2str(sessionType) ' Session' num2str(session) ' CS+']);
xlabel('Time (/10 ms)');
ylabel('Trials');
figure(2); %CS-
imagesc(csMinus);
title([mouseName ' SessionType ' num2str(sessionType) ' Session' num2str(session) ' CS-']);
xlabel('Time (/10 ms)');
ylabel('Trials');