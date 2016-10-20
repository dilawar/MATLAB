%Plots variables from a dataset
clear all
close all

R1 = 0;
C1 = 0;
R2 = 48;
C2 = 1000;
nAnimals = 1;

x = zeros(nAnimals, R2+1, C2+1);
y = zeros(nAnimals, R2+1, C2+1);

for i = 1:1
    i_str = num2str(i+8);
    mouse = ['MouseK' i_str]; %string
    sessionType = '1'; %string
    session = '7'; %string
    
    dataset = [mouse, '_SessionType', sessionType, '_Session' session];
    %saveDirec = ('/home/giulia/Desktop/Work/Analysis/eyeBlinkMatlab/');
    saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis');
    saveFolder = [saveDirec, '/', mouse '/' dataset '/'];
    x(i,:,:) = csvread([saveFolder, dataset, '_csPlus_butterFiltered.csv'], R1, C1, [R1,C1,R2,C2]);
    y(i,:,:) = csvread([saveFolder, dataset, '_csMinus_butterFiltered.csv'], R1, C1, [R1,C1,R2,C2]);
end

%g = [10 20 30 40];
%g = probeTrials; % only works if probeTrials is in memory
g = csvread([saveFolder, dataset, '_probeTrials.csv']);
z = x(:,g,:);

clims = [-30 30];

subplot(2,2,1);
imagesc(squeeze(x), clims);
colorbar;
title([mouse ' SessionType' sessionType ' Session' session ' CS+ Trials']);
xlabel('Time (/10 ms)');
ylabel('CS- Trials');

subplot(2,2,2);
imagesc(squeeze(y), clims);
colorbar;
title([mouse ' SessionType' sessionType ' Session' session ' CS- Trials']);
xlabel('Time (/10 ms)');
ylabel('CS- Trials');

a = mean(x,1); %Trial average of CS+
b = mean(y,1); %Trial average of CS-

subplot(2,2,3);
plot(squeeze(squeeze((median(median(x,1),2)))),'blue');
hold on;
plot(squeeze(squeeze((median(median(y,1),2)))),'red');
hold on;
plot(squeeze(squeeze((median(median(z,1),2)))), 'green');
title([mouse ' SessionType' sessionType ' Session' session ' Trial Medians']);

subplot(2,2,4);
plot(squeeze(squeeze((mean(mean(x,1),2)))),'blue');
hold on;
plot(squeeze(squeeze((mean(mean(y,1),2)))),'red');
hold on;
plot(squeeze(squeeze((mean(mean(z,1),2)))), 'green');
title([mouse ' SessionType' sessionType ' Session' session ' Trial Averaged']);
print([saveFolder, dataset, '_CSplus_plots'], '-dpng');