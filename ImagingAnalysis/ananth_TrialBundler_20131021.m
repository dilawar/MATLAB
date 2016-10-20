%bundler for type specific trials
%Ananthamurthy

ntrials=50;

ncells=118;
nframes=191;
%Remember: indexing is as "(cells, trials, frames)"

a=load('/Users/ananthamurthy/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/calbdf.mat');
b=load('/Users/ananthamurthy/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/caltr.mat');

tone_trials=[1,2,3,4,5,11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45];
puff_trials=[6,7,8,9,10,16,17,18,19,20,26,27,28,29,30,36,37,38,39,40,46,47,48,49,50];

%for Tone
for tone=1:length(tone_trials)
    %tone_session_mouse4_tonenpuff((1:ncells),tone,(1:nframes))=a.calbdf((1:ncells),tone_trials(tone),(1:nframes));
    tone_session_mouse4_tonenpuff((1:ncells),tone,(1:nframes))=b.caltr((1:ncells),tone_trials(tone),(1:nframes));
end

%for Puff
for puff=1:length(puff_trials)
    %puff_session_mouse4_tonenpuff((1:ncells),puff,(1:nframes))=a.calbdf((1:ncells),puff_trials(puff),(1:nframes));
    puff_session_mouse4_tonenpuff((1:ncells),puff,(1:nframes))=b.caltr((1:ncells),puff_trials(puff),(1:nframes));
end

tone_avg=mean(tone_session_mouse4_tonenpuff,2);
puff_avg=mean(puff_session_mouse4_tonenpuff,2);
session_avg=mean(calbdf,2);

%Plotting
for j=1:25
    tone_trials(j)
    imagesc(squeeze(squeeze(tone_session_mouse4_tonenpuff(:,j,:))))
    pause
end
clear j
for j=1:25
    puff_trials(j)
    imagesc(squeeze(squeeze(puff_session_mouse4_tonenpuff(:,j,:))))
    pause
end
clear j

%Plotting - Mean (Trial Average)

imagesc(squeeze(squeeze(100*(tone_avg))));
pause

imagesc(squeeze(squeeze(100*(puff_avg))));
pause

%Plotting - 5 trials at a time

tone5_1=mean(tone_session_mouse4_tonenpuff(:,(1:5),:),2);
imagesc(squeeze(squeeze(tone5_1)))
pause
tone5_2=mean(tone_session_mouse4_tonenpuff(:,(6:10),:),2);
imagesc(squeeze(squeeze(tone5_2)))
pause
tone5_3=mean(tone_session_mouse4_tonenpuff(:,(11:15),:),2);
imagesc(squeeze(squeeze(tone5_3)))
pause
tone5_4=mean(tone_session_mouse4_tonenpuff(:,(16:20),:),2);
imagesc(squeeze(squeeze(tone5_4)))
pause
tone5_5=mean(tone_session_mouse4_tonenpuff(:,(21:25),:),2);
imagesc(squeeze(squeeze(tone5_5)))
pause

puff5_1=mean(puff_session_mouse4_tonenpuff(:,(1:5),:),2);
imagesc(squeeze(squeeze(puff5_1)))
pause
puff5_2=mean(puff_session_mouse4_tonenpuff(:,(6:10),:),2);
imagesc(squeeze(squeeze(puff5_2)))
pause
puff5_3=mean(puff_session_mouse4_tonenpuff(:,(11:15),:),2);
imagesc(squeeze(squeeze(puff5_3)))
pause
puff5_4=mean(puff_session_mouse4_tonenpuff(:,(16:20),:),2);
imagesc(squeeze(squeeze(puff5_4)))
pause
puff5_5=mean(puff_session_mouse4_tonenpuff(:,(21:25),:),2);
imagesc(squeeze(squeeze(puff5_5)))

%Plotting - sets
tone_set1=mean(tone_session_mouse4_tonenpuff(:,(1:5:21),:),2);
imagesc(squeeze(squeeze(tone_set1)))
pause
tone_set2=mean(tone_session_mouse4_tonenpuff(:,(2:5:21),:),2);
imagesc(squeeze(squeeze(tone_set2)))
pause
tone_set3=mean(tone_session_mouse4_tonenpuff(:,(3:5:23),:),2);
imagesc(squeeze(squeeze(tone_set3)))
pause
tone_set4=mean(tone_session_mouse4_tonenpuff(:,(4:5:24),:),2);
imagesc(squeeze(squeeze(tone_set4)))
pause
tone_set5=mean(tone_session_mouse4_tonenpuff(:,(5:5:25),:),2);
imagesc(squeeze(squeeze(tone_set5)))
pause

puff_set1=mean(puff_session_mouse4_tonenpuff(:,(1:5:21),:),2);
imagesc(squeeze(squeeze(puff_set1)))
pause
puff_set2=mean(puff_session_mouse4_tonenpuff(:,(2:5:21),:),2);
imagesc(squeeze(squeeze(puff_set2)))
pause
puff_set3=mean(puff_session_mouse4_tonenpuff(:,(3:5:23),:),2);
imagesc(squeeze(squeeze(puff_set3)))
pause
puff_set4=mean(puff_session_mouse4_tonenpuff(:,(4:5:24),:),2);
imagesc(squeeze(squeeze(puff_set4)))
pause
puff_set5=mean(puff_session_mouse4_tonenpuff(:,(5:5:25),:),2);
imagesc(squeeze(squeeze(puff_set5)))
pause