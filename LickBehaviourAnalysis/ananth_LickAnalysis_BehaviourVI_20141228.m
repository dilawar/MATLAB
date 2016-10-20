% Purpose: To analyse the Lick Behaviour Training data from Behaviour VI
% Author: Kambadur Ananthamurthy

clear all
close all

%NOTE: Perform Analysis with one block of data at a time
save_direc='/Users/ananth/Desktop/Work/Analysis/LickBehaviourAnalysis';

%The line below makes the program cycle through all datasets. To analyse
%just one, paste its path into a file called folder_list.txt and edit the
%path below to reflect this

list_path='/Users/ananth/Desktop/Work/LickBehaviour/Training/folder_list.txt'; %Check path ID
fid=fopen(list_path);


toneduration=1000;
step=100; %ms

animal=0;
while 1 %You might see a variable "ans" = 0, after the program runs to completion
    tline=fgetl(fid);
    if ~ischar(tline), break, end
    
    direc=[tline '/'];
    
    %General Information
    slashi=strfind(direc, '/');
    dataset=direc(1, (slashi(1, (length(slashi)-1))+1:((length(direc)-1)))); % Dataset name, in the format of MouseNN_BlockXX_sessionY
    
    uscorei=strfind(dataset, '_');
    mouse=dataset(1:uscorei(1,1)-1);
    
    task=dataset((uscorei(1,1)+1):(uscorei(1,2)-1));
    
    block=dataset((uscorei(1,2)+1):(uscorei(1,3)-1));
    
    %block='Block3'; %string
    %blocknumber=3; %number
    %sessioni=findstr(dataset, 'ession'); % just to avoid confusion between "Session" and "session"
    session=dataset((uscorei(1,3)+1):length(dataset));
    sessionnumber_string=session(8:length(session));
    sessionnumber=str2num(sessionnumber_string);
    
    %Read the results.m file for the training session results and save into a matrix "results"
    results_path=[direc mouse '_' task '_' block '_' session '_result.m'];
    
    results=load(results_path);
    %The indexing for the matrix "results" is defined as
    %1. Trial Number
    %2. Pre-Tone loops
    %3. Tone 1 loops
    %4. Tone 2 lick
    %5. Actual Tone 1 duration (ms)
    %6. Actual Tone 2 duration (ms)
    %7. Success Rate (%)
    %8. Lick Dependency
    %9. Correct Trials
    [ntrials,nparams]=size(results); % Number of trials that actually occurred with the number of parameters measured/saved
    %Don't use the length function as this simply calculates the length of the largest dimension (not necessarily the number of trials)
    
    %     TotalCorrectTrials=results(ntrials,nparams);
    %SuccessRate=(TotalCorrectTrials/ntrials)*100; % for all animals and sessions considered
    %NoGoLoops=mean(results(:,3));
    %PreToneLoops=mean(results(:,2));
    
    %     figure(1);
    %     plot(results(:,3),'red');
    %     hold on; plot(results(:,2),'blue');
    %
    %SR=(TotalCorrectTrials/ntrials)*100;
    %Performance=((SR/100)*(1-(NoGoLoops/1000))*(1-(PreToneLoops/1000)))*100;
    
    % Finding perfect trials
    %PerfectTrials(animal,sessionnumber,(1:ntrials))=ismember(NoGoLoops(animal,sessionnumber,(1:ntrials)),1);
    %PerfectSuccessRate(totalanimals,totalsessions)=(sum(PerfectTrials(animal,sessionnumber))/ntrials)*100;
    
    %     %Save Data
    %     if isdir([save_direc '/' mouse])~=1
    %         mkdir([save_direc '/' mouse]);
    %     else
    %     end
    
    %disp(sessionnumber);
    %disp('done!');
    
    %
    %Plots - one dataset at a time!
    %figure(1); plot(results(:,2),'red'); %Tone 1 loops
    %hold on; plot((results(:,6)+300),'black'); % Current Success Rate
    %hold on; plot((results(:,3)*100),'green*'); % Tone 2 licks
    
    %     figure(2);
    %     bar(GoToneLickHist_norm);
    %     for i=1:ntrials
    %         PreToneLickTimes(i)=mean(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']));
    %         NoGoToneLickTimes(i)=mean(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']));
    %     end
    %     clear i
    
    %The problem is that every trial has a different length
    for i=1:ntrials
        PreToneLickTimes_length(i)=size(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']),1);
        %PreToneLickTimes(i,PreToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']));
        
        NoGoToneLickTimes_length(i)=size(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']),1);
        %NoGoToneLickTimes(i,NoGoToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']));
    end
    clear i
    
    %Is MATLAB stupid or am I?
    for i=1:ntrials
        PreToneLickTimes_relative(i,1:PreToneLickTimes_length(i))=dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']);
        PreToneLickTimes_relative(PreToneLickTimes_relative>(max(PreToneLickTimes_relative(i,:))-500))=0; %house-keeping: get rid of timestamps representing the whole duration
        PreToneLickTimes_relative(PreToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
        
        NoGoToneLickTimes_relative(i,1:NoGoToneLickTimes_length(i))=dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']);
        NoGoToneLickTimes_relative(NoGoToneLickTimes_relative>(max(NoGoToneLickTimes_relative(i,:))-500))=0; %house-keeping: get rid of timestamps representing the whole duration
        NoGoToneLickTimes_relative(NoGoToneLickTimes_relative<30)=0; %house-keeping: avoid reading the same lick again
    end
    clear i
    
    PreToneLickTimes=cumsum(PreToneLickTimes_relative,2);
    NoGoToneLickTimes=cumsum(NoGoToneLickTimes_relative,2);
    
    GoToneLickTimes_relative=results(:,6);
    %GoToneLickTimes_relative(GoToneLickTimes_relative>950)=0; %house-keeping: get rid of timestamps representing the whole duration
    %GoToneLickTimes_relative(GoToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
    %GoToneLickTimes=cumsum(GoToneLickTimes_relative);
    
    %histograms
    j=5;
    bin=1;
    for k=step:step:max(max(PreToneLickTimes)) %since this is a matrix
        [m n]=size(PreToneLickTimes(PreToneLickTimes>j & PreToneLickTimes<k));
        PreToneLickHist(bin)=m*n;
        j=k;
        bin=bin+1;
    end
    %disp(['check' num2str(loop)]);
    clear i
    clear k
    clear m
    clear n
    clear bin
    
    j=30;
    bin=1;
    for k=step:step:max(max(NoGoToneLickTimes)) %since this is a matrix
        [m n]=size(NoGoToneLickTimes(NoGoToneLickTimes>j & NoGoToneLickTimes<k));
        NoGoToneLickHist(bin)=m*n;
        j=k;
        bin=bin+1;
    end
    clear i
    clear k
    clear m
    clear n
    clear bin
    
    j=5;
    bin=1;
    for k=step:step:max(GoToneLickTimes_relative) %since this is a vector
        [m n]=size(GoToneLickTimes_relative(GoToneLickTimes_relative>j & GoToneLickTimes_relative<k)); %Since there is no looping for the Go Tone
        GoToneLickHist(bin)=m*n;
        j=k;
        bin=bin+1;
    end
    
    clear i
    clear k
    clear j
    clear m
    clear n
    clear bin
    
    %Normalization
    normPreToneLickHist=PreToneLickHist/max(PreToneLickHist);
    normNoGoToneLickHist=NoGoToneLickHist/max(NoGoToneLickHist);
    normGoToneLickHist=GoToneLickHist/max(GoToneLickHist);
    
    animal=animal+1;
    PreTone_norm=max(max(PreToneLickHist));
    PreToneLickHist_norm(animal,:)=PreToneLickHist/PreTone_norm;
    
    NoGoTone_norm=max(max(NoGoToneLickHist));
    NoGoToneLickHist_norm(animal,:)=NoGoToneLickHist/NoGoTone_norm;
    
    GoTone_norm=max(max(GoToneLickHist));
    GoToneLickHist_norm(animal,:)=GoToneLickHist/GoTone_norm;
    
    all_PreToneLickHist(animal,:)=PreToneLickHist;
    all_NoGoToneLickHist(animal,:)=NoGoToneLickHist;
    all_GoToneLickHist(animal,:)=GoToneLickHist;
    
    %     figure(2);
    %     bar(PreToneLickHist_norm,'blue');
    %     figure(3);
    %     bar(NoGoToneLickHist_norm,'red');
    %     figure(4);
    %     bar(GoToneLickHist_norm,'green');
    %
    %     %Next Dataset
    %     reply = input('Proceed to the next Dataset? y/n ', 's');
    %     if isempty(reply)
    %         reply = 'y';
    %     end
    
end
fclose(fid);
close all

%Sums
% sum_all_PreToneLickHist=sum(all_PreToneLickHist);
% sum_all_NoGoToneLickHist=sum(all_NoGoToneLickHist);
% sum_all_GoToneLickHist=sum(all_GoToneLickHist);
%
% figure(2);
% bar(sum_all_PreToneLickHist,'blue');
% figure(3);
% bar(sum_all_NoGoToneLickHist,'red');
% figure(4);
% bar(sum_all_GoToneLickHist,'green');
%
%
% figure(5);
% bar(PreToneLickHist,'blue');
% figure(6);
% bar(NoGoToneLickHist,'red');
% figure(7);
% bar(GoToneLickHist,'green');

sumPreToneLickHist_norm=sum(PreToneLickHist_norm)/3;
sumNoGoToneLickHist_norm=sum(NoGoToneLickHist_norm)/3;
sumGoToneLickHist_norm=sum(GoToneLickHist_norm)/3;
figure(5);
bar(sumPreToneLickHist_norm,'blue');
figure(6);
bar(sumNoGoToneLickHist_norm,'red');
figure(7);
bar(sumGoToneLickHist_norm,'green');