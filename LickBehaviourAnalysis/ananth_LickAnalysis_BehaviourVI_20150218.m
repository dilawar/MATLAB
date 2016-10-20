% Purpose: To analyse the Lick Behaviour Training data from Behaviour VI
% Author: Kambadur Ananthamurthy

clear all
close all

%NOTE: Perform Analysis with one block of data at a time
save_direc='/home/rai/Desktop/NCBS/Behaviour/Analyze';

%The line below makes the program cycle through all datasets. To analyse
%just one, paste its path into a file called folder_list.txt and edit the
%path below to reflect this

list_path='/home/rai/Desktop/NCBS/Behaviour/Data/folder_list.txt'; %Check path ID
fid=fopen(list_path);

totalanimals=1;  %*! Must be edited for every experiment
totalsessions=20; %*! Must be edited for every animal
totaltrials=202; %*! Must be edited for every animal
TotalCorrectTrials=zeros(totalanimals,totalsessions,totaltrials);
SuccessRate=zeros(totalanimals,totalsessions,totaltrials);
CorrectTrials=zeros(totalanimals,totalsessions,totaltrials);
Tone1Loops=zeros(totalanimals,totalsessions,totaltrials);
PerfectTrials=zeros(totalanimals,totalsessions,totaltrials);
PerfectSuccessRate=zeros(totalanimals,totalsessions,totaltrials);

maxExpectedIncorrectLicks=1000;
PTlicks=0;
NGlicks=0;
Glicks=0;

toneduration=1000; % in ms
binsize=50; %in ms
GoToneLickHist=zeros((toneduration/binsize),1);
ExpectedLickDuration=2; %in ms

mouse_check='MouseS2'; %First animal's name
animal=0;

%Program Modules
multiDataset=0;
readProtocol=0;

while 1 %You might see a variable "ans" = 0, after the program runs to completion
    tline= fgetl(fid);
    if ~ischar(tline), break, end
    
    direc=[tline '/'];
    
    %General Information
    slashi=strfind(direc, '/');
    dataset=direc(1, (slashi(1, (length(slashi)-1))+1:((length(direc)-1)))); % Dataset name, in the format of MouseNN_BlockXX_sessionY
    
    uscorei=strfind(dataset, '_');
    mouse=dataset(1:uscorei(1,1)-1);
    
    if multiDataset==1;
        if ismember(mouse==mouse_check,0)
            animal=animal+1;
        else
            animal=1;
        end
        mouse_check=mouse; %assigning for a new animal
    else
        animal=1;
    end
    
    task=dataset((uscorei(1,1)+1):(uscorei(1,2)-1));
    
    block=dataset((uscorei(1,2)+1):(uscorei(1,3)-1));
    
    %block='Block3'; %string
    %blocknumber=3; %number
    
    %sessioni=findstr(dataset, 'ession'); % just to avoid confusion between "Session" and "session"
    session=dataset((uscorei(1,3)+1):length(dataset));
    sessionnumber_string=session(8:length(session));
    sessionnumber=str2num(sessionnumber_string);
    
    %Read Protocol information from a text file labelled protocol.txt
    protocol_path=[direc mouse '_' task '_' block '_' session '_protocol.txt'];
    fid_protocol=fopen(protocol_path);
    protocol=zeros(1,14); %I have used 14 since there are 14 parameters stored in the text file protocol.txt, at this time
    for i=1:14
        temp=fgetl(fid_protocol); %Reads the first/next line in the .txt file; replaces the text from previous line
        if temp==-1
            
        else
            a=strfind(temp, ' - '); %Finds " - "
            
            if isempty(a)
                %skipping headings
            else
                value=temp((a+3):length(temp));
                if i>2
                   protocol(1,i)=str2num(value); % everything from after " - " to the end of the line
                end
                
                % The indexing for the array "protocol" is defined as
%                1. Task 
%                2. Block
%                3. Session 
%                4. Criterion
%                5. Pre-Tone Duration (ms) 
%                6. CS A (Hz) 
%                7. CS B (Hz) 
%                8. CS 1 Duration Max (s) 
%                9. CS 1 Duration Min (s) 
%                10. Delay Interval Duration (ms) 
%                11. CS 2 Duration (ms) 
%                12. Water Time (ms) 
%                13. Max ITI (s) 
%                14. Min ITI (s) 
               
                clear temp
                clear a
            end
        end
    end
    clear i
    fclose(fid_protocol);
    
    %Read the results.m file for the training session results and save into a matrix "results"
    results_path=[direc mouse '_' task '_' block '_' session '_result.m'];
    
    results=load(results_path);
    
    %The indexing for the matrix "results" is defined as
    %     1. Trial Number
    T=1;
    %     2. Pre-Tone Loops
    PTL=2;
    %     3. Tone 1 loops
    T1L=3;
    %     4. Tone 2 lick
    T2L=4;
    %     5. Actual Tone 1 duration (ms)
    AT1D=5;
    %     6. Actual Tone 2 duration (ms)
    AT2D=6;
    %     7. Success Rate (%)
    CSR=7; % Current SR
    %     8. Lick Dependency
    LD=8;
    %     9. Correct Trials
    CT=9;
    
    [ntrials,nparams]=size(results); % Number of trials that actually occurred with the number of parameters measured/saved
    %Don't use the length function as this simply calculates the length of the largest dimension (not necessarily the number of trials)
    
    %     if ntrials~=600;
    %         disp('Check ntrials');
    %     end

    %     %Save Data
    %     if isdir([save_direc '/' mouse])~=1
    %         mkdir([save_direc '/' mouse]);
    %     else
    %     end
    
    %disp(sessionnumber);
    %disp('done!');
    %The problem is that every trial has a different length
    for i=1:ntrials
        PreToneLickTimes_length(i)=size(dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']),1);
        %PreToneLickTimes(i,PreToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']));
        
        NoGoToneLickTimes_length(i)=size(dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']),1);
        %NoGoToneLickTimes(i,NoGoToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']));
    end
    clear i
    
    %distinctPreToneLicks=zeros(ntrials,maxExpectedIncorrectLicks);
    %distinctNoGoToneLicks=zeros(ntrials,maxExpectedIncorrectLicks);
    PreToneLicks=zeros(ntrials,1);
    NoGoToneLicks=zeros(ntrials,1);
    
    %Is MATLAB stupid or am I?
    for i=1:ntrials
        PreToneLickTimes_relative(i,1:PreToneLickTimes_length(i))=dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']);
        PreToneLickTimes_relative(PreToneLickTimes_relative>(max(PreToneLickTimes_relative(i,:))-500))=0; %house-keeping: get rid of timestamps representing the whole duration
        PreToneLickTimes_relative(PreToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
        %distinctPreToneLicks=zeros(size(PreToneLickTimes_relative));
              
        
        count=0;
        for gj=1:PreToneLickTimes_length(i)
            if PreToneLickTimes_relative(i,gj)>0
                count=count+1;
            else
            end
            if count>=1
               PreToneLicks(i)=PreToneLicks(i)+1;
            else
            end
            count=0;
            
        end
        PTlicks= PTlicks+PreToneLicks(i);
        
        clear gj
        clear count
%         distinctPreToneLicks(PreToneLickTimes_relative>ExpectedLickDuration)=1;
%         nPreToneLicks(i)=sum(distinctPreToneLicks(i));
        
        NoGoToneLickTimes_relative(i,1:NoGoToneLickTimes_length(i))=dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']);
        NoGoToneLickTimes_relative(NoGoToneLickTimes_relative>(max(NoGoToneLickTimes_relative(i,:))-500))=0; %house-keeping: get rid of timestamps representing the whole duration
        NoGoToneLickTimes_relative(NoGoToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
        %distinctNoGoToneLicks=zeros(size(NoGoToneLickTimes_relative));
%         distinctNoGoToneLicks(NoGoToneLickTimes_relative>ExpectedLickDuration)=1;
%         nNoGoToneLicks(i)=sum(distinctNoGoToneLicks(i));
        count=0;
        for gj=1:NoGoToneLickTimes_length(i)
            if NoGoToneLickTimes_relative(i,gj)>0
                count=count+1;
            else
                count=0;
            end
            if count>=1
                NoGoToneLicks(i)=PreToneLicks(i)+1;
            else
            end
            count=0;
        end
        
        if results(i,T2L)>0
            factor=1;
        else
            factor=0;
        end
        
        %Finding Performance for each trial
        Performance(i,1) = (100/(1+NoGoToneLicks(i)+PreToneLicks(i)))*factor;
        
        NGlicks= NGlicks+NoGoToneLicks(i);
        Glicks=results(ntrials,nparams);
        AllLicks= [PTlicks,NGlicks,Glicks];
       
    end
    clear i
    
    %Finding Performance for whole session
    SessionPerformance(protocol(1,3),1) = sum(Performance)/ntrials
    
    % "The cumulative sum of relative time is absolute time" - Anonymous
    PreToneLickTimes=cumsum(PreToneLickTimes_relative,2);
    NoGoToneLickTimes=cumsum(NoGoToneLickTimes_relative,2);
    
    GoToneLickTimes_relative=results(:,AT2D);
    %GoToneLickTimes_relative(GoToneLickTimes_relative>950)=0; %house-keeping: get rid of timestamps representing the whole duration
    %GoToneLickTimes_relative(GoToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
    %GoToneLickTimes=cumsum(GoToneLickTimes_relative);
    
  
    
    
    TotalCorrectTrials(animal,sessionnumber)=results(ntrials,nparams);
    SuccessRate(animal,sessionnumber)=(TotalCorrectTrials(animal,sessionnumber)/ntrials)*100; % for all animals and sessions considered
    
    %TotalCorrectTrials=results(ntrials,nparams); %Since the last column is saved with the sum of the correct trials
    %SuccessRate=(TotalCorrectTrials/ntrials)*100; % for all animals and sessions considered
    CorrectTrials(animal,sessionnumber,(1:ntrials))=results(:,CT);
    
    Tone1Loops(animal,sessionnumber,(1:ntrials))=results(:,3); %Look at indexing
    Tone2Licks_binary(animal,sessionnumber,(1:ntrials))=results(:,T2L);
    
    % Finding perfect trials
    for i=1:ntrials
        if ismember(Tone1Loops(animal,sessionnumber,i),1) && ismember(Tone2Licks_binary(animal,sessionnumber,i),1)
            PerfectTrials(animal,sessionnumber,i)=1;
        else
            PerfectTrials(animal,sessionnumber,i)=0;
        end
    end
    clear i
    PerfectSuccessRate(animal,sessionnumber)=(sum(PerfectTrials(animal,sessionnumber,:))/ntrials)*100;
    
   
    %histograms
    gj=5;
    bin=1;
    for k=binsize:binsize:max(max(PreToneLickTimes)) %since this is a matrix
        [m n]=size(PreToneLickTimes(PreToneLickTimes>gj & PreToneLickTimes<k));
        PreToneLickHist(bin)=m*n;
        gj=k;
        bin=bin+1;
    end
    %disp(['check' num2str(loop)]);
    clear i
    clear k
    clear m
    clear n
    clear bin
    
    gj=30;
    bin=1;
    for k=binsize:binsize:max(max(NoGoToneLickTimes)) %since this is a matrix
        [m n]=size(NoGoToneLickTimes(NoGoToneLickTimes>gj & NoGoToneLickTimes<k));
        NoGoToneLickHist(bin)=m*n;
        gj=k;
        bin=bin+1;
    end
    clear i
    clear k
    clear m
    clear n
    clear bin
    
    gj=5;
    bin=1;
    for k=binsize:binsize:max(GoToneLickTimes_relative) %since this is a vector
        [m n]=size(GoToneLickTimes_relative(GoToneLickTimes_relative>gj & GoToneLickTimes_relative<k)); %Since there is no looping for the Go Tone
        GoToneLickHist(bin)=m*n;
        gj=k;
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

    PreTone_norm=max(max(PreToneLickHist));
    PreToneLickHist_norm(animal,:)=PreToneLickHist/PreTone_norm;
    
    NoGoTone_norm=max(max(NoGoToneLickHist));
    NoGoToneLickHist_norm(animal,:)=NoGoToneLickHist/NoGoTone_norm;
    
    GoTone_norm=max(max(GoToneLickHist));
    GoToneLickHist_norm(animal,:)=GoToneLickHist/GoTone_norm;
    
    all_PreToneLickHist(animal,:)=PreToneLickHist;
    all_NoGoToneLickHist(animal,:)=NoGoToneLickHist;
    all_GoToneLickHist(animal,:)=GoToneLickHist;
    
     %Plots - one dataset at a time!
%     figure(1); plot(results(:,T1L),'blue'); %Tone 1 loops
%     %hold on; plot((results(:,CSR)+300),'black'); % Current Success Rate
%     %%%%%
%     hold on; plot((results(:,T2L)*100),'green*'); % Tone 2 licks -yes/no
%     
    figure(2);
    bar(PreToneLickHist_norm,'blue');
    title('Pre-Tone Lick Histogram');
    figure(3);
    bar(NoGoToneLickHist_norm,'red');
    title('No Go Tone Lick Histogram');
    figure(4);
    bar(GoToneLickHist_norm,'green');
    title('Go Tone Lick Histogram');
    figure(5);
    bar(AllLicks);
    title('Comparison over three phases');
    ylabel('No. of licks');
    xlabel('Phases');
    
%     figure(6)
%     plot(SessionPerformance)
    
    %Next Dataset
    reply = input('Proceed to the next Dataset? y/n ', 's');
    if isempty(reply)
        reply = 'y';
    end
    
end
fclose(fid);

% %Save
% save([save_direc '/' mouse '/TotalCorrectTrials.txt'], 'TotalCorrectTrials', '-ASCII'); %.txt
% save([save_direc '/' mouse '/TotalCorrectTrials'], 'TotalCorrectTrials'); %.mat
% xlswrite([save_direc '/' mouse '/TotalCorrectTrials.xls'], TotalCorrectTrials); %.xls
% save([save_direc '/' mouse '/SuccessRate.txt'], 'SuccessRate', '-ASCII'); %.txt
% save([save_direc '/' mouse '/SuccessRate'], 'SuccessRate'); %.mat
% xlswrite([save_direc '/' mouse '/SuccessRate.xls'], SuccessRate); %.xls
