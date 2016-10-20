clear all
close all

% Execution Modules
readProtocol=0; %to read protocol informaiton (1), or not (0)
makeHistograms=0; %to make lick histograms (1), or not (0)
saveData=0; %to save data (1), or not (0)
plotData=0; %to plot data (1), or not (0)
%userControlNext=1; %to require user input to go to next dataset (1), or not (0)
newDatasets=1; %Ture (1); False (0): in the old datasets (newDatsets=0), the Pre-Tone and No Go lick times are saved without the task and block details in the filename

save_direc='/Users/ananth/Desktop/Work/Analysis/LickBehaviourAnalysis';

list_path='/Users/ananth/Desktop/Work/LickBehaviour/Training/folder_list.txt'; %Check path ID
fid=fopen(list_path);

% weight_path='/Users/ananth/Desktop/Work/LickBehaviour/Training/weight_list.txt';
% fid_weight=fopen(weight_path);

totalanimals=1;  %*! Must be edited for every experiment
totalsessions=1; %*! Must be edited for every animal
totaltrials=300; %*! Must be edited for every animal
%TotalCorrectTrials=zeros(totalanimals,totalsessions,totaltrials);
%SuccessRate=zeros(totalanimals,totalsessions);
%CorrectTrials=zeros(totalanimals,totalsessions,totaltrials);
%Tone1Loops=zeros(totalanimals,totalsessions,totaltrials);
%PerfectTrials=zeros(totalanimals,totalsessions,totaltrials);
%PerfectSuccessRate=zeros(totalanimals,totalsessions);
PTlicks=zeros(totaltrials,1);
NGlicks=zeros(totaltrials,1);
Glicks=zeros(totaltrials,1);

%SessionPerformance=zeros(totalanimals,totalsessions);

PreToneDuration=1000; % in ms
MinNoGoToneDuration=1000; %in ms
ToneDuration=1000; % in ms
binsize=50; %in ms
GoToneLickHist=zeros((ToneDuration/binsize),1);

firstmouse='MouseS2'; %!@ edit for every run of this code
animal=1;

while 1 %the loop basically runs to infinity unless asked to break and end
    tline= fgetl(fid);
    if ~ischar(tline), break, end
    
    direc=[tline '/'];
    
    %General Information
    slashi=strfind(direc, '/');
    dataset=direc(1, (slashi(1, (length(slashi)-1))+1:((length(direc)-1)))); % Dataset name, in the format of MouseNN_BlockXX_sessionY
    
    uscorei=strfind(dataset, '_');
    mouse=dataset(1:uscorei(1,1)-1);
    
    %criteria to check if the current mouse count is the same or different
    %     if sum(mouse==firstmouse)/(length(mouse))==1
    %         %do nothing
    %     else
    %         animal=animal+1;
    %         firstmouse=mouse;
    %     end
    
    task=dataset((uscorei(1,1)+1):(uscorei(1,2)-1));
    block=dataset((uscorei(1,2)+1):(uscorei(1,3)-1));
    
    %sessioni=findstr(dataset, 'ession'); % just to avomouse_check='MouseS2'; %First animal's nameid confusion between "Session" and "session"
    session=dataset((uscorei(1,3)+1):length(dataset));
    sessionnumber_string=session(8:length(session));
    sessionnumber=str2num(sessionnumber_string);
    
    disp(['Animal:', num2str(animal) ' Session:', int2str(sessionnumber)])
    
    %Read Protocol information from a text file labelled protocol.txt
    if readProtocol==1
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
    else
    end
    
    %         %Reading the weights
    %         weight=fgetl(fid_weight)
    %         WeightChange(sessionnumber,1)=str2double(weight)
    
    %Read the results.m file for the training session results and save into a matrix "results"
    results_path=[direc mouse '_' task '_' block '_' session '_result.m'];
    results=load(results_path);
    
    %The indexing for the matrix "results" is defined as
    T=1;       %     1. Trial Number
    PTL=2;     %     2. Pre-Tone Loops
    T1L=3;     %     3. Tone 1 loops
    T2L=4;     %     4. Tone 2 lick
    AT1D=5;    %     5. Actual Tone 1 duration (ms)
    AT2D=6;    %     6. Actual Tone 2 duration (ms)
    CSR=7;     %     7. Current Success Rate (%)
    LD=8;      %     8. Lick Dependency
    CT=9;      %     9. Correct Trials
    
    [ntrials,nparams]=size(results); % Number of trials that actually occurred with the number of parameters measured/saved
    %Don't use the length function as this simply calculates the length of the largest dimension (not necessarily the number of trials)
    
    if ntrials<200;
        disp(['Check Animal:', num2str(animal) ' Session:', int2str(sessionnumber) ' "ntrials"'])
    else
    end
    
    %To make the lick frequency histograms
    
    %The problem is that every trial has a different length
    if newDatasets==1 %normal case
        for i=1:ntrials
            PreToneLickTimes_length(i)=size(dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']),1);
            %PreToneLickTimes(i,PreToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']));
            
            NoGoToneLickTimes_length(i)=size(dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']),1);
            %NoGoToneLickTimes(i,NoGoToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']));
        end
    else %the filenames for Pre-Tone and No Go lick times do not have the task and block details
        for i=1:ntrials
            PreToneLickTimes_length(i)=size(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']),1);
            %PreToneLickTimes(i,PreToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']));
            
            NoGoToneLickTimes_length(i)=size(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']),1);
            %NoGoToneLickTimes(i,NoGoToneLickTimes_length(i))=cumsum(dlmread([direc mouse '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']));
        end
    end
    clear i
    
    PreToneLicks=zeros(ntrials,1);
    NoGoToneLicks=zeros(ntrials,1);
    TrialPerformance=zeros(ntrials,1);
    
    %Is MATLAB stupid or am I?
    
    for i=1:ntrials
        
        PreToneLickTimes_relative(i,1:PreToneLickTimes_length(i))=dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_Pre-ToneTimes.xls']);
        PreToneLickTimes_relative(PreToneLickTimes_relative>(max(PreToneLickTimes_relative(i,:))-500))=0; %house-keeping: get rid of timestamps representing the whole duration
        PreToneLickTimes_relative(PreToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
        
        for gj=1:PreToneLickTimes_length(i)
            if PreToneLickTimes_relative(i,gj)>10
                PreToneLicks(i)=PreToneLicks(i)+1;
            end
        end
        clear gj
        PTlicks(i)= PTlicks(i)+PreToneLicks(i);
        
        NoGoToneLickTimes_relative(i,1:NoGoToneLickTimes_length(i))=dlmread([direc mouse '_' task '_' block '_' session '_' 'Trial' num2str(i) '_No-GoTimes.xls']);
        NoGoToneLickTimes_relative(NoGoToneLickTimes_relative>(max(NoGoToneLickTimes_relative(i,:))-500))=0; %house-keeping: get rid of timestamps representing the whole duration
        NoGoToneLickTimes_relative(NoGoToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
        
        for gj=1:NoGoToneLickTimes_length(i)
            if NoGoToneLickTimes_relative(i,gj)>10
                NoGoToneLicks(i)=NoGoToneLicks(i)+1;
            end
        end
        clear gj
        NGlicks(i)= NGlicks(i)+NoGoToneLicks(i);
        
        %Glicks=results(ntrials,nparams);
        %AllLicks(sessionnumber,(1:3))= [PTlicks,NGlicks,Glicks];
        %Establishing Performance for each trial
        
        %TrialPerformance(i) = (1/(1+NoGoToneLicks(i)+PreToneLicks(i)))*results(i,T2L)*100; %in %
        %TrialPerformance(i)=(1-(PreToneLoops(animal,sessionnumber,i)/pretoneduration))*(1-(NoGoToneLoops(animal,sessionnumber,i)/minnogotoneduration))*(GoToneLoops(animal,sessionnumber,i))*100;
    end
    
    clear i
    
    %Number of Licks
    PreToneLicks=sum(PTlicks)
    NoGoToneLicks=sum(NGlicks)
    
    %Finding Performance for whole session
    %SessionPerformance(animal,sessionnumber) = sum(TrialPerformance)/ntrials;
    PreToneLoops=mean(results(:,PTL)) %to get an average
    NoGoToneLoops=mean(results(:,T1L)) %to get an average
    GoToneLicks=sum(results(:,T2L)) %to get the total hit trials (Lick during Go Tones)
    
    TotalCorrectTrials=results(ntrials,CT); %Since the last column is saved with the sum of the correct trials
    SuccessRate=TotalCorrectTrials(animal,sessionnumber)/ntrials;
    SessionPerformance=(1-(PreToneLoops/PreToneDuration))*(1-(NoGoToneLoops/MinNoGoToneDuration))*SuccessRate*100
    
    
    % "The cumulative sum of relative time is absolute time" - Anonymous
    % PreToneLickTimes=cumsum(PreToneLickTimes_relative,2);
    % NoGoToneLickTimes=cumsum(NoGoToneLickTimes_relative,2);
    % GoToneLickTimes_relative=results(:,AT2D);
    %GoToneLickTimes_relative(GoToneLickTimes_relative>950)=0; %house-keeping: get rid of timestamps representing the whole duration
    %GoToneLickTimes_relative(GoToneLickTimes_relative<3)=0; %house-keeping: avoid reading the same lick again
    %GoToneLickTimes=cumsum(GoToneLickTimes_relative);
    
    %Finding perfect trials
    PerfectTrials=0;
    for i=1:ntrials
        if ismember(results(i,T1L),1) && ismember(results(i,T2L),1)
            PerfectTrials=PerfectTrials+1;
        else
        end
    end
    clear i
    PerfectTrials
    % PerfectSuccessRate(animal,sessionnumber)=(sum(PerfectTrials(animal,sessionnumber,:))/ntrials)*100;
    
    %histograms
    if makeHistograms==1
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
        
        
        %Plots - one dataset at a time!
        %     figure(1); plot(results(:,T1L),'blue'); %Tone 1 loops
        %     hold on; plot((results(:,CSR)+300),'black'); % Current Success Rate
        %     hold on; plot((results(:,T2L)*100),'green*'); % Tone 2 licks -yes/no
        
        %         figure(2);
        %         bar(PreToneLickHist,'blue');
        %         title('Pre-Tone Lick Histogram');
        %         figure(3);
        %         bar(NoGoToneLickHist,'red');
        %         title('No Go Tone Lick Histogram');
        %         figure(4);
        %         bar(GoToneLickHist,'green');
        %         title('Go Tone Lick Histogram');
        
    else
    end
    
    if plotData==1
        if sessionnumber==1
        else
            figure(2);
            bar(PreToneLickHist,'blue');
            title('Pre-Tone Lick Histogram');
            figure(3);
            bar(NoGoToneLickHist,'red');
            title('No Go Tone Lick Histogram');
            figure(4);
            bar(GoToneLickHist,'green');
            title('Go Tone Lick Histogram');
            
            %         figure(5);
            %         plot(PreToneLoops,'blue', 'LineWidth', 2);
            %         hold on; plot(NoGoToneLoops,'red', 'LineWidth', 2);
            %         hold on; plot(GoToneLoops,'green', 'LineWidth', 2);
            %         title('Comparison over three phases');
            %         ylabel('No. of licks');
            %         xlabel('Phases');
            %
            figure(6)
            plot(mean(SessionPerformance,1), 'LineWidth', 2) %mean across animals
            axis([1 30 0 100])
            xlabel ('Sessions')
            ylabel('Performance')
            title('Session Performance')
            
        end
        
        %Next Dataset
        reply = input('Proceed to the next Dataset? y/n ', 's');
        if isempty(reply)
            reply = 'y';
        end
    end
end
fclose(fid);

% figure(7)
% plot(WeightChange, SessionPerformance,'h')
% axis([0 15 1 100])
% xlabel('Percent Weight Change')
% ylabel('Session Performance')

%Save
if saveData==1
    %check for directory
    if isdir([save_direc '/' date num2str(newDatasets) '_' block])~=1
        mkdir([save_direc '/' date num2str(newDatasets) '_' block]);
    else
    end
    
    %save data
    save([save_direc '/' date num2str(newDatasets) '_' block '/SessionPerformance.txt'], 'SessionPerformance', '-ASCII'); %.txt
    %save([save_direc '/' date num2str(newDatasets) '_' block' '/SessionPerformance'], 'SessionPerformance'); %.mat
    xlswrite([save_direc '/' date num2str(newDatasets) '_' block '/SessionPerformance.xls'], SessionPerformance); %.xls
    
    save([save_direc '/' date num2str(newDatasets) '_' block '/PTlicks.txt'], 'PTlicks', '-ASCII'); %.txt
    save([save_direc '/' date num2str(newDatasets) '_' block '/PTlicks'], 'PTlicks'); %.mat
    xlswrite([save_direc '/' date num2str(newDatasets) '_' block '/PTlicks.xls'], PTlicks); %.xls
    
    save([save_direc '/' date num2str(newDatasets) '_' block '/NGlicks.txt'], 'NGlicks', '-ASCII'); %.txt
    save([save_direc '/' date num2str(newDatasets) '_' block '/NGlicks'], 'NGlicks'); %.mat
    xlswrite([save_direc '/' date num2str(newDatasets) '_' block 'NGlicks.xls'], NGlicks); %.xls
    
    save([save_direc '/' date num2str(newDatasets) '_' block '/Glicks.txt'], 'Glicks', '-ASCII'); %.txt
    save([save_direc '/' date num2str(newDatasets) '_' block '/Glicks'], 'Glicks'); %.mat
    xlswrite([save_direc '/' date num2str(newDatasets) '_' block '/Glicks.xls'], SessionPerformance); %.xls
    
    save([save_direc '/' date num2str(newDatasets) '_' block '/PreToneLoops.txt'], 'PreToneLoops', '-ASCII'); %.txt
    save([save_direc '/' date num2str(newDatasets) '_' block '/PreToneLoops'], 'PreToneLoops'); %.mat
    xlswrite([save_direc '/' date num2str(newDatasets) '_' block '/PreToneLoops.xls'], PreToneLoops); %.xls
    
    save([save_direc '/' date num2str(newDatasets) '_' block '/NoGoToneLoops.txt'], 'NoGoToneLoops', '-ASCII'); %.txt
    save([save_direc '/' date num2str(newDatasets) '_' block '/NoGoToneLoops'], 'NoGoToneLoops'); %.mat
    xlswrite([save_direc '/' date num2str(newDatasets) '_' block 'NoGoToneLoops.xls'], NoGoToneLoops); %.xls
    
    save([save_direc '/' date num2str(newDatasets) '_' block '/GoToneLoops.txt'], 'GoToneLoops', '-ASCII'); %.txt
    save([save_direc '/' date num2str(newDatasets) '_' block '/GoToneLoops'], 'GoToneLoops'); %.mat
    xlswrite([save_direc '/' date num2str(newDatasets) '_' block '/GoToneLoops.xls'], GoToneLicks); %.xls
    
else
end
