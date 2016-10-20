% Purpose: To analyse the Delayed Non-Match to Sample (DNMTS) Training data
% Author: Kambadur Ananthamurthy

%clear all
%close all

%NOTE: Perform Analysis with one block from one animal's worth of data at a time. Else,
%comment out the calculations for SEM, etc.

save_direc='/Users/ananth/Desktop/Work/Analysis/DNMTSAnalysis';

%The line below makes the program cycle through all datasets. To analyse
%just one, paste its path into a file called folder_list.txt and edit the
%path below to reflect this.

list_path='/Users/ananth/Desktop/Work/DNMTS/DNMTS_folder_list.txt'; %Check path ID
fid=fopen(list_path);

totalsessions=1; %!@ Specify every time
PercentLicks=zeros(4,totalsessions);
% pLicks=zeros(4,totalsessions);

SuccessRate=zeros(totalsessions); %in "%"
FailureRate=zeros(totalsessions); %in "%"
Performance=zeros(totalsessions); %in "%"

while 1 %You might see a variable "ans" = 0, after the program runs to completion
    tline=fgetl(fid);
    if ~ischar(tline), break, end
    
    direc=[tline '/'];
    
    %General Information
    slashi=strfind(direc, '/');
    dataset=direc(1, (slashi(1, ((length(slashi)-1) ))+1):(length(direc)-1)); % Dataset name, in the format of MouseX_DNMTS_sessionY_extrainfo
    
    uscorei=strfind(dataset, '_'); %index of all the "_"s in the dataset name
    mouse=dataset(1:uscorei(1,1)-1); %Mouse name
    task=dataset((uscorei(1,1)+1):(uscorei(1,2)-1)); %Task
    block=dataset((uscorei(1,2)+1):(uscorei(1,3)-1)); %Session
    session=dataset((uscorei(1,3)+1):length(dataset)); %Session
    sessionnumber=str2num(session(8:length(session)));
    
    %Read Protocol information from a text file labelled protocol.txt
%     protocol_path=[direc mouse '_' task '_' session '_protocol.txt'];
%     fid_protocol=fopen(protocol_path);
%     protocol=zeros(1,11); %I have used 11 since there are 11 parameters stored in the text file protocol.txt, at this time
%     for i=1:11
%         temp=fgetl(fid_protocol); %Reads the first/next line in the .txt file; replaces the text from previous line
%         if temp==-1
%             
%         else
%             a=strfind(temp, ' - '); %Finds " - "
%             
%             if isempty(a)
%                 %skipping headings
%             else
%                 protocol(i)=str2num(temp((a+2):length(temp))); % everything from after " - " to the end of the line
%                 
%                 % The indexing for the array "protocol" is defined as
%                 %         1. Number of Trials
%                 %         2. Pre-tone Interval (ms)
%                 %         3. Tone 1 Frequency (Hz)
%                 %         4. CS 1 Duration (ms)
%                 %         5. Delay Interval (ms)
%                 %         6. Tone 2 Frequency (Hz)
%                 %         7. CS 2 Duration (ms)
%                 %         8. Max ITI (s)
%                 %         9. Min ITI (s)
%                 %         10. Water Time (ms)
%                 %         11. Min Timeout Wait (s)
%                 clear temp
%                 clear a
%             end
%         end
%     end
%     
%     clear i
%     fclose(fid_protocol);
    
    %Read the results.m file for the training session results and save into a matrix "results"
    results_path=[direc mouse '_' task '_' block '_' session '_result.m'];
    
    results=load(results_path);
    
%     The indexing for the matrix "results" is defined as
%     1. Trial Number
%     2. Tone Lick
%     3. Trial Condition
%     4. Actual CS duration
%     5. Lick Dependency
%     6. Correct trials
%     7. Incorrect trials
%     8. Failure Rate
%     9. Performance
%     10. Success Rate
    

    [ntrials,nparams]=size(results); % Number of trials that actually occurred with the number of parameters measured/saved
    %Don't use the length function as this simply calculates the length of the largest dimension (not necessarily the number of trials)
    
    %     if ntrials~=600;
    %         disp('Check ntrials');
    %     end
    %     if nparams~=17;
    %         disp('Check nparams');
    %     end
    
    % 
    
%     %Plot frequency histograms
%     %Numbers (Licks)
%     nPreToneLicks=sum(results(:,2));
%     nTone1Licks=sum(results(:,3));
%     nDelayIntervalLicks=sum(results(:,4));
%     nTone2Licks=sum(results(:,5));
%     
%     TotalLicks=(nPreToneLicks+nTone1Licks+nDelayIntervalLicks+nTone2Licks);
%     
%     %Fractions (Licks)
%     fPreToneLick=nPreToneLicks/TotalLicks;
%     fTone1Lick=nTone1Licks/TotalLicks;
%     fDelayIntervalLick=nDelayIntervalLicks/TotalLicks;
%     fTone2Lick=nTone2Licks/TotalLicks;
%     
%     %Percentages (Licks)
%     PercentLicks((1:4),sessionnumber)=100*[fPreToneLick,fTone1Lick,fDelayIntervalLick,fTone2Lick]; %builds the matrix
%     
%     %Probabilities (Licks)
%     %     pPreToneLick=nPreToneLicks/ntrials;
%     %     pTone1Lick=nToneLicks/ntrials;
%     %     pDelayIntervalLick=nRewardLicks/ntrials;
%     %     pTone2Lick=nITILicks/ntrials;
%     %
%     %pLicks(blocknumber,(1:4),sessionnumber)=[pPreToneLick,pToneLick,pCorrectLick,pITILick]; %builds the matrix
%     
%     %Xaxis=('Pre-Tone' 'Tone1' 'Delay' 'Tone2');
%     
    %Success Rates (%), Failure Rates (%) and Performance (%)
    CorrectTrials=sum(results(ntrials, 6));
    nNMTs=sum(results(:,3)); %number of Non-Match Trials (sums all the "1"s in the column designated to "Trial Condition")
    SuccessRate(sessionnumber)=(CorrectTrials/nNMTs)*100; %in "%" for all sessions considered
    
    IncorrectTrials=sum(results(:,8));
    nMTs=ntrials-nNMTs; %number of Match Trials (Effectively, sums all the "0"s in the column designated to "Trial Condition")
    FailureRate(sessionnumber)=(IncorrectTrials/nMTs)*100; %in "%" for all sessions considered
    
    %AbortedTrials=(results(ntrials,17));
    
    Performance(sessionnumber)=((nNMTs/ntrials)*SuccessRate(sessionnumber))+((nMTs/ntrials)*(100-FailureRate(sessionnumber))); %in "%" for all sessions considered
    
    %Lick Profile Histogram (per binsize ms) - TONE
%     binsize=10; %in ms
%     nbins=length(ToneLickProfile(:,1))/binsize;
%     
%     tone_lickhist=zeros(nbins,totalsessions);
%     for i=1:nbins
%         tone_lickhist(i,sessionnumber)=sum(ToneLickProfile((((i-1)*binsize)+1):(((i-1)*binsize)+binsize)));
%     end
%     normalizingfactor_tone_lickhist=max(max(tone_lickhist));
%     tone_lickhist=tone_lickhist/normalizingfactor_tone_lickhist;
%     mean_tone_lickhist=mean(tone_lickhist,2);
%     normalizingfactor_mean_tone_lickhist=max(max(mean_tone_lickhist));
%     mean_tone_lickhist=mean_tone_lickhist/normalizingfactor_mean_tone_lickhist;
%     
%     %Save Data
%     if isdir([save_direc '/' mouse '_' task '_' session])~=1
%     mkdir([save_direc '/' mouse '_' task '_' session]);
%     else
%     end
    
    %save([save_direc '/' mouse '/' block '/PercentLicks.txt'], 'PercentLicks', '-ASCII'); %.txt
    %save([save_direc '/' mouse '/' block '/PercentLicks'], 'PercentLicks'); %.mat
    %     save([save_direc '/' mouse '/' block '/pLicks.txt'], 'pLicks', '-ASCII'); %.txt
    %     save([save_direc '/' mouse '/' block '/pLicks'], 'pLicks'); %.mat
    %     save([save_direc '/' mouse '/' block '/lickprofile.txt'], 'lickprofile', '-ASCII'); %.txt
    %     save([save_direc '/' mouse '/' block '/lickprofile'], 'lickprofile'); %.mat
    %     save([save_direc '/' mouse '/' block '/tone_lickhist.txt'], 'tone_lickhist', '-ASCII'); %.txt
    %     save([save_direc '/' mouse '/' block '/tone_lickhist'], 'tone_lickhist'); %.mat
    
    %     if blocknumber_string==1;
    %         %do nothing
    %     else
    %         %Session Success Rate
    %         SessionSuccessRate(blocknumber,sessionnumber)=results(ntrials,11); %refer above
    %
    %         %Instantaneous Success Rate
    %         Inst_SuccessRate(blocknumber,(1:ntrials),sessionnumber)=results(:,11); %refer above
    %
    %Next Dataset
    %     reply = input('Proceed to the next Dataset? y/n ', 's');
    %     if isempty(reply)
    %         reply = 'y';
    %     end
    disp(sessionnumber);
    disp('done!');
end
%end
fclose(fid);

% %MeanInstSuccessRate
% mean_Inst_SuccessRate=mean(Inst_SuccessRate,2);
% sem_Inst_SuccessRate=std(Inst_SuccessRate,1,2)/sqrt(sessionnumber); %considers the latest session number (last session)
%
% save([save_direc '/' mouse '/' block '/mean_Inst_SuccessRate'], 'mean_Inst_SuccessRate'); %.mat
% save([save_direc '/' mouse '/' block '/sem_Inst_SuccessRate'], 'sem_Inst_SuccessRate'); %.mat
%
% %SessionSuccessRate
% save([save_direc '/' mouse '/' block '/SessionSuccessRate'], 'SessionSuccessRate'); %.mat

%Plots
figure(1); bar(PercentLicks(:,:));
title('Percentage of Total Licks');
figure(2); bar(SuccessRate);
title('Success Rate (%)');
figure(3); bar(FailureRate);
title('Failure Rate (%)');
figure(4); bar(Performance);
title('Performance (%)');
% bar(ToneLickProfile);
% bar(mean(ToneLickProfile,2));
% errorbar(mean_Inst_SuccessRate,sem_Inst_SuccessRate);

%Save Data
if isdir([save_direc '/' mouse])~=1
    mkdir([save_direc '/' mouse]);
else
end

%Save
save([save_direc '/' mouse '/PercentLicks.txt'], 'PercentLicks', '-ASCII'); %.txt
save([save_direc '/' mouse '/PercentLicks'], 'PercentLicks'); %.mat
csvwrite([save_direc '/' mouse '/PercentLicks.csv'], PercentLicks);

save([save_direc '/' mouse '/SuccessRate.txt'], 'SuccessRate', '-ASCII'); %.txt
save([save_direc '/' mouse '/SuccessRate'], 'SuccessRate'); %.mat
csvwrite([save_direc '/' mouse '/SuccessRate.csv'], SuccessRate);

save([save_direc '/' mouse '/FailureRate.txt'], 'FailureRate', '-ASCII'); %.txt
save([save_direc '/' mouse '/FailureRate'], 'FailureRate'); %.mat
csvwrite([save_direc '/' mouse '/FailureRate.csv'], FailureRate);

save([save_direc '/' mouse '/Performance.txt'], 'Performance', '-ASCII'); %.txt
save([save_direc '/' mouse '/Performance'], 'Performance'); %.mat
csvwrite([save_direc '/' mouse '/Performance.csv'], Performance);
% save([save_direc '/' mouse '/tone_lickhist.txt'], 'tone_lickhist', '-ASCII'); %.txt
% save([save_direc '/' mouse '/tone_lickhist'], 'tone_lickhist'); %.mat
% save([save_direc '/' mouse '/mean_tone_lickhist.txt'], 'mean_tone_lickhist', '-ASCII'); %.txt
% save([save_direc '/' mouse '/mean_tone_lickhist'], 'mean_tone_lickhist'); %.mat

%csvwrite('csvlist.dat',m,0,2) %filename, variable, row, column %% values of row and column represent offsets
%type csvlist.dat

% xlswrite([save_direc '/' mouse '/PercentLicks.xls'], 'PercentLicks');
