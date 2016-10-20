% Purpose: To analyse the Lick Behaviour Training data
% Author: Kambadur Ananthamurthy

%clear all
%close all

%NOTE: Perform Analysis with one block from one animal's worth of data at a time. Else,
%comment out the calculations for SEM, etc.

save_direc='/Users/ananth/Desktop/Work/Analysis/LickBehaviourAnalysis';

%The line below makes the program cycle through all datasets. To analyse
%just one, paste its path into a file called folder_list.txt and edit the
%path below to reflect this.

list_path='/Users/ananth/Desktop/Work/LickBehaviour/Training/folder_list.txt'; %Check path ID
fid=fopen(list_path);

totalsessions=7;
toneduration=300; % in ms
ToneLickProfile=zeros(300,totalsessions);

while 1 %You might see a variable "ans" = 0, after the program runs to completion
    tline=fgetl(fid);
    if ~ischar(tline), break, end
    
    direc=[tline '/'];
    
    %Read Protocol information from a text file labelled protocol.txt
    protocol_path=[direc 'protocol.txt'];
    fid_protocol=fopen(protocol_path);
    protocol=zeros(1,10); %I have used 10 since there are 10 parameters stored in the text file protocol.txt, at this time
    for i=1:11
        temp=fgetl(fid_protocol); %Reads the first/next line in the .txt file; replaces the text from previous line
        if temp==-1
            
        else
            a=findstr(temp, ' - '); %Finds " - "
            
            if isempty(a)
                %skipping headings
            else
                protocol(i)=str2num(temp((a+2):length(temp))); % everything from after " - " to the end of the line
                
                % The indexing for the array "protocol" is defined as
                %         1. Number of Trials
                %         2. Pre-tone Interval (ms)
                %         3. Tone Frequency (Hz)
                %         4. CS Duration (ms)
                %         5. Reward Duration (ms)
                %         6. Max ITI (s)
                %         7. Min ITI (s)
                %         8. Water Time (ms)
                %         9. Wait Before Reward (ms)
                %         10. Min Timeout Wait (s)
                %         11. Direct Water Till Trial
                clear temp
                clear a
            end
        end
    end
    
    clear i
    fclose(fid_protocol);
    
    %General Information
    slashi=findstr(direc, '/');
    dataset=direc(1, (slashi(1, ((length(slashi)-1) ))+1):(length(direc)-1)); % Dataset name, in the format of MouseX_lick_sessionY_extrainfo
    
    uscorei=findstr(dataset, '_');
    mouse=dataset(1:uscorei(1,1)-1);
    
    %     block=dataset((uscorei(1,1)+1):(uscorei(1,2)-1));
    %     blocknumber_string=block(6:length(block));
    %     blocknumber=str2num(blocknumber_string);
    block='Block3'; %string
    blocknumber=3; %number
    %     sessioni=findstr(dataset, 'ession'); % just to avoid confusion between "Session" and "session"
    session=dataset((uscorei(1,2)+1):length(dataset));
    sessionnumber_string=session(8:length(session));
    sessionnumber=str2num(sessionnumber_string);
    
    %Read the results.m file for the training session results and save into a matrix "results"
    results_path=[direc block '_' session '_results.m'];
    
    results=load(results_path);
    
    %The indexing for the matrix "results" is defined as
    %     1. Trial Number
    %     2. Pre-Tone Lick
    %     3. Tone Lick
    %     4. Reward Lick
    %     5. ITI Lick
    %     6. Pre-Tone Lick Time (ms)
    %     7. Tone Lick Time (ms)
    %     8. Reward Lick Time (ms)
    %     9. Current ITI (s)
    %     10. ITI Lick Time (ms)
    %     11. Current Success Rate (%)
    
    [ntrials,nparams]=size(results); % Number of trials that actually occurred with the number of parameters measured/saved
    %Don't use the length function as this simply calculates the length of the largest dimension (not necessarily the number of trials)
    
    %     if ntrials~=600;
    %         disp('Check ntrials');
    %     end
    
    %toneduration=protocol(4);
    
    %Plot frequency histograms
    %Numbers
    nPreToneLicks=sum(results(:,2));
    nToneLicks=sum(results(:,3));
    nRewardLicks=sum(results(:,4));
    nITILicks=sum(results(:,5));
    
    TotalLicks=(nPreToneLicks+nToneLicks+nRewardLicks+nITILicks);
    
    %Fractions
    fPreToneLick=nPreToneLicks/TotalLicks;
    fToneLick=nToneLicks/TotalLicks;
    fRewardLick=nRewardLicks/TotalLicks;
    fITILick=nITILicks/TotalLicks;
    
    %Percentages % Think about preallocating
    %PercentLicks(blocknumber,(1:4),sessionnumber)=100*[fPreToneLick,fToneLick,fRewardLick,fITILick]; %builds the matrix
    PercentLicks((1:4),sessionnumber)=100*[fPreToneLick,fToneLick,fRewardLick,fITILick]; %builds the matrix
    
    %Probabilities
    %     pPreToneLick=nPreToneLicks/ntrials;
    %     pToneLick=nToneLicks/ntrials;
    %     pCorrectLick=nRewardLicks/ntrials;
    %     pITILick=nITILicks/ntrials;
    %
    %pLicks(blocknumber,(1:4),sessionnumber)=[pPreToneLick,pToneLick,pCorrectLick,pITILick]; %builds the matrix
    
    %Xaxis=['Pre-Tone' 'Tone' 'Reward' 'ITI'];
    
    
    %Lick Profile (per ms)
    % lickprofile=zeros(1,(protocol(1,2)+protocol(1,4)+protocol(1,9)+protocol(1,5)+((protocol(1,6)+protocol(1,7))*1000/2)));
    CorrectTrials=0; % only for initialization
    for i=1:ntrials
        %         if (results(i,6))==0
        %         else
        %             lickprofile(1,(results(i,6)))=1; %adding pretone licks
        %         end
        
        if results(i,3)==1
            if results(i,7)>300
                ToneLickProfile(300,sessionnumber)=ToneLickProfile(300,sessionnumber)+1; %adding tone licks (considers all sessions specified)
            else
                ToneLickProfile(results(i,7),sessionnumber)=ToneLickProfile(results(i,7),sessionnumber)+1; %adding tone licks (considers all sessions specified)
            end
        end
        
        %         if (results(i,8)+protocol(1,2)+protocol(1,4)+protocol(1,9))==0
        %         else
        %             lickprofile(1,(results(i,8)+protocol(1,2)+protocol(1,4)+protocol(1,9)))=1; %adding reward phase licks
        %         end
        %
        %         if (results(i,10)+protocol(1,2)+protocol(1,4)+protocol(1,9)+protocol(1,5))==0
        %         else
        %             lickprofile(1,(results(i,10)+protocol(1,2)+protocol(1,4)+protocol(1,9)+protocol(1,5)))=1; %adding ITI licks
        %         end
        
        if results(i,4)==1
            CorrectTrials=CorrectTrials+1;
        end
    end
    clear i
    
    SessionSuccessRate(sessionnumber)=(CorrectTrials/ntrials)*100; % for all sessions considered
    
    %Lick Profile Histogram (per binsize ms) - TONE
    binsize=10; %in ms
    nbins=length(ToneLickProfile(:,1))/binsize;
    
    tone_lickhist=zeros(nbins,totalsessions);
    for i=1:nbins
        tone_lickhist(i,sessionnumber)=sum(ToneLickProfile((((i-1)*binsize)+1):(((i-1)*binsize)+binsize)));
    end
    normalizingfactor_tone_lickhist=max(max(tone_lickhist));
    tone_lickhist=tone_lickhist/normalizingfactor_tone_lickhist;
    mean_tone_lickhist=mean(tone_lickhist,2);
    normalizingfactor_mean_tone_lickhist=max(max(mean_tone_lickhist));
    mean_tone_lickhist=mean_tone_lickhist/normalizingfactor_mean_tone_lickhist;
    
    %Save Data
    if isdir([save_direc '/' mouse '/' block])~=1
    mkdir([save_direc '/' mouse '/' block]);
    else
    end
    
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
figure(2); bar(SessionSuccessRate);
figure(3); bar(tone_lickhist);
figure(4); bar(mean_tone_lickhist);
% bar(ToneLickProfile);
% bar(mean(ToneLickProfile,2);
%errorbar(mean_Inst_SuccessRate,sem_Inst_SuccessRate);

%Save
save([save_direc '/' mouse '/PercentLicks.txt'], 'PercentLicks', '-ASCII'); %.txt
save([save_direc '/' mouse '/PercentLicks'], 'PercentLicks'); %.mat
save([save_direc '/' mouse '/SessionSuccessRate.txt'], 'SessionSuccessRate', '-ASCII'); %.txt
save([save_direc '/' mouse '/SessionSuccessRate'], 'SessionSuccessRate'); %.mat
save([save_direc '/' mouse '/tone_lickhist.txt'], 'tone_lickhist', '-ASCII'); %.txt
save([save_direc '/' mouse '/tone_lickhist'], 'tone_lickhist'); %.mat
save([save_direc '/' mouse '/mean_tone_lickhist.txt'], 'mean_tone_lickhist', '-ASCII'); %.txt
save([save_direc '/' mouse '/mean_tone_lickhist'], 'mean_tone_lickhist'); %.mat