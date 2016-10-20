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

totalanimals=3;  %*! Must be edited for every experiment
totalsessions=6; %*! Must be edited for every animal
totaltrials=400; %*! Must be edited for every animal
TotalCorrectTrials=zeros(totalanimals,totalsessions,totaltrials);
SuccessRate=zeros(totalanimals,totalsessions,totaltrials);
CorrectTrials=zeros(totalanimals,totalsessions,totaltrials);
Tone1Loops=zeros(totalanimals,totalsessions,totaltrials);
PerfectTrials=zeros(totalanimals,totalsessions,totaltrials);
PerfectSuccessRate=zeros(totalanimals,totalsessions,totaltrials);

toneduration=500;
step=30;
GoToneLickHist=zeros((toneduration/step),1);

mouse_check='MouseK1'; %First animal's name
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
    
    if mouse==mouse_check
        animal=1;
    else
        animal=animal+1;
    end
    
    mouse_check=mouse; %assigning for a new animal
    
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
    protocol=zeros(1,8); %I have used 10 since there are 10 parameters stored in the text file protocol.txt, at this time
    for i=1:8
        temp=fgetl(fid_protocol); %Reads the first/next line in the .txt file; replaces the text from previous line
        if temp==-1
            
        else
            a=strfind(temp, ' - '); %Finds " - "
            
            if isempty(a)
                %skipping headings
            else
                value=temp((a+2):length(temp));
                if isnumeric(value)
                    protocol(1,i)=str2num(value); % everything from after " - " to the end of the line
                else
                    %protocol(1,i)=value;
                    protocol(1,i)=0;
                end
                
                % The indexing for the array "protocol" is defined as
                %         1. Task
                %         2. Block
                %         3. Session
                %         4. CS A (Hz)
                %         5. CS B (Hz)
                %         6. CS 1 Duration Min (ms)
                %         7. CS 2 Duration (ms)
                %         8. Water Time (ms)
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
    %1. Trial Number
    %2. Tone 1 loops
    %3. Tone 2 lick
    %4. Actual Tone 1 duration (ms)
    %5. Actual Tone 2 duration (ms)
    %6. Success Rate (%)
    %7. Lick Dependency
    %8. Correct Trials
    [ntrials,nparams]=size(results); % Number of trials that actually occurred with the number of parameters measured/saved
    %Don't use the length function as this simply calculates the length of the largest dimension (not necessarily the number of trials)
    
    %     if ntrials~=600;
    %         disp('Check ntrials');
    %     end
    
    TotalCorrectTrials(animal,sessionnumber)=results(ntrials,nparams);
    SuccessRate(animal,sessionnumber)=(TotalCorrectTrials(animal,sessionnumber)/ntrials)*100; % for all animals and sessions considered
    
    %TotalCorrectTrials=results(ntrials,nparams); %Since the last column is saved with the sum of the correct trials
    %SuccessRate=(TotalCorrectTrials/ntrials)*100; % for all animals and sessions considered
    CorrectTrials(animal,sessionnumber,(1:ntrials))=results(:,nparams);
    
    Tone1Loops(animal,sessionnumber,(1:ntrials))=results(:,2);
    
    % Finding perfect trials
    PerfectTrials(animal,sessionnumber,(1:ntrials))=ismember(Tone1Loops(animal,sessionnumber,(1:ntrials)),1);
    PerfectSuccessRate(totalanimals,totalsessions)=(sum(PerfectTrials(animal,sessionnumber))/ntrials)*100;
    
    %     %Save Data
    %     if isdir([save_direc '/' mouse])~=1
    %         mkdir([save_direc '/' mouse]);
    %     else
    %     end
    
    %disp(sessionnumber);
    %disp('done!');
    
    GoToneLickTimes=sort(results(:,5),'ascend');
    j=0;
    count=1;
    for i=step:step:toneduration
        GoToneLickHist(count)=length(GoToneLickTimes(GoToneLickTimes>j & GoToneLickTimes<i));
        j=i;
        count=count+1;
    end
    clear i
    clear j
    clear count
    
    norm=max(max(GoToneLickHist));
    GoToneLickHist_norm=GoToneLickHist/norm;
    
    %Plots - one dataset at a time!
    figure(1); plot(results(:,2),'blue'); %Tone 1 loops
    hold on; plot((results(:,6)+300),'black'); % Current Success Rate
    hold on; plot((results(:,3)*100),'green*'); % Tone 2 licks
    
    figure(2);
    bar(GoToneLickHist_norm);
    
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