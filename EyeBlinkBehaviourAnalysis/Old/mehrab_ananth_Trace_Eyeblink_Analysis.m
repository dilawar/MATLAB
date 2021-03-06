clear all
close all

save_direc = '/Users/kambadurananthamurthy/Desktop/Work/Analysed Data';


%the line below makes the program cycle through all datasets. to analyse
%just one, paste its path into a file called folder_list.txt and edit the
%path below to reflect this.
list_path = '/Users/kambadurananthamurthy/Desktop/Work - NCBS/All_Directories.txt'; %Check path ID
fid = fopen(list_path);
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    
    direc = tline;
    direc = [direc '/'];
    
    %reading in Protocol informatino in tiff tag
    info = imfinfo([direc 'Trial no - 1.tif']);
    temp = info.ImageDescription;
    
    a = findstr(temp, 'Trial Length');
    trial_time = str2num(temp( (a+19):(a+23)));                      %in ms
    
    a = findstr(temp, 'Time of CS');
    time_of_beep = str2num(temp( (a+12):(a+15)));                    %in ms
    
    a = findstr(temp, 'CS-US lag');
    time_of_puff = str2num(temp( (a+12):(a+15))) + time_of_beep;     %in ms
    
    clear a
    
    
    %threshold formulae given in the next two lines
    %thresh = (t_mean.*baseline_mean + t_SD.*basline_SD)
    %pre_beep_thresh = (t_mean.*baseline_mean + t_SD.*basline_SD).*pre_thresh
    t_mean = .6;        %mean multiplier for post beep thresh
    t_SD = 1;           %SD multiplier to determine thresh
    pre_thresh = .7;    %multiplier for pre-beep thresh
    startle_time = 50;  %time after tone onset for which CRs not counted (startle)
    
    
    
    %parsing direc for dataset name and mouse-name (SENSITIVE TO POSITION OF _ )
    %date = findstr(direc, '201');
    
    slashi = findstr(direc, '/');
    
    dataset_name = direc(1, (slashi(1, ((length(slashi)-1) ))+1):(length(direc)-1));
    uscorei = findstr(dataset_name, '_');
    mouse_name = dataset_name(1:uscorei(1,1)-1);
    secl_slash = slashi(1, (length(slashi)-1) );
    control_direc = direc(1, 1:(secl_slash + uscorei(1,1)-1) );
    control_direc = [control_direc '_control_no-puff'];
    direci = direc;
    
    
    
    %--------------
    %loop to first analyse corressponding control_no-puff trials and then loop
    %through the trials of the current dataset
    
    for control_direc_counter = 1:2
        if control_direc_counter == 1
            direc = control_direc;
        else
            direc = direci;
        end
        
        info = load([direc '/info.xls']);
        info_times = load([direc '/times.xls']);
        
        dir_contents = dir(direc);
        
        no_trials = size(dir_contents, 1) -6;
        no_frames = size(info_times, 2);
        frame_time = round(trial_time./no_frames);
        no_frames = size(info_times, 2);
        beep_frame = round(time_of_beep./frame_time);
        puff_frame = round(time_of_puff./frame_time);
        startle_count = round(startle_time./frame_time);
        
        
        %checking for skipped frames
        bad_trials = zeros(no_trials, 1);
        ideal_times = 0:frame_time:( (no_frames-1).*frame_time);
        for trial_no = 1:((no_trials)-1)
            if  max(info_times(trial_no, :) - ideal_times) > 5
                bad_trials(trial_no, 1) = 1;
            elseif max(info_times(trial_no, :) - ideal_times) < 5
            end
            
        end
        
        clear ideal_times
        
        
        %READING IN DATA
        score_mat = zeros(no_trials, no_frames);
        for trial_no = 1: ((no_trials)-1)
            %skipping probe trials
            %if rem(trial_no, 5) == 0
            %    continue
            %    keyboard
            %else
            %end
            
            path = [direc '/Trial no - ' int2str(trial_no) '.tif'];
            im = imread(path, 1);
            frame_orig = zeros(size(im));
            
            %averaging first 8 frames to get initial frame for comparison (averaging reduces differences due to noise)
            for i = 1:8
                frame = imread(path, i);
                frame = double(frame);
                frame_orig = frame_orig + frame;
                
            end
            
            frame_orig = frame_orig./5;
            pixels_orig = reshape(frame_orig, 1, []);
            
            
            
            
            for frame_no = 1:no_frames
                
                frame = imread(path, frame_no);
                frame = double(frame);
                pixels_curr = reshape(frame, 1, []);
                corr_score = corrcoef(pixels_orig, pixels_curr);
                score_mat(trial_no, frame_no) = corr_score(1,2);
                
            end
            
        end
        
        %assigning a higher score for lower correlation (from 0 to 1)
        score_mat = 1-score_mat;
        
        %Calculating threshold for CR detection as a multiple of score SD above mean score during first 8 trials
        baseline_score_traces = score_mat(:, 1:8); %HAHAHA
        baseline_points = reshape(baseline_score_traces, 1, []);
        baseline_mean = mean(baseline_points);
        baseline_SD = std(baseline_points);
        
        clear baseline_score_traces
        clear baseline_points
        
        
        %apportioning score traces for CR and UR periods
        CR_score_traces = score_mat(:, (beep_frame + startle_count):(puff_frame-1) );
        UR_score_traces = score_mat(:, puff_frame:round(puff_frame + 500./frame_time) );
        preCR_score_traces = score_mat(:, 1:(beep_frame-1) );
        
        %classifying as response 1 or 0 on comparing to thresh
        CR_score_vec = zeros(no_trials, 1);
        UR_score_vec = zeros(no_trials, 1);
        record = zeros(no_trials, 1);
        
        for trial_no = 1:((no_trials)-1)
            %positive CR condition
            if sum(CR_score_traces(trial_no, :)) > (t_mean.*baseline_mean + t_SD.*baseline_SD) .*size(CR_score_traces, 2)
                CR_score_vec (trial_no, 1) = sum(CR_score_traces(trial_no, :));          
                %area under the curve as score
            else
                %trial_no %HERE!!
            end
            
            %negative pre-beep condition
            if sum(preCR_score_traces(trial_no, :)) > (t_mean.*baseline_mean + t_SD.*baseline_SD).*size(preCR_score_traces,2).*pre_thresh
                CR_score_vec(trial_no, 1) = 0;
                record(trial_no, 1) = 1;
            else
            end
            
            
            if sum(UR_score_traces(trial_no, :)) > (t_mean.*baseline_mean + t_SD.*baseline_SD).*size(UR_score_traces, 2)
                UR_score_vec (trial_no, 1) = sum(UR_score_traces(trial_no, :));
                %area under the curve as score
            else
            end
            if control_direc_counter == 1
                %clear UR_score_vec
                %clear UR_score_traces
                control_score_vec = CR_score_vec;
                
            elseif control_direc_counter == 2
                pad = zeros((no_trials - length(control_score_vec) ), 1 );
                control_score_vec = [control_score_vec; pad];
                
            end
        end
    end
    
    %keyboard
    
    %Calculating Percentage response with a moving window
    window_width = 6;
    window_step = 1;
    CR_percentage_vec = zeros(((no_trials./window_step) - (window_width-window_step)./window_step ) - 1, 1);
    UR_percentage_vec = zeros(((no_trials./window_step) - (window_width-window_step)./window_step ) - 1, 1);
    control_percentage_vec = zeros(((no_trials./window_step) - (window_width-window_step)./window_step ) - 1, 1);
    for i = 0:((no_trials./window_step) - (window_width-window_step)./window_step ) - 1;
        CR_percentage_vec(i+1, 1) = mean(CR_score_vec( (i.*window_step + 1):(i.*window_step) + window_width) ).*100;
        UR_percentage_vec(i+1, 1) = mean(UR_score_vec( (i.*window_step + 1):(i.*window_step) + window_width) ).*100;
        control_percentage_vec(i+1, 1) = mean(control_score_vec( (i.*window_step + 1):(i.*window_step) + window_width) ).*100;
        
    end
    
    
    
    
    
    %--------------------
    
    
    %normalising CR and UR scores and percentages to UR mean
    UR_score_mean = mean(UR_score_vec);
    UR_percentage_mean = mean(UR_percentage_vec);
    
    UR_score_vec = UR_score_vec./UR_score_mean;
    CR_score_vec = CR_score_vec./UR_score_mean;
    control_score_vec = control_score_vec./UR_score_mean;
    
    UR_percentage_vec = UR_percentage_vec./UR_percentage_mean;
    CR_percentage_vec = CR_percentage_vec./UR_percentage_mean;
    control_percentage_vec = control_percentage_vec./UR_percentage_mean;
    
    %Kolmogorov - Smirnov Test
    %[h,p] = kstest2(CR_score_vec,control_score_vec)
    
    %criterion for mouse having learned
    
    learn = [];
    
    
    for i = 1:(length (CR_score_vec)-1)
        b = CR_score_vec(i);
        c = [];
        for k = 1:29
            c = control_score_vec(k,:);
        end
        
        if b > mode(c)
            learn(i) = 1;
        else learn(i) = 0;
        end
    end
    
    sum_learn = sum(learn);
    per_acc = (sum_learn/length(CR_percentage_vec))*100;
    
    
    %plotting
%     figure(3)
%     plot(UR_percentage_vec, '*')
%     hold on
%     plot(CR_percentage_vec, 'r*')
%     title(dataset_name)
%     plot(control_percentage_vec, 'g*')
%     hold off
%     
    
    
%     figure(4)
%     plot(UR_score_vec, '*')
%     hold on
%     plot(CR_score_vec, 'r*')
%     title(dataset_name)
%     plot(control_score_vec, 'g*')
%     hold off
%     input('So, do we go to the next session?')
%     
    save_data = [UR_percentage_vec, CR_percentage_vec];
%     save([save_direc '/' dataset_name '.txt'], 'save_data', '-ASCII');
%     
%     updating central list with 'learned or not' status
%     all_data = load('/Users/kambadurananthamurthy/Desktop/Work - NCBS/Users.txt');
%     
%     padding dataset_name
    pad = zeros(1, (50 - size(dataset_name, 2)) ) + 32;
    pad = char(pad);
    dataset_name = [dataset_name, pad];
    
    
    if size(int2str(max(CR_percentage_vec)*100), 2)<2
        CR_max = ['0' int2str(max(CR_percentage_vec)*100)];
    else
        CR_max = [int2str(max(CR_percentage_vec)*100)];
    end
    
     all_data = [UR_score_vec, CR_score_vec, control_score_vec];
     
     mkdir([save_direc '/' dataset_name '/']);
     save([save_direc '/' dataset_name '/all_data.txt'], 'all_data', '-ASCII');
     
     save([save_direc '/' dataset_name '/percentage_accuracy.txt'],'per_acc', '-ASCII');
     
    %Next Session
    %input('Hit the Return key to Go to the next session/end')
    
end
fclose(fid);