%incorporating meta-k-means clustering as per Ashesh Dhawale's code -
%12/1/13

clear all
close all

direc_list = 'C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\trace_folder_list_20120516.txt';
data_direc = 'C:\Data\data\BlinkView\';
an_direc = 'C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\';

dataset_type = 1;        %1 - trace data; 2 - pseudo data; 3 - non learners

rho_control = 1;         %0 - analyses control no puff datasets, 1 - analyses training (trace or pseudo) datasets

Ca_width = 250;          %in ms

if dataset_type == 1 || dataset_type == 3
    stimOI = 1;              %stimulus of interest - 1 = tone, 2 = puff (useful with pseudorand datasets)
    learned_only = 1;        %1 - use only datasets for animals that have learned the task 0 - use all datasets in list, 2 - analyses only non-learner datasets
elseif dataset_type == 2
    stimOI = 2;               %stimulus of interest - 1 = tone, 2 = puff (useful with pseudorand datasets)
    learned_only = 0;        %1 - use only datasets for animals that have learned the task 0 - use all datasets in list, 2 - analyses only non-learner datasets
elseif dataset_type == 4
    learned_only = 2;
    stimOI = 1;
else
end

filtering_on = 0;        %1 - indiv. trial traces filtered. 0 - un-filtered.
blink_trials_only = 0;   %1 - uses only blink trials to classify cells using RMSE, 0 - uses all trials after pk_behav_trial
negs_nan = 1;            %0 - uses un-changed Ca amlit values, 1 - replaces negative values with nans  
pk_b_trial_only = 0;     %1 - uses only trials after pk behav trial for calculations; 0 - uses all trials   
r_iters = 100;           %number of iterations of randomisation used to find averaged r-shifted rb ratio - might have to go as high as 3000.
non_ov_trials = 0;       %1 - non-overlapping trial sets used for kernel estimation and rb ratio calculation, 0 - all trials used for both
saving = 0;
saved_corr_scores_inter = [];
saved_corr_scores_intra = [];
all_saved_corr_scores = [];
un_clust = [];
%DONT CHANGE THIS - ONLY SPONT ACTIVITY DATA TO BE USED
bk_period_control = 1;   %0 - analyses tone-puff period, 1 - analyses background period


%loading in lists of cell time-fields
if dataset_type == 1
    cell_field_lists = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\trace.mat');
    cell_field_lists = cell_field_lists.pk_times_saved;
    orig_data = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\trace_data.mat');
    orig_data = orig_data.cell_rem_lists;
elseif dataset_type == 2
    cell_field_lists = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\pseudo.mat');
    cell_field_lists = cell_field_lists.pk_times_saved;
    orig_data = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\pseudo_data.mat');
    orig_data = orig_data.cell_rem_lists;
elseif dataset_type == 3
    cell_field_lists = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\bk.mat');
    cell_field_lists = cell_field_lists.pk_times_saved;
    orig_data = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\bk_data.mat');
    orig_data = orig_data.cell_rem_lists;
elseif dataset_type == 4
    cell_field_lists = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\nonl.mat');
    cell_field_lists = cell_field_lists.pk_times_saved;
    orig_data = load('C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\pk_time_saver\nonl_data.mat');
    orig_data = orig_data.cell_rem_lists;
else
end


% %only trace, training datasets are split into three parts, based on trial numbers
% if dataset_type == 1
%     no_s_parts = 3;
% else
%     no_s_parts = 1;
% end

no_s_parts = 1;
corrcut_vec = [];


for session_part = 1:no_s_parts
    if no_s_parts == 3
        if rho_control == 1
            if session_part == 1
                t_list_orig = 1:10;
            elseif session_part == 2
                t_list_orig = 21:30;
            elseif session_part == 3
                t_list_orig = 40:49;
            end
        elseif rho_control == 0
            if session_part == 1
                t_list_orig = 1:10;
            elseif session_part == 2
                t_list_orig = 11:20;
            elseif session_part == 3
                t_list_orig = 21:30;
            end
        end
    elseif no_s_parts == 1
        if rho_control == 1
            t_list_orig = 1:50;
        elseif rho_control == 0
            t_list_orig = 1:30;
        else
        end
    else
    end

    dir_counter = 0;

    saved_scores = [];
    saved_scores_i = [];
    learned_list = [];
    all_scores = [];
    iter = 0;
    saved_pk_times = [];
    
    gp_lists = [];
    cell_rem_lists = [];
    cell_rem_lists2 = [];
    frame_times = [];
    
    set_counter = 0;
    fid = fopen(direc_list);
    while 1
        dir_counter = dir_counter + 1;

            list_direc = fgetl(fid);
            if ~ischar(list_direc),   break,   end


            %loading raw data and setting up initial variables
            [raw_data_mat direc set_type rand_times_list trial_time no_frames no_cells no_trials frame_time...
                CS_onset_time CS_onset_frame CS_duration CS_US_delay US_onset_frame US_duration trial_type_vec learned blink_list] = load_raw_data_mat(list_direc, an_direc, data_direc, rho_control);



            learned_list = [learned_list; learned];
            blink_trials = find(blink_list == 1);

            %condition to use only animals that learned the task
            if learned_only == 1
                if learned == 0
                    continue
                elseif learned == 1

                end
            elseif learned_only == 0
            elseif learned_only == 2
                if learned == 1
                    continue
                elseif learned == 0
                end
            end


            set_counter = set_counter + 1;
            
            %only used for trace data - for pseudo data, no of t_fields
            %calculated from frame_time
            if dataset_type == 1
                a = cell_field_lists{set_counter, 1};
                no_fields = max(a);
            else
            end

            %CLIPPING data near point of interest (tone or puff)
            if bk_period_control == 0        
                pre_clip = 2000;
                post_clip = 2000;
            elseif bk_period_control == 1
                pre_clip = -2000;
                post_clip = 8000; 
            end


            [raw_data_mat no_frames CS_onset_frame US_onset_frame] = rand_raw_data_clipper(raw_data_mat, set_type,...
                stimOI, rand_times_list, pre_clip, post_clip, frame_time, CS_onset_frame, US_onset_frame, CS_US_delay, bk_period_control);


            %calculating dF/F and suppressing peaks of width lesser than Ca_width
            [dff_data_mat dffdata_peaks] = dff_maker(raw_data_mat, frame_time, Ca_width);


            %sanitising data based on no. of infs/nans
            [dff_data_mat dff_data_mat_orig dffdata_peaks dffdata_peaks_orig trial_type_vec bad_cellsi, bad_trialsi] = imaging_data_sanitiser(dff_data_mat, dffdata_peaks, trial_type_vec);
            sanit_ratio = ((size(dff_data_mat, 1).* size(dff_data_mat, 2).* size(dff_data_mat, 1))./(size(dff_data_mat_orig, 1).* size(dff_data_mat_orig, 2).* size(dff_data_mat_orig, 1)));
            no_cells = size(dff_data_mat, 2);
            no_trials_orig = no_trials;
            no_trials = size(dff_data_mat, 3);

            cell_list_x = bad_cellsi;

            if no_cells < 5

                continue
            elseif no_trials_orig - no_trials > (no_trials_orig./2.5)
            %if no_trials_orig - no_trials > (no_trials_orig./5)     
                continue
            else

            end
            
            %loading ROI matrices
            ROI_mat = load([direc 'ROIs.txt']);
            ROI_mat = bwlabel(ROI_mat);
            
            %putting all trials back into dff_data_mat and removing bad
            %trials in next para of code
            dff_data_mat = dff_data_mat_orig;
            dff_data_mat(:, bad_cellsi, :) = [];
            
            %removing bad cell ROIs from ROI matrix
            bad_cellsi_x = sort(bad_cellsi, 'descend');
            for b_cell_no = 1:length(bad_cellsi)
                b_cell = bad_cellsi_x(b_cell_no);
                temp = find(ROI_mat == b_cell);
                ROI_mat(temp) = 0;                          %removing bad cell
                temp = find(ROI_mat > b_cell);
                ROI_mat(temp) = ROI_mat(temp) - 1;          %re-numbering cells numbered higher than b_cell
            end
            clear cell_no
            clear bad_cellsi_x
            
           
            if dataset_type == 2 || rho_control == 0 || dataset_type == 4
                t_list_orig = 1:no_trials;
            else
            end                
            
            
            %identifying bad trials in t_list and removing them.
            [del rem_i] = intersect(t_list_orig, bad_trialsi);
            t_list = t_list_orig;
            t_list(rem_i) = [];      %removing the trials in t_list that also fall in bad_trialsi
           
            
            %replacing neg/huge values with nans
            if negs_nan == 1
                std_val = nanstd(reshape(dff_data_mat, [], 1));
                temp = find(dff_data_mat < -25.*std_val);
                dff_data_mat(temp) = nan;
                
            else
            end
            
            
            bk_data = dff_data_mat(:, :, t_list);
            no_trials = size(bk_data, 3);
            
            
            
            %stitching all trials together to make one long vector for each cell
            bk_data_mat = zeros( (size(bk_data, 1).*no_trials), no_cells);
            for cell_no = 1:no_cells
                cell_dat = squeeze(bk_data(:, cell_no, :));
                for trial_no = 0:(no_trials - 1)
                    bk_data_mat( ((trial_no.*no_frames) + 1):(trial_no + 1).*no_frames, cell_no) = cell_dat(:, (trial_no + 1) );
                end
            end
            clear bk_data
           

            X = bk_data_mat';       %activity data matrix - cells x frames, with trials concatenated
            
            thr = 0.8;              %threshold of kmeans++ runs a pair of cells must co-occur in to be considered for clustering
            no_iters = 1000;         %no of times meta_k_means_tank calls k_means_plusplus; set to 100 for preliminary corrcut scan
            clustat = [];
            
            
            
            
             
            ccorr = corrcoef(bk_data_mat, 'rows', 'pairwise');
            
            %loop to try out various values of corrcut
            ct_range = 0.65:0.05:0.9;
            for corrcut = 0.65:0.05:0.9
                %corrcut                             %corrcut is the threshold inter-cluster correlation coefficient above which they are merged
                [a b c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);
                
                
                clustcorr = [];
                clustnum = [];
                for i = 1 : size(a, 1)
                    temp = a{i,1};
                    clustnum = [clustnum length(temp)];
                    tempcorr = 0;
                    count = 1;
                    for x = 1 : length(temp)-1
                        for y = x + 1 : length(temp)
                            tempcorr = tempcorr + ccorr(temp(x), temp(y));
                            count = count + 1;
                        end
                    end
                    clustcorr = [clustcorr, tempcorr/count];
                end
                clustat = [clustat; [mean(clustnum) (sum(clustcorr.*clustnum)/sum(clustnum))]];
            end
            
            %identifying best value of corrcut
            x = clustat(:, 1);
            x = (x - min(x))/(max(x)-min(x));
            y = clustat(:, 2);
            y = (y - min(y))/(max(y)-min(y));
            [null, corrcuti] = max(x.*y);
            
            corrcut = ct_range(corrcuti);
            
            corrcut_vec = [corrcut_vec; corrcut, session_part];

            clear null
            
            %running meta-k-means to identify final clusters with
            %well-chosen value of corrcut
            no_iters = 1000;            
            [a b c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);
            
            %building list of cells and their group numbers
            no_clusts = size(a, 1);
            cell_gp_vec = zeros(no_cells, 1);
            for c_num = 1:no_clusts
                c_list = a{c_num, 1};
                cell_gp_vec(c_list, 1) = c_num;
            end
            
            clear no_clusts
            clear c_num
            clear c_list
            
            cell_gp_vec = [cell_gp_vec, (1:1:no_cells)'];
            cell_gp_vec_orig = cell_gp_vec;

            %sorting cells within groups as per their corrcoeffs with
            %group-averaged trace
            c_lengths = [];
            c_list_f = cell_gp_vec_orig(:, 1);                  %full list of all cells' cluster numbers
            
            
            bk_data_mat_o = [];
            ROI_mat_t = ROI_mat;
            %accounting for un-clustered cells
            c_list = find(c_list_f == 0);
            un_clust = [un_clust; length(c_list), no_cells];
            for c_no = 1:length(c_list)
                c_noi = c_list(c_no);
                temp = find(ROI_mat == c_noi); 
                ROI_mat_t(temp) = 0;                            %not counting un-clustered cells
            end
            
            for clust_no = 1:max(cell_gp_vec_orig(:, 1));
                c_list = find(c_list_f == clust_no);            %cell numbers in old list that belong to current cluster
                %condition to skip clusters with < 5 cells in them
                if length(c_list) < 5
                    for c_no = 1:length(c_list)
                        c_noi = c_list(c_no);
                        tempi = find(ROI_mat == c_noi);
                        ROI_mat_t(tempi) = 0;                   %cells belonging to very small clusters not to be marked
                    end
                    
                    continue
                else
                end                
                c_lengths = [c_lengths; length(c_list)]; 
                cr_list = c(c_list, clust_no);                  %list of corrcoeffs of cells in c_list, with their own cluster's averaged trace 
                curr_trace_mat = bk_data_mat(:, c_list)';       %matrix of activity data for each cell that belongs to this cluster
                
                %assigning cluster numbers to ROIs
                for c_no = 1:length(c_list)
                    c_noi = c_list(c_no);
                    tempi = find(ROI_mat == c_noi);
                    ROI_mat_t(tempi) = clust_no;
                end
                
                curr_trace_mat = [cr_list, curr_trace_mat];     %concatenating corrcoefs to activity traces to sort both together
                curr_trace_mat = sortrows(curr_trace_mat, -1);      %sorting rows by first column ie corrcoefs
                bk_data_mat_o = [bk_data_mat_o; curr_trace_mat(:, 2:end)];
            end
            
            bk_data_mat_o = bk_data_mat_o';

            %calculating new corrcoeff mat with re-arranged cells
            figure(3)
            imagesc(bk_data_mat_o');
           
            c_o = corrcoef(bk_data_mat_o, 'rows', 'pairwise');
            c_oi = c_o;
            for c_no = 1:length(c_lengths)
                
                if c_no > 1
                    temp = c_oi( (sum(c_lengths(1:(c_no-1))) + 1), ((sum(c_lengths(1:(c_no-1))) + 1)) );
                    c_oi( (sum(c_lengths(1:(c_no-1))) + 1), ((sum(c_lengths(1:(c_no-1))) + 1)) ) = nan;
                else
                    temp = c_oi( (1:(c_lengths(1))), (1:(c_lengths(1))) );
                    c_oi( (1:(c_lengths(1))):(1:(c_lengths(1))) ) = nan;
                end
                saved_corr_scores_intra = [saved_corr_scores_intra; reshape(temp, [], 1)];      %saving intra group corr coeffs
            end
            a = reshape(c_oi, [], 1);
            temp = find(isnan(a) == 1);
            a(temp) = [];
            saved_corr_scores_inter = [saved_corr_scores_inter; a];                             %%saving inter group corr coeffs
            clear temp
            clear a
            
            %saving within group and across group corrcoefs
            
            
            
            
            drawnow
            imagesc(c_o)
            drawnow
            figure(4)
            imagesc(c_o);
            title(['no fields - ' int2str(max(cell_gp_vec_orig(:, 1)))])
            figure(5)
            imagesc(bk_data_mat_o')
            figure(6)
            imagesc(ROI_mat_t)
            keyboard
            
            cell_list_x2 = [];



            gp_lists = [gp_lists; {cell_gp_vec_orig(:, 1)}];
            cell_rem_lists = [cell_rem_lists; {cell_list_x}];
            frame_times = [frame_times; frame_time];
            cell_rem_lists2 = [cell_rem_lists2; {cell_list_x2}];
            beep

            
            
            disp(int2str(set_counter))
            all_saved_corr_scores = [all_saved_corr_scores; {saved_corr_scores_intra}, {saved_corr_scores_inter}];
    end
    fclose(fid);

    %saving group lists to file
    if rho_control == 1 && dataset_type == 1
        if session_part == 1
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_early.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_early.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_early2.mat'], 'cell_rem_lists2');
        elseif session_part == 2
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_mid.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_mid.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_mid2.mat'], 'cell_rem_lists2');
        elseif session_part == 3
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_late.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_late.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_late2.mat'], 'cell_rem_lists2');
        else
        end

    elseif rho_control == 1 && dataset_type == 2
        if session_part == 1
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_r_early.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_r_early.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_r2_early.mat'], 'cell_rem_lists2');
        elseif session_part == 2
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_r_mid.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_r_mid.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_r2_mid.mat'], 'cell_rem_lists2');
        elseif session_part == 3
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_r_late.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_r_late.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_r2_late.mat'], 'cell_rem_lists2');
        end
            
    elseif rho_control == 0 && dataset_type == 1
        if session_part == 1
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_ctr_early.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_early.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr2_early.mat'], 'cell_rem_lists2');
        elseif session_part == 2
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_ctr_mid.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_mid.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr2_mid.mat'], 'cell_rem_lists2');
        elseif session_part == 3
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_ctr_late.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_late.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr2_late.mat'], 'cell_rem_lists2');
        else
        end
    elseif rho_control == 0 && dataset_type == 2
        if session_part == 1
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_ctr_r_early.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_r_early.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_r2_early.mat'], 'cell_rem_lists2');
        elseif session_part == 2
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_ctr_r_mid.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_r_mid.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_r2_mid.mat'], 'cell_rem_lists2');
        elseif session_part == 3
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_ctr_r_late.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_r_late.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_ctr_r2_late.mat'], 'cell_rem_lists2');
        else
        end
        
        
    elseif rho_control == 1 && dataset_type == 4            %non learners, training session
        if session_part == 1
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_nl_early.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_nl_early.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_nl2_early.mat'], 'cell_rem_lists2');
        elseif session_part == 2
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_nl_mid.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_nl_mid.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_nl2_mid.mat'], 'cell_rem_lists2');
        elseif session_part == 3
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\gp_lists_nl_late.mat'], 'gp_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_nl_late.mat'], 'cell_rem_lists');
            save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\orig_data_nl2_late.mat'], 'cell_rem_lists2');
        else
        end
    else
    end

    save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_kmeans\frame_times.mat'], 'frame_times');
    save(['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\bk_corr_mkmeans\saved_corr_scores_' int2str(set_type) '_'...
        int2str(rho_control) '.mat'], 'all_saved_corr_scores');
end 
beep
