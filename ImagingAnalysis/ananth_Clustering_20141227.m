% Everything after developing the non-baseline corrected raw intensity time
% series for all cells, trials and frames; till clustering and analysis of
% the clusters

% Kambadur Ananthamurthy

dataset='20140312/Mouse6/Mouse6-ROI';

% load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/' dataset '/cellmask.mat'])
% load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/' dataset '/calbdf.mat'])

% ncells=144;
% nsessiontrials=30;
% nframes=266;

tau = 0.5;  % time-course of calcium transient in seconds
k = 3;  % number of clusters

% Meta k means
%first get the optimal corrcut (cutoff)

calbdf1((1:ncells),(1:nsessiontrials),(1:85))=calbdf(:,:,(1:85)); %calbdf_spont or pre-stim
calbdf2((1:ncells),(1:nsessiontrials),(1:(length(89:178))))=calbdf(:,:,(89:178)); %calbdf_stimn9s or 1s stimulus with 9 seconds of afters
calbdf3((1:ncells),(1:nsessiontrials),(1:(length(179:nframes))))=calbdf(:,:,(179:nframes)); %%calbdf_poststimn10s or from 10s after the stim onset

%To get rid of some trials
for i = 1:ncells
  tone_mouse6_final = tone_mouse6;
  tone_mouse6_final(:,[6,11,13,14,19],:) = [];
end
l=size(tone_mouse6_final,2);
tone1((1:ncells),(1:l),(1:85))=tone_mouse6_final(:,:,(1:85)); %calbdf_spont or pre-stim
tone2((1:ncells),(1:l),(1:(length(89:178))))=tone_mouse6_final(:,:,(89:178)); %calbdf_stimn9s or 1s stimulus with 9 seconds of afters
tone3((1:ncells),(1:l),(1:(length(179:nframes))))=tone_mouse6_final(:,:,(179:nframes)); %%calbdf_poststimn10s or from 10s after the stim onset

% Trial-shuffled calbdf dataset

% rand('seed', 0);
% calbdfr_spont = calbdf_spont * 0;
% for c = 1 : max(max(cellmask))
%     calbdfr_spont(c,(1:49), (1:95)) = calbdf_spont(c, randperm(49), :);
% end

Z=tone1;

X = reshape(permute(Z(:, :, :), [1 3 2]), [size(calbdf, 1) numel(squeeze(Z(1, :, :)))]);

ccorr = tcorr; %(2:end, 2:end);
thr = 0.8;
clustat = [];

niterations=1000;
for corrcut = 0.6:0.05:0.9
    corrcut
    [a b c] = ananth_meta_k_means_tank(X, 'correlation', corrcut, thr, niterations); %1000 iterations
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
        clustcorr = [clustcorr tempcorr/count];
    end
    clustat = [clustat; [mean(clustnum) (sum(clustcorr.*clustnum)/sum(clustnum))]];
end

% Dombeck et al's heuristic to determine correct correlation cutoff

x = clustat(:, 1);
x = (x - min(x))/(max(x)-min(x));
y = clustat(:, 2);
y = (y - min(y))/(max(y)-min(y));
x.*y
maxvalue = max(x.*y)

maxvalue_indices=find(x.*y == maxvalue);
corrcuti=max(maxvalue_indices);
ct_range = 0.6:0.05:0.9;
corrcut = ct_range(corrcuti)

clear a
clear b
clear c

%Finally, with a chosen corrcut,
niterations=2500;
[a b c] = ananth_meta_k_means_tank(X, 'correlation', corrcut, thr, niterations); %2500 iterations


no_clusts = size(a, 1);
% no_cells=max(max(max(cellmask)));
no_cells=118;

cell_gp_vec = zeros(no_cells, 1);

for clustnum = 1:no_clusts
    c_list = a{clustnum, 1};
    cell_gp_vec(c_list, 1) = clustnum;
end

cell_gp_vec = [cell_gp_vec, (1:1:no_cells)'];
cell_gp_vec_orig = cell_gp_vec;

%sorting cells within groups as per their corrcoeffs with
%group-averaged trace
c_lengths = [];
c_list_f = cell_gp_vec_orig(:, 1);                  %full list of all cells' cluster numbers

bk_data_mat_o = [];
curr_trace_mat = [];

cellmask_t=cellmask;

un_clust = [];
saved_corr_scores_inter = [];
saved_corr_scores_intra = [];

%accounting for un-clustered cells
c_list = find(c_list_f == 0);
un_clust = [un_clust; length(c_list), no_cells];
for c_no = 1:length(c_list)
    c_noi = c_list(c_no);
    temp = find(cellmask == c_noi);
    cellmask_t(temp) = 0;                            %not counting un-clustered cells
end

for clust_no = 1:max(cell_gp_vec_orig(:, 1));
    c_list = find(c_list_f == clust_no);            %cell numbers in old list that belong to current cluster
    
    %condition to skip clusters with < 5 cells in them
    if length(c_list) < 5
        for c_no = 1:length(c_list)
            c_noi = c_list(c_no);
            tempi = find(cellmask == c_noi);
            cellmask_t(tempi) = 0;                   %cells belonging to very small clusters not to be marked
        end
        
        
        continue
    else
    end
    
    c_lengths = [c_lengths; length(c_list)];
    cr_list = c(c_list, clust_no);                  %list of corrcoeffs of cells in c_list, with their own cluster's averaged trace
    curr_trace_mat = X(c_list,:);       %matrix of activity data for each cell that belongs to this cluster
    
    %assigning cluster numbers to ROIs
    for c_no = 1:length(c_list)
        c_noi = c_list(c_no);
        tempi = find(cellmask == c_noi);
        cellmask_t(tempi) = clust_no;
    end
    
    curr_trace_mat = [cr_list,curr_trace_mat];     %concatenating corrcoefs to activity traces to sort both together
    curr_trace_mat_sorted = sortrows(curr_trace_mat, -1);      %sorting rows by first column ie corrcoefs
    bk_data_mat_o = [bk_data_mat_o; curr_trace_mat_sorted(:, 2:end)];
end

%calculating new corrcoeff mat with re-arranged cells
%figure(3)
%imagesc(bk_data_mat_o');

c_o = corrcoef(bk_data_mat_o', 'rows', 'pairwise'); %sorted
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
reshaped_c_oi = reshape(c_oi, [], 1);
temp = find(isnan(reshaped_c_oi) == 1);
reshaped_c_oi(temp) = [];
saved_corr_scores_inter = [saved_corr_scores_inter; reshaped_c_oi];

corr_unsorted=corrcoef(X','rows','pairwise');

tempi = find(cellmask == 0); % so that only background is selected and not unclustered cells
cellmask_t(tempi) = (max(cell_gp_vec_orig(:, 1))+2);

%drawnow
%imagesc(c_o)
drawnow
figure(2)
imagesc(c_o);
title(['number of fields - ' int2str(max(cell_gp_vec_orig(:, 1)))])
%figure(5)
%imagesc(bk_data_mat_o')
figure(4)
imagesc(cellmask_t)
colormap hot
figure(5)
imagesc(corr_unsorted)
title('Unsorted')
%
% beep

% save('/Users/ananth/Desktop/trialshuffled_spont_sorted_corrcoef.mat', 'c_o')
% save('/Users/ananth/Desktop/trialshuffled_spont_unsorted_corrcoef.mat', 'corr_unsorted')
% save('/Users/ananth/Desktop/trialshuffled_spont_cellmask.mat', 'cellmask_t')


% load('/Users/ananth/Desktop/spont_sorted_corrcoef.mat')
% load('/Users/ananth/Desktop/spont_unsorted_corrcoef.mat')
% load('/Users/ananth/Desktop/spont_cellmask.mat')

% save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat', 'a');
% save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat', 'b');
% save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat', 'c');