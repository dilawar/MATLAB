% Everything after developing the non-baseline corrected raw intensity time
% series for all cells, trials and frames; till clustering and analysis of
% the clusters

% Kambadur Ananthamurthy

dataset='20140312/Mouse6/Mouse6-ROI';

load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/' dataset '/cal.mat'])
load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/' dataset '/cellmask.mat'])
load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/' dataset '/calbdf.mat'])
load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/' dataset '/caltr.mat'])

tau = 0.5;  % time-course of calcium transient in seconds
k = 3;  % number of clusters

caltr = zeros(size(cal));
calbaseline = squeeze(mean(cal, 3));
nstd = 2;
nstdlow = 0.5;
minoncount = 3; % min number of consecutive frames for a calcium transient
trlength=20; % trial length in seconds
nframes=191;
trials=50;

samp = trlength/nframes; % sampling rate
winfilt = ceil(15/samp); % Filtering half window
dcount = 0;

% Response times (||| spike times)

thresh = 2;
timepc = cell(max(max(cellmask)), trials);
timep = zeros(max(max(cellmask)), trials, nframes);

dec = exp(-(0 : ceil(2 * tau/samp))/(tau/samp));
pad = zeros(length(dec)-1, 1);

for c = 1 : max(max(cellmask))
    for t = 1 : trials
        temp1 = deconv([squeeze(caltr(c, t, :)); pad], dec);
        %         temp2 = [0; diff(squeeze(caltr(c, t, :)))];
        temp2 = deconv([squeeze(calbdf(c, t, :)); pad], dec);
        
        % find only peaks which are > 2 std above mean
        maxtab = peakfinder(temp2, thresh * nanstd(temp2));
        maxsig = find(temp1 > (thresh * nanstd(temp2)));
        
        if ~isempty(maxsig)
            timep(c, t, intersect(maxtab, maxsig)) = 1;
            timepc {c, t} = intersect(maxtab, maxsig);
        end
    end
end

% Correlation - spontaneous and odour evoked

%afr = [(1:floor(20/samp)), floor(80/samp:90/samp), floor(150/samp:160/samp), floor(220/samp:230/samp)];
afr = (1:floor(20/samp));
acorr = zeros(max(max(cellmask)));
% ocorr = zeros(max(max(cellmask)));
tcorr = zeros(max(max(cellmask)));
nancorr = tcorr;
spcorr = tcorr;

X = calbdf; %change made

for t = 1 : 3
    atcorr = corr(squeeze(X(:, t, (1:95)))', squeeze(X(:, t, (1:95)))');
    %     otcorr = corr(squeeze(caltr(:, t, odstart:odend))', squeeze(caltr(:, t, odstart:odend))');
    ttcorr = corr(squeeze(X(:, t, :))', squeeze(calbdf(:, t, :))');
    sptcorr = corr(squeeze(timep(:, t, :))', squeeze(timep(:, t, :))');
    nancorr = nancorr + (~isnan(ttcorr));
    
    atcorr(isnan(atcorr)) = 0;
    %     otcorr(isnan(otcorr)) = 0;
    ttcorr(isnan(ttcorr)) = 0;
    sptcorr(isnan(sptcorr)) = 0;
    
    acorr = acorr + atcorr;
    %     ocorr = ocorr + otcorr;
    tcorr = tcorr + ttcorr;
    spcorr = spcorr + sptcorr;
    
end

acorr = acorr./nancorr;
% ocorr = ocorr./nancorr;
tcorr = tcorr./nancorr;
sptcorr = sptcorr./nancorr;


% Distance matrix

centroids = regionprops(cellmask, 'centroid');
celld = zeros(max(max(cellmask)));

for c1 = 1 : max(max(cellmask))
    for c2 = 1 : max(max(cellmask))
        cent1 = centroids(c1).Centroid;
        cent2 = centroids(c2).Centroid;
        celld(c1, c2) = sqrt((cent1(1) - cent2(1))^2 + (cent1(2) - cent2(2))^2);
    end
end

% centroid_matrix=cell2mat(struct2cell(centroids)); %row vector of all x and y
% centroid_x=centroid_matrix(1,[1:2:end]); % x coordinates
% centroid_y=centroid_matrix(1,[2:2:end]); % y coordinates
% centroid_xy=[centroid_x' centroid_y'];
% celld=pdist(centroid_xy, 'euclidean'); % distance matrix of all cells

% Clustering

X = reshape(permute(calbdf(:,:,(1:95)), [1 3 2]), [size(calbdf, 1) numel(squeeze(calbdf(1, :, (1:95))))]);

IDX = [];
it = 100;  % iterations for k-means

for i = 1 : it
    %     [I ce] = kmeans(X, k, 'distance', 'correlation');
    [I ce] = kmeans_plusplus(X, k);
    if i == 1        % Template centroids
        CE = ce;
    else            % Centroid matching to templates to preserve cluster identity
        cemap = corr(ce', CE');
        [m ind] = max(cemap);
        tempI = I;
        for x = 1 : k
            tempI(I == x) = ind(x);
        end
        I = tempI;
    end
    IDX = [IDX I];
end

cellmaskc = zeros(size(cellmask, 1), size(cellmask, 2), 3);

for x = 1:max(max(cellmask))
    for y = 1 : k
        temp = cellmask * 0;
        temp (find(cellmask == x)) = sum(squeeze(IDX(x, :) == y))/it;
        cellmaskc(:, :, y) = cellmaskc(:, :, y) + temp;       % Spatial layout
    end
end

% Trial-shuffled calbdf dataset

rand('seed', 0);
calbdfr = calbdf * 0;
for c = 1 : max(max(cellmask))
    calbdfr(c, :, :) = calbdf(c, randperm(trials), :);
end


% Correlation versus distance

corrd = [];
count = 1;

for c1 = 1 : max(max(cellmask)) - 1
    for c2 = (c1 + 1) : max(max(cellmask))
        corrd (count, :) = [tcorr(c1, c2) celld(c1, c2)];
        count = count + 1;
    end
end

% Meta k means
%first get the optimal corrcut (cutoff)
calbdf_spont((1:118),(1:30),(1:95))=calbdf(:,:,(1:95));
%calbdfr_spont((1:118),(1:50),(1:95))=calbdfr(:,:,(1:95));

calbdf_poststim((1:118),(1:30),(1:86))=calbdf(:,:,(106:191));
%calbdf_stimonwards((1:118),(1:50),(1:96))=calbdf(:,:,(96:191));

%tone_trials=[1,2,3,4,5,11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45];
%puff_trials=[6,7,8,9,10,16,17,18,19,20,26,27,28,29,30,36,37,38,39,40,46,47,48,49,50];

ncells=118;

% %for Tone
% for tone=1:length(tone_trials)
%     %tone_session_mouse4_tonenpuff((1:ncells),tone,(1:nframes))=a.calbdf((1:ncells),tone_trials(tone),(1:nframes));
%     tone_session_mouse4_tonenpuff((1:ncells),tone,(1:nframes))=calbdf((1:ncells),tone_trials(tone),(1:nframes));
% end

% %for Puff
% for puff=1:length(puff_trials)
%     %puff_session_mouse4_tonenpuff((1:ncells),puff,(1:nframes))=a.calbdf((1:ncells),puff_trials(puff),(1:nframes));
%     puff_session_mouse4_tonenpuff((1:ncells),puff,(1:nframes))=calbdf((1:ncells),puff_trials(puff),(1:nframes));
% end

%calbdf_pretone((1:118),(1:25),(1:95))=tone_session_mouse4_tonenpuff(:,:,(1:95));
%calbdf_prepuff((1:118),(1:25),(1:95))=puff_session_mouse4_tonenpuff(:,:,(1:95));

%calbdf_posttone((1:118),(1:25),(1:86))=tone_session_mouse4_tonenpuff(:,:,(106:191));
%calbdf_postpuff((1:118),(1:25),(1:86))=puff_session_mouse4_tonenpuff(:,:,(106:191));


calbdf_spont_stitched(1:118,(1:15),(1:95))=calbdf_spont(:,(1:15),:);
calbdf_spont_stitched(1:118,(16:49),(1:95))=calbdf_spont(:,(17:50),:);

% Trial-shuffled calbdf dataset

rand('seed', 0);
calbdfr_spont_stitched = calbdf_spont_stitched * 0;
for c = 1 : max(max(cellmask))
    calbdfr_spont_stitched(c,(1:49), (1:95)) = calbdf_spont_stitched(c, randperm(49), :);
end


calbdf_prepuff_stitched(1:118,(1:5),(1:95))=calbdf_prepuff(:,(1:5),:);
calbdf_prepuff_stitched(1:118,(6:24),(1:95))=calbdf_prepuff(:,(7:25),:);

calbdf_stimonwards_stitched(1:118,(1:15),(1:96))=calbdf_stimonwards(:,(1:15),:);
calbdf_stimonwards_stitched(1:118,(16:49),(1:96))=calbdf_stimonwards(:,(17:50),:);

calbdf_poststim_stitched(1:118,(1:15),(1:86))=calbdf_poststim(:,(1:15),:);
calbdf_poststim_stitched(1:118,(16:49),(1:86))=calbdf_poststim(:,(17:50),:);

calbdf_postpuff_stitched(1:118,(1:5),(1:86))=calbdf_postpuff(:,(1:5),:);
calbdf_postpuff_stitched(1:118,(6:24),(1:86))=calbdf_postpuff(:,(7:25),:);

calbdf_postpuff_stitched_end((1:118),(1:24),(1:43))=calbdf_postpuff_stitched((1:118),:,(44:86));
calbdf_postpuff_stitched_begin((1:118),(1:24),(1:43))=calbdf_postpuff_stitched((1:118),:,(1:43));
calbdf_postpuff_stitched_fagend((1:118),(1:24),(1:19))=calbdf_postpuff_stitched((1:118),:,(68:86));
calbdf_postpuff_stitched_secondlastend((1:118),(1:24),(1:19))=calbdf_postpuff_stitched((1:118),:,(48:66));

Z=calbdf_spont_stitched;

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

% save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat','a')
% save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/b.mat','b')
% save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/c.mat','c')

%load('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat')
%load('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/b.mat')
%load('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/c.mat')

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
% drawnow
% figure(4)
% imagesc(c_o);
% title(['number of fields - ' int2str(max(cell_gp_vec_orig(:, 1)))])
% %figure(5)
% %imagesc(bk_data_mat_o')
% figure(6)
% imagesc(cellmask_t)
% colormap hot
% figure(7)
% imagesc(corr_unsorted)
% title('Unsorted')
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