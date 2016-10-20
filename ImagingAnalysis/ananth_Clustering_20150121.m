% 2. Clustering
% Use only after generating the dF/F calcium traces
% - Kambadur Ananthamurthy

%clear all
close all

% Execution modules
shuffleTrials=0; % boolean
getDataset=0; % boolean

nphases=1; %the number of phases in the trials to analyze, e.g. - pre-stim, post-stim, etc.
tau = 0.5;  % time-course of calcium transient in seconds
k = 3;  % number of clusters
stimOn=10 ; % in s, when stimulus onset occurs
stimduration=1; %in s, the duration of the stimulus

fid = fopen('/Users/ananth/Desktop/Work/Imaging/analyze.txt');  % in analyze.txt, add url till the ROIindex folder name, ending, without "/"
while 1 %fid is valid when 1
    direc = fgetl(fid);
    
    if ~ischar(direc)
        break
    else
    end
    
    if getDataset==1
        %dataset_date=('20131218');
        date = findstr(direc, '201');
        dataset_date = direc(date:date+7);
        
        %dataset_mousename=('Mouse4');
        %uscorei=findstr(direc, '_');
        fslashi=findstr(direc,'/');
        dataset_mousename=direc(fslashi(length(fslashi)):(direc-4));
        
        dataset_ROIindex=('-ROI');
        %dataset_ROIindex=direc((uscorei+1):length(direc));
        %dataset = [dataset_date '/' dataset_mousename '/' dataset_ROIindex];
        dataset=[dataset_date '/' dataset_mousename];
        mkdir([save_direc '/' dataset]);
        
        tiftag = imfinfo([direc '/1_Trial-1_ROI-1.tif']); %use the tiff file to be used to make cellmask
        img_desc = tiftag.ImageDescription;
        
        framei = findstr(img_desc, 'Frames');
        nframes = str2num(img_desc((framei+9):(framei+12)));
        %nframes=266;
        
        samplingrate = nframes/triallength; % Sampling rate
        
        save_direc = ('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis');
        dataset_='20140312/Mouse6/Mouse6-ROI';
        
        load([save_direc '/' dataset_date '/toneTrials.mat'], 'toneTrials');
        load([save_direc '/' dataset_date '/puffTrials.mat'], 'puffTrials');
        load([save_direc '/' dataset_date '/cellmask.mat'], 'cellmask');
    else
    end
    
    % Meta k means
    %first get the optimal corrcut (cutoff)
    
    nstd = 2;
    nstdlow = 0.5;
    minoncount = 3; % min number of consecutive frames for a calcium transient
    trlength=20; % trial length in seconds
    nframes=191;
    trials=50;
    
    samp = samplingrate; % sampling rate
    winfilt = ceil(15/samp); % Filtering half window
    dcount = 0;
    
    % Response times (||| spike times)
    
    thresh = 2;
    timepc = cell(max(max(cellmask)), trials);
    timep = zeros(max(max(cellmask)), trials, nframes);
    
    % Correlation
    
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
    
    spontaneousPeriod=1:(floor(stimOn*samplingrate));
    stimPeriod=floor(stimOn*samplingrate):ceil((stimOn*samplingrate)+samplingrate);
    allstimPeriod=floor(stimOn*samplingrate):size(toneTrials,3);
    onlystimPeriod=(floor((stimOn*samplingrate):((stimOn*samplingrate)+(5*ceil(samplingrate))))); % everything from stimulus onset to 5s after stimulus
    onlypoststimPeriod=(floor((stimOn*samplingrate)+(stimduration*samplingrate)+samplingrate):size(toneTrials,3)); % everything from 1 s after the offset of the stimulus
    
    tone1=toneTrials(:,:,spontaneousPeriod); %spontaneous periods in tone trials
    tone2=toneTrials(:,:,allstimPeriod); %all post-stim periods in tone trials
    
    puff1=puffTrials(:,:,spontaneousPeriod); %spontaneous periods in puff trials
    puff2=puffTrials(:,:,allstimPeriod); %post-stim periods in puff trials
    puff3=puffTrials(:,:,onlystimPeriod); %post-stim periods, skipping the stim
    puff4=puffTrials(:,:,onlypoststimPeriod); %post-stim periods, skipping the stim
    
    spont=[tone1 puff1];
    for i=1:nphases
        if i==1
            %Z=tone1;
            Z=puff3;
        elseif i==2
            Z=tone2;
        elseif i==3
            Z=puff1;
        else
            Z=puff2;
        end
        
        % Trial-shuffled calbdf dataset
        if shuffleTrials==1
            rand('seed', 0);
            %Z = puff2 * 0;
            for c = 1 : max(max(cellmask))
                Z(c,(1:(nsessiontrials-1)), (1:size(Z,3))) = Z(c, randperm((nsessiontrials-1)), :);
            end
        else
        end
        
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
        
        beep
        
        % save('/Users/ananth/Desktop/trialshuffled_spont_sorted_corrcoef.mat', 'c_o')
        % save('/Users/ananth/Desktop/trialshuffled_spont_unsorted_corrcoef.mat', 'corr_unsorted')
        % save('/Users/ananth/Desktop/trialshuffled_spont_cellmask.mat', 'cellmask_t')
        
        
        % load('/Users/ananth/Desktop/spont_sorted_corrcoef.mat')
        % load('/Users/ananth/Desktop/spont_unsorted_corrcoef.mat')
        % load('/Users/ananth/Desktop/spont_cellmask.mat')
        
        % save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat', 'a');
        % save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat', 'b');
        % save('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/a.mat', 'c');
        
        if nphases==1
        else
            reply = input('Quit looking at datasets? ', 's'); %  Y to quit; Enter to continue
            if numel(reply)>0
                break
            else
            end
        end
    end
end