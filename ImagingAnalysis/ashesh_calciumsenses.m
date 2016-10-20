% File: calciumsenses.m
% Function: Generate a df/f calcium trace for every cell in a data set
% Author: Ashesh Dhawale 
% Date: 11/01/2011

clear all; 

javaaddpath 'C:\Program Files\MATLAB\R2008b\java\ij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2008b\java\mij.jar'

tau = 0.5;  % time-course of calcium transient in seconds
trlength = 320;  % trial length in seconds

trials = 4;

loadcal = 1; % load 'cal' instead of computing it

nroi = 1;
frstart = 4;

k = 2;  % number of clusters

mskfile = 1; % mask made from which file - trial #

% Segment and label cell mask
msk = imread('map.tif');
msk1 = watershed(-bwdist(~msk));
msk(msk1 == 0) = 0;
cellmsk = bwlabel(msk);

linker = '1_Trial-';
tiftag = imfinfo([linker '1_ROI-1.tif']);
nframes = numel(tiftag);

ncell = max(max(cellmsk)); % # of cells
samp = trlength/nframes; % Sampling rate
            
rframes = 500;

% Load ref image
img = zeros(rframes, tiftag(1).Height, tiftag(1).Width);

for count = 1 : rframes
    img (count, :, :) = imread([linker num2str(mskfile) '_ROI-1.tif']);
%     img (count, :, :) = imread('F:\NCBS\Data\20110814\tone-puff_2-ROI\Trial - 1-1.tif');
end
refimg = squeeze(median(img, 1));

% refimg = imread('refimg.tif');

cal = zeros(max(max(cellmsk)), trials, nframes);
caldf = cal;
caldfsig = cal;
calbdf = cal;
calbdfsig = cal;

if loadcal
    load('cal.mat');
    cal = cal(:, 1:trials, :);
else
    
    MIJ.start ('F:\CSHL\Code\Fiji.app\');  % Start ImageJ
    MIJ.createImage('Refimg', int16(refimg), 1);

    for tr = 1 : trials
        tr
        for r = 1 : nroi

            tiftag = imfinfo([linker num2str(tr) '_ROI-' num2str(r) '.tif']);

            % Read in image sequence. 
%             img = zeros(nframes, tiftag(1).Height, tiftag(1).Width);            
             
        
            for count = 1 : nframes
                tempimg = imread([linker num2str(tr) '_ROI-' num2str(r) '.tif'], count);
                  
    % % Register img to refimg frame by frame by cross corr
    %             cc = normxcorr2(tempimg, refimg);
    %             [max_cc, imax] = max(abs(cc(:)));
    %             [ypeak, xpeak] = ind2sub(size(cc),imax(1));
    %             corr_offset = [(xpeak-size(refimg,2)) (ypeak-size(refimg,1))];
    % 
    %             tform = maketform('affine',[1 0 0; 0 1 0; corr_offset(1) corr_offset(2) 1]);
    %             img (count, :, :) = imtransform(tempimg, tform, 'XData', [1 size(refimg, 2)], 'YData', [1 size(refimg, 1)]);

    % Register img to refimg using Turboreg in ImageJ

                MIJ.createImage(['Img' num2str(count)], int16(tempimg), 1);

                MIJ.run('TurboReg ', ['-align -window Img', num2str(count), ' 0 0 ', num2str(tiftag(1).Width-1),...
                        ' ', num2str(tiftag(1).Height-1), ' -window Refimg 0 0 ', ...
                        num2str(tiftag(1).Width-1), ' ', num2str(tiftag(1).Height-1),...
                        ' -translation ', num2str(tiftag(1).Width/2), ' ', ...
                        num2str(tiftag(1).Height/2), ' ', num2str(tiftag(1).Width/2), ...
                        ' ', num2str(tiftag(1).Height/2), ' -showOutput']);
                     
                temp = MIJ.getImage ('Output');
                tempimg2 = squeeze(uint16(temp(:, :, 1)));
                MIJ.run('Close', 'Output');
                MIJ.run('Close', ['Img' num2str(count)]);
                
                % Generate intensity traces for all cells
                for c = 1 : max(max(cellmsk))
                  cal(c, tr, count) = squeeze(mean(mean(tempimg2(cellmsk == c), 2), 1));
                end
            end


%             % Generate intensity traces for all cells
%             for c = 1 : max(max(cellmsk))
%                 if min(min(min(img(:, cellmsk == c)))) > 0
%                     cal(c, tr, :) = squeeze(mean(mean(img(:, cellmsk == c), 3), 2));
%                 else
%                     cal(c, tr, :) = zeros(nframes, 1) * NaN;
%                 end
%             end

        end

    end

    MIJ.exit;
    
    for c = 1 : ncell
        for tr = 1 : trials
            if min(squeeze(cal(c, tr, :))) == 0
                cal(c, tr, :) = zeros(nframes, 1) * NaN;
            end
        end
    end
end

caltr = zeros(size(cal));
calbaseline = squeeze(mean(cal, 3));
nstd = 2;
nstdlow = 0.5;
minoncount = 3; % min number of consecutive frames for a calcium transient
winfilt = ceil(15/samp); % Filtering half window
dcount = 0;

for c = 1 : max(max(cellmsk))
    c;
    for t = 1 : trials
        caldf(c, t, :) = (cal(c, t, :) - calbaseline(c, t))/calbaseline(c, t);
        caldfsig(c, t, :) = caldf(c, t, :) > (nstd * std(squeeze(caldf(c, t, :))));

        % Baseline correction
        [maxtab mintab] = peakdet(squeeze(caldf(c, t, 2:end)), 0.2);
        temp = squeeze(cal(c, t, 1:end));
        for m = 1 : size(maxtab, 1)
            temp(max(uint32(maxtab(m,1)-5)+1, 1): min(maxtab(m,1)+30, nframes)) = NaN;
        end
        
        % Window filtering
        cabase = zeros(nframes, 1);
        
        for count = winfilt + 1 : nframes - winfilt - 1
%             cabase (count) = nanmean(temp(count-winfilt:count+winfilt));
            
            % 10 percentile value correction
            cabase (count) = prctile(squeeze(cal(c, t, count-winfilt:count+winfilt)), 10);
        end
        
        cabase(1:winfilt) = cabase(winfilt+1);
        cabase(end-winfilt:end) = cabase(end-winfilt-1);
        
        calbdf(c, t, :) = (squeeze(cal(c, t, :)) - cabase)./cabase;
        
%         % Exponential fit
%         x = 1:length(temp);
%         xf = x(~isnan(temp));
%         tempf = temp(~isnan(temp));
%         
%         if numel(xf) > 1
%             tfit = fit(xf', tempf, 'exp1');
%             cabase = feval(tfit, x);
%         else
%             cabase = zeros(nframes, 1) * NaN;
%         end
%         
        castd = nanstd(temp - cabase);


        % Calculation of SD from negative going transients
%         
%         temp2 = squeeze(calbdf(c, t, :));
%         castd = nanstd([(temp2(temp2<=0)); (-temp2(temp2<0))]);
                
%         calbdfsig(c, t, :) = abs(calbdf(c, t, :)) > (nstd * castd/cabase);
        
        % Generate signifcant traces only in caltr
        
        count = 1;
        on = 0;
        oncount = 0;
        
        for count = 1 : nframes
            if on
                tempstd = nstdlow;
            else
                tempstd = nstd;
            end
            if calbdf(c, t, count) < -(tempstd * castd/cabase(count))
                if ~on
                    dcount = dcount + 1;
                end
                on = 1;
                caltr(c, t, count) = calbdf(c, t, count);
                oncount = oncount + 1;
            else
                % Impose min trlength
                if on && oncount < minoncount
                    caltr(c, t, count-oncount:count-1) = 0;
                end
                oncount = 0;
                
                on = 0;
                
            end
        end
        
    end
end



% Response times (||| spike times)

thresh = 2;
timepc = cell(max(max(cellmsk)), trials);
timep = zeros(max(max(cellmsk)), trials, nframes);

dec = exp(-(0 : ceil(2 * tau/samp))/(tau/samp));
pad = zeros(length(dec)-1, 1);

for c = 1 : max(max(cellmsk))
    for t = 1 : trials
        temp1 = deconv([squeeze(caltr(c, t, :)); pad], dec);
%         temp2 = [0; diff(squeeze(caltr(c, t, :)))]; 
        temp2 = deconv([squeeze(calbdf(c, t, :)); pad], dec);

        % find only peaks which are > 2 std above mean
        maxtab = peakfinder(temp1, thresh * nanstd(temp2));
        maxsig = find(temp1 > (thresh * nanstd(temp2)));
        
        if ~isempty(maxsig)
            timep(c, t, intersect(maxtab, maxsig)) = 1;
            timepc {c, t} = intersect(maxtab, maxsig);
        end
    end
end


% Frame to frame movement computation

movt = zeros(trials, nframes, 2);

if loadcal
    load('movt.mat');
    movt = movt(1:trials,:, :);
else
    for tr = 1 : trials
        tr
        for count = 1 : nframes
                tempimg = imread([linker num2str(tr) '_ROI-1.tif'], count);
                cc = normxcorr2(tempimg, refimg);
                movt(tr, count, 1) = cc(size(refimg,1), size(refimg,2));
                [max_cc, imax] = max(abs(cc(:)));
                [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                movt(tr, count, 2) = sqrt((xpeak-size(refimg,2))^2 + (ypeak-size(refimg,1))^2);
        end
    end
end

% Sensory responses, and movt

stime = 1; % stimulus duration in sec
siti = 5; % inter stimulus interval in sec
sstart = [20 90 160 230]; % odour tone puff light

wind = 3; % half window size in seconds
windfr = ceil(wind/samp);
calsens = zeros(4, 10, ncell, trials, 2*windfr+1);
calsensmov = zeros(4, 10, trials, 2*windfr+1);

for stim = 1 : 4
    for p = 1 : 10
        temp = round((((p-1) * (siti+stime)) + sstart(stim))/samp);
        calsens(stim, p, :, :, :) = timep(:, :, temp-windfr:temp+windfr);
        calsensmov(stim, p, :, :) = squeeze(movt(:, temp-windfr:temp+windfr, 1));
    end
end



% Response Statistics

% Mean fluorescence

sareap = zeros(5, max(max(cellmsk)), 2);
sareap2 = zeros(4, max(max(cellmsk)), 2);
% sfr = ceil(stime/samp);
sfr = 7;        % ~500 ms

for c = 1 : max(max(cellmsk))
    for stim = 1 : 4
        temp = reshape(permute(calsens(stim, :, c, :, :), [1 3 2 4 5]), trials*10, 2*windfr+1);
        
        sareap(stim, c, 1) = squeeze(mean(mean(temp(:, windfr+1:windfr+sfr+1))) - mean(mean(temp(:, 1:1+sfr))));
        [h sareap(stim, c, 2)] = kstest2(squeeze(mean(temp(:, windfr+1:windfr+sfr+1), 2)), squeeze(mean(temp(:, 1:1+sfr), 2)));
        
        sareap2(stim, c, 1) = squeeze(mean(mean(temp(:, 1+sfr+(windfr+1:windfr+sfr+1)))) - mean(mean(temp(:, 1:1+sfr))));
        [h sareap2(stim, c, 2)] = kstest2(squeeze(mean(temp(:, 1+sfr+(windfr+1:windfr+sfr+1)), 2)), squeeze(mean(temp(:, 1:1+sfr), 2)));
    end
    sareap(5, c, 1) = squeeze(nanmean(mean(temp(:, 1:1+sfr))) - mean(mean(temp(:, 1+sfr+(1:1+sfr)))));
    [h sareap(5, c, 2)] = ttest(squeeze(mean(temp(:, 1:1+sfr), 2)), squeeze(mean(temp(:, 1+sfr+(1:1+sfr)), 2)));
end

slist = sareap(1:4, :, 2)<0.05;
alist = sareap(5, :, 2)<0.05;


% Using multiple air-periods

% sfr = ceil(stime/samp);
sfr = 7;        % ~500 ms
sareap = zeros(5, max(max(cellmsk)), length(1:sfr:windfr-sfr-1), 2);
sareap2 = zeros(4, max(max(cellmsk)), length(1:sfr:windfr-sfr-1), 2);

for c = 1 : max(max(cellmsk))
    for stim = 1 : 4
        count = 1;
        for fr = 1:sfr: windfr-sfr-1
            temp = reshape(permute(calsens(stim, :, c, :, :), [1 3 2 4 5]), trials*10, 2*windfr+1);

            sareap(stim, c, count, 1) = squeeze(mean(mean(temp(:, windfr+1:windfr+sfr+1))) - mean(mean(temp(:, fr:fr+sfr))));
            [h sareap(stim, c, count, 2)] = kstest2(squeeze(mean(temp(:, windfr+1:windfr+sfr+1), 2)), squeeze(mean(temp(:, fr:fr+sfr), 2)));

            sareap2(stim, c, count, 1) = squeeze(mean(mean(temp(:, 1+sfr+(windfr+1:windfr+sfr+1)))) - mean(mean(temp(:, fr:fr+sfr))));
            [h sareap2(stim, c, count, 2)] = kstest2(squeeze(mean(temp(:, 1+sfr+(windfr+1:windfr+sfr+1)), 2)), squeeze(mean(temp(:, fr:fr+sfr), 2)));
    
            sareap(5, c, count, 1) = squeeze(nanmean(mean(temp(:, 1:1+sfr))) - mean(mean(temp(:, fr:fr+sfr))));
            [h sareap(5, c, count, 2)] = ttest(squeeze(mean(temp(:, 1:1+sfr), 2)), squeeze(mean(temp(:, fr:fr+sfr), 2)));
    
            count = count + 1;
        end
    end
    
end

slist = ((median(sareap(1:4, :, :, 2),3)<0.05) + (median(sareap2(1:4, :, :, 2),3)<0.05))>0;
alist = sareap(5, :, 2)<0.05;


% Correlation - spontaneous and odour evoked

afr = [(1:floor(20/samp)), floor(80/samp:90/samp), floor(150/samp:160/samp), floor(220/samp:230/samp)]; 
acorr = zeros(max(max(cellmsk)));
% ocorr = zeros(max(max(cellmsk)));
tcorr = zeros(max(max(cellmsk)));
nancorr = tcorr;
spcorr = tcorr;

X = caltr;

for t = 1 : 3
    atcorr = corr(squeeze(X(:, t, afr))', squeeze(X(:, t, afr))');
%     otcorr = corr(squeeze(caltr(:, t, odstart:odend))', squeeze(caltr(:, t, odstart:odend))');
    ttcorr = corr(squeeze(X(:, t, :))', squeeze(caltr(:, t, :))');
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


% Event based and continuous trace cross correlation matrices

wind = 4;  % half-window size in seconds
windfr = ceil(wind/samp); % half-window size in frames
crcorr = zeros(ncell, ncell, 2*windfr+1);

for c1 = 1 : ncell
    for c2 = 1 : ncell
        sph = zeros(1, 2*windfr+1);
        for t = 1 : trials
            % Event based
            spikes1 = timepc{c1, t};
            spikes2 = timepc{c2, t};
            for sp = 1 : length(spikes1)
                if (spikes1(sp) > windfr) && (spikes1(sp) < (nframes-windfr))
                    sptemp = spikes2(spikes2 < (spikes1(sp)+windfr) & spikes2 > (spikes1(sp)-windfr)) - spikes1(sp);
                                        %                     sphist = hist(sptemp, [-windfr:windfr]);
                    sph(sptemp + windfr + 1) = sph(sptemp + windfr + 1) + 1;
                end
            end
        end
        crcorr(c1, c2, :) = sph;
    end
end

% Avalanche size and prob

dt = 1; % frame bin

aval = zeros(1, 5000);

for t = 1 : trials
    avsize = zeros(1, floor(nframes/dt));
    count = 1;
    for bin = 1 : dt : (floor(nframes/dt)-1)*dt
        avsize(count) = sum(sum(squeeze(timep(:, t, [bin:bin+dt]))));
        count = count + 1;
    end
    zav = find(avsize == 0);
    for z = 1 : (numel(zav)-1)
        temp = sum(avsize(zav(z):zav(z+1)));
        aval(temp + 1) = aval(temp + 1) + 1;
    end
end


% Distance matrix

centroids = regionprops(cellmsk, 'centroid');
celld = zeros(max(max(cellmsk)));

for c1 = 1 : max(max(cellmsk))
    for c2 = 1 : max(max(cellmsk))
        cent1 = centroids(c1).Centroid;
        cent2 = centroids(c2).Centroid;
        celld(c1, c2) = sqrt((cent1(1) - cent2(1))^2 + (cent1(2) - cent2(2))^2);
    end
end
    
% Read in Behaviour data

beh = zeros(trials, trlength*1000);

for tr = 1 : trials
    temp = dlmread([linker num2str(tr) '.xls']);
    beh (tr, :) = temp (:, 1)';
end


% Clustering

X = reshape(permute(caltr(:,1:2,:), [1 3 2]), [size(caltr, 1) numel(squeeze(caltr(1, 1:2, :)))]);

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

cellmskc = zeros(size(cellmsk, 1), size(cellmsk, 2), 3);

for x = 1:max(max(cellmsk))
    for y = 1 : k
        temp = cellmsk * 0;
        temp (find(cellmsk == x)) = sum(squeeze(IDX(x, :) == y))/it;
        cellmskc(:, :, y) = cellmskc(:, :, y) + temp;       % Spatial layout
    end
end

% Clustering a segment of the data

caltrp = permute(timep, [1 3 2]);
X = reshape(caltrp(:, 100:end, :), [size(caltr, 1) numel(squeeze(caltr(1, :, 100:end)))]);

IDX = [];
it = 100;  % iterations for k-means

for i = 1 : it
    [I ce] = kmeans(X, k, 'distance', 'correlation');
    if i ~= 1
        
        cemap = zeros(k, k);        % Matching to existing clusters
        for p = 1 : k
            tI = I; tI(I==p) = 1; tI(I~=p) = 0;
            for q = 1 : k
                cI = squeeze(IDX(:,1));
                tId = squeeze(IDX(:,1)); tId(cI==q) = 1; tId(cI~=q) = 0;
                cemap(p, q) = corr(tI, tId);
            end
        end
        
        [m ind] = max(cemap);
        tempI = I;
        for x = 1 : k
            tempI(I == x) = ind(x);
        end
        I = tempI;
        
    end
    IDX = [IDX I];
end

cellmskc = zeros(size(cellmsk, 1), size(cellmsk, 2), 3);

for x = 1:max(max(cellmsk))
    for y = 1 : k
        temp = cellmsk * 0;
        temp (cellmsk == x) = sum(squeeze(IDX(x, :) == y))/it;
        cellmskc(:, :, y) = cellmskc(:, :, y) + temp;       % Spatial layout
    end
end

% Single run of kmeans

ind = [];
stim = 3;

for i = 1 : 10
    ind = [ind (round(((i-1)*6 + sstart(stim))/samp):round(((i-1)*6 + sstart(stim) + 1)/samp))];
end
        
X = reshape(permute(caltr(:, 1, ind), [1 3 2]), [size(caltr, 1) numel(squeeze(caltr(1, 1, ind)))]);

[IDX ce] = kmeans(X, k, 'distance', 'correlation');
[temp I] = sort(IDX);

cellmskc = cellmsk * 0;
for x = 1 : max(max(cellmsk))
    cellmskc(cellmsk == x) = IDX(x);
end

% Trial-shuffled calbdf dataset

rand('seed', 0);
calbdfr = calbdf * 0;
for c = 1 : max(max(cellmsk))
    calbdfr(c, :, :) = calbdf(c, randperm(trials), :);
end


% Correlation versus distance

corrd = [];
count = 1;

for c1 = 1 : max(max(cellmsk)) - 1
    for c2 = (c1 + 1) : max(max(cellmsk))
        corrd (count, :) = [tcorr(c1, c2) celld(c1, c2)];
        count = count + 1;
    end
end
        

% Meta k means

X = reshape(permute(caltr(:, 1:3, :), [1 3 2]), [size(caltr, 1) numel(squeeze(caltr(1, 1:3, :)))]);

ccorr = tcorr; %(2:end, 2:end);
corrcut = 0.7;
thr = 0.8;
clustat = [];

for corrcut = 0.4:0.05:0.9
    corrcut
    [a b c] = meta_k_means_tank(X, 'correlation', corrcut, thr);
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
    
    
cellmskc = cellmsk*0;
for i = 1 : size(a, 1)
    temp = a{i, 1};
    for x = 1 : length(temp)
        cellmskc(cellmsk == (temp(x))) = i;
    end
end

% Dombeck et al's heuristic to determine correct correlation cutoff

x = clustat(:, 1);
x = (x - min(x))/(max(x)-min(x));
y = clustat(:, 2);
y = (y - min(y))/(max(y)-min(y));
x.*y
