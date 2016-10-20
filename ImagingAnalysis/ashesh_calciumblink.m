% File: calciumblink.m
% Function: Generate a df/f calcium trace for every cell in a data set
% Author: Ashesh Dhawale 
% Date: 11/01/2011

clear all; 

javaaddpath '/Users/kambadurananthamurthy/Documents/MATLAB/java/ij.jar'
javaaddpath '/Users/kambadurananthamurthy/Documents/MATLAB/java/mij.jar'

tau = 0.5;  % time-course of calcium transient in seconds
trlength = 15;  % trial length in seconds

trials = 50;

% exlist = [64 70 80 125];
exlist = [];

loadcal = 1; % load cal instead of computing it

nroi = 1;
frstart = 4;
cstart = 30;
ustart = 38;
sdur = 7;
astart = 20;

k = 3;  % number of clusters

mskfile = 1; % mask made from which file - trial #

% Segment and label cell mask
msk = imread('map.tif');
msk1 = watershed(-bwdist(~msk));
msk(msk1 == 0) = 0;
cellmsk = bwlabel(msk);

ncellold = max(max(cellmsk));

for x = 1 : length(exlist)
    cellmsk(cellmsk == exlist(x)) = 0;
end

cellmsk = bwlabel(cellmsk > 0);

linker = 'Trial - ';
tiftag = imfinfo([linker '0-1.tif']);
nframes = numel(tiftag);

ncell = max(max(cellmsk)); % # of cells
samp = trlength/nframes; % Sampling rate
            
% Load ref image
img = zeros(nframes, tiftag(1).Height, tiftag(1).Width);
for count = 1 : nframes
    img (count, :, :) = imread([linker num2str(mskfile) '-1.tif']);
%     img (count, :, :) = imread('F:\NCBS\Data\20110814\tone-puff_2-ROI\Trial - 1-1.tif');
end
refimg = squeeze(median(img, 1));

cal = zeros(max(max(cellmsk)), trials, nframes);
caldf = cal;
caldfsig = cal;
calbdf = cal;
calbdfsig = cal;

if loadcal
    load('cal.mat');
    cal = cal(~ismember(1:ncellold, exlist), 1:trials, :);

else
    
    MIJ.start ('F:\CSHL\Code\Fiji.app\');  % Start ImageJ
    MIJ.createImage('Refimg', int16(refimg), 1);

    for tr = 1 : trials
        tr
        for r = 1 : nroi

            tiftag = imfinfo([linker num2str(tr-1) '-' num2str(r) '.tif']);

            % Read in image sequence. 
            img = zeros(nframes, tiftag(1).Height, tiftag(1).Width);            
            for count = 1 : nframes
                tempimg = imread([linker num2str(tr-1) '-' num2str(r) '.tif'], count);

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
                img (count, :, :) = squeeze(uint16(temp(:, :, 1)));
                MIJ.run('Close', 'Output');
                MIJ.run('Close', ['Img' num2str(count)]);
    %      
            end


            % Generate intensity traces for all cells
            for c = 1 : max(max(cellmsk))
                if min(min(min(img(:, cellmsk == c)))) > 0
                    cal(c, tr, :) = squeeze(mean(mean(img(:, cellmsk == c), 3), 2));
                else
                    cal(c, tr, :) = zeros(nframes, 1) * NaN;
                end
            end

        end

    end

    MIJ.exit;
end

caltr = zeros(size(cal));
calbaseline = squeeze(mean(cal, 3));
nstd = 2;
nstdlow = 0.5;
minoncount = 3; % min number of consecutive frames for a calcium transient
dcount = 0;

for c = 1 : max(max(cellmsk))
    for t = 1 : trials
        caldf(c, t, :) = (cal(c, t, :) - calbaseline(c, t))/calbaseline(c, t);
        caldfsig(c, t, :) = caldf(c, t, :) > (nstd * std(squeeze(caldf(c, t, :))));

        % Baseline correction
        [maxtab mintab] = peakdet(squeeze(caldf(c, t, :)), 0.2);
        temp = squeeze(cal(c, t, :));
        for m = 1 : size(maxtab, 1)
            temp(uint32(maxtab(m,1)-5)+1: maxtab(m,1)+30) = NaN;
        end
%         cabase = nanmean(temp(2:end));
        
        % 10 percentile value correction
        cabase = prctile(squeeze(cal(c, t, :)), 10)';
        
%         castd = nanstd(temp(2:end));
        castd = nanstd(temp-cabase);
        
        calbdf(c, t, :) = (cal(c, t, :) - cabase)/cabase;
        calbdfsig(c, t, :) = abs(calbdf(c, t, :)) > (nstd * castd/cabase);
        
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
            if calbdf(c, t, count) > (tempstd * castd/cabase)
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



% for c1 = 1 : max(max(cellmsk))
%     for c2 = 1 : max(max(cellmsk))
%         
%         atcorr = 0; otcorr = 0; ttcorr = 0;
%         t = 0;
%         
%         for t = 1 : trials
%             atemp = corr(squeeze(calbdf(c1, t, frstart:(odstart-1))), squeeze(calbdf(c2, t, frstart:(odstart-1))));
%             otemp = corr(squeeze(calbdf(c1, t, odstart:odend)), squeeze(calbdf(c2, t, odstart:odend)));
%             ttemp = corr(squeeze(calbdf(c1, t, :)), squeeze(calbdf(c2, t, :)));
%             if ~isnan(atemp)
%                 atcorr = atcorr + atemp;
%                 otcorr = otcorr + otemp;
%                 ttcorr = ttcorr + ttemp;
%                 t = t + 1;
%             end
%         end
%         
%         acorr(c1, c2) = atcorr/t;
%         ocorr(c1, c2) = otcorr/t;
%         tcorr(c1, c2) = ttcorr/t;
%     end
% end


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

        % find only peaks which are > (thresh) std above mean
        maxtab = peakfinder(temp1, thresh * nanstd(temp2));
        maxsig = find(temp1 > (thresh * nanstd(temp2)));
        
        if ~isempty(maxsig)
            timep(c, t, intersect(maxtab, maxsig)) = 1;
            timepc {c, t} = intersect(maxtab, maxsig);
        end
    end
end


% Correlation - spontaneous and odour evoked

acorr = zeros(max(max(cellmsk)));
cscorr = zeros(max(max(cellmsk)));
uscorr = zeros(max(max(cellmsk)));
tcorr = cscorr;
nanccorr = cscorr; nanucorr = cscorr; nantcorr = cscorr; nanscorr = cscorr; 
nanacorr = cscorr;
spcorr = cscorr;

X = caltr;

for t = 1 : trials
    atcorr = corr(squeeze(X(:, t, astart:astart+sdur))', squeeze(X(:, t, astart:astart+sdur))');
    cstcorr = corr(squeeze(X(:, t, cstart:(cstart+sdur)))', squeeze(X(:, t, cstart:(cstart+sdur)))');
    ustcorr = corr(squeeze(X(:, t, ustart:ustart+sdur))', squeeze(X(:, t, ustart:ustart+sdur))');
    ttcorr = corr(squeeze(X(:, t, :))', squeeze(X(:, t, :))');
    sptcorr = corr(squeeze(timep(:, t, :))', squeeze(timep(:, t, :))');
    nanacorr = nanacorr + (~isnan(atcorr));
    nanccorr = nanccorr + (~isnan(cstcorr));
    nanucorr = nanucorr + (~isnan(ustcorr));
    nantcorr = nantcorr + (~isnan(ttcorr));
    nanscorr = nanscorr + (~isnan(sptcorr));
    atcorr(isnan(atcorr)) = 0;
    cstcorr(isnan(cstcorr)) = 0;
    ustcorr(isnan(ustcorr)) = 0;
    ttcorr(isnan(ttcorr)) = 0;
    sptcorr(isnan(sptcorr)) = 0;
    
    acorr = acorr + atcorr;
    cscorr = cscorr + cstcorr;
    uscorr = uscorr + ustcorr;
    tcorr = tcorr + ttcorr;
    spcorr = spcorr + sptcorr;
    
end

acorr = acorr./nanacorr;
cscorr = cscorr./nanccorr;
uscorr = uscorr./nanucorr;
tcorr = tcorr./nantcorr;
spcorr = spcorr./nanscorr;


% --------------------
% Response statistics

% Area

cex = [6:10:trials];
uex = [1:5:trials];
X = caltr;
cX = caltr(:, ~ismember([1:trials], cex), :);
uX = caltr(:, ~ismember([1:trials], uex), :);

careap = zeros(max(max(cellmsk)), 2);
uareap = careap;
aareap = careap;

for c = 1 : max(max(cellmsk))
    careap(c, 1) = squeeze(nanmean(mean(cX(c, :, cstart:cstart+sdur), 3), 2) - mean(mean(X(c, :, astart:astart+sdur), 3), 2));
    uareap(c, 1) = squeeze(nanmean(mean(uX(c, :, ustart:ustart+sdur), 3), 2) - mean(mean(X(c, :, astart:astart+sdur), 3), 2));
    aareap(c, 1) = squeeze(nanmean(mean(X(c, :, astart:astart+sdur), 3), 2) - mean(mean(X(c, :, 1+sdur+(astart:astart+sdur)), 3), 2));
    
    [h careap(c, 2)] = kstest2(squeeze(mean(cX(c, :, cstart:cstart+sdur), 3)), squeeze(mean(X(c, :, astart:astart+sdur), 3)));
    [h uareap(c, 2)] = kstest2(squeeze(mean(uX(c, :, ustart:ustart+sdur), 3)), squeeze(mean(X(c, :, astart:astart+sdur), 3)));
    [h aareap(c, 2)] = kstest2(squeeze(mean(X(c, :, astart:astart+sdur), 3)), squeeze(mean(X(c, :, 1+sdur+(astart:astart+sdur)), 3)));
end

clist = careap(:,2)<0.05;
ulist = uareap(:,2)<0.05;
alist = aareap(:,2)<0.05;


% Averaged p value for multiple instances of the test

cex = [6:10:trials];
uex = [[1:10] [11:5:trials]];
X = caltr;
cX = caltr(:, ~ismember([1:trials], cex), :);
uX = caltr(:, ~ismember([1:trials], uex), :);

careap = zeros(max(max(cellmsk)), length(frstart:sdur:(cstart-3*sdur-1)), 2);
uareap = careap;
aareap = careap;

for c = 1 : max(max(cellmsk))
    count = 1;
    for fr = frstart:sdur:(cstart-3*sdur-1)
        careap(c, count, 1) = squeeze(nanmean(mean(abs(cX(c, :, cstart:cstart+sdur)), 3), 2) - mean(mean(X(c, :, fr:fr+sdur), 3), 2));
        uareap(c, count, 1) = squeeze(nanmean(mean(abs(uX(c, :, ustart:ustart+sdur)), 3), 2) - mean(mean(X(c, :, fr:fr+sdur), 3), 2));
        aareap(c, count, 1) = squeeze(nanmean(mean(abs(X(c, :, astart:astart+sdur)), 3), 2) - mean(mean(X(c, :, fr:fr+sdur), 3), 2));

        [h careap(c, count, 2)] = kstest2(squeeze(mean(cX(c, :, cstart:cstart+sdur), 3)), squeeze(mean(abs(X(c, :, fr:fr+sdur)), 3)));
        [h uareap(c, count, 2)] = kstest2(squeeze(mean(uX(c, :, ustart:ustart+sdur), 3)), squeeze(mean(abs(X(c, :, fr:fr+sdur)), 3)));
        [h aareap(c, count, 2)] = kstest2(squeeze(mean(X(c, :, astart:astart+sdur), 3)), squeeze(mean(abs(X(c, :, fr:fr+sdur)), 3)));
        
        count = count + 1;
    end
end

clist = median(careap(:,:,2), 2) < 0.05;
ulist = median(uareap(:,:,2), 2) < 0.05;
alist = median(aareap(:,:,2), 2) < 0.05;

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

aval = zeros(1, 1000);

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

beh = zeros(trials, trlength*10000);

for tr = 1 : trials
    temp = dlmread([linker num2str(tr-1) '-1.xls']);
    beh (tr, :) = temp (:, 1)';
end


% Clustering

X = reshape(permute(caltr, [1 3 2]), [size(caltr, 1) numel(squeeze(caltr(1, :, :)))]);

IDX = [];
it = 100;  % iterations for k-means

for i = 1 : it
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

cellmskc = zeros(size(cellmsk, 1), size(cellmsk, 2), k);

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
    [I ce] = kmeans_plusplus(X, k);
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

X = reshape(permute(caltr, [1 3 2]), [size(caltr, 1) numel(squeeze(caltr(1, :, :)))]);

[IDX ce] = kmeans_plusplus(X, k);
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

X = reshape(permute(caltr, [1 3 2]), [size(caltr, 1) numel(squeeze(caltr(1, :, :)))]);
ccorr = tcorr;
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
    
x = clustat(:, 1);
x = (x - min(x))/(max(x)-min(x));
y = clustat(:, 2);
y = (y - min(y))/(max(y)-min(y));
x.*y

cellmskc = cellmsk*0;
for i = 1 : size(a, 1)
    temp = a{i, 1};
    for x = 1 : length(temp)
        cellmskc(cellmsk == (temp(x))) = i;
    end
end
