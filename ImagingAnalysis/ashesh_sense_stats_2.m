% Author: Ashesh Dhawale
% Date: 6/11/11

cp = [];
up = [];
ccorrp = [];
ucorrp = [];
ov = [];
prop = [];
reliab = cell(1,3);
reliaban = cell(1,3);
combpresall = cell(1,3);
combrall = cell(2,3);
bcopt = zeros(3,2);
buopt = zeros(3,2);
mnr = cell(1,3);


for i = 1 : 3
    % Load data
    load([num2str(i) '_careap.mat']);
    load([num2str(i) '_aareap.mat']);
    load([num2str(i) '_acorr.mat']);
    load([num2str(i) '_tcorr.mat']);
    load([num2str(i) '_clusters.mat']);
    load([num2str(i) '_cellmsk.mat']);
    if i < 3
        load([num2str(i) '_uareap.mat']);
    end
        
    IDX = zeros(size(tcorr,1), 1);      % another cluster representation
    for x = 1 : size(a,1)
        IDX(a{x, 1}) = x;
    end
    rIDX = IDX(randperm(length(IDX)));  % shuffled clusters
    
    
    % Expected distribution
    
    ecl = [];
    for p = 1 : size(a, 1)
        ecl(p) = sum(IDX == p)/length(IDX);
    end
    
    
    % Real distribution
    if i < 3
        clist = median(careap(:,:,2), 2) < 0.05;
        cres = clist.*squeeze(mean(careap(:, :, 1), 2));
        ulist = median(uareap(:,:,2), 2) < 0.05;
        ures = ulist.*squeeze(mean(uareap(:, :, 1), 2));
    else
        clist1 = squeeze((median(careap(1, 2, :, :, 2),4)<0.05));
        clist2 = squeeze((median(careap(2, 2, :, :, 2),4)<0.05));
        clist = or(clist1, clist2);
        cres1 = clist1.*squeeze(mean(careap(1, 2, :, :, 1), 4));
        cres2 = clist2.*squeeze(mean(careap(2, 2, :, :, 1), 4));
        cres1(cres1 == 0) = NaN;
        cres2(cres2 == 0) = NaN;
        cres = nanmean([cres1 cres2], 2);
        cres(isnan(cres)) = 0;
        
        ulist1 = squeeze((median(careap(1, 3, :, :, 2),4)<0.05));
        ulist2 = squeeze((median(careap(2, 3, :, :, 2),4)<0.05));
        ulist = or(ulist1, ulist2);
        ures1 = ulist1.*squeeze(mean(careap(1, 3, :, :, 1), 4));
        ures2 = ulist2.*squeeze(mean(careap(2, 3, :, :, 1), 4));
        ures1(ures1 == 0) = NaN;
        ures2(ures2 == 0) = NaN;
        ures = nanmean([ures1 ures2], 2);
        ures(isnan(ures)) = 0;
    end
    
    ccl = [];
    ucl = [];
    for p = 1 : size(a, 1)
       ccl(p) = sum(clist & (IDX == p));
       ucl(p) = sum(ulist & (IDX == p));
    end
    
    ccl = ccl/sum(clist);
    ucl = ucl/sum(ulist);
    
    cdist = sqrt(sum((ccl - ecl).^2));
    udist = sqrt(sum((ucl - ecl).^2));
    
    % Shuffled distribution
    
    rand('seed', 0);
    rcdist = zeros(10000, 1);
    rudist = zeros(10000, 1);
    for t = 1 : 10000
        rclist = clist(randperm(length(clist)));
        rulist = ulist(randperm(length(ulist)));
        rcl = [];
        rul = [];
        for p = 1 : size(a, 1)
            rcl(p) = sum(rclist & (IDX==p));
            rul(p) = sum(rulist & (IDX==p));
        end
        rcl = rcl/sum(rclist);
        rul = rul/sum(rulist);
        
        rcdist(t) = sqrt(sum((rcl - ecl).^2));
        rudist(t) = sqrt(sum((rul - ecl).^2));
    end

    % p value
    cp(i) = sum(rcdist > cdist)/length(rcdist);
    up(i) = sum(rudist > udist)/length(rudist);
    
    
    
    % Correlations of tone and puff sensitive neurons
    X = acorr;
    ccorr = [];
    ucorr = [];
    
    ctemp = X(clist, clist);
    for x = 1 : sum(clist)-1
        for y = x+1 : sum(clist)
            ccorr = [ccorr ctemp(x, y)];
        end
    end
    
    utemp = X(ulist, ulist);
    for x = 1 : sum(ulist)-1
        for y = x+1 : sum(ulist)
            ucorr = [ucorr utemp(x, y)];
        end
    end
    
    % Correlations between other neurons
    nccorr = [];
    nucorr = [];
    
    ctemp = X(~clist, ~clist);
    for x = 1 : sum(~clist)-1
        for y = x+1 : sum(~clist)
            nccorr = [nccorr ctemp(x, y)];
        end
    end
    
    utemp = X(~ulist, ~ulist);
    for x = 1 : sum(~ulist)-1
        for y = x+1 : sum(~ulist)
            nucorr = [nucorr utemp(x, y)];
        end
    end
    
    
    % Correlations between all neurons
    allcorr = [];
    
    for x = 1 : size(X, 1)-1
        for y = x+1 : size(X, 1)
            allcorr = [allcorr X(x, y)];
        end
    end
    
    [h ccorrp(i, 1)] = kstest2(ccorr, allcorr, 0.05);
    [h ucorrp(i, 1)] = kstest2(ucorr, allcorr, 0.05);
    [h ccorrp(i, 2)] = kstest2(ccorr, nccorr, 0.05);
    [h ucorrp(i, 2)] = kstest2(ucorr, nucorr, 0.05);
    
    
    % Both responsive cells proportion
    ov(i,1) = sum(clist & ulist)/length(clist);
    ov(i,2) = sum(clist)*sum(ulist)/(length(clist)^2);
    
    rov = [];
    rand('seed', 0);
    for x = 1 : 1000
        rctemp = clist(randperm(length(clist)));
        rutemp = ulist(randperm(length(ulist)));
        rov = [rov sum(rctemp & rutemp)/length(clist)];
    end
    
    ov(i,3) = sum(rov > ov(i,1))/length(rov);   % significance of overlap
    
    prop(i,1) = sum(clist)/length(clist);   % proportion tone responsive
    prop(i,2) = sum(cres<0)/sum(clist);     % proportion inhibitory
    prop(i,3) = sum(ulist)/length(ulist);   % proportion puff responsive
    prop(i,4) = sum(ures<0)/sum(ulist);     % proportion inhibitory
  
  
    % Reliability of responses
    load([num2str(i) '_timep.mat']);
    load([num2str(i) '_caltr.mat']);
    trials = size(timep, 2);
    ncell = max(max(cellmsk));
    sdur = 2;
    
    switch i
        case 1
            samp = 15/size(timep, 3);
            cex = 6:10:trials;
            uex = 1:5:trials;
            cX = timep(:, ~ismember(1:trials, cex), :);
            uX = timep(:, ~ismember(1:trials, uex), :);
            cC = caltr(:, ~ismember(1:trials, cex), :);
            uC = caltr(:, ~ismember(1:trials, uex), :);
            frstart = 4;
            cstart = 78;
            ustart = 86;
            
        case 2
            samp = 6/size(timep, 3);
            cex = 6:10:trials;
            uex = 1:5:trials;
            cX = timep(:, ~ismember(1:trials, cex), :);
            uX = timep(:, ~ismember(1:trials, uex), :);
            cC = caltr(:, ~ismember(1:trials, cex), :);
            uC = caltr(:, ~ismember(1:trials, uex), :);
            frstart = 4;
            cstart = 32;
            ustart = 39;
            
        case 3
            samp = 320/size(timep, 3);
            stime = 1; % stimulus duration in sec
            siti = 5; % inter stimulus interval in sec
            sstart = [20 90 160 230]; % odour tone puff light
            wind = 3; % half window size in seconds
            windfr = ceil(wind/samp);
            calsens = zeros(4, 10, ncell, trials, 2*windfr+1);
            for stim = 1 : 4
                for p = 1 : 10
                    temp = round((((p-1) * (siti+stime)) + sstart(stim))/samp);
                    calsens(stim, p, :, :, :) = timep(:, :, temp-windfr:temp+windfr);
                end
            end
            cX = reshape(permute(calsens(2, :, :, :, :), [1 3 4 2 5]), ncell, trials*10, 2*windfr+1);
            uX = reshape(permute(calsens(3, :, :, :, :), [1 3 4 2 5]), ncell, trials*10, 2*windfr+1);
            
            calsensC = zeros(4, 10, ncell, trials, 2*windfr+1);
            for stim = 1 : 4
                for p = 1 : 10
                    temp = round((((p-1) * (siti+stime)) + sstart(stim))/samp);
                    calsensC(stim, p, :, :, :) = caltr(:, :, temp-windfr:temp+windfr);
                end
            end
            cC = reshape(permute(calsensC(2, :, :, :, :), [1 3 4 2 5]), ncell, trials*10, 2*windfr+1);
            uC = reshape(permute(calsensC(3, :, :, :, :), [1 3 4 2 5]), ncell, trials*10, 2*windfr+1);
            
            cstart = windfr+3; ustart = windfr+3; frstart = 1;
            trials = size(cX, 2);
    end
    
    ctr = size(cX, 2);
    utr = size(uX, 2);
    
    presa = squeeze(sum(sum(cX(:, :, frstart:cstart-1), 3), 2)/(trials*length(frstart:cstart-1)));
    presad = 1 - ((1-presa).^(sdur+1));
    
    presaa = zeros(ncell, 1);
    for air = frstart : sdur+1 : cstart-sdur-1
        presaa = presaa + (squeeze(sum(sum(cX(:, :, air:air+sdur), 3) > 0, 2)))/ctr;
    end 
    presaa = presaa/length(frstart:sdur+1:cstart-sdur-1);

    presc = (squeeze(sum(sum(cX(:, :, cstart:cstart+sdur), 3) > 0, 2)))/(ctr);
    presu = (squeeze(sum(sum(uX(:, :, ustart:ustart+sdur), 3) > 0, 2)))/(utr);
    
    presaac = presaa; presaac(presc < presaa) = 1 - presaa(presc < presaa); % To account for inhibitory responses
    presc(presc < presaa) = 1 - presc(presc < presaa);
    
    presaau = presaa; presaau(presu < presaa) = 1 - presaa(presu < presaa);
    presu(presu < presaa) = 1 - presu(presu < presaa);
    
    reliab{i} = [presc presaac presu presaau];
    
    
    % ROC curve based on analog values
    
    cset = sum(cC(:, :, cstart:cstart+sdur), 3);
    uset = sum(uC(:, :, ustart:ustart+sdur), 3);
    
    aset = zeros(ncell, ctr, length(frstart:sdur+1:cstart-sdur-1));
    count = 1;
    for air = frstart : sdur+1 : cstart-sdur-1
        aset(:, :, count) = sum(cC(:, :, air:air+sdur), 3);
        count = count + 1;
    end
    aset = reshape(aset, ncell, ctr*length(frstart:sdur+1:cstart-sdur-1));
    
    presanc = zeros(ncell, 2);
    presanu = zeros(ncell, 2);
    cthresh = zeros(ncell, 2);
    uthresh = zeros(ncell, 2);
    
    for n = 1 : ncell
        atemp = sort(aset(n, :), 'descend');
        ctemp = sort(cset(n, :), 'descend');
        utemp = sort(uset(n, :), 'descend');
       
        croc = zeros(length(atemp), 2);
        uroc = zeros(length(atemp), 2);
        for x = 1 : length(atemp)
            croc(x, 1) = sum(atemp >= atemp(x))/length(atemp);
            croc(x, 2) = sum(ctemp >= atemp(x))/length(ctemp);
            uroc(x, 1) = sum(atemp >= atemp(x))/length(atemp);
            uroc(x, 2) = sum(utemp >= atemp(x))/length(utemp);
        end
        [temp1 ind1] = min(sqrt(croc(:,1).^2 + (croc(:,2)-1).^2));
        [temp2 ind2] = min(sqrt(croc(:,2).^2 + (croc(:,1)-1).^2));
        if temp1 < temp2
           presanc(n, 1) = croc(ind1, 1);
           presanc(n, 2) = croc(ind1, 2);
           cthresh(n, 1) = atemp(ind1);
           cthresh(n, 2) = 1;
        else
           presanc(n, 1) = croc(ind2, 2);
           presanc(n, 2) = croc(ind2, 1);
           cthresh(n, 1) = atemp(ind2);
           cthresh(n, 2) = 0;
        end
        
        
        [temp1 ind1] = min(sqrt(uroc(:,1).^2 + (uroc(:,2)-1).^2));
        [temp2 ind2] = min(sqrt(uroc(:,2).^2 + (uroc(:,1)-1).^2));
        if temp1 < temp2
           presanu(n, 1) = uroc(ind1, 1);
           presanu(n, 2) = uroc(ind1, 2);
           uthresh(n, 1) = atemp(ind1);
           uthresh(n, 2) = 1;
        else
           presanu(n, 1) = uroc(ind2, 2);
           presanu(n, 2) = uroc(ind2, 1);
           uthresh(n, 1) = atemp(ind2);
           uthresh(n, 2) = 0;
        end
    end
    
    reliaban{i} = [presanc presanu];
    
    
        
    % Cell combinations
    
%     maxnum = 1;
%     combpres = cell(1, maxnum);
%     for num = 1 : maxnum
%         num
%         C = nchoosek(1:ncell, num);
%         tempresc = zeros(num, size(C, 1));
%         tempresu = tempresc;
%         tempresa = tempresc;
%         for comb = 1 : size(C, 1)
%             for thresh = 0 : num-1
%                 tempresc(thresh+1,comb) = (squeeze(sum(sum(sum(cX(C(comb, :), :, cstart:cstart+sdur), 3) > 0, 1) > thresh, 2)))/trials;
%                 tempresu(thresh+1,comb) = (squeeze(sum(sum(sum(uX(C(comb, :), :, ustart:ustart+sdur), 3) > 0, 1) > thresh, 2)))/size(uX, 2);
%                 temp = 0;
%                 for air = frstart : sdur+1 : cstart-sdur-1
%                     temp = temp + (squeeze(sum(sum(sum(cX(C(comb, :), :, air:air+sdur), 3) > 0, 1) > thresh, 2)))/trials;
%                 end
%                 tempresa(thresh+1,comb) = temp/length(frstart:sdur+1:cstart-sdur-1);
%             end
%         end
%         combpres{1, num} = [tempresa; tempresc; tempresu];
%     end
%     
%     combpresall{i} = combpres;
%     
%     
    % Cell combinations - top 10% of subset of all combinations
    % Digital as well as analog representations
    
    tries = 1000;
    combr = cell(ncell, 1);
    combran = cell(ncell, 1);
    
    aX = zeros(ncell, ctr, length(frstart:sdur+1:cstart-sdur-1));
    count = 1;
    for air = frstart : sdur+1 : cstart-sdur-1
        aX(:, :, count) = sum(cX(:, :, air:air+sdur), 3) > 0;
        count = count + 1;
    end
    aX = reshape(aX, ncell, ctr*length(frstart:sdur+1:cstart-sdur-1));
    
    cXan = zeros(size(cset)); uXan = zeros(size(uset)); 
    acXan = zeros(size(aset)); auXan = acXan;
    for n = 1 : ncell
        if cthresh(n, 2)
            cXan(n, :) = cset(n, :) > cthresh(n, 1);
            acXan(n, :) = aset(n, :) > cthresh(n, 1);
        else
            cXan(n, :) = cset(n, :) < cthresh(n, 1);
            acXan(n, :) = aset(n, :) < cthresh(n, 1);
        end
        
        if uthresh(n, 2)
            uXan(n, :) = uset(n, :) > uthresh(n, 1);
            auXan(n, :) = aset(n, :) > uthresh(n, 1);
        else
            uXan(n, :) = uset(n, :) < uthresh(n, 1);
            auXan(n, :) = aset(n, :) < uthresh(n, 1);
        end
    end
    
    rand('seed', 0);
%     
%     for num = 1 : ncell
%         num
%         if nchoosek(ncell, num) < 10*tries
%             combs = nchoosek(1:ncell, num);
%             if nchoosek(ncell, num) > tries
%                 temp = randperm(nchoosek(ncell, num));
%                 combs = combs(temp(1:tries), :);
%             end
%         else
%             combs = zeros(tries, num);
%             for t = 1 : tries
%                 temp = randperm(ncell);
%                 combs(t, :) = temp(1:num);
%             end
%         end
        
%         tempresc = zeros(num, size(combs, 1));
%         tempresu = tempresc;
%         tempresa = tempresc;
%         
%         prescan = zeros(num, size(combs, 1));
%         presuan = prescan;
%         presacan = prescan;
%         presauan = prescan;
%         
%         for t = 1 : size(combs, 1)
%             for thresh = 0 : num-1
%                 % Digital
%                 tempresc(thresh+1, t) = (squeeze(sum(sum(sum(cX(combs(t, :), :, cstart:cstart+sdur), 3) > 0, 1) > thresh, 2)))/ctr;
%                 tempresu(thresh+1, t) = (squeeze(sum(sum(sum(uX(combs(t, :), :, ustart:ustart+sdur), 3) > 0, 1) > thresh, 2)))/utr;
%                 tempresa(thresh+1, t) = squeeze(sum(sum(aX(combs(t, :), :), 1) > thresh, 2))/size(aX, 2);
%                 
%                 % Analog
%                 prescan(thresh+1, t) = squeeze(sum(sum(cXan(combs(t, :), :), 1) > thresh, 2))/ctr;
%                 presuan(thresh+1, t) = squeeze(sum(sum(uXan(combs(t, :), :), 1) > thresh, 2))/utr;
%                 presacan(thresh+1, t) = squeeze(sum(sum(acXan(combs(t, :), :), 1) > thresh, 2))/size(acXan, 2);
%                 presauan(thresh+1, t) = squeeze(sum(sum(auXan(combs(t, :), :), 1) > thresh, 2))/size(auXan, 2);
%             end
%         end
%         
%         combr{num} = [tempresa; tempresc; tempresu];
%         combran{num} = [prescan; presacan; presuan; presauan];
%         
%     end
%     
%     combrall{1, i} = combr;
%     combrall{2, i} = combran;
    
    
    % Burst size
    
    burstc = sum(sum(cX(:, :, cstart:cstart+sdur), 3), 1);
    burstu = sum(sum(uX(:, :, ustart:ustart+sdur), 3), 1);
    bursta = [];
    for air = frstart : sdur+1 : cstart-sdur-1
        bursta = [bursta sum(sum(cX(:, :, air:air+sdur), 3), 1)];
    end
    bursta = sort(bursta, 'descend');
    
    bcroc = zeros(length(bursta), 2);
    buroc = zeros(length(bursta), 2);
    for x = 1 : length(bursta)
        bcroc(x, 1) = sum(bursta >= bursta(x))/length(bursta);
        bcroc(x, 2) = sum(burstc >= bursta(x))/length(burstc);
        buroc(x, 1) = sum(bursta >= bursta(x))/length(bursta);
        buroc(x, 2) = sum(burstu >= bursta(x))/length(burstu);
    end
    
    [temp1 ind1] = min(sqrt(bcroc(:,1).^2 + (bcroc(:,2)-1).^2));
    [temp2 ind2] = min(sqrt(bcroc(:,2).^2 + (bcroc(:,1)-1).^2));
    if temp1 < temp2
       bcopt(i, :) = [bcroc(ind1, 1) bcroc(ind1, 2)];
    else
       bcopt(i, :) = [bcroc(ind2, 2) bcroc(ind2, 1)];
    end
    
    [temp1 ind1] = min(sqrt(buroc(:,1).^2 + (buroc(:,2)-1).^2));
    [temp2 ind2] = min(sqrt(buroc(:,2).^2 + (buroc(:,1)-1).^2));
    if temp1 < temp2
       buopt(i, :) = [buroc(ind1, 1) buroc(ind1, 2)];
    else
       buopt(i, :) = [buroc(ind2, 2) buroc(ind2, 1)];
    end
    
    
    % Cluster based burst responses
    
    clburstc = zeros(size(a, 1), size(cX, 2));
    clburstu = zeros(size(a, 1), size(uX, 2)); 
    clbursta = zeros(size(a, 1), size(cX, 2)*length(frstart : sdur+1 : cstart-sdur-1));
    
    for x = 1 : size(a,1)
        clburstc(x, :) = sum(sum(cX(a{x, 1}, :, cstart:cstart+sdur), 3), 1);
        clburstu(x, :) = sum(sum(uX(a{x, 1}, :, ustart:ustart+sdur), 3), 1);
        temp = [];
        for air = frstart : sdur+1 : cstart-sdur-1
            temp = [temp sum(sum(cX(a{x, 1}, :, air:air+sdur), 3), 1)];
        end
        clbursta(x, :) = temp;
        
        [rocurve tpr fpr] = doroc(clburstc(x, :), clbursta(x, :));
%         figure(1); 
%         plot(rocurve(:,1), rocurve(:,2), 'b-');
%         hold on;
%         plot(fpr, tpr, 'bo');
        
%         figure(2); hold on
%         plot(sqrt(bcopt(i,1).^2 + (bcopt(i,2)-1).^2), sqrt(fpr.^2 + (tpr-1).^2), 'bo')
        
        [rocurve tpr fpr] = doroc(clburstu(x, :), clbursta(x, :));
        
%         figure(1);
%         plot(rocurve(:,1), rocurve(:,2), 'r-');
%         plot(fpr, tpr, 'ro');
%         
%         figure(2); 
%         plot(sqrt(buopt(i,1).^2 + (buopt(i,2)-1).^2), sqrt(fpr.^2 + (tpr-1).^2), 'ro')
        
    end
    
    
    % Multinominal logistic regression
    
    cR = (sum(cX(:, :, cstart:cstart+sdur), 3) > 0)';
    uR = (sum(uX(:, :, ustart:ustart+sdur), 3) > 0)';
    
    nair = round(ncell*10/size(cX,2)) - 1;
    
    aR = [];
    for air = frstart : sdur+1 : min([(frstart + nair*(sdur+1)) (cstart-sdur-1)])
        aR = [aR (sum(cX(:, :, air:air+sdur), 3) > 0)];
    end
    aR = aR';
    
    cR = [cR; aR];
    uR = [uR; aR];
    
    
    aset = zeros(ncell, ctr, length(frstart : sdur+1 : min([(frstart + nair*(sdur+1)) (cstart-sdur-1)])));
    count = 1;
    for air = frstart : sdur+1 : min([(frstart + nair*(sdur+1)) (cstart-sdur-1)])
        aset(:, :, count) = sum(cC(:, :, air:air+sdur), 3);
        count = count + 1;
    end
    aset = reshape(aset, ncell, numel(aset(1,:,:)));
    
    acXan = zeros(size(aset)); auXan = acXan;
    for n = 1 : ncell
        if cthresh(n, 2)
            acXan(n, :) = aset(n, :) > cthresh(n, 1);
        else
            acXan(n, :) = aset(n, :) < cthresh(n, 1);
        end
        
        if uthresh(n, 2)
            auXan(n, :) = aset(n, :) > uthresh(n, 1);
        else
            auXan(n, :) = aset(n, :) < uthresh(n, 1);
        end
    end
    
    cRan = [cXan'; acXan'];
    uRan = [uXan'; auXan'];
    
    cY = zeros(size(cR,1), 2);
    cY(1:size(cX,2), 1) = 1; cY(size(cX,2)+1:end, 2) = 1;
    uY = zeros(size(uR,1), 2);
    uY(1:size(uX,2), 1) = 1; uY(size(uX,2)+1:end, 2) = 1;
    
    cB = mnrfit(cR, cY); cP = mnrval(cB, cR);
    uB = mnrfit(uR, uY); uP = mnrval(uB, uR);
    
    % Trial shuffled dataset
    
    rand('seed', 0);
    
    rcset = zeros(size(cset));
    ruset = zeros(size(uset));
    raset = zeros(size(aset));
    
    for n = 1 : ncell
        rcset(n,:) = cset(n, randperm(size(cset,2)));
        ruset(n,:) = uset(n, randperm(size(uset,2)));
        raset(n,:) = aset(n, randperm(size(aset,2)));
    end
    
    
    % MNR decoding performance versus cell number
    
    % Ordered by coefficient weight in complete dataset
    
%     [crap cI] = sort(abs(cB(1:end-1)), 'descend');
%     [crap uI] = sort(abs(uB(1:end-1)), 'descend');
%     
%     cPc = zeros(ncell, 2);
%     uPc = cPc;
%     
%     for x = 1 : ncell
%         ctemp = mnrval(mnrfit(cR(:,cI(1:x)), cY), cR(:,cI(1:x)));
%         utemp = mnrval(mnrfit(uR(:,uI(1:x)), uY), uR(:,uI(1:x)));
%         ctemp = ctemp > 0.5;
%         utemp = utemp > 0.5;
%         cPc(x, :) = [(sum(ctemp(1:size(cX, 2),1))/size(cX, 2)) (sum(ctemp(size(cX, 2)+1:end,1))/size(aR, 1))];
%         uPc(x, :) = [(sum(utemp(1:size(uX, 2),1))/size(uX, 2)) (sum(utemp(size(uX, 2)+1:end,1))/size(aR, 1))];
%     end
%     
    
    
    % Top 10 percentile of multiple random combinations
    
    tries = 100;
    cPc = zeros(ncell, tries, 2);
    uPc = cPc;
    cPcan = cPc;
    uPcan = cPc;
    Pall = cell(ncell, 1);
    
    rand('seed', 0);
    cYr = cY(randperm(size(cY, 1)), :);     % Shuffled stimulus vector
    rPc = cPc;
    rPcan = cPc;
    
    cRan = [cset aset]';
    uRan = [uset aset]';
    
    cR = [rcset raset]';
    uR = [ruset raset]';
    
    warning off all;
%     
%     rand('seed', 0);
%     for x = 1 : ncell
%         x
%         P = zeros(6, tries, x+1);
%         
%         if nchoosek(ncell, x) < 10*tries
%             combs = nchoosek(1:ncell, x);
%             if nchoosek(ncell, x) > tries
%                 temp = randperm(nchoosek(ncell, x));
%                 combs = combs(temp(1:tries), :);
%             end
%         else
%             combs = zeros(tries, x);
%             for t = 1 : tries
%                 temp = randperm(ncell);
%                 combs(t, :) = temp(1:x);
%             end
%         end
%         
%         for t = 1 : size(combs, 1)
%             P(1, t, :) = mnrfit(cR(:,combs(t, :)), cY);
%             P(2, t, :) = mnrfit(uR(:,combs(t, :)), uY);
%             P(3, t, :) = mnrfit(cR(:,combs(t, :)), cYr);
%             ctemp = mnrval(squeeze(P(1, t, :)), cR(:,combs(t, :)));
%             utemp = mnrval(squeeze(P(2, t, :)), uR(:,combs(t, :)));
%             rtemp = mnrval(squeeze(P(3, t, :)), cR(:,combs(t, :)));
%             ctemp = ctemp > 0.5;
%             utemp = utemp > 0.5;
%             rtemp = rtemp > 0.5;
%             cPc(x, t, :) = [(sum(ctemp(1:size(cX, 2),1))/size(cX, 2)) (sum(ctemp(size(cX, 2)+1:end,1))/size(aR, 1))];
%             uPc(x, t, :) = [(sum(utemp(1:size(uX, 2),1))/size(uX, 2)) (sum(utemp(size(uX, 2)+1:end,1))/size(aR, 1))];
%             rPc(x, t, :) = [(sum(rtemp(:,1)&cYr(:,1))/size(cX, 2)) (sum(xor(rtemp(:,1), cYr(:,1))&rtemp(:,1))/size(aR, 1))];
%             
%             P(4, t, :) = mnrfit(cRan(:,combs(t, :)), cY);
%             P(5, t, :) = mnrfit(uRan(:,combs(t, :)), uY);
%             P(6, t, :) = mnrfit(cRan(:,combs(t, :)), cYr);
%             ctemp = mnrvalll(squeeze(P(4, t, :)), cRan(:,combs(t, :)));
%             utemp = mnrval(squeeze(P(5, t, :)), uRan(:,combs(t, :)));
%             rtemp = mnrval(squeeze(P(6, t, :)), cRan(:,combs(t, :)));
%             ctemp = ctemp > 0.5;
%             utemp = utemp > 0.5;
%             rtemp = rtemp > 0.5;
%             cPcan(x, t, :) = [(sum(ctemp(1:size(cX, 2),1))/size(cX, 2)) (sum(ctemp(size(cX, 2)+1:end,1))/size(aR, 1))];
%             uPcan(x, t, :) = [(sum(utemp(1:size(uX, 2),1))/size(uX, 2)) (sum(utemp(size(uX, 2)+1:end,1))/size(aR, 1))];
%             rPcan(x, t, :) = [(sum(rtemp(:,1)&cYr(:,1))/size(cX, 2)) (sum(xor(rtemp(:,1), cYr(:,1))&rtemp(:,1))/size(aR, 1))];
%         end
%         
%         if t < tries
%             cPc(x, t+1:end, :) = NaN;
%             uPc(x, t+1:end, :) = NaN;
%             rPc(x, t+1:end, :) = NaN;
%             cPcan(x, t+1:end, :) = NaN;
%             uPcan(x, t+1:end, :) = NaN;
%             rPcan(x, t+1:end, :) = NaN;
%         end
%         Pall{x} = P;
%     end 
%     
%     warning on all
%     
%     mnrall = [];
%     mnrall(1, :, :, :) = cPc;
%     mnrall(2, :, :, :) = uPc;
%     mnrall(3, :, :, :) = rPc;
%     mnrall(4, :, :, :) = cPcan;
%     mnrall(5, :, :, :) = uPcan;
%     mnrall(6, :, :, :) = rPcan;
%     
%     
%     mnr{i} = mnrall;
%     
%     save([num2str(i) '_mnr.mat'], 'mnrall');
%     save([num2str(i) '_Pall.mat'], 'Pall');

    % TPR and FPR with MLR of shuffled data as a function of training sample size
    
    iter = 100;
    rPsc = zeros(length(size(cX,2)*2 : 50 : size(cRan, 1)), iter, 2);
    rand('seed', 0);
    
    for t = 1:iter
        t
        cYr = cY(randperm(size(cY,1)), :);
        [crap I] = sort(cYr(:,1), 'descend');
        cYr = cYr(I,:); 
        cRI = cRan(I,:);
        
        count = 1;
        for ssize = size(cX,2)*2 : 50 : size(cR, 1) 
            rtemp = mnrval(mnrfit(cRI(1:ssize,:), cYr(1:ssize,:)), cRI(1:ssize,:));
            rtemp = rtemp > 0.5;
            rPsc(count,t,:) = [(sum(rtemp(:,1)&cYr(1:ssize,1))/size(cX, 2)) (sum(xor(rtemp(:,1), cYr(1:ssize,1))&rtemp(:,1))/(ssize-size(cX,2)))];
            count = count + 1;
        end
    end
    
end
% save('mnr.mat', 'mnr');
return
for i = 1
        
    
    % Are clusters better at decoding? Within cluster MLR decoding
    
    cld = cell(size(a,1), 4);
    iter = 100;
    rand('seed', 0);
    
    for cl = 1 : size(a,1)
        alist = a{cl, 1};
        ctemp = mnrval(mnrfit(cR(:,alist), cY), cR(:,alist));
        utemp = mnrval(mnrfit(uR(:,alist), uY), uR(:,alist));
        ctemp = ctemp > 0.5;
        utemp = utemp > 0.5;
        cldc = [(sum(ctemp(1:size(cX, 2),1))/size(cX, 2)) (sum(ctemp(size(cX, 2)+1:end,1))/size(aR, 1))];
        cldu = [(sum(utemp(1:size(uX, 2),1))/size(uX, 2)) (sum(utemp(size(uX, 2)+1:end,1))/size(aR, 1))];
        
        rcldc = zeros(iter, 2);
        rcldu = zeros(iter, 2);
        for t = 1 : iter
            t
            temp = randperm(ncell);
            rlist = temp(1:length(alist));
            
            rctemp = mnrval(mnrfit(cR(:,rlist), cY), cR(:,rlist));
            rctemp = rctemp > 0.5;
            rutemp = mnrval(mnrfit(uR(:,rlist), uY), uR(:,rlist));
            rutemp = rutemp > 0.5;
            
            rcldc(t, :) = [(sum(rctemp(1:size(cX, 2),1))/size(cX, 2)) (sum(rctemp(size(cX, 2)+1:end,1))/size(aR, 1))];
            rcldu(t, :) = [(sum(rutemp(1:size(uX, 2),1))/size(uX, 2)) (sum(rutemp(size(uX, 2)+1:end,1))/size(aR, 1))]; 
        end
        
        cld{cl, 1} = cldc;
        cld{cl, 2} = cldu;
        cld{cl, 3} = rcldc;
        cld{cl, 4} = rcldu;
    end
    
    
    
    % Are clusters better at decoding? Pooled cluster MLR decoding
    
    IDX = zeros(ncell, 1);      % another cluster representation
    
    PcR = zeros(size(cR, 1), size(a,1));
    PuR = zeros(size(uR, 1), size(a,1));
    iter = 100;
    
    for x = 1 : size(a,1)
        IDX(a{x, 1}) = x;
        PcR(:, x) = sum(cR(:, a{x,1}), 2);
        PuR(:, x) = sum(uR(:, a{x,1}), 2);
    end
    
    ctemp = mnrval(mnrfit(PcR, cY), PcR);
    utemp = mnrval(mnrfit(PuR, uY), PuR);
    ctemp = ctemp > 0.5;
    utemp = utemp > 0.5;
    PcP = [(sum(ctemp(1:size(cX, 2),1))/size(cX, 2)) (sum(ctemp(size(cX, 2)+1:end,1))/size(aR, 1))];
    PuP = [(sum(utemp(1:size(uX, 2),1))/size(uX, 2)) (sum(utemp(size(uX, 2)+1:end,1))/size(aR, 1))];
    
    rand('seed', 0);
    rPcP = zeros(iter, 2);
    rPuP = zeros(iter, 2);
    
    for t = 1 : iter
        t
        rIDX = IDX(randperm(length(IDX)));  % shuffled clusters
        rPcR = zeros(size(cR, 1), size(a,1));
        rPuR = zeros(size(uR, 1), size(a,1));
        for x = 1 : size(a,1)
            rPcR(:, x) = sum(cR(:, rIDX==x), 2);
            rPuR(:, x) = sum(uR(:, rIDX==x), 2);    
        end
        
        ctemp = mnrval(mnrfit(rPcR, cY), rPcR);
        utemp = mnrval(mnrfit(rPuR, uY), rPuR);
        ctemp = ctemp > 0.5;
        utemp = utemp > 0.5;
        rPcP(t, :) = [(sum(ctemp(1:size(cX, 2),1))/size(cX, 2)) (sum(ctemp(size(cX, 2)+1:end,1))/size(aR, 1))];
        rPuP(t, :) = [(sum(utemp(1:size(uX, 2),1))/size(uX, 2)) (sum(utemp(size(uX, 2)+1:end,1))/size(aR, 1))];
    end
    
    
    % Averaged coefficient vector through resampling
%     
%     cBx = [];
%     uBx = [];
%     cvratio = 3;
%     
%     cRa = [(sum(cX(:, :, cstart:cstart+sdur), 3) > 0) (sum(cX(:, :, frstart:frstart+sdur), 3) > 0)]';
%     uRa = [(sum(uX(:, :, ustart:ustart+sdur), 3) > 0) (sum(uX(:, :, frstart:frstart+sdur), 3) > 0)]';
%     cYa = zeros(size(cRa,1), 2);
%     cYa(1:size(cX,2), 1) = 1; cYa(size(cX,2)+1:end, 2) = 1;
%     uYa = zeros(size(uRa,1), 2);
%     uYa(1:size(uX,2), 1) = 1; uYa(size(uX,2)+1:end, 2) = 1;
%     
%     rand('seed', 10000);
%     
%     warning off all
%     
%     cR = cRa(1:2:end, :); cY = cYa(1:2:end, :);
%     uR = uRa(1:2:end, :); uY = uYa(1:2:end, :);
%     
%     for x = 1 : 100
%         x
%         crlist = randperm(size(cR, 1));
%         urlist = randperm(size(uR, 1));
%         
%         cBx = [cBx glmfit(cR(crlist(round(size(cR, 1)/cvratio):end), :), cY(crlist(round(size(cR, 1)/cvratio):end), 1), 'binomial')];
%         uBx = [uBx glmfit(uR(urlist(round(size(uR, 1)/cvratio):end), :), uY(urlist(round(size(uR, 1)/cvratio):end), 1), 'binomial')];
%     end
%     
%     cP = mnrval(median(cBx, 2), cRa(1:2:end,:));
%     uP = mnrval(median(uBx, 2), uRa(1:2:end,:));
%     
%     warning on all
    
        
    % Decoding ability of analog df/f values
    
    
end



% Plotting combination data

tries = 1000;

figure(1); hold on;
figure(2); hold on;
figure(3); hold on;
figure(4); hold on;
figure(5); hold on;

c2dhist = zeros(3, length((0:0.01:1)), length((0:0.01:1)));
u2dhist = c2dhist;

for i = 1 : 3
    
    cscat = [];
    uscat = [];
    
    combr = combrall{1, i};
    combran = combrall{2, i};
    
    ncell = length(combr);
    
    load([num2str(i) '_clusters.mat']);
    
    IDX = zeros(ncell, 1);      % another cluster representation
    for x = 1 : size(a,1)
        IDX(a{x, 1}) = x;
    end
    
    % Digital
    cmean = zeros(ncell, 2, 2);
    ctop10 = zeros(ncell, 2, 2);
    umean = zeros(ncell, 2, 2);
    utop10 = zeros(ncell, 2, 2);
    
    % Analog
    ancmean = zeros(ncell, 2, 2);
    anctop10 = zeros(ncell, 2, 2);
    anumean = zeros(ncell, 2, 2);
    anutop10 = zeros(ncell, 2, 2);
    
    rand('seed', 0);
    for n = 1 : ncell
        warning off all
        % Cluster homogeniety
        if nchoosek(ncell, n) < 10*tries
            combs = nchoosek(1:ncell, n);
            if nchoosek(ncell, n) > tries
                temp = randperm(nchoosek(ncell, n));
                combs = combs(temp(1:tries), :);
            end
        else
            combs = zeros(tries, n);
            for t = 1 : tries
                temp = randperm(ncell);
                combs(t, :) = temp(1:n);
            end
        end
        warning on all
        
        % Cluster homogeneity index
        clhom = zeros(1, size(combs, 1));
        for x = 1 : size(combs, 1)
            clsum = [];
            for c = 1 : size(a, 1)
                clsum = [clsum sum(IDX(combs(x, :))== c)/n];
            end
            clhom(x) = max(clsum);
        end
        
        
        % Combination data
        
        % Digital
        temp1 = combr{n}; 
        tempresa = temp1(1:n, :);
        tempresc = temp1(n+1:2*n, :);
        tempresu = temp1(2*n+1:end, :);
%         tempresa(tempresu == 0) = NaN; tempresc(tempresu == 0) = NaN; tempresu(tempresu == 0) = NaN;
        
        cdist1 = sqrt((tempresc-1).^2 + tempresa.^2);
        cdist2 = sqrt((tempresa-1).^2 + tempresc.^2);
        [cdist cI] = min(min(cdist1, cdist2), [], 1);
        cI = cI + (0:n:numel(tempresc)-1);
        [crap cord] = sort(cdist);
        cI10 = cI(cord(1:round(length(cI)/10)));  % top 10 percentile
        
        udist1 = sqrt((tempresu-1).^2 + tempresa.^2);
        udist2 = sqrt((tempresa-1).^2 + tempresu.^2);
        [udist uI] = min(min(udist1, udist2), [], 1);
        uI = uI + (0:n:numel(tempresu)-1);
        [crap uord] = sort(udist);
        uI10 = uI(uord(1:round(length(uI)/10)));  % top 10 percentile
        
%         figure(4*i); hold on;
%         if n>20 && n<60
%             plot(clhom, cdist, 'b.')
%             plot(clhom, udist, 'r.')
%         end
        
        cmean(n, 1, 1) = nanmean(tempresc(cI));
        cmean(n, 1, 2) = nanstd(tempresc(cI));
        cmean(n, 2, 1) = nanmean(tempresa(cI));
        cmean(n, 2, 2) = nanstd(tempresa(cI));
        
        ctop10(n, 1, 1) = nanmean(tempresc(cI10));
        ctop10(n, 1, 2) = nanstd(tempresc(cI10));
        ctop10(n, 2, 1) = nanmean(tempresa(cI10));
        ctop10(n, 2, 2) = nanstd(tempresa(cI10));
        
        umean(n, 1, 1) = nanmean(tempresu(uI));
        umean(n, 1, 2) = nanstd(tempresu(uI));
        umean(n, 2, 1) = nanmean(tempresa(uI));
        umean(n, 2, 2) = nanstd(tempresa(uI));
        
        utop10(n, 1, 1) = nanmean(tempresu(uI10));
        utop10(n, 1, 2) = nanstd(tempresu(uI10));
        utop10(n, 2, 1) = nanmean(tempresa(uI10));
        utop10(n, 2, 2) = nanstd(tempresa(uI10));
        
        % Analog
        temp2 = combran{n};
        prescan = temp2(1:n, :);
        presacan = temp2(n+1:2*n, :);
        presuan = temp2(2*n+1:3*n, :);
        presauan = temp2(3*n+1:end, :);
        
        cdist1 = sqrt((prescan-1).^2 + presacan.^2);
        cdist2 = sqrt((presacan-1).^2 + prescan.^2);
        [cdist cI] = min(min(cdist1, cdist2), [], 1);
        cI = cI + (0:n:numel(prescan)-1);
        [crap cord] = sort(cdist);
        cI10 = cI(cord(1:round(length(cI)/10)));  % top 10 percentile
        
        udist1 = sqrt((presuan-1).^2 + presauan.^2);
        udist2 = sqrt((presauan-1).^2 + presuan.^2);
        [udist uI] = min(min(udist1, udist2), [], 1);
        uI = uI + (0:n:numel(presuan)-1);
        [crap uord] = sort(udist);
        uI10 = uI(uord(1:round(length(uI)/10)));  % top 10 percentile
        
%         figure(5*i); hold on
%         if n>20 && n<60
%             plot(clhom, cdist, 'b.')
%             plot(clhom, udist, 'r.')
%         end

        chtemp = hist3([clhom; cdist]', 'Edges',  {(0:0.01:1), (0:0.01:1)});
        c2dhist(i, :, :) = squeeze(c2dhist(i, :, :)) + chtemp;
        uhtemp = hist3([clhom; udist]', 'Edges',  {(0:0.01:1), (0:0.01:1)});
        u2dhist(i, :, :) = squeeze(u2dhist(i, :, :)) + uhtemp;
        
        cscat = [cscat [clhom; cdist]];
        uscat = [uscat [clhom; udist]];

        ancmean(n, 1, 1) = nanmean(prescan(cI));
        ancmean(n, 1, 2) = nanstd(prescan(cI));
        ancmean(n, 2, 1) = nanmean(presacan(cI));
        ancmean(n, 2, 2) = nanstd(presacan(cI));
        
        anctop10(n, 1, 1) = nanmean(prescan(cI10));
        anctop10(n, 1, 2) = nanstd(prescan(cI10));
        anctop10(n, 2, 1) = nanmean(presacan(cI10));
        anctop10(n, 2, 2) = nanstd(presacan(cI10));
        
        anumean(n, 1, 1) = nanmean(presuan(uI));
        anumean(n, 1, 2) = nanstd(presuan(uI));
        anumean(n, 2, 1) = nanmean(presauan(uI));
        anumean(n, 2, 2) = nanstd(presauan(uI));
        
        anutop10(n, 1, 1) = nanmean(presuan(uI10));
        anutop10(n, 1, 2) = nanstd(presuan(uI10));
        anutop10(n, 2, 1) = nanmean(presauan(uI10));
        anutop10(n, 2, 2) = nanstd(presauan(uI10));
        
    end
    
    figure(1);
    plot(squeeze(anctop10(:, 2, 1)), squeeze(anctop10(:, 1, 1)), 'b.');
    plot(squeeze(anutop10(:, 2, 1)), squeeze(anutop10(:, 1, 1)), 'r.');
    
    figure(3*i);
    plot(1:ncell, squeeze(anctop10(:, 1, 1)), 'b-'); hold on;
    plot(1:ncell, squeeze(anctop10(:, 2, 1)), 'b--');
    plot(1:ncell, squeeze(anutop10(:, 1, 1)), 'r-');
    plot(1:ncell, squeeze(anutop10(:, 2, 1)), 'r--');
    
    figure(2);
    plot(squeeze(anctop10(:, 2, 1)), squeeze(anctop10(:, 1, 1)), 'b.');
    plot(squeeze(anutop10(:, 2, 1)), squeeze(anutop10(:, 1, 1)), 'r.');
    
    figure(4*i); hold on
    plot(cscat(1,:), cscat(2,:), 'b.');
    
    figure(5*i); hold on
    plot(uscat(1,:), uscat(2,:), 'r.');
    
end


