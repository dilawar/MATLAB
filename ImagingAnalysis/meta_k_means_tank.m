% Modified by Ashesh Dhawale to include kmeans ++ seeding 
% and using correlation to merge clusters.


function [allclusters, centroidcorr, dendmem]=meta_k_means_tank(eventsmat, distmetric, corrcut, thr)

if nargin==1;
    distmetric='correlation';
end

expr = eventsmat;
numden = size(expr,1);
numtimes = size(expr, 2);
counts = zeros(numden);
centroidmat = [];

dendevents=sum(eventsmat');
silentdends=find(dendevents<=1);
silentdendcnt=zeros(numden,1);
silentdendcnt(silentdends)=1;

if length(silentdends)>=numden-4;
    allclusters=[];
    dendmem=[];
    centroidcorr=[];
    dunnsinitial=[];
    
else
    expr=zeros(numden-silentdends,numtimes);
    cnt=1;
    for i=1:numden;
        
        if silentdendcnt(i)~=1;
            expr(cnt,:)=eventsmat(i,:);
            cnt=cnt+1;
        end
    end
    numden = size(expr,1);
    numtimes = size(expr, 2);
    counts = zeros(numden);
    centroidmat = [];

    
    % k clusters
    for k = 3
       
        % repeat t times
        for t = 1:2500
            try
                
                [IDX, C, SUMD, D] = kmeans_plusplus(expr, k);

                % update counts
                % returns a matrix showing the number of times each pair of
                % dendrites were clustered together
                for x = 1:numden
                    for y = 1:numden
                        if IDX(x)==IDX(y)
                            counts(x,y) = counts(x,y) + 1;
                        end
                    end
                end

            catch
                IDX=[];
            end

        end
   
        if ~isempty(IDX);
            
            clusters = [];
            for xx = 1:numden
                for yy = 1:xx-1
                    if counts(xx, yy) >= thr*t
                        clusters = [clusters xx];
                        clusters = [clusters yy];
                    end
                end
            end

            % sort list of dendrites
            clusters = sort(clusters);

            % delete the repeated dendrites from the list
            check = clusters(1);

            list = check;

            for c = 2:length(clusters)
                if clusters(c) ~= check
                    list = [list clusters(c)];
                    check = clusters(c);
                end
            end

            % vector of average point to centroid distances for each cluster
            numvect = [];

            cnt=1;
            dendmem=[];
            dendmem=dendmem';
            % find the clusters and print out graded memberships for the dendrites
            for n = 1:length(list)
                item = list(n);
                if item ~= 0
                    members = [item];
                    for nn = 1:numden
                        if nn ~= item
                            if counts(nn, item) >= thr*t
                                members = [members nn];
                                list(find(list == nn)) = 0;
                            end
                        end
                    end
                    eval(['cluster' num2str(cnt) '=members;']);
                    cnt=cnt+1;
                    % find centroid for cluster
                    centroid = zeros(1, numtimes);
                    for m = 1:length(members)
                        centroid = centroid + expr(members(m), :);
                    end
                    centroid = centroid/length(members);

                    centroidmat = [centroidmat; centroid];

                    numerator = 0;
                    % calculate the validity numerator
                    for m = 1:length(members)
                        numerator = numerator + (sum((expr(members(m),:)-centroid).^2))^0.5;
                    end
                    numerator = numerator/length(members);
                    numvect = [numvect numerator];

                    %graded membership using correlations
                    grmem = zeros(numden, 1);

                    for g = 1:numden
                        grmem(g, 1) = corr2(centroid, expr(g, :));
                    end
                    dendmem=[dendmem grmem];
                    grmem = [(1:numden); grmem'];

                end
            end

            allclusters=cell(cnt-1,2);
            for i=1:cnt-1;
                eval(['tempcluster=cluster' num2str(i) ';']);
                allclusters{i,1}=tempcluster;
                allclusters{i,2}=centroidmat(i,:);
            end


            centroidcorr = corrcoef(centroidmat');



%             combclusters=0;

%             [newclusters, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, combclusters);

            endclustering=0;

%             numclusters=10;
            if size(allclusters, 1) < 3
                endclustering = 1;
            end

            while endclustering==0;
                numclusters=size(allclusters,1);
                combineclusters=[];
                for i=1:numclusters-1;
                    for j=(i+1):numclusters;
                        if centroidcorr(i, j)>=corrcut;
                            combineclusters=[i j];
                        end
                    end
                end
                if ~isempty(combineclusters);
%                     [allclusters, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, [combineclusters(1) combineclusters(2)]);
                    newclusters = cell(numclusters-1, 2);
                    lst = (1:numclusters);
                    clist = lst(lst ~= combineclusters(1) & lst ~= combineclusters(2));
                    newclusters(1:numclusters-2, 1) = allclusters(clist, 1);
                    mergecluster = [allclusters{combineclusters(1), 1} allclusters{combineclusters(2), 1}];
                    newclusters{numclusters-1, 1} = mergecluster;
                    
                    centroidmat = [];
                    dendmem = [];
                    dendmem = dendmem';
                    for c = 1 : size(newclusters, 1)
                        temp = newclusters{c, 1};
                        centroid = mean(expr(temp, :), 1);
                        centroidmat = [centroidmat; centroid];
                        newclusters {c, 2} = centroid;

                        grmem = zeros(numden, 1);
                        for g = 1:numden
                            grmem(g, 1) = corr2(centroid, expr(g, :));
                        end
                        dendmem=[dendmem grmem];
                    end
                    
                    centroidcorr = corrcoef(centroidmat');
                    
                    allclusters = newclusters;
                    
                    
                else
                    endclustering=1;
                end
                if numclusters==3;
                    endclustering=1;
                end
            end

            numclusters=size(allclusters,1);

            for zz=1:numclusters;
                clusters=allclusters{zz,1};
                for zzz=1:length(clusters);
                    clusters(zzz)=clusters(zzz)+sum(silentdendcnt(1:clusters(zzz)));
                end
                allclusters{zz,1}=clusters;
            end
            
        else
            allclusters=[];
            centroidcorr=[];
            dendmem=[]; 
            dunnsinitial=0;
        end
        
    end
    
end

