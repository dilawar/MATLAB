% Author: Ashesh Dhawale
% Date: 2/11/11

icld = [];      % intra cluster distances
ocld = [];      % inter cluster distances
rcld = [];      % intra cluster distances for rand clusters
acld = [];      % all inter cell distances

iclc = [];      % intra cluster correlations
oclc = [];      % inter cluster correlations
rclc = [];      % intra cluster correlations for rand clusters
aclc = [];      % all inter cell correlations

nn = 8;         % # of neighbours to consider
ncli = [];      % cluster homogeneity index of neigbours aka neighbour identity
rcli = [];      % cluster homogeneity index of neighbours for rand clusters
ccorr = [];     % Correlation of cell to others in its cluster
ccorrn = [];    % Correlation of cell to cells in all other clusters

scale = 0.842105263; % From spatial scale.xls

for i = 1 : 7
    % Load data
    load([num2str(i) '_clusters.mat']);
    load([num2str(i) '_centroids.mat']);
    load([num2str(i) '_corr.mat']);
    load([num2str(i) '_cellmsk.mat']);
    
    IDX = zeros(size(scorr,1), 1);      % another cluster representation
    
    for x = 1 : size(a,1)
        IDX(a{x, 1}) = x;
    end
    
    rIDX = IDX(randperm(length(IDX)));  % shuffled clusters
    
    % Distance matrix
    sdist = zeros(length(IDX));
    for x = 1 : length(IDX)
        for y = 1 : length(IDX)
            cent1 = centroids(x).Centroid;
            cent2 = centroids(y).Centroid;
            sdist(x,y) = sqrt((cent1(1) - cent2(1))^2 + (cent1(2) - cent2(2))^2);
        end
    end
    
    for x = 1 : length(IDX)-1
        for y = x+1 : length(IDX)
            acld = [acld sdist(x, y)];
            aclc = [aclc scorr(x, y)];
        end
    end
    
    ciclc = [];
    coclc = [];
    
    for p = size(a, 1)
        temp = find(IDX == p);
        rtemp = find(rIDX == p);
        otemp = find(IDX ~= p);
        
        % Intracluster distances and correlations
        for q = 1 : length(temp)-1
            for r = q + 1 : length(temp)
                icld = [icld sdist(temp(q), temp(r))];  % Distance
                iclc = [iclc scorr(temp(q), temp(r))];  % Correlation
                ciclc = [ciclc scorr(temp(q), temp(r))];  % Correlation
            end
        end
        
        % Shuffled data
        for q = 1 : length(rtemp)-1
            for r = q + 1 : length(rtemp)
                rcld = [rcld sdist(rtemp(q), rtemp(r))];  % Distance
                rclc = [rclc scorr(rtemp(q), rtemp(r))];  % Correlation
            end
        end
        
        % Intercluster distances and correlations
        for q = 1 : length(temp)
            for r = 1 : length(otemp)
                ocld = [ocld sdist(temp(q), otemp(r))];  % Distance
                oclc = [oclc scorr(temp(q), otemp(r))];  % Correlation
                coclc = [coclc scorr(temp(q), otemp(r))];  % Correlation
            end
        end
        
        % Nearest neighbour cluster identities and correlations
        for q = 1 : length(temp)
            [crap qdist] = sort(sdist(temp(q), :));
            qcli = mean(IDX(qdist(2:nn+1)) == p);
            qrcli = mean(rIDX(qdist(2:nn+1)) == p);
            ncli = [ncli qcli];
            rcli = [rcli qrcli];
            
            qcl = IDX == p; qcl(q) = 0;
            ccorr = [ccorr squeeze(mean(scorr(q, qcl)))];
            qcln = IDX ~= p;
            ccorrn = [ccorrn squeeze(mean(scorr(q, qcln)))];
        end
        
        
        
        
        % Cluster maps
        cellmskc = cellmsk*0;
        for i = 1 : size(a, 1)
            temp = a{i, 1};
            for x = 1 : length(temp)
                cellmskc(cellmsk == (temp(x))) = i;
            end
        end
%         figure; imagesc(cellmskc)       
    end
    
%     figure(1); hold on
%     errorbar([0 1], [mean(ciclc) mean(coclc)], [(std(ciclc)/sqrt(length(ciclc))) (std(ciclc)/sqrt(length(ciclc)))], 'k');
    
end


% Some stats

aclds = acld * scale;
step = 20;

stat = [];
for x = 0 : step : 250
    temp = aclc(aclds>x & aclds<=(x+step));
    stat = [stat; [(x+step/2) mean(temp) (std(temp)/sqrt(length(temp)))]];
end

nstat1 = [];
nstat2 = [];
for x = -1:.25:1
    temp1 = ccorr(ncli == x);
    temp2 = ccorrn(ncli == x);
    nstat1 = [nstat1; [x mean(temp1) (std(temp1)/sqrt(length(temp1)))]];
    nstat2 = [nstat2; [x mean(temp2) (std(temp2)/sqrt(length(temp2)))]];
end