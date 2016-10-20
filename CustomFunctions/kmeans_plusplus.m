% Implementation of kmeans ++ based on Arthur and Vassilvitskii 2007
% Ashesh Dhawale 23/10/11

function [IDX, C, SUMD, D] = kmeans_plusplus(datmat, k)

% rand('seed', 0);

n = size(datmat, 1);
nlist = 1 : n;

% Pick first centre randomly
centres = ceil(rand * n);

% Calculate distances to centre
datdist = 1 - corrcoef(datmat');

for x = 2 : k
%     calculating centre selection prob based on
%     distance of each point to closest centre
    distvect = zeros(n, 1);
    for num = 1 : n
        distvect(num) = min(datdist(centres, num).^2);
    end
    distvect = distvect/sum(distvect);
    
    % Pick next centre based on weighed probabilites
    newcentre = find(rand < cumsum(distvect), 1);
    centres = [centres newcentre];
end

[IDX, C, SUMD, D] = kmeans(datmat, k, 'distance', 'correlation', 'start', datmat(centres,:));

    