% Filename: doroc.m
% Author: Ashesh Dhawale
% Date: 5/12/11

function [rocurve, tpr, fpr] = doroc(svect, avect)

rocurve = zeros(length(avect), 2);

avect = sort(avect, 'descend');

for x = 1 : length(avect)
    rocurve(x, 1) = sum(avect >= avect(x))/length(avect);
    rocurve(x, 2) = sum(svect >= avect(x))/length(svect);
end

[temp1 ind1] = min(sqrt(rocurve(:,1).^2 + (rocurve(:,2)-1).^2));
[temp2 ind2] = min(sqrt(rocurve(:,2).^2 + (rocurve(:,1)-1).^2));
if temp1 < temp2
   fpr = rocurve(ind1, 1);
   tpr = rocurve(ind1, 2);
else
   fpr = rocurve(ind2, 2); 
   tpr = rocurve(ind2, 1);
   rocurve = rocurve(:, [2 1]);
end